//! Per-base methylation calls and Bismark-style edit-distance summaries.
//!
//! Walks a single record's alignment in genomic orientation and emits four
//! BAM tags:
//!
//! * `XM:Z` -- Bismark methylation call string, derived from the observed
//!   `SEQ` vs the unconverted reference. Same length and orientation as
//!   `SEQ` (genomic). Alphabet: `Z`/`z` (methylated/unmethylated CpG),
//!   `X`/`x` (CHG), `H`/`h` (CHH), `.` (non-cytosine on the read's source
//!   strand, indels, soft-clips, or contexts that fall off the reference).
//!
//! * `YM:Z` -- holodeck-specific, truth-derived methylation call string.
//!   Same shape as `XM:Z` but the methylated/unmethylated decision comes
//!   from the per-haplotype methylation bitmap rather than from the
//!   observed read base. `YM == XM` for runs with zero errors; they
//!   diverge when a sequencing error perturbs a cytosine -- e.g. a
//!   `C → A` error at a methylated CpG yields `XM = '.'` (mismatch, not a
//!   recognised methylation event) but `YM = 'Z'` (truth says
//!   methylated). When invoked via [`populate_pair_call_tags`] (the
//!   standard simulator path), Holodeck only stores CpG methylation
//!   truth, so the closure short-circuits to `false` for CHG/CHH and
//!   those contexts are always lowercase in `YM`. Direct callers of
//!   [`compute_call_tags`] that pass a closure returning `true` for
//!   non-CpG contexts will see uppercase `X`/`H` for those positions —
//!   that's a library-level escape hatch, not the simulator default.
//!
//! * `NM:i` -- Bismark-style edit distance: mismatches against the
//!   *unconverted* reference, with bisulfite-allowed events suppressed
//!   (ref `C` → seq `T` for top-strand reads; ref `G` → seq `A` for
//!   bottom-strand reads). Insertions and deletions count one per base.
//!
//! * `MD:Z` -- Bismark-style match/mismatch description against the
//!   unconverted reference. Bisulfite-allowed mismatches are folded into
//!   the surrounding match runs. Insertions are not represented (per the
//!   MD spec); deletions are emitted as `^<bases>`.

use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record_buf::Cigar;

use crate::haplotype::Haplotype;
use crate::meth::{ConversionType, MethylationAnnotation, MethylationMode, MethylationTable};
use crate::read::SimulatedRead;
use crate::read_naming::TruthAlignment;

/// Methylation context of a cytosine on the read's source strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CpgContext {
    /// `CpG`: cytosine followed by guanine on the same strand.
    Cpg,
    /// `CHG`: cytosine followed by H (A/C/T) followed by guanine.
    Chg,
    /// `CHH`: cytosine followed by H followed by H.
    Chh,
}

impl CpgContext {
    /// Methylated character: uppercase `Z`/`X`/`H`.
    const fn methylated(self) -> u8 {
        match self {
            Self::Cpg => b'Z',
            Self::Chg => b'X',
            Self::Chh => b'H',
        }
    }

    /// Unmethylated character: lowercase `z`/`x`/`h`.
    const fn unmethylated(self) -> u8 {
        match self {
            Self::Cpg => b'z',
            Self::Chg => b'x',
            Self::Chh => b'h',
        }
    }
}

/// Classify the methylation context of the cytosine at `ref_idx` in
/// `ref_bytes`, on the read's source strand. Returns `None` when the base
/// at `ref_idx` is not a cytosine on the requested strand, or when there
/// is not enough flanking reference to classify (e.g. the C is the last
/// base of the contig).
///
/// On the top strand, a cytosine is the literal base `C` and context
/// reads forward (`ref[i+1]`, `ref[i+2]`).
///
/// On the bottom strand, a cytosine appears as `G` on the top strand
/// (`G` is the complement of `C`) and context reads *backward* on the top
/// strand because bottom-strand 5'→3' order is right-to-left in
/// top-strand coordinates: bottom[i+1] = comp(top[i-1]), etc. So a
/// CpG on the bottom strand at top-strand position `i` requires
/// `top[i] == 'G'` (the bottom C) and `top[i-1] == 'C'` (which is the
/// bottom-strand `G` that completes the dinucleotide).
fn classify_context(ref_bytes: &[u8], ref_idx: usize, is_top_strand: bool) -> Option<CpgContext> {
    if is_top_strand {
        let here = ref_bytes.get(ref_idx)?;
        if !here.eq_ignore_ascii_case(&b'C') {
            return None;
        }
        let next1 = *ref_bytes.get(ref_idx + 1)?;
        if next1.eq_ignore_ascii_case(&b'G') {
            return Some(CpgContext::Cpg);
        }
        let next2 = *ref_bytes.get(ref_idx + 2)?;
        if next2.eq_ignore_ascii_case(&b'G') {
            return Some(CpgContext::Chg);
        }
        Some(CpgContext::Chh)
    } else {
        let here = ref_bytes.get(ref_idx)?;
        if !here.eq_ignore_ascii_case(&b'G') {
            return None;
        }
        if ref_idx == 0 {
            return None;
        }
        let prev1 = ref_bytes[ref_idx - 1];
        if prev1.eq_ignore_ascii_case(&b'C') {
            return Some(CpgContext::Cpg);
        }
        if ref_idx < 2 {
            return None;
        }
        let prev2 = ref_bytes[ref_idx - 2];
        if prev2.eq_ignore_ascii_case(&b'C') {
            return Some(CpgContext::Chg);
        }
        Some(CpgContext::Chh)
    }
}

/// Output of [`compute_call_tags`] for a single mate.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CallTags {
    /// Bismark `XM:Z` -- observation-derived per-base methylation call.
    /// Genomic orientation; same length as `SEQ`.
    pub xm: Vec<u8>,
    /// holodeck `YM:Z` -- truth-derived per-base methylation call.
    /// Genomic orientation; same length as `SEQ`.
    pub ym: Vec<u8>,
    /// Bismark `NM:i` -- edit distance against the unconverted reference,
    /// with bisulfite-allowed events suppressed.
    pub nm: u32,
    /// Bismark `MD:Z` -- match/mismatch description against the
    /// unconverted reference, with bisulfite-allowed events suppressed.
    pub md: String,
}

/// Compute the four Bismark-style methylation tags for one BAM record.
///
/// `seq_genomic` is the read sequence in *genomic* orientation (i.e. how
/// it will appear in the BAM `SEQ` field -- reverse-complemented from the
/// FASTQ read for reverse-strand alignments).
///
/// `cigar` is the alignment CIGAR in genomic orientation, which already
/// reflects the leftmost-aligned-position layout of the BAM record.
///
/// `ref_contig` is the full contig sequence as a slice. `ref_aln_start_0`
/// is the 0-based start of the alignment within the contig (i.e.
/// `truth.position - 1`). The walker reads at most a few bases past the
/// alignment end for CHG/CHH lookahead; if the alignment runs against the
/// contig boundary the unclassifiable bases emit `.`.
///
/// `is_top_strand` controls both the cytosine class (top: literal `C`;
/// bottom: literal `G`) and the bisulfite-allowed-mismatch direction
/// (top: ref `C` → seq `T`; bottom: ref `G` → seq `A`).
///
/// `methylation_mode` controls the meaning of an observed `C→T` (top) or
/// `G→A` (bottom) event when computing `XM:Z`:
///
/// * [`MethylationMode::EmSeq`] -- preserved cytosines (seq `C`/`G`) are
///   methylated; converted cytosines (seq `T`/`A`) are unmethylated.
/// * [`MethylationMode::Taps`] -- the inverse: preserved cytosines are
///   unmethylated; converted cytosines are methylated.
///
/// `is_methylated_at_ref` is a closure that, given a 0-based reference
/// position, returns whether that position's read-strand cytosine is
/// methylated according to the haplotype's truth bitmap. Always called
/// with positions where [`classify_context`] returned `Some`. Holodeck's
/// methylation model only stores state for CpG sites, so the closure
/// can short-circuit to `false` for CHG/CHH inputs.
///
/// # Panics
///
/// Panics if the CIGAR's read-consumed length disagrees with
/// `seq_genomic.len()`: every `Match`/`SequenceMatch`/`SequenceMismatch`/
/// `Insertion`/`SoftClip` op increments `seq_idx` and the walker then
/// indexes `seq_genomic[seq_idx + k]` directly, so an over-consuming CIGAR
/// is an unconditional out-of-bounds in both debug and release builds.
/// This is an internal-invariant panic: the simulator constructs CIGARs
/// from its own `Fragment` walker, so a mismatch indicates a holodeck bug
/// rather than malformed external input.
#[must_use]
#[allow(clippy::too_many_lines, reason = "single-pass walker handles every CIGAR op")]
pub fn compute_call_tags(
    seq_genomic: &[u8],
    cigar: &Cigar,
    ref_contig: &[u8],
    ref_aln_start_0: u32,
    is_top_strand: bool,
    methylation_mode: MethylationMode,
    is_methylated_at_ref: impl Fn(u32) -> bool,
) -> CallTags {
    let mut xm: Vec<u8> = Vec::with_capacity(seq_genomic.len());
    let mut ym: Vec<u8> = Vec::with_capacity(seq_genomic.len());
    let mut nm: u32 = 0;
    let mut md = String::new();
    let mut md_match_run: u32 = 0;

    let mut seq_idx: usize = 0;
    let mut ref_idx: usize = ref_aln_start_0 as usize;

    let target_ref_base: u8 = if is_top_strand { b'C' } else { b'G' };
    let unmeth_seq_base: u8 = if is_top_strand { b'T' } else { b'A' };

    for op in cigar.as_ref() {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for k in 0..len {
                    let s = seq_genomic[seq_idx + k];
                    let r = ref_contig.get(ref_idx + k).copied().unwrap_or(b'N');

                    let context = classify_context(ref_contig, ref_idx + k, is_top_strand);

                    // XM (observation-based). Under em-seq, a preserved
                    // cytosine (seq matches `target_ref_base`) is methylated
                    // and a converted cytosine (seq matches `unmeth_seq_base`)
                    // is unmethylated. TAPS inverts this: a converted
                    // cytosine signals methylation.
                    let xm_char = match context {
                        Some(ctx) => {
                            let (preserved_call, converted_call) = match methylation_mode {
                                MethylationMode::EmSeq => (ctx.methylated(), ctx.unmethylated()),
                                MethylationMode::Taps => (ctx.unmethylated(), ctx.methylated()),
                            };
                            let s_up = s.to_ascii_uppercase();
                            if s_up == target_ref_base {
                                preserved_call
                            } else if s_up == unmeth_seq_base {
                                converted_call
                            } else {
                                b'.'
                            }
                        }
                        None => b'.',
                    };
                    xm.push(xm_char);

                    // YM (truth-based)
                    let ym_char = match context {
                        Some(ctx) => {
                            #[expect(
                                clippy::cast_possible_truncation,
                                reason = "ref_idx + k bounded by contig length"
                            )]
                            let abs_ref_pos = (ref_idx + k) as u32;
                            if is_methylated_at_ref(abs_ref_pos) {
                                ctx.methylated()
                            } else {
                                ctx.unmethylated()
                            }
                        }
                        None => b'.',
                    };
                    ym.push(ym_char);

                    // NM/MD: bisulfite-allowed mismatch?
                    let r_up = r.to_ascii_uppercase();
                    let s_up = s.to_ascii_uppercase();
                    let is_match = r_up == s_up;
                    let is_bs_allowed = r_up == target_ref_base && s_up == unmeth_seq_base;

                    if is_match || is_bs_allowed {
                        md_match_run += 1;
                    } else {
                        md.push_str(&md_match_run.to_string());
                        md.push(r_up as char);
                        md_match_run = 0;
                        nm += 1;
                    }
                }
                seq_idx += len;
                ref_idx += len;
            }
            Kind::Insertion => {
                for _ in 0..len {
                    xm.push(b'.');
                    ym.push(b'.');
                }
                #[expect(clippy::cast_possible_truncation, reason = "CIGAR op lens fit in u32")]
                let len_u32 = len as u32;
                nm += len_u32;
                seq_idx += len;
                // No MD emission for insertions per spec.
            }
            Kind::Deletion => {
                md.push_str(&md_match_run.to_string());
                md.push('^');
                for k in 0..len {
                    let r = ref_contig.get(ref_idx + k).copied().unwrap_or(b'N');
                    md.push(r.to_ascii_uppercase() as char);
                }
                md_match_run = 0;
                #[expect(clippy::cast_possible_truncation, reason = "CIGAR op lens fit in u32")]
                let len_u32 = len as u32;
                nm += len_u32;
                ref_idx += len;
            }
            Kind::SoftClip => {
                for _ in 0..len {
                    xm.push(b'.');
                    ym.push(b'.');
                }
                seq_idx += len;
                // Soft clip does not consume reference and does not
                // contribute to NM/MD.
            }
            Kind::HardClip | Kind::Pad => {
                // Hard clip and pad consume neither read nor reference.
            }
            Kind::Skip => {
                // N (ref skip / intron): consumes reference but neither
                // read nor edit distance.
                ref_idx += len;
            }
        }
    }

    md.push_str(&md_match_run.to_string());

    debug_assert_eq!(seq_idx, seq_genomic.len(), "CIGAR did not consume exactly the read length");

    CallTags { xm, ym, nm, md }
}

/// Compute call tags for one mate and stash them on the annotation.
///
/// The walker needs `seq` in *genomic* orientation. The pipeline keeps
/// `read.bases` in FASTQ (read 5'→3') orientation, so reverse-complement
/// when the truth alignment is on the reverse strand.
///
/// Methylation truth lookup goes ref-pos → hap-pos via the haplotype's
/// coordinate map, then queries the per-strand methylation bitmap with
/// `is_negative_strand = !is_top_strand` (the table API takes "negative
/// strand" rather than "top/bottom"). The genome-strand of the read,
/// expressed as `is_top_strand`, comes from the fragment's
/// [`ConversionType`] (`Ct → top`, `Ga → bottom`).
#[allow(clippy::too_many_arguments)]
fn compute_one_mate(
    read: &SimulatedRead,
    truth: &TruthAlignment,
    cigar: &Cigar,
    ref_contig: &[u8],
    methylation_table: &MethylationTable,
    haplotype: &Haplotype,
    conversion_type: ConversionType,
    methylation_mode: MethylationMode,
) -> CallTags {
    let is_top_strand = matches!(conversion_type, ConversionType::Ct);
    let is_negative_strand = !is_top_strand;

    // Convert FASTQ-orientation bases to genomic orientation if the
    // record is on the reverse strand.
    let seq_genomic: Vec<u8> = if truth.is_forward {
        read.bases.clone()
    } else {
        let mut tmp = read.bases.clone();
        crate::fragment::reverse_complement(&mut tmp);
        tmp
    };

    // The truth alignment stores position as 1-based; the walker takes
    // 0-based. Saturating-sub guards against the (impossible) zero.
    let ref_aln_start_0 = truth.position.saturating_sub(1);

    compute_call_tags(
        &seq_genomic,
        cigar,
        ref_contig,
        ref_aln_start_0,
        is_top_strand,
        methylation_mode,
        |ref_pos| {
            let hap_pos = haplotype.hap_position_for(ref_pos);
            methylation_table.is_methylated(hap_pos, is_negative_strand)
        },
    )
}

/// Populate the `r1_call_tags` and `r2_call_tags` fields of the given
/// [`MethylationAnnotation`] for one read pair.
///
/// `methylation_mode` selects the chemistry semantics threaded into
/// `XM:Z`; see [`compute_call_tags`].
///
/// Idempotent for a given pair — overwrites whatever was there.
#[allow(clippy::too_many_arguments)]
pub fn populate_pair_call_tags(
    read1: &SimulatedRead,
    read2: Option<&SimulatedRead>,
    r1_truth: &TruthAlignment,
    r2_truth: Option<&TruthAlignment>,
    r1_cigar: &Cigar,
    r2_cigar: Option<&Cigar>,
    ref_contig: &[u8],
    methylation_table: &MethylationTable,
    haplotype: &Haplotype,
    methylation_mode: MethylationMode,
    annotation: &mut MethylationAnnotation,
) {
    let conv = annotation.conversion_type;

    annotation.r1_call_tags = Some(compute_one_mate(
        read1,
        r1_truth,
        r1_cigar,
        ref_contig,
        methylation_table,
        haplotype,
        conv,
        methylation_mode,
    ));

    annotation.r2_call_tags = if let (Some(r2), Some(r2t), Some(r2c)) = (read2, r2_truth, r2_cigar)
    {
        Some(compute_one_mate(
            r2,
            r2t,
            r2c,
            ref_contig,
            methylation_table,
            haplotype,
            conv,
            methylation_mode,
        ))
    } else {
        None
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::record::cigar::op::Op;

    fn cigar(ops: &[(Kind, usize)]) -> Cigar {
        Cigar::from(ops.iter().map(|&(k, n)| Op::new(k, n)).collect::<Vec<_>>())
    }

    fn meth_none(_: u32) -> bool {
        false
    }

    fn meth_all(_: u32) -> bool {
        true
    }

    #[test]
    fn test_classify_context_top_cpg() {
        // ACGT: 'C' at index 1 followed by 'G' at index 2 → CpG.
        assert_eq!(classify_context(b"ACGT", 1, true), Some(CpgContext::Cpg));
    }

    #[test]
    fn test_classify_context_top_chg() {
        // CAG: 'C' at 0, 'A' at 1, 'G' at 2 → CHG.
        assert_eq!(classify_context(b"CAG", 0, true), Some(CpgContext::Chg));
    }

    #[test]
    fn test_classify_context_top_chh() {
        // CAA: 'C' at 0, 'A' at 1, 'A' at 2 → CHH.
        assert_eq!(classify_context(b"CAA", 0, true), Some(CpgContext::Chh));
    }

    #[test]
    fn test_classify_context_top_not_c() {
        assert_eq!(classify_context(b"AGTA", 0, true), None);
        assert_eq!(classify_context(b"AGTA", 1, true), None);
    }

    #[test]
    fn test_classify_context_top_no_lookahead() {
        // 'C' at last position → can't classify.
        assert_eq!(classify_context(b"AAC", 2, true), None);
    }

    #[test]
    fn test_classify_context_bottom_cpg() {
        // Top "ACGT": bottom-strand C is at top-pos 2 ('G'); top[1]='C'
        // satisfies bottom-strand CpG context.
        assert_eq!(classify_context(b"ACGT", 2, false), Some(CpgContext::Cpg));
    }

    #[test]
    fn test_classify_context_bottom_chg() {
        // Top "CAG": bottom-strand C is at top-pos 2 ('G'); top[1]='A' ≠ 'C',
        // top[0]='C' → CHG (bottom 5'→3': G[at top 2], A→T, C→G = CHG... wait
        // let me re-derive: bottom-strand reading 5'→3' from top-pos 2:
        // bottom[0]=comp(top[2])=comp(G)=C, bottom[1]=comp(top[1])=comp(A)=T,
        // bottom[2]=comp(top[0])=comp(C)=G → "CTG" which is C-T-G = CHG ✓
        assert_eq!(classify_context(b"CAG", 2, false), Some(CpgContext::Chg));
    }

    #[test]
    fn test_classify_context_bottom_chh() {
        // Top "TAG": bottom-strand C at top-pos 2 ('G'); reading bottom 5'→3':
        // C-T-A → CHH.
        assert_eq!(classify_context(b"TAG", 2, false), Some(CpgContext::Chh));
    }

    #[test]
    fn test_classify_context_bottom_at_pos_zero_returns_none() {
        // 'G' at index 0 has no top[i-1] to read, so context is unclassifiable.
        assert_eq!(classify_context(b"GAT", 0, false), None);
    }

    #[test]
    fn test_simple_match_full_methylation_top() {
        // Top-strand read against reference "ACGT": 4 bases, perfect match.
        // Position 1 is a CpG-context C (truth: methylated → 'Z').
        // Position 2 is a 'G' (top-strand context: not C → '.').
        let seq = b"ACGT";
        let cig = cigar(&[(Kind::Match, 4)]);
        let tags = compute_call_tags(seq, &cig, b"ACGT", 0, true, MethylationMode::EmSeq, meth_all);
        assert_eq!(tags.xm, b".Z..".to_vec());
        assert_eq!(tags.ym, b".Z..".to_vec());
        assert_eq!(tags.nm, 0);
        assert_eq!(tags.md, "4");
    }

    #[test]
    fn test_simple_match_zero_methylation_with_conversion_top() {
        // Top-strand read against ref "ACGT" with the CpG C unmethylated and
        // converted to T (em-seq-style): seq = "ATGT".
        // - position 0: A, not a C in any orientation → '.'
        // - position 1: ref C (CpG), seq T → unmethylated 'z'.
        // - position 2: ref G, no top-strand context → '.'
        // - position 3: T, not a C → '.'
        // The C→T at position 1 is bisulfite-allowed → not counted in NM,
        // and the MD walk treats it as a match.
        let seq = b"ATGT";
        let cig = cigar(&[(Kind::Match, 4)]);
        let tags =
            compute_call_tags(seq, &cig, b"ACGT", 0, true, MethylationMode::EmSeq, meth_none);
        assert_eq!(tags.xm, b".z..".to_vec());
        assert_eq!(tags.ym, b".z..".to_vec());
        assert_eq!(tags.nm, 0);
        assert_eq!(tags.md, "4");
    }

    #[test]
    fn test_taps_methylated_converted_top() {
        // TAPS top-strand: ref "ACGT" with the CpG C methylated, chemistry
        // converts methylated cytosines to thymine: seq = "ATGT".
        // - position 1: ref C (CpG), seq T → methylated 'Z' under TAPS
        //   (converted base signals methylation).
        // - YM at position 1 is 'Z' (truth: methylated).
        // The C→T at position 1 is the chemistry event, not a real mismatch.
        let seq = b"ATGT";
        let cig = cigar(&[(Kind::Match, 4)]);
        let tags = compute_call_tags(seq, &cig, b"ACGT", 0, true, MethylationMode::Taps, meth_all);
        assert_eq!(tags.xm, b".Z..".to_vec());
        assert_eq!(tags.ym, b".Z..".to_vec());
        assert_eq!(tags.nm, 0);
        assert_eq!(tags.md, "4");
    }

    #[test]
    fn test_taps_unmethylated_preserved_top() {
        // TAPS top-strand: ref "ACGT" with the CpG C unmethylated, chemistry
        // leaves unmethylated cytosines untouched: seq = "ACGT".
        // - position 1: ref C (CpG), seq C → unmethylated 'z' under TAPS
        //   (preserved base signals unmethylation).
        let seq = b"ACGT";
        let cig = cigar(&[(Kind::Match, 4)]);
        let tags = compute_call_tags(seq, &cig, b"ACGT", 0, true, MethylationMode::Taps, meth_none);
        assert_eq!(tags.xm, b".z..".to_vec());
        assert_eq!(tags.ym, b".z..".to_vec());
        assert_eq!(tags.nm, 0);
        assert_eq!(tags.md, "4");
    }

    #[test]
    fn test_taps_xm_eq_ym_under_zero_errors() {
        // Invariant: under zero sequencing errors XM == YM for both
        // chemistries. Mix methylated and unmethylated CpG sites (closure
        // says: pos 1 methylated, pos 4 unmethylated). Ref "ACGTACGT" has
        // CpGs at top-pos 1 and 5.
        let cig = cigar(&[(Kind::Match, 8)]);
        let truth = |p: u32| p == 1; // CpG@1 methylated, CpG@5 unmethylated.

        // em-seq: meth → C preserved, unmeth → T converted.
        let em_seq = b"ACGTATGT";
        let em_tags =
            compute_call_tags(em_seq, &cig, b"ACGTACGT", 0, true, MethylationMode::EmSeq, truth);
        assert_eq!(em_tags.xm, em_tags.ym, "XM must equal YM under zero errors (em-seq)");

        // TAPS: meth → T converted, unmeth → C preserved.
        let taps_seq = b"ATGTACGT";
        let taps_tags =
            compute_call_tags(taps_seq, &cig, b"ACGTACGT", 0, true, MethylationMode::Taps, truth);
        assert_eq!(taps_tags.xm, taps_tags.ym, "XM must equal YM under zero errors (TAPS)");
    }

    #[test]
    fn test_real_mismatch_increments_nm() {
        // Top-strand: ref "ACGT", seq "AAGT" — a real A→A is matched at 0,
        // but ref[1]='C' vs seq[1]='A' is a non-bisulfite mismatch.
        let seq = b"AAGT";
        let cig = cigar(&[(Kind::Match, 4)]);
        let tags =
            compute_call_tags(seq, &cig, b"ACGT", 0, true, MethylationMode::EmSeq, meth_none);
        assert_eq!(tags.nm, 1);
        assert_eq!(tags.md, "1C2");
        // XM at position 1: ref is C (CpG context) but seq is A — neither
        // C (methylated) nor T (unmethylated) → '.'.
        assert_eq!(tags.xm, b"....".to_vec());
    }

    #[test]
    fn test_ym_diverges_from_xm_under_error_at_methylated_cpg() {
        // Ref "ACGT", CpG C at position 1 is methylated (truth). An error
        // converts it to A (post-error SEQ = "AAGT"). XM reports '.' at
        // position 1 (mismatch); YM reports 'Z' (truth says methylated).
        let seq = b"AAGT";
        let cig = cigar(&[(Kind::Match, 4)]);
        let tags =
            compute_call_tags(seq, &cig, b"ACGT", 0, true, MethylationMode::EmSeq, |p| p == 1);
        assert_eq!(tags.xm, b"....".to_vec());
        assert_eq!(tags.ym, b".Z..".to_vec());
        assert_eq!(tags.nm, 1);
    }

    #[test]
    fn test_bottom_strand_ga_is_bisulfite_allowed() {
        // Bottom-strand read aligned to ref "ACGT": bottom-strand C is at
        // top-pos 2 ('G'). Under em-seq the bottom C unmethylated would
        // bisulfite-convert: the seq base in genomic orientation becomes A
        // (G → A). So seq = "ACAT".
        // - position 1: ref C, but on bottom strand we don't classify
        //   top-strand C as a target → XM='.'.
        // - position 2: ref G, bottom strand CpG context → 'z' (unmeth).
        // The G→A at position 2 is bisulfite-allowed for is_top=false.
        let seq = b"ACAT";
        let cig = cigar(&[(Kind::Match, 4)]);
        let tags =
            compute_call_tags(seq, &cig, b"ACGT", 0, false, MethylationMode::EmSeq, meth_none);
        assert_eq!(tags.xm, b"..z.".to_vec());
        assert_eq!(tags.nm, 0);
        assert_eq!(tags.md, "4");
    }

    #[test]
    fn test_taps_bottom_strand_ga_signals_methylation() {
        // TAPS bottom-strand: ref "ACGT", bottom CpG C at top-pos 2 ('G')
        // methylated → chemistry converts G→A. seq = "ACAT".
        // - position 2: ref G, bottom-strand CpG context, seq A → 'Z'
        //   under TAPS (converted base signals methylation).
        let seq = b"ACAT";
        let cig = cigar(&[(Kind::Match, 4)]);
        let tags = compute_call_tags(seq, &cig, b"ACGT", 0, false, MethylationMode::Taps, meth_all);
        assert_eq!(tags.xm, b"..Z.".to_vec());
        assert_eq!(tags.ym, b"..Z.".to_vec());
        assert_eq!(tags.nm, 0);
        assert_eq!(tags.md, "4");
    }

    #[test]
    fn test_insertion_increments_nm_and_emits_dots() {
        // Ref "ACGT", read "ACXGT" with an inserted X between positions 1
        // and 2. CIGAR: 2M 1I 2M.
        let seq = b"ACXGT";
        let cig = cigar(&[(Kind::Match, 2), (Kind::Insertion, 1), (Kind::Match, 2)]);
        let tags = compute_call_tags(seq, &cig, b"ACGT", 0, true, MethylationMode::EmSeq, meth_all);
        // XM/YM: position 0 '.', position 1 'Z' (CpG meth), position 2 '.'
        // (insertion), positions 3,4 '.','.'.
        assert_eq!(tags.xm, b".Z...".to_vec());
        assert_eq!(tags.ym, b".Z...".to_vec());
        // NM: 1 inserted base.
        assert_eq!(tags.nm, 1);
        // MD: 4 matches, no representation of the insertion.
        assert_eq!(tags.md, "4");
    }

    #[test]
    fn test_deletion_increments_nm_and_emits_caret() {
        // Ref "ACGTA", read "ACTA" with the G at ref-pos 2 deleted. CIGAR:
        // 2M 1D 2M.
        let seq = b"ACTA";
        let cig = cigar(&[(Kind::Match, 2), (Kind::Deletion, 1), (Kind::Match, 2)]);
        let tags =
            compute_call_tags(seq, &cig, b"ACGTA", 0, true, MethylationMode::EmSeq, meth_none);
        // XM positions: 0 '.', 1 'z' (CpG, unmeth — context "CG" still
        // intact in ref slice even though read is missing the G), 2 '.', 3
        // '.'. (Note: position 1 in the ref is C; ref[2]='G' so context is
        // CpG.)
        assert_eq!(tags.xm.len(), 4);
        // NM: 1 deleted base.
        assert_eq!(tags.nm, 1);
        // MD: "2^G2".
        assert_eq!(tags.md, "2^G2");
    }

    #[test]
    fn test_soft_clip_does_not_count() {
        // Ref "ACGT", read "XACGTY" with 1S clip on each end. CIGAR:
        // 1S 4M 1S.
        let seq = b"XACGTY";
        let cig = cigar(&[(Kind::SoftClip, 1), (Kind::Match, 4), (Kind::SoftClip, 1)]);
        let tags = compute_call_tags(seq, &cig, b"ACGT", 0, true, MethylationMode::EmSeq, meth_all);
        // XM: '.', '.', 'Z', '.', '.', '.' (soft-clipped ends + the 4M block
        // with CpG at pos 1 of ref).
        assert_eq!(tags.xm, b"..Z...".to_vec());
        assert_eq!(tags.nm, 0);
        assert_eq!(tags.md, "4");
    }

    #[test]
    fn test_chg_chh_always_lowercase_in_ym() {
        // Ref "CAG" — CHG context at position 0. Truth says methylated; YM
        // would still be 'X' (uppercase methylated CHG) -- but holodeck's
        // model means the "truth" closure should never return true for
        // non-CpG sites. We still verify the *walker* honours its closure
        // (the policy of returning false for non-CpG lives upstream).
        // This test pins the walker's behaviour: when given a methylated
        // CHG, it emits 'X'. The CpG-only invariant is enforced by the
        // upstream closure builder.
        let seq = b"CAG";
        let cig = cigar(&[(Kind::Match, 3)]);
        let tags = compute_call_tags(seq, &cig, b"CAG", 0, true, MethylationMode::EmSeq, meth_all);
        assert_eq!(tags.xm, b"X..".to_vec());
        assert_eq!(tags.ym, b"X..".to_vec());
    }

    #[test]
    fn test_multiple_mismatches_in_md() {
        // Ref "AAAA", seq "ATTA" — two real mismatches at positions 1 and 2.
        let seq = b"ATTA";
        let cig = cigar(&[(Kind::Match, 4)]);
        let tags =
            compute_call_tags(seq, &cig, b"AAAA", 0, true, MethylationMode::EmSeq, meth_none);
        assert_eq!(tags.nm, 2);
        assert_eq!(tags.md, "1A0A1");
    }

    #[test]
    fn test_populate_pair_call_tags_clears_stale_r2_when_se() {
        // Regression: `populate_pair_call_tags` documents itself as
        // "overwrites whatever was there", but the SE branch (read2 is None)
        // previously left a stale r2_call_tags untouched.
        use crate::meth::{ConversionType, MethylationAnnotation, MethylationTable};
        use crate::read::SimulatedRead;
        use crate::read_naming::TruthAlignment;
        use rand::SeedableRng;
        use rand::rngs::SmallRng;

        let reference: &[u8] = b"ACGT";
        let hap = crate::haplotype::build_haplotypes(&[], 1, &mut SmallRng::seed_from_u64(0))
            .into_iter()
            .next()
            .unwrap();
        let table = MethylationTable::empty(reference.len());

        let read1 = SimulatedRead {
            name: "r1".to_string(),
            bases: b"ACGT".to_vec(),
            qualities: vec![b'I'; 4],
        };
        let r1_truth = TruthAlignment {
            contig: "chr1".to_string(),
            position: 1,
            is_forward: true,
            haplotype: 0,
            fragment_length: 4,
            n_errors: 0,
        };
        let r1_cigar = cigar(&[(Kind::Match, 4)]);

        let stale = CallTags {
            xm: b"STALE".to_vec(),
            ym: b"STALE".to_vec(),
            nm: 99,
            md: "stale".to_string(),
        };
        let mut anno = MethylationAnnotation {
            conversion_type: ConversionType::Ct,
            r1_pre_conversion_bases: None,
            r2_pre_conversion_bases: None,
            r1_call_tags: None,
            r2_call_tags: Some(stale),
        };

        populate_pair_call_tags(
            &read1,
            None,
            &r1_truth,
            None,
            &r1_cigar,
            None,
            reference,
            &table,
            &hap,
            crate::meth::MethylationMode::EmSeq,
            &mut anno,
        );

        assert!(anno.r1_call_tags.is_some(), "R1 tags must be populated");
        assert!(anno.r2_call_tags.is_none(), "R2 tags must be cleared in SE branch");
    }

    #[test]
    fn test_alignment_offset_into_contig() {
        // Place the alignment at ref-pos 5 in a longer contig. Verifies
        // that ref_aln_start_0 is honoured for both reference walks and
        // YM truth lookups.
        let contig: &[u8] = b"NNNNNACGTNNNN";
        // Alignment starts at index 5 (the 'A'), 4 bases long.
        let seq = b"ACGT";
        let cig = cigar(&[(Kind::Match, 4)]);
        let tags =
            compute_call_tags(seq, &cig, contig, 5, true, MethylationMode::EmSeq, |p| p == 6);
        // CpG C is at ref-pos 6 → truth methylated.
        assert_eq!(tags.xm, b".Z..".to_vec());
        assert_eq!(tags.ym, b".Z..".to_vec());
        assert_eq!(tags.nm, 0);
        assert_eq!(tags.md, "4");
    }
}
