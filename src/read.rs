//! Simulated read pair generation.
//!
//! Combines fragment extraction, error model application, and read naming
//! into complete [`ReadPair`] objects ready for FASTQ output. Computes
//! proper CIGARs from the haplotype-to-reference coordinate mapping for
//! golden BAM output.

use noodles::sam::alignment::record::cigar::op::{Kind, Op};
use noodles::sam::alignment::record_buf::Cigar;
use rand::Rng;

use crate::error_model::{self, ErrorModel, ReadEnd};
use crate::fragment::{
    Fragment, extract_read_bases, lowercase_fraction, reverse_complement, uppercase_in_place,
};
use crate::meth::{
    ConversionType, MethylationAnnotation, MethylationConfig, apply_methylation_conversion,
};
use crate::read_naming::{
    ReadNaming, TruthAlignment, encoded_pe_name, encoded_se_name, simple_name,
};

/// A single simulated read with bases, quality scores, and metadata.
#[derive(Debug, Clone)]
pub struct SimulatedRead {
    /// Read name (shared between R1 and R2 of a pair).
    pub name: String,
    /// Base sequence (possibly with errors applied).
    pub bases: Vec<u8>,
    /// Quality scores (Phred+33 encoded).
    pub qualities: Vec<u8>,
}

/// A simulated read pair (or single read for SE mode).
#[derive(Debug)]
pub struct ReadPair {
    /// First read (always present).
    pub read1: SimulatedRead,
    /// Second read (present only for paired-end mode).
    pub read2: Option<SimulatedRead>,
    /// Truth alignment for R1.
    pub r1_truth: TruthAlignment,
    /// Truth alignment for R2 (present only for paired-end mode).
    pub r2_truth: Option<TruthAlignment>,
    /// Truth CIGAR for R1, reflecting haplotype variants and adapter
    /// soft-clipping.
    pub r1_cigar: Cigar,
    /// Truth CIGAR for R2. `None` for single-end reads.
    pub r2_cigar: Option<Cigar>,
    /// Methylation annotation, populated iff methylation simulation was
    /// enabled for this pair. Carries the conversion type and
    /// pre-conversion bases for the Bismark-compatible golden-BAM tag set
    /// (`XG`, `XR`, `XM`, `YM`, `YS`, `NM`, `MD`).
    pub methylation: Option<crate::meth::MethylationAnnotation>,
}

/// Output of building one mate of a read pair.
struct MateOutput {
    /// Final post-conversion, post-error bases.
    bases: Vec<u8>,
    /// Phred+33 qualities, length equal to `bases`.
    qualities: Vec<u8>,
    /// Pre-conversion bases (5'→3' read orientation). `None` when methylation
    /// chemistry was not applied or capture was not requested.
    pre_conversion: Option<Vec<u8>>,
    /// CIGAR including any adapter soft-clip.
    cigar: Cigar,
    /// Truth alignment for this mate.
    truth: TruthAlignment,
}

/// Apply BS chemistry to a fragment in source-strand orientation, returning
/// the chemistry-applied bases in TOP-strand orientation (to match
/// `Fragment::bases`).
///
/// In a directional library both R1 and R2 derive from the same source
/// strand: R1 reads the BS-converted source 5'→3', R2 reads the
/// PCR-synthesized complement (= revcomp of source 5'→3'). So chemistry
/// must be applied once, at fragment scale, on the source strand — not
/// independently per mate.
///
/// `is_forward` selects the source: top for CT (top-strand-derived
/// fragments) or bottom for GA. The methylation table is indexed by
/// top-strand haplotype coordinates; for GA fragments we revcomp the
/// fragment into bottom orientation, run [`apply_methylation_conversion`]
/// with `is_negative_strand = true` (which selects the bottom bitmap and
/// reverses the index), then revcomp back to top orientation before
/// returning.
///
/// Returns the post-chemistry top-strand bases together with whether the
/// molecule was drawn as a conversion failure (a molecule property, the same
/// for both mates).
fn apply_fragment_chemistry(
    pre_chem_top: &[u8],
    hap_start: u32,
    is_forward: bool,
    haplotype_index: usize,
    config: &MethylationConfig<'_>,
    rng: &mut impl Rng,
) -> (Vec<u8>, bool) {
    let mut bases = pre_chem_top.to_vec();
    if !is_forward {
        reverse_complement(&mut bases);
    }
    let n = bases.len();
    let conversion_failed = apply_methylation_conversion(
        &mut bases,
        n,
        !is_forward,
        hap_start,
        haplotype_index,
        config,
        rng,
    );
    if !is_forward {
        reverse_complement(&mut bases);
    }
    (bases, conversion_failed)
}

/// Build one mate (R1 or R2) from already-extracted bases.
///
/// `bases` is post-chemistry, post-uppercase, in read 5'→3' orientation
/// (the same orientation [`extract_read_bases`] produces). The caller is
/// responsible for running the ambiguity-fraction filter before invoking
/// this function: errors mutate `bases` and advance `rng`, so rejection
/// must happen first to avoid desynchronising the RNG stream.
///
/// Pipeline (per-mate): apply_errors → CIGAR + truth.
/// Chemistry runs once at fragment scale upstream — see
/// [`apply_fragment_chemistry`].
#[allow(clippy::too_many_arguments)]
fn build_mate(
    fragment: &Fragment,
    contig_name: &str,
    end: ReadEnd,
    is_negative_strand: bool,
    mut bases: Vec<u8>,
    pre_conversion: Option<Vec<u8>>,
    adapter_bases: usize,
    model: &impl ErrorModel,
    rng: &mut impl Rng,
) -> MateOutput {
    let frag_len = fragment.bases.len();
    let genomic = frag_len.min(bases.len());
    let right_start = frag_len.saturating_sub(genomic);

    // Fragment ref_positions are always in ascending (forward) reference
    // order; negative-strand reads take from the right end of the fragment.
    let positions = if is_negative_strand {
        &fragment.ref_positions[right_start..frag_len]
    } else {
        &fragment.ref_positions[..genomic]
    };

    let (n_errors, qualities) = error_model::apply_errors(model, &mut bases, end, rng);

    let cigar = cigar_from_ref_positions(positions, adapter_bases, is_negative_strand);

    let ref_pos = if fragment.ref_positions.is_empty() {
        0
    } else if is_negative_strand {
        fragment.ref_positions[right_start] + 1
    } else {
        fragment.ref_positions[0] + 1
    };

    #[expect(clippy::cast_possible_truncation, reason = "fragment length fits in u32")]
    let fragment_length = frag_len as u32;

    let truth = TruthAlignment {
        contig: contig_name.to_string(),
        position: ref_pos,
        is_forward: match end {
            ReadEnd::Read1 => fragment.is_forward,
            ReadEnd::Read2 => !fragment.is_forward,
        },
        haplotype: fragment.haplotype_index,
        fragment_length,
        n_errors,
    };

    MateOutput { bases, qualities, pre_conversion, cigar, truth }
}

/// Generate a simulated read pair from a fragment.
///
/// Extracts R1 and optionally R2 bases from the fragment, applies the error
/// model, constructs read names with truth information, and computes CIGARs
/// from the fragment's reference coordinate mapping.
///
/// Returns `None` if either R1 or R2 has a lowercase-base fraction greater
/// than `max_n_frac` (i.e. too many bases came from ambiguity-resolved
/// reference positions); the caller should resample. See the [`fragment`]
/// module documentation for the lowercase-marker convention.
///
/// [`fragment`]: crate::fragment
///
/// # Arguments
/// * `fragment` — The source fragment with bases and reference positions.
/// * `contig_name` — Contig name for read naming.
/// * `read_num` — 1-based read pair number.
/// * `read_length` — Desired read length.
/// * `paired` — Whether to generate both R1 and R2.
/// * `adapter_r1` — Adapter sequence for R1.
/// * `adapter_r2` — Adapter sequence for R2.
/// * `max_n_frac` — Reject the pair if R1 or R2 has a lowercase fraction
///   exceeding this threshold. Use `1.0` to disable.
/// * `error_model` — Error model to apply.
/// * `naming` — Read-naming strategy (encoded truth, simple sequential,
///   or Illumina-style with a shared header).
/// * `methylation` — Optional methylation chemistry configuration. When
///   `Some`, the qualifying class of cytosines (unmethylated under em-seq,
///   methylated under TAPS) in the genomic portion of each read is
///   converted to T with probability `conversion_rate` before the error
///   model is applied.
/// * `capture_pre_conversion` — When `true` AND `methylation` is `Some`,
///   capture the pre-conversion bases of each mate for the `YS:Z` golden-BAM
///   tag. Has no effect when `methylation` is `None`.
/// * `rng` — Random number generator.
///
#[allow(clippy::too_many_arguments, clippy::too_many_lines)] // Orchestrator for the read-pair pipeline
pub fn generate_read_pair(
    fragment: &Fragment,
    contig_name: &str,
    read_num: u64,
    read_length: usize,
    paired: bool,
    adapter_r1: &[u8],
    adapter_r2: &[u8],
    max_n_frac: f64,
    model: &impl ErrorModel,
    naming: &ReadNaming,
    methylation: Option<&MethylationConfig>,
    capture_pre_conversion: bool,
    rng: &mut impl Rng,
) -> Option<ReadPair> {
    let frag_len = fragment.bases.len();
    let adapter_bases = read_length.saturating_sub(frag_len.min(read_length));

    let r1_negative_strand = !fragment.is_forward;
    let r2_negative_strand = fragment.is_forward;

    // Step 1 — ambiguity filter on raw extracted bases (case preserved so
    // lowercase_fraction can detect ambiguity-resolved positions). Reject
    // BEFORE drawing any error-model RNG so a rejection doesn't desync the
    // stream.
    let r1_raw = extract_read_bases(&fragment.bases, read_length, adapter_r1, r1_negative_strand);
    if lowercase_fraction(&r1_raw) > max_n_frac {
        return None;
    }
    if paired {
        let r2_raw =
            extract_read_bases(&fragment.bases, read_length, adapter_r2, r2_negative_strand);
        if lowercase_fraction(&r2_raw) > max_n_frac {
            return None;
        }
    }

    // Step 2 — uppercased pre-chemistry top-strand fragment. This is the
    // single source from which both R1 and R2 are derived (post-chemistry).
    // Capturing pre-conversion bases for YS:Z is done by extracting from
    // this same buffer before chemistry runs.
    let mut pre_chem_top = fragment.bases.clone();
    uppercase_in_place(&mut pre_chem_top);

    // Step 3 — apply chemistry once, at fragment scale, on the source
    // strand. This produces top-strand-oriented bases reflecting the
    // appropriate strand's chemistry (em-seq / TAPS, with CpG context
    // resolved per-haplotype).
    let (post_chem_top, conversion_failed) = match methylation {
        Some(mc) => apply_fragment_chemistry(
            &pre_chem_top,
            fragment.hap_start,
            fragment.is_forward,
            fragment.haplotype_index,
            mc,
            rng,
        ),
        None => (pre_chem_top.clone(), false),
    };

    // Step 4 — derive per-mate read bases from the chemistry-applied
    // fragment. extract_read_bases handles orientation (revcomp when
    // is_negative_strand) and adapter padding for short fragments.
    let r1_bases = extract_read_bases(&post_chem_top, read_length, adapter_r1, r1_negative_strand);

    // Capture pre-conversion mate bases for YS:Z when requested. Read these
    // from `pre_chem_top` so the orientation matches the post-chemistry
    // mate bases (extract_read_bases revcomps the same way for both).
    let want_pre = capture_pre_conversion && methylation.is_some();
    let r1_pre_conversion = want_pre
        .then(|| extract_read_bases(&pre_chem_top, read_length, adapter_r1, r1_negative_strand));

    let r1 = build_mate(
        fragment,
        contig_name,
        ReadEnd::Read1,
        r1_negative_strand,
        r1_bases,
        r1_pre_conversion,
        adapter_bases,
        model,
        rng,
    );

    if !paired {
        let name = match naming {
            ReadNaming::Encoded => encoded_se_name(read_num, &r1.truth),
            ReadNaming::Simple => simple_name(read_num),
            ReadNaming::Illumina(header) => header.format_name(read_num),
        };

        let methylation_annotation = methylation.map(|_| MethylationAnnotation {
            conversion_type: ConversionType::from_strand(fragment.is_forward),
            conversion_failed,
            r1_pre_conversion_bases: r1.pre_conversion,
            r2_pre_conversion_bases: None,
            r1_call_tags: None,
            r2_call_tags: None,
        });

        return Some(ReadPair {
            read1: SimulatedRead { name, bases: r1.bases, qualities: r1.qualities },
            read2: None,
            r1_truth: r1.truth,
            r2_truth: None,
            r1_cigar: r1.cigar,
            r2_cigar: None,
            methylation: methylation_annotation,
        });
    }

    let r2_bases = extract_read_bases(&post_chem_top, read_length, adapter_r2, r2_negative_strand);
    let r2_pre_conversion = want_pre
        .then(|| extract_read_bases(&pre_chem_top, read_length, adapter_r2, r2_negative_strand));
    let r2 = build_mate(
        fragment,
        contig_name,
        ReadEnd::Read2,
        r2_negative_strand,
        r2_bases,
        r2_pre_conversion,
        adapter_bases,
        model,
        rng,
    );

    let name = match naming {
        ReadNaming::Encoded => encoded_pe_name(read_num, &r1.truth, &r2.truth),
        ReadNaming::Simple => simple_name(read_num),
        ReadNaming::Illumina(header) => header.format_name(read_num),
    };

    let methylation_annotation = methylation.map(|_| MethylationAnnotation {
        conversion_type: ConversionType::from_strand(fragment.is_forward),
        conversion_failed,
        r1_pre_conversion_bases: r1.pre_conversion,
        r2_pre_conversion_bases: r2.pre_conversion,
        r1_call_tags: None,
        r2_call_tags: None,
    });

    Some(ReadPair {
        read1: SimulatedRead { name: name.clone(), bases: r1.bases, qualities: r1.qualities },
        read2: Some(SimulatedRead { name, bases: r2.bases, qualities: r2.qualities }),
        r1_truth: r1.truth,
        r2_truth: Some(r2.truth),
        r1_cigar: r1.cigar,
        r2_cigar: Some(r2.cigar),
        methylation: methylation_annotation,
    })
}

/// Compute a CIGAR from a slice of ascending reference positions, plus
/// optional adapter soft-clipping.
///
/// Positions must be in ascending order (forward strand). Consecutive
/// positions incrementing by 1 produce M ops, same position produces I ops,
/// and gaps produce D ops.
///
/// This works for both R1 and R2: BAM CIGARs are always expressed in forward
/// reference order from the leftmost aligned position, so even negative-strand
/// reads use ascending positions.
///
/// Adapter bases are appended as a soft-clip (S) operation. For
/// negative-strand reads (`negative_strand = true`) the adapter is sequenced
/// at the 3' end of the read but sits at the left (5') end of the stored BAM
/// record after reverse-complementing, so the S op is placed at the *start*
/// of the CIGAR. For forward-strand reads the S op is placed at the *end*.
#[must_use]
pub fn cigar_from_ref_positions(
    positions: &[u32],
    adapter_bases: usize,
    negative_strand: bool,
) -> Cigar {
    let mut ops: Vec<Op> = Vec::new();

    if positions.is_empty() {
        if adapter_bases > 0 {
            ops.push(Op::new(Kind::SoftClip, adapter_bases));
        }
        return Cigar::from(ops);
    }

    let mut match_run: usize = 1; // First base is always a match.
    let mut ins_run: usize = 0;

    for i in 1..positions.len() {
        let prev = positions[i - 1];
        let curr = positions[i];

        if curr == prev + 1 {
            // Sequential ascending position: flush any insertion, extend match.
            if ins_run > 0 {
                ops.push(Op::new(Kind::Insertion, ins_run));
                ins_run = 0;
            }
            match_run += 1;
        } else if curr == prev {
            // Same position: insertion base. Flush match run, extend insertion.
            if match_run > 0 {
                ops.push(Op::new(Kind::Match, match_run));
                match_run = 0;
            }
            ins_run += 1;
        } else {
            // Gap in positions: deletion. Flush current runs.
            if ins_run > 0 {
                ops.push(Op::new(Kind::Insertion, ins_run));
                ins_run = 0;
            }
            if match_run > 0 {
                ops.push(Op::new(Kind::Match, match_run));
            }
            let gap = curr.saturating_sub(prev).saturating_sub(1);
            if gap > 0 {
                ops.push(Op::new(Kind::Deletion, gap as usize));
            }
            match_run = 1; // Current base starts a new match.
        }
    }

    // Flush remaining runs.
    if ins_run > 0 {
        ops.push(Op::new(Kind::Insertion, ins_run));
    }
    if match_run > 0 {
        ops.push(Op::new(Kind::Match, match_run));
    }

    // Adapter bases as soft-clip. Negative-strand reads are stored
    // reverse-complemented in BAM, so the adapter (at the 3' end in read
    // order) moves to the left (5') end of the record — a leading S op.
    if adapter_bases > 0 {
        if negative_strand {
            ops.insert(0, Op::new(Kind::SoftClip, adapter_bases));
        } else {
            ops.push(Op::new(Kind::SoftClip, adapter_bases));
        }
    }

    Cigar::from(ops)
}

/// Format a CIGAR as a human-readable string (e.g. "50M3I20M2S").
#[must_use]
pub fn cigar_to_string(cigar: &Cigar) -> String {
    use std::fmt::Write;
    cigar.as_ref().iter().fold(String::new(), |mut s, op| {
        let kind_char = match op.kind() {
            Kind::Match => 'M',
            Kind::Insertion => 'I',
            Kind::Deletion => 'D',
            Kind::SoftClip => 'S',
            Kind::HardClip => 'H',
            Kind::Skip => 'N',
            Kind::Pad => 'P',
            Kind::SequenceMatch => '=',
            Kind::SequenceMismatch => 'X',
        };
        let _ = write!(s, "{}{kind_char}", op.len());
        s
    })
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand::rngs::SmallRng;

    use super::*;
    use crate::error_model::illumina::IlluminaErrorModel;

    /// Build a simple fragment for testing.
    fn test_fragment(bases: &[u8], ref_start: u32) -> Fragment {
        #[expect(clippy::cast_possible_truncation, reason = "test data is small")]
        let ref_positions: Vec<u32> = (ref_start..ref_start + bases.len() as u32).collect();
        Fragment {
            bases: bases.to_vec(),
            ref_positions,
            ref_start,
            hap_start: ref_start,
            is_forward: true,
            haplotype_index: 0,
        }
    }

    // --- CIGAR generation tests ---

    #[test]
    fn test_cigar_all_match() {
        let cigar = cigar_from_ref_positions(&[0, 1, 2, 3, 4], 0, false);
        assert_eq!(cigar_to_string(&cigar), "5M");
    }

    #[test]
    fn test_cigar_with_insertion() {
        // Positions 0,1,2,2,2,3,4: two inserted bases at ref pos 2.
        let cigar = cigar_from_ref_positions(&[0, 1, 2, 2, 2, 3, 4], 0, false);
        assert_eq!(cigar_to_string(&cigar), "3M2I2M");
    }

    #[test]
    fn test_cigar_with_deletion() {
        // Gap from 2 to 5: 2 deleted ref bases.
        let cigar = cigar_from_ref_positions(&[0, 1, 2, 5, 6], 0, false);
        assert_eq!(cigar_to_string(&cigar), "3M2D2M");
    }

    #[test]
    fn test_cigar_with_adapter_softclip_forward() {
        // Forward-strand: adapter is a trailing soft-clip.
        let cigar = cigar_from_ref_positions(&[0, 1, 2], 2, false);
        assert_eq!(cigar_to_string(&cigar), "3M2S");
    }

    #[test]
    fn test_cigar_with_adapter_softclip_negative_strand() {
        // Negative-strand: adapter moves to a leading soft-clip after RC.
        let cigar = cigar_from_ref_positions(&[0, 1, 2], 2, true);
        assert_eq!(cigar_to_string(&cigar), "2S3M");
    }

    #[test]
    fn test_cigar_all_adapter() {
        // All-adapter read: placement doesn't matter, but it still round-trips.
        let cigar = cigar_from_ref_positions(&[], 5, false);
        assert_eq!(cigar_to_string(&cigar), "5S");
    }

    #[test]
    fn test_cigar_with_insertion_and_deletion() {
        // Insertion at pos 2 (two extra bases), then deletion of 2 ref bases.
        let cigar = cigar_from_ref_positions(&[0, 1, 2, 2, 5, 6], 0, false);
        assert_eq!(cigar_to_string(&cigar), "3M1I2D2M");
    }

    #[test]
    fn test_cigar_with_adapter_and_deletion() {
        let cigar = cigar_from_ref_positions(&[0, 1, 4, 5], 3, false);
        assert_eq!(cigar_to_string(&cigar), "2M2D2M3S");
    }

    #[test]
    fn test_cigar_single_base() {
        let cigar = cigar_from_ref_positions(&[42], 0, false);
        assert_eq!(cigar_to_string(&cigar), "1M");
    }

    #[test]
    fn test_cigar_high_positions() {
        // Negative-strand R2: positions start from a high offset (ascending),
        // no adapter — CIGAR is identical to forward strand.
        let cigar = cigar_from_ref_positions(&[100, 101, 102, 103, 104], 0, true);
        assert_eq!(cigar_to_string(&cigar), "5M");
    }

    // --- Read pair generation tests ---

    #[test]
    fn test_generate_pe_read_pair() {
        let fragment = test_fragment(b"ACGTACGTACGTACGTACGT", 100);
        let model = IlluminaErrorModel::new(10, 0.0, 0.0);
        let mut rng = SmallRng::seed_from_u64(42);

        let pair = generate_read_pair(
            &fragment,
            "chr1",
            1,
            10,
            true,
            b"ADAPTER",
            b"ADAPTER",
            1.0,
            &model,
            &ReadNaming::Encoded,
            None,
            false,
            &mut rng,
        )
        .expect("no ambiguous bases — should not reject");

        assert_eq!(pair.read1.bases, b"ACGTACGTAC");
        assert!(pair.read2.is_some());
        assert_eq!(pair.r1_truth.position, 101);
        assert!(pair.r2_truth.is_some());
        assert_eq!(cigar_to_string(&pair.r1_cigar), "10M");
        assert_eq!(cigar_to_string(pair.r2_cigar.as_ref().unwrap()), "10M");
    }

    #[test]
    fn test_generate_se_read() {
        let fragment = test_fragment(b"ACGTACGTAC", 100);
        let model = IlluminaErrorModel::new(10, 0.0, 0.0);
        let mut rng = SmallRng::seed_from_u64(42);

        let pair = generate_read_pair(
            &fragment,
            "chr1",
            5,
            10,
            false,
            b"ADAPTER",
            b"ADAPTER",
            1.0,
            &model,
            &ReadNaming::Encoded,
            None,
            false,
            &mut rng,
        )
        .unwrap();

        assert!(pair.read2.is_none());
        assert!(pair.r2_cigar.is_none());
        assert_eq!(cigar_to_string(&pair.r1_cigar), "10M");
    }

    #[test]
    fn test_adapter_cigar_softclip() {
        // Forward-strand fragment: R1 is positive-strand (trailing S), R2 is
        // negative-strand (leading S after BAM reverse-complement).
        let fragment = test_fragment(b"AC", 0);
        let model = IlluminaErrorModel::new(5, 0.0, 0.0);
        let mut rng = SmallRng::seed_from_u64(42);

        let pair = generate_read_pair(
            &fragment,
            "chr1",
            1,
            5,
            true,
            b"TTTTT",
            b"GGGGG",
            1.0,
            &model,
            &ReadNaming::Encoded,
            None,
            false,
            &mut rng,
        )
        .unwrap();

        assert_eq!(cigar_to_string(&pair.r1_cigar), "2M3S");
        assert_eq!(cigar_to_string(pair.r2_cigar.as_ref().unwrap()), "3S2M");
    }

    #[test]
    fn test_simple_name_mode() {
        let fragment = test_fragment(b"ACGT", 0);
        let model = IlluminaErrorModel::new(4, 0.0, 0.0);
        let mut rng = SmallRng::seed_from_u64(42);

        let pair = generate_read_pair(
            &fragment,
            "chr1",
            42,
            4,
            true,
            b"A",
            b"A",
            1.0,
            &model,
            &ReadNaming::Simple,
            None,
            false,
            &mut rng,
        )
        .unwrap();

        assert_eq!(pair.read1.name, "holodeck::42");
    }

    #[test]
    fn test_quality_scores_correct_length() {
        let fragment = test_fragment(b"ACGTACGTAC", 0);
        let model = IlluminaErrorModel::new(10, 0.001, 0.01);
        let mut rng = SmallRng::seed_from_u64(42);

        let pair = generate_read_pair(
            &fragment,
            "chr1",
            1,
            10,
            true,
            b"A",
            b"A",
            1.0,
            &model,
            &ReadNaming::Encoded,
            None,
            false,
            &mut rng,
        )
        .unwrap();

        assert_eq!(pair.read1.qualities.len(), 10);
        assert_eq!(pair.read2.as_ref().unwrap().qualities.len(), 10);
    }

    #[test]
    fn test_rejects_when_r1_exceeds_max_n_frac() {
        // Fragment entirely lowercase: every base on R1 (forward) is
        // flagged, lowercase_fraction == 1.0, which exceeds any threshold < 1.
        let mut fragment = test_fragment(b"acgtacgtac", 100);
        // Mark as forward so R1 is the positive-strand, genomic-bases read.
        fragment.is_forward = true;
        let model = IlluminaErrorModel::new(10, 0.0, 0.0);
        let mut rng = SmallRng::seed_from_u64(42);

        let pair = generate_read_pair(
            &fragment,
            "chr1",
            1,
            10,
            true,
            b"ADAPTER",
            b"ADAPTER",
            0.5,
            &model,
            &ReadNaming::Encoded,
            None,
            false,
            &mut rng,
        );

        assert!(pair.is_none(), "all-lowercase fragment should be rejected at threshold 0.5");
    }

    #[test]
    fn test_accepts_when_lowercase_below_threshold() {
        // 3 lowercase out of 10 = 0.3, below threshold 0.5 — should accept
        // and also emit uppercase bases.
        let fragment = test_fragment(b"ACaGcTAtCA", 0);
        let model = IlluminaErrorModel::new(10, 0.0, 0.0);
        let mut rng = SmallRng::seed_from_u64(42);

        let pair = generate_read_pair(
            &fragment,
            "chr1",
            1,
            10,
            true,
            b"ADAPTER",
            b"ADAPTER",
            0.5,
            &model,
            &ReadNaming::Encoded,
            None,
            false,
            &mut rng,
        )
        .expect("0.3 < 0.5 — should accept");

        // Emitted bases must be uppercase ACGT only.
        for &b in &pair.read1.bases {
            assert!(matches!(b, b'A' | b'C' | b'G' | b'T' | b'N'), "r1 got {b:?}");
        }
        for &b in &pair.read2.as_ref().unwrap().bases {
            assert!(matches!(b, b'A' | b'C' | b'G' | b'T' | b'N'), "r2 got {b:?}");
        }
    }

    #[test]
    fn test_generate_read_pair_with_em_seq_converts_forward_r1() {
        use crate::meth::{
            ContigMethylation, MethylationConfig, MethylationMode, MethylationTable,
        };

        let cm = ContigMethylation::from_tables(vec![MethylationTable::empty(1000)]);
        let mc = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::EmSeq,
            conversion_rate: 1.0,
            failure_rate: 0.0,
        };

        // Fragment with C's at known positions, forward strand → R1 sees them
        // directly; with 0% methylation and 100% conversion, every C in the
        // genomic portion of R1 must become T.
        let fragment = test_fragment(b"ACGTACGTAC", 100);
        let model = IlluminaErrorModel::new(10, 0.0, 0.0);
        let mut rng = SmallRng::seed_from_u64(42);

        let pair = generate_read_pair(
            &fragment,
            "chr1",
            1,
            10,
            true,
            b"ADAPTER",
            b"ADAPTER",
            1.0,
            &model,
            &ReadNaming::Encoded,
            Some(&mc),
            true,
            &mut rng,
        )
        .unwrap();

        // No C's should remain in R1 (genomic portion). Adapter is ADAPTER —
        // contains a 'C' but no genomic C should remain.
        #[expect(clippy::naive_bytecount, reason = "tiny test slice; clarity over speed")]
        let n_c_in_r1 = pair.read1.bases.iter().filter(|&&b| b == b'C').count();
        assert_eq!(
            n_c_in_r1, 0,
            "expected all C's in R1 converted, got bases {:?}",
            pair.read1.bases
        );

        let ann = pair.methylation.as_ref().expect("methylation annotation must be set");
        assert_eq!(ann.conversion_type, crate::meth::ConversionType::Ct);
        // Pre-conversion R1 retains C's; post-conversion R1 does not.
        let r1_pre = ann
            .r1_pre_conversion_bases
            .as_ref()
            .expect("capture_pre_conversion=true → R1 pre-conversion bases must be Some");
        #[expect(clippy::naive_bytecount, reason = "tiny test slice; clarity over speed")]
        let n_c_pre = r1_pre.iter().filter(|&&b| b == b'C').count();
        assert!(n_c_pre > 0, "pre-conversion R1 should still have C's, got {n_c_pre}");
        assert_eq!(
            r1_pre.len(),
            pair.read1.bases.len(),
            "pre-conversion length must match post-conversion read length"
        );
        // PE pair → R2 pre-conversion must also be present (this is a paired test).
        assert!(
            ann.r2_pre_conversion_bases.is_some(),
            "PE methylation annotation must include R2 pre-conversion bases"
        );
    }

    #[test]
    fn test_generate_read_pair_with_em_seq_reverse_strand_yields_ga_conversion_type() {
        use crate::meth::{
            ContigMethylation, ConversionType, MethylationConfig, MethylationMode, MethylationTable,
        };

        let cm = ContigMethylation::from_tables(vec![MethylationTable::empty(1000)]);
        let mc = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::EmSeq,
            conversion_rate: 1.0,
            failure_rate: 0.0,
        };

        // Reverse-strand fragment → reads come from the bottom strand → XG=GA.
        let mut fragment = test_fragment(b"ACGTACGTAC", 100);
        fragment.is_forward = false;
        let model = IlluminaErrorModel::new(10, 0.0, 0.0);
        let mut rng = SmallRng::seed_from_u64(42);

        let pair = generate_read_pair(
            &fragment,
            "chr1",
            1,
            10,
            true,
            b"ADAPTER",
            b"ADAPTER",
            1.0,
            &model,
            &ReadNaming::Encoded,
            Some(&mc),
            true,
            &mut rng,
        )
        .unwrap();

        let ann = pair.methylation.as_ref().expect("methylation annotation must be set");
        assert_eq!(ann.conversion_type, ConversionType::Ga);
        assert!(ann.r2_pre_conversion_bases.is_some());

        // Directional GA fragment: source strand is bottom of genome.
        //   R1 5'→3' = c2t(bottom). All bottom C's → T. R1 has no C's.
        //   R2 5'→3' = revcomp(c2t(bottom)). c2t(bottom) has no C's (all
        //     converted) → revcomp has no G's. C content of R2 reflects
        //     original-bottom G positions (= original-top C positions).
        // So R2 has no G's (not no C's).
        #[expect(clippy::naive_bytecount, reason = "tiny test slice; clarity over speed")]
        let cytosines_in_r1 = pair.read1.bases.iter().filter(|&&b| b == b'C').count();
        assert_eq!(
            cytosines_in_r1, 0,
            "R1 should have no C's after full conversion of source strand"
        );
        let r2 = pair.read2.as_ref().unwrap();
        #[expect(clippy::naive_bytecount, reason = "tiny test slice; clarity over speed")]
        let guanines_in_r2 = r2.bases.iter().filter(|&&b| b == b'G').count();
        assert_eq!(
            guanines_in_r2, 0,
            "R2 5'→3' = revcomp(c2t(bottom)) must have no G's after full conversion"
        );
    }

    #[test]
    fn test_directional_r2_of_ct_fragment_is_revcomp_of_c2t_top() {
        use crate::meth::{
            ContigMethylation, MethylationConfig, MethylationMode, MethylationTable,
        };

        // CpG-free top "ACAGACAGACAG" (12 bp). Under em-seq full conversion:
        //   c2t(top) = "ATAGATAGATAG"
        // Real directional libraries: R1 reads the BS-converted source strand;
        // R2 reads its PCR-synthesized complement. So for a CT fragment:
        //   R1 5'→3' = c2t(top)
        //   R2 5'→3' = revcomp(c2t(top))
        // — and crucially NOT c2t(revcomp(top)) (= c2t(bottom)).
        //
        // Verified against real Twist EM-seq base composition: R2 has C
        // content that mirrors R1 G content (within ~1pp), with R2 G content
        // ~3% (only methylated CpGs survive). The c2t(bottom) model would
        // give R2 ~25% G content.
        let cm = ContigMethylation::from_tables(vec![MethylationTable::empty(100)]);
        let mc = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::EmSeq,
            conversion_rate: 1.0,
            failure_rate: 0.0,
        };

        let fragment = test_fragment(b"ACAGACAGACAG", 0);
        let model = IlluminaErrorModel::new(12, 0.0, 0.0);
        let mut rng = SmallRng::seed_from_u64(42);

        let pair = generate_read_pair(
            &fragment,
            "chr1",
            1,
            12,
            true,
            b"ADAPTER",
            b"ADAPTER",
            1.0,
            &model,
            &ReadNaming::Encoded,
            Some(&mc),
            false,
            &mut rng,
        )
        .unwrap();

        assert_eq!(&pair.read1.bases, b"ATAGATAGATAG", "R1 must equal c2t(top)");
        let r2 = pair.read2.unwrap();
        assert_eq!(
            &r2.bases, b"CTATCTATCTAT",
            "R2 must equal revcomp(c2t(top)) for directional behavior; \
             c2t(revcomp(top)) = TTGTTTGTTTGT would indicate the wrong (per-mate) chemistry model"
        );
    }

    #[test]
    fn test_directional_r2_of_ga_fragment_is_revcomp_of_c2t_bottom() {
        use crate::meth::{
            ContigMethylation, MethylationConfig, MethylationMode, MethylationTable,
        };

        // For a GA (bottom-strand-derived) fragment, source strand is the
        // genome's bottom strand. fragment.bases is still in TOP orientation;
        // is_forward=false signals GA.
        //   bottom = revcomp(top) = revcomp("ACAGACAGACAG") = "CTGTCTGTCTGT"
        //   c2t(bottom) = "TTGTTTGTTTGT"
        //   R1 5'→3' = c2t(bottom) = "TTGTTTGTTTGT"
        //   R2 5'→3' = revcomp(c2t(bottom)) = "ACAAACAAACAA"
        // Under the prior per-mate model, R2 of GA was c2t(top) =
        // "ATAGATAGATAG" — wrong direction of chemistry for a GA fragment.
        let cm = ContigMethylation::from_tables(vec![MethylationTable::empty(100)]);
        let mc = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::EmSeq,
            conversion_rate: 1.0,
            failure_rate: 0.0,
        };

        let mut fragment = test_fragment(b"ACAGACAGACAG", 0);
        fragment.is_forward = false;
        let model = IlluminaErrorModel::new(12, 0.0, 0.0);
        let mut rng = SmallRng::seed_from_u64(42);

        let pair = generate_read_pair(
            &fragment,
            "chr1",
            1,
            12,
            true,
            b"ADAPTER",
            b"ADAPTER",
            1.0,
            &model,
            &ReadNaming::Encoded,
            Some(&mc),
            false,
            &mut rng,
        )
        .unwrap();

        assert_eq!(&pair.read1.bases, b"TTGTTTGTTTGT", "R1 of GA must equal c2t(bottom)");
        let r2 = pair.read2.unwrap();
        assert_eq!(
            &r2.bases, b"ACAAACAAACAA",
            "R2 of GA must equal revcomp(c2t(bottom)) for directional behavior"
        );
    }

    #[test]
    fn test_rejects_before_applying_errors_to_r1() {
        // Clean forward half (becomes R1), all-lowercase reverse half (becomes
        // R2 via reverse-complement). R1 passes the filter, but the pair
        // must be rejected *before* R1's errors advance the RNG — so a
        // subsequent call on a clean fragment produces identical output to
        // skipping the rejected call entirely.
        let mut fragment_half_bad = test_fragment(&[b'A'; 20], 0);
        // Second half lowercase: becomes R2 when is_forward=true.
        for b in &mut fragment_half_bad.bases[10..] {
            *b = b'a';
        }
        fragment_half_bad.is_forward = true;

        let clean_fragment = test_fragment(&[b'A'; 20], 0);
        let model = IlluminaErrorModel::new(10, 0.5, 0.5); // High rate so errors are observable.

        // Run A: rejected pair, then clean pair.
        let mut rng_a = SmallRng::seed_from_u64(123);
        let rejected = generate_read_pair(
            &fragment_half_bad,
            "chr1",
            1,
            10,
            true,
            b"TTTTTTTTTT",
            b"TTTTTTTTTT",
            0.5,
            &model,
            &ReadNaming::Encoded,
            None,
            false,
            &mut rng_a,
        );
        assert!(rejected.is_none(), "expected rejection when R2 is all-lowercase");
        let after_reject = generate_read_pair(
            &clean_fragment,
            "chr1",
            2,
            10,
            true,
            b"TTTTTTTTTT",
            b"TTTTTTTTTT",
            0.5,
            &model,
            &ReadNaming::Encoded,
            None,
            false,
            &mut rng_a,
        )
        .unwrap();

        // Run B: skip the rejected call entirely on an identical RNG seed.
        let mut rng_b = SmallRng::seed_from_u64(123);
        let direct = generate_read_pair(
            &clean_fragment,
            "chr1",
            2,
            10,
            true,
            b"TTTTTTTTTT",
            b"TTTTTTTTTT",
            0.5,
            &model,
            &ReadNaming::Encoded,
            None,
            false,
            &mut rng_b,
        )
        .unwrap();

        // If rejection had consumed any RNG draws (for R1's apply_errors or
        // qualities), the clean pair's bases/qualities would diverge.
        assert_eq!(after_reject.read1.bases, direct.read1.bases);
        assert_eq!(after_reject.read1.qualities, direct.read1.qualities);
        assert_eq!(
            after_reject.read2.as_ref().unwrap().bases,
            direct.read2.as_ref().unwrap().bases
        );
    }

    // --- Illumina naming integration tests ---

    #[test]
    fn test_illumina_name_pe() {
        use crate::read_naming::IlluminaHeader;

        let fragment = test_fragment(b"ACGTACGTACGTACGTACGT", 100);
        let model = IlluminaErrorModel::new(10, 0.0, 0.0);
        let mut rng = SmallRng::seed_from_u64(42);
        let header = IlluminaHeader::from_seed(99);
        let naming = ReadNaming::Illumina(&header);

        let pair = generate_read_pair(
            &fragment, "chr1", 1, 10, true, b"ADAPTER", b"ADAPTER", 1.0, &model, &naming, None,
            false, &mut rng,
        )
        .unwrap();

        let fields: Vec<&str> = pair.read1.name.split(':').collect();
        assert_eq!(fields.len(), 7, "Illumina PE name must have 7 colon-separated fields");
        assert_eq!(
            pair.read1.name,
            pair.read2.as_ref().unwrap().name,
            "R1 and R2 must share the same name"
        );
    }

    #[test]
    fn test_illumina_name_se() {
        use crate::read_naming::IlluminaHeader;

        let fragment = test_fragment(b"ACGTACGTAC", 100);
        let model = IlluminaErrorModel::new(10, 0.0, 0.0);
        let mut rng = SmallRng::seed_from_u64(42);
        let header = IlluminaHeader::from_seed(99);
        let naming = ReadNaming::Illumina(&header);

        let pair = generate_read_pair(
            &fragment, "chr1", 5, 10, false, b"ADAPTER", b"ADAPTER", 1.0, &model, &naming, None,
            false, &mut rng,
        )
        .unwrap();

        let fields: Vec<&str> = pair.read1.name.split(':').collect();
        assert_eq!(fields.len(), 7, "Illumina SE name must have 7 colon-separated fields");
        assert!(pair.read2.is_none());
    }
}
