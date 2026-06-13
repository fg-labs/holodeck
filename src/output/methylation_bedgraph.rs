//! Population-fraction methylation bedGraph in MethylDackel `extract` format.
//!
//! Derived closed-form from the per-haplotype methylation bitmap — no
//! simulated reads. For each reference CpG, aggregates over
//! (haplotype × strand) presence/absence and emits
//! `rate = round(100 * n_methylated / (n_methylated + n_unmethylated))`
//! as an integer percentage (matching MethylDackel `extract`).
//!
//! This is intentionally distinct from `simulate`'s `--cpg-truth-bedgraph`,
//! which is coverage-weighted (depends on simulated read counts).

use std::io::Write;

use anyhow::{Result, ensure};

use crate::haplotype::Haplotype;
use crate::meth::ContigMethylation;
use crate::sequence_dict::SequenceDictionary;

/// Write the `track` header line for a MethylDackel-format population-fraction
/// bedGraph.
///
/// Call once before streaming per-contig records with
/// [`write_bedgraph_records`].
///
/// # Errors
///
/// Returns an error if writing to `writer` fails.
pub fn write_bedgraph_header<W: Write>(writer: &mut W) -> Result<()> {
    writeln!(
        writer,
        "track type=\"bedGraph\" description=\"holodeck methylate population fraction\""
    )?;
    Ok(())
}

/// Write per-CpG bedGraph records for a single contig, without the track
/// header line.
///
/// For each reference CpG (top-strand `CG` dinucleotide), aggregates over all
/// haplotype × strand combinations whose haplotype actually *preserves* that
/// CpG:
/// - `n_methylated` — count of (haplotype, strand) entries whose methylation
///   bit is set at this site.
/// - `n_unmethylated` — preserving haplotype entries that are not methylated.
/// - `rate` — `round(100 * n_methylated / (n_methylated + n_unmethylated))`,
///   an integer percentage matching MethylDackel `extract`.
///
/// Reference CpGs destroyed on every haplotype (SNP at the C or G, indel
/// straddling either base) are omitted entirely — `n_meth = n_unmeth = 0`
/// makes the rate undefined, not zero. Variants that only shift downstream
/// coordinates are handled by mapping each ref position through that
/// haplotype's [`Haplotype::extract_fragment`] / `hap_start` before indexing
/// the haplotype-coordinate methylation bitmap.
///
/// `haplotypes` must be either empty (variant-free contig — every ref CpG
/// is queried at its reference coordinates directly) or carry exactly one
/// entry per haplotype tracked by `methylation`, in the same order. A
/// length mismatch is a programming error from a caller that built the
/// methylation table and the haplotypes slice from inconsistent inputs;
/// returning an error here surfaces that mismatch instead of silently
/// reading reference offsets out of haplotype-coordinate bitmaps.
///
/// The `chrom` start and end columns use 0-based half-open coordinates
/// matching MethylDackel's `extract` output: `start = top-C ref position`,
/// `end = start + 1`.
///
/// # Errors
///
/// Returns an error if writing to `writer` fails, or if `haplotypes` is
/// non-empty and its length differs from `methylation.len()`.
pub fn write_bedgraph_records<W: Write>(
    writer: &mut W,
    chrom: &str,
    reference: &[u8],
    haplotypes: &[Haplotype],
    methylation: &ContigMethylation,
) -> Result<()> {
    ensure!(
        haplotypes.is_empty() || haplotypes.len() == methylation.len(),
        "haplotypes length ({}) must equal methylation haplotype count ({}) \
         on contig {chrom:?}",
        haplotypes.len(),
        methylation.len(),
    );
    if reference.len() < 2 {
        return Ok(());
    }

    // Aggregate per reference CpG. `BTreeMap` so emit order is the sorted
    // top-C reference position, matching MethylDackel's output and the
    // previous loop's natural order. Each value is `(n_meth, n_unmeth)`.
    let mut by_ref_pos: std::collections::BTreeMap<u32, (u32, u32)> =
        std::collections::BTreeMap::new();

    if haplotypes.is_empty() {
        // No variants → ref == hap. Walk the reference once.
        for i in 0..reference.len() - 1 {
            if !reference[i].eq_ignore_ascii_case(&b'C')
                || !reference[i + 1].eq_ignore_ascii_case(&b'G')
            {
                continue;
            }
            #[expect(clippy::cast_possible_truncation, reason = "ref pos fits u32")]
            let top_c = i as u32;
            let entry = by_ref_pos.entry(top_c).or_insert((0, 0));
            for hap_idx in 0..methylation.len() {
                let table = methylation.table_for(hap_idx);
                if table.is_methylated(top_c, false) {
                    entry.0 += 1;
                } else {
                    entry.1 += 1;
                }
                if table.is_methylated(top_c + 1, true) {
                    entry.0 += 1;
                } else {
                    entry.1 += 1;
                }
            }
        }
    } else {
        // Variants present → materialize each haplotype once, walk its
        // bases looking for CpGs, map each haplotype CpG back to a ref
        // position via the per-base `ref_positions` returned by
        // `extract_fragment`. The per-CpG `extract_fragment(..., 2)` call
        // we used to do per ref CpG × haplotype is gone; this is O(L + V)
        // per haplotype instead of O(C × log V).
        for (hap_idx, hap) in haplotypes.iter().enumerate() {
            #[expect(clippy::cast_possible_truncation, reason = "ref length fits u32")]
            let hap_len = hap.hap_position_for(reference.len() as u32) as usize;
            let (hap_bases, ref_positions, _hap_start) =
                hap.extract_fragment(reference, 0, hap_len);
            if hap_bases.len() < 2 {
                continue;
            }
            let table = methylation.table_for(hap_idx);
            for h in 0..hap_bases.len() - 1 {
                if !hap_bases[h].eq_ignore_ascii_case(&b'C')
                    || !hap_bases[h + 1].eq_ignore_ascii_case(&b'G')
                {
                    continue;
                }
                // Only count haplotype CpGs that correspond to a true ref
                // CpG at adjacent ref positions. A CpG formed across an
                // insertion (`ref_positions[h+1] == ref_positions[h]`) has
                // no reference coordinate to attribute it to — MethylDackel
                // wouldn't see it either, and we want output parity.
                let ref_top_c = ref_positions[h];
                let ref_bot_c = ref_positions[h + 1];
                if ref_bot_c != ref_top_c + 1 {
                    continue;
                }
                let ref_top_idx = ref_top_c as usize;
                let ref_bot_idx = ref_bot_c as usize;
                if ref_bot_idx >= reference.len()
                    || !reference[ref_top_idx].eq_ignore_ascii_case(&b'C')
                    || !reference[ref_bot_idx].eq_ignore_ascii_case(&b'G')
                {
                    continue;
                }
                #[expect(clippy::cast_possible_truncation, reason = "hap pos fits u32")]
                let hap_top = h as u32;
                let entry = by_ref_pos.entry(ref_top_c).or_insert((0, 0));
                if table.is_methylated(hap_top, false) {
                    entry.0 += 1;
                } else {
                    entry.1 += 1;
                }
                if table.is_methylated(hap_top + 1, true) {
                    entry.0 += 1;
                } else {
                    entry.1 += 1;
                }
            }
        }
    }

    for (top_c, (n_meth, n_unmeth)) in by_ref_pos {
        let denom = n_meth + n_unmeth;
        if denom == 0 {
            continue;
        }
        // MethylDackel's `extract` reports the rate as an integer percentage;
        // round to match so this output and `simulate --cpg-truth-bedgraph`
        // (which rounds identically) can be compared by the same downstream
        // tooling.
        let rate = (f64::from(n_meth) / f64::from(denom) * 100.0).round();
        #[expect(clippy::cast_possible_truncation, reason = "rate is in [0, 100]")]
        #[expect(clippy::cast_sign_loss, reason = "rate is non-negative")]
        let rate = rate as u32;
        writeln!(writer, "{chrom}\t{top_c}\t{}\t{rate}\t{n_meth}\t{n_unmeth}", top_c + 1)?;
    }
    Ok(())
}

/// Write a complete population-fraction bedGraph — header line followed by
/// per-CpG records for every contig in `per_contig`.
///
/// Each `per_contig` entry is `(contig_name, ContigMethylation, reference,
/// haplotypes_for_this_contig)`. The `haplotypes` slice must be the same
/// per-haplotype layout used to build `ContigMethylation` so per-haplotype
/// ref→hap coordinate mappings line up; pass an empty slice for variant-free
/// contigs.
///
/// `dict` is accepted for future contig-order validation but is not used yet.
///
/// # Errors
///
/// Returns an error if writing to `writer` fails.
pub fn write_bedgraph<W: Write>(
    writer: &mut W,
    dict: &SequenceDictionary,
    per_contig: &[(String, ContigMethylation, Vec<u8>, Vec<Haplotype>)],
) -> Result<()> {
    let _ = dict; // reserved for future contig-order validation
    write_bedgraph_header(writer)?;
    for (chrom, cm, reference, haplotypes) in per_contig {
        write_bedgraph_records(writer, chrom, reference, haplotypes, cm)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::haplotype::{Haplotype, build_haplotypes};
    use crate::meth::{ContigMethylation, MethylationTable};
    use crate::sequence_dict::{SequenceDictionary, SequenceMetadata};
    use crate::vcf::genotype::{Genotype, VariantRecord};

    /// Methylation state for one haplotype at the single CpG in the `ACGT`
    /// test reference (top-C at pos 1, bottom-C at pos 2).
    #[derive(Clone, Copy)]
    struct HapMethState {
        top: bool,
        bottom: bool,
    }

    /// Build a diploid [`ContigMethylation`] for the reference `ACGT`, which
    /// has a single CpG: top-C at position 1, bottom-C at position 2.
    fn build_diploid_at_one_cpg(
        hap0: HapMethState,
        hap1: HapMethState,
    ) -> (Vec<u8>, ContigMethylation) {
        // Reference ACGT: CpG at top-C pos 1 and bottom-C pos 2.
        let reference = b"ACGT".to_vec();
        let mut h0 = MethylationTable::with_len(4);
        if hap0.top {
            h0.set_top(1, true);
        }
        if hap0.bottom {
            h0.set_bottom(2, true);
        }
        let mut h1 = MethylationTable::with_len(4);
        if hap1.top {
            h1.set_top(1, true);
        }
        if hap1.bottom {
            h1.set_bottom(2, true);
        }
        let cm = ContigMethylation::from_tables(vec![h0, h1]);
        (reference, cm)
    }

    /// Build a minimal [`SequenceDictionary`] for a single contig.
    fn single_contig_dict(name: &str, len: usize) -> SequenceDictionary {
        let entry = SequenceMetadata::new(0, name.to_string(), len);
        SequenceDictionary::from_entries(vec![entry])
    }

    #[test]
    fn population_fraction_bedgraph_diploid_all_methylated() {
        let (reference, cm) = build_diploid_at_one_cpg(
            HapMethState { top: true, bottom: true },
            HapMethState { top: true, bottom: true },
        );
        let dict = single_contig_dict("chr1", reference.len());
        let mut buf = Vec::new();
        write_bedgraph(&mut buf, &dict, &[("chr1".to_string(), cm, reference, Vec::new())])
            .unwrap();
        let s = String::from_utf8(buf).unwrap();
        // Track header must be present.
        assert!(s.starts_with("track type="), "missing track header: {s}");
        // One CpG record: start=1, end=2, rate=100, n_meth=4, n_unmeth=0.
        assert!(s.contains("chr1\t1\t2\t100"), "expected rate 100: {s}");
        assert!(s.contains("4\t0"), "expected n_meth=4 n_unmeth=0: {s}");
    }

    #[test]
    fn population_fraction_bedgraph_diploid_hemi() {
        // Only haplotype-0 top strand is methylated → n_meth=1 out of 4.
        let (reference, cm) = build_diploid_at_one_cpg(
            HapMethState { top: true, bottom: false },
            HapMethState { top: false, bottom: false },
        );
        let dict = single_contig_dict("chr1", reference.len());
        let mut buf = Vec::new();
        write_bedgraph(&mut buf, &dict, &[("chr1".to_string(), cm, reference, Vec::new())])
            .unwrap();
        let s = String::from_utf8(buf).unwrap();
        // rate = 100 * 1/4 = 25.
        assert!(s.contains("chr1\t1\t2\t25"), "expected rate 25: {s}");
        assert!(s.contains("1\t3"), "expected n_meth=1 n_unmeth=3: {s}");
    }

    #[test]
    fn population_fraction_bedgraph_no_cpg_emits_only_header() {
        // Reference with no CpG: "AATT".
        let reference = b"AATT".to_vec();
        let h0 = MethylationTable::with_len(4);
        let cm = ContigMethylation::from_tables(vec![h0]);
        let dict = single_contig_dict("chr1", reference.len());
        let mut buf = Vec::new();
        write_bedgraph(&mut buf, &dict, &[("chr1".to_string(), cm, reference, Vec::new())])
            .unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("track type="), "missing track header: {s}");
        // No data lines.
        let data_lines: Vec<&str> = s.lines().filter(|l| !l.starts_with("track")).collect();
        assert!(data_lines.is_empty(), "unexpected data lines: {s}");
    }

    #[test]
    fn population_fraction_bedgraph_multiple_cpgs() {
        // Reference ACGTACG has two CpGs: top-C at pos 1 and pos 5.
        // Diploid, all four (hap × strand) bits set for both sites.
        let reference = b"ACGTACG".to_vec();
        let mut h0 = MethylationTable::with_len(7);
        h0.set_top(1, true);
        h0.set_bottom(2, true);
        h0.set_top(5, true);
        h0.set_bottom(6, true);
        let mut h1 = MethylationTable::with_len(7);
        h1.set_top(1, true);
        h1.set_bottom(2, true);
        h1.set_top(5, true);
        h1.set_bottom(6, true);
        let cm = ContigMethylation::from_tables(vec![h0, h1]);
        let dict = single_contig_dict("chr1", reference.len());
        let mut buf = Vec::new();
        write_bedgraph(&mut buf, &dict, &[("chr1".to_string(), cm, reference, Vec::new())])
            .unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("chr1\t1\t2\t100"), "missing first CpG: {s}");
        assert!(s.contains("chr1\t5\t6\t100"), "missing second CpG: {s}");
    }

    #[test]
    fn population_fraction_bedgraph_rounds_non_integer_rate() {
        // Triploid reference ACGT (one CpG, top-C at 1, bottom-C at 2) gives
        // 6 (haplotype × strand) entries. Set 2 of them methylated → 2/6 =
        // 33.33%, which MethylDackel reports as the integer `33`. This pins
        // the rounding: an unrounded f64 would print `33.33333333333333`.
        let reference = b"ACGT".to_vec();
        let mut h0 = MethylationTable::with_len(4);
        h0.set_top(1, true);
        let mut h1 = MethylationTable::with_len(4);
        h1.set_top(1, true);
        let h2 = MethylationTable::with_len(4);
        let cm = ContigMethylation::from_tables(vec![h0, h1, h2]);
        let dict = single_contig_dict("chr1", reference.len());
        let mut buf = Vec::new();
        write_bedgraph(&mut buf, &dict, &[("chr1".to_string(), cm, reference, Vec::new())])
            .unwrap();
        let s = String::from_utf8(buf).unwrap();
        // rate = round(100 * 2/6) = 33; n_meth=2, n_unmeth=4.
        assert!(s.contains("chr1\t1\t2\t33\t2\t4"), "expected rounded rate 33: {s}");
        assert!(!s.contains("33.3"), "rate must be an integer, not a float: {s}");
    }

    #[test]
    fn write_bedgraph_records_streaming_matches_write_bedgraph() {
        // Verify that write_bedgraph_header + write_bedgraph_records produces
        // the same output as write_bedgraph.
        let (reference, cm) = build_diploid_at_one_cpg(
            HapMethState { top: true, bottom: false },
            HapMethState { top: true, bottom: true },
        );
        let dict = single_contig_dict("chr1", reference.len());

        let mut all_at_once = Vec::new();
        write_bedgraph(
            &mut all_at_once,
            &dict,
            &[("chr1".to_string(), cm.clone(), reference.clone(), Vec::new())],
        )
        .unwrap();

        let mut streamed = Vec::new();
        write_bedgraph_header(&mut streamed).unwrap();
        let no_haps: Vec<Haplotype> = Vec::new();
        write_bedgraph_records(&mut streamed, "chr1", &reference, &no_haps, &cm).unwrap();

        assert_eq!(all_at_once, streamed, "streaming API must match batch API");
    }

    /// Insertion variant: a CpG downstream of an indel must be read from the
    /// haplotype-shifted coordinate, not from its raw reference offset.
    ///
    /// Regression: before this fix, `write_bedgraph_records` indexed
    /// `ContigMethylation` (haplotype-coord) with reference indices, so the
    /// methylated CpG on the inserted haplotype was missed entirely and the
    /// emitted rate was 0% instead of 50%.
    #[test]
    fn population_fraction_bedgraph_respects_haplotype_coordinates_for_indel() {
        // Reference: AACGAA. The single CpG is at ref pos 2 (top-C) / 3 (G).
        // Heterozygous insertion at pos 1 on hap1 turns the A at pos 1 into
        // "AAA", shifting downstream positions by +2. On hap1 the CpG lives
        // at hap-coords 4 / 5.
        let reference = b"AACGAA".to_vec();
        let variants = vec![VariantRecord {
            position: 1,
            ref_allele: b"A".to_vec(),
            alt_alleles: vec![b"AAA".to_vec()],
            genotype: Genotype::parse("0|1").unwrap(),
        }];
        let haplotypes = build_haplotypes(&variants, 2, &mut rand::rng());
        assert_eq!(haplotypes.len(), 2);
        // Hap0 has no variants → length matches reference. Its CpG is
        // unmethylated (default-zero bits).
        // Hap1 has the +2 insertion → length = 8. Its CpG is methylated on
        // both strands at the shifted haplotype coordinates.
        let h0 = MethylationTable::with_len(reference.len());
        let mut h1 = MethylationTable::with_len(reference.len() + 2);
        h1.set_top(4, true);
        h1.set_bottom(5, true);
        let cm = ContigMethylation::from_tables(vec![h0, h1]);

        let mut buf = Vec::new();
        write_bedgraph_header(&mut buf).unwrap();
        write_bedgraph_records(&mut buf, "chr1", &reference, &haplotypes, &cm).unwrap();
        let s = String::from_utf8(buf).unwrap();

        // Expected: hap0 contributes 0+2 unmethylated, hap1 contributes 2+0
        // methylated. Rate = 100 * 2 / 4 = 50.
        assert!(s.contains("chr1\t2\t3\t50\t2\t2"), "expected rate 50, counts 2/2: {s}");
    }

    /// SNP that destroys the CpG on a haplotype: that haplotype is excluded
    /// from the denominator entirely instead of being miscounted as
    /// unmethylated.
    #[test]
    fn population_fraction_bedgraph_drops_haplotypes_with_destroyed_cpg() {
        // Reference: ACGT, single CpG at pos 1/2. Het SNP at pos 1 (C→T)
        // destroys the CpG on hap1.
        let reference = b"ACGT".to_vec();
        let variants = vec![VariantRecord {
            position: 1,
            ref_allele: b"C".to_vec(),
            alt_alleles: vec![b"T".to_vec()],
            genotype: Genotype::parse("0|1").unwrap(),
        }];
        let haplotypes = build_haplotypes(&variants, 2, &mut rand::rng());
        let mut h0 = MethylationTable::with_len(reference.len());
        let h1 = MethylationTable::with_len(reference.len());
        // Hap0 carries the reference CpG and is fully methylated.
        h0.set_top(1, true);
        h0.set_bottom(2, true);
        let cm = ContigMethylation::from_tables(vec![h0, h1]);

        let mut buf = Vec::new();
        write_bedgraph_header(&mut buf).unwrap();
        write_bedgraph_records(&mut buf, "chr1", &reference, &haplotypes, &cm).unwrap();
        let s = String::from_utf8(buf).unwrap();

        // Only hap0 preserves the CpG and is fully methylated → rate=100,
        // counts 2/0 (NOT 2/2 as the old reference-index-only writer would
        // have produced by counting hap1's destroyed CpG as unmethylated).
        assert!(s.contains("chr1\t1\t2\t100\t2\t0"), "expected rate 100, counts 2/0: {s}");
    }

    /// Length mismatch between `haplotypes` and the methylation table is a
    /// caller bug — surface it as an error rather than silently reading
    /// haplotype-coordinate bitmaps with reference offsets.
    #[test]
    fn population_fraction_bedgraph_errors_when_haplotypes_len_mismatches_methylation() {
        let reference = b"ACGT".to_vec();
        // methylation says 2 haplotypes…
        let cm = ContigMethylation::from_tables(vec![
            MethylationTable::with_len(reference.len()),
            MethylationTable::with_len(reference.len()),
        ]);
        // …but caller passes only 1 haplotype.
        let haplotypes = build_haplotypes(&[], 1, &mut rand::rng());
        assert_eq!(haplotypes.len(), 1);

        let mut buf = Vec::new();
        let err = write_bedgraph_records(&mut buf, "chr1", &reference, &haplotypes, &cm)
            .expect_err("length mismatch should error");
        let msg = format!("{err}");
        assert!(msg.contains("haplotypes length"), "unexpected error message: {msg}");
        assert!(msg.contains("methylation haplotype count"), "unexpected error: {msg}");
    }

    /// CpG destroyed on every haplotype: emit no row at all (denominator 0
    /// means the population rate is undefined, not zero).
    #[test]
    fn population_fraction_bedgraph_skips_cpg_destroyed_on_every_haplotype() {
        // Reference: ACGT, single CpG. Hom-alt SNP destroys it on both haps.
        let reference = b"ACGT".to_vec();
        let variants = vec![VariantRecord {
            position: 1,
            ref_allele: b"C".to_vec(),
            alt_alleles: vec![b"T".to_vec()],
            genotype: Genotype::parse("1|1").unwrap(),
        }];
        let haplotypes = build_haplotypes(&variants, 2, &mut rand::rng());
        let h0 = MethylationTable::with_len(reference.len());
        let h1 = MethylationTable::with_len(reference.len());
        let cm = ContigMethylation::from_tables(vec![h0, h1]);

        let mut buf = Vec::new();
        write_bedgraph_header(&mut buf).unwrap();
        write_bedgraph_records(&mut buf, "chr1", &reference, &haplotypes, &cm).unwrap();
        let s = String::from_utf8(buf).unwrap();

        let data_lines: Vec<&str> = s.lines().filter(|l| !l.starts_with("track")).collect();
        assert!(data_lines.is_empty(), "expected no data lines, got: {s}");
    }
}
