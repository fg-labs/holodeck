//! Per-CpG ground-truth methylation tally and bedGraph writer.
//!
//! Counts, per reference CpG site, how many simulated read mates "called"
//! the site as methylated vs unmethylated based on the simulator's
//! per-haplotype, per-strand methylation bitmap. The output format matches
//! [MethylDackel]'s `extract` CpG `.bedGraph` exactly:
//!
//! ```text
//! track type="bedGraph" description="..."
//! chrom  start  end  rate(0-100)  n_methylated  n_unmethylated
//! ```
//!
//! so the same downstream concordance scripts that consume MethylDackel
//! output (e.g. comparing aligner-derived methylation calls to truth) can
//! be pointed at the truth bedGraph without modification.
//!
//! [MethylDackel]: https://github.com/dpryan79/MethylDackel
//!
//! # Counting model
//!
//! For each simulated read mate (R1 and R2 counted independently — matching
//! MethylDackel's default per-read tallying):
//!
//! - For every reference CpG (top-strand `CG` dinucleotide at top positions
//!   `(p, p+1)`) whose strand-specific C is covered by the mate's genomic
//!   span:
//!     - **CT (top-strand-derived) fragments** carry top-strand chemistry;
//!       the call is `top[hap_pos_for(p)]`.
//!     - **GA (bottom-strand-derived) fragments** carry bottom-strand
//!       chemistry; the call is `bottom[hap_pos_for(p+1)]`.
//! - The call increments `n_methylated` (bit set) or `n_unmethylated` (bit
//!   clear).
//!
//! Both R1 and R2 of a directional pair derive from the same source-strand
//! chemistry. When both mates cover the same CpG site (the genomic overlap
//! of the pair), they make the same call (since chemistry is now applied
//! once at fragment scale) and both contribute.
//!
//! Calls are made against the per-haplotype methylation bitmap, not the
//! post-conversion sequenced base. So the truth reflects the simulator's
//! biological intent and is unaffected by `--methylation-conversion-rate`
//! (chemistry inefficiency) or sequencing errors — useful precisely because
//! it lets a downstream evaluator measure how much of the alignment-tool
//! disagreement is due to chemistry/error noise vs alignment failure.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};

use crate::haplotype::Haplotype;
use crate::meth::ContigMethylation;
use crate::sequence_dict::SequenceDictionary;

/// Per-(contig, top-C reference position) tally of methylated vs unmethylated
/// calls. Indexed by `(contig_index, ref_pos_of_top_C)`; the implied bedGraph
/// `end` is `ref_pos + 1` (matches MethylDackel).
#[derive(Debug, Default)]
pub struct CpgTruthTally {
    counts: HashMap<(usize, u32), (u32, u32)>,
}

impl CpgTruthTally {
    /// Empty tally.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of distinct CpG sites with at least one call. Test-only.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn len(&self) -> usize {
        self.counts.len()
    }

    /// Look up the (n_methylated, n_unmethylated) tuple for a given site.
    /// Test-only.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn get(&self, contig_idx: usize, ref_pos: u32) -> Option<(u32, u32)> {
        self.counts.get(&(contig_idx, ref_pos)).copied()
    }

    /// Record one mate's contribution. `mate_ref_positions` is the per-base
    /// reference position list for the mate's GENOMIC portion (no adapter),
    /// in ascending reference order. `ref_cpgs` is the precomputed sorted
    /// list of top-C reference CpG positions for this contig (caller pre-
    /// computes once per contig). `is_forward_fragment` distinguishes CT
    /// (top-strand-derived) from GA (bottom-strand-derived) fragments.
    ///
    /// For CT fragments the call is read off `top[hap_position_for(p)]`;
    /// for GA fragments off `bottom[hap_position_for(p+1)]`. Both look up
    /// against the haplotype the fragment was sampled from. The bitmap is
    /// indexed by haplotype coordinates, so we map each reference CpG
    /// position through `Haplotype::hap_position_for`.
    ///
    /// # Panics
    ///
    /// Panics if `mate_ref_positions.last()` is unreachable after the
    /// non-empty check (it cannot be — the early return covers that case).
    #[allow(clippy::too_many_arguments)]
    pub fn record_mate(
        &mut self,
        contig_idx: usize,
        mate_ref_positions: &[u32],
        ref_cpgs: &[u32],
        methylation: &ContigMethylation,
        haplotype: &Haplotype,
        haplotype_index: usize,
        is_forward_fragment: bool,
    ) {
        if mate_ref_positions.is_empty() || ref_cpgs.is_empty() {
            return;
        }
        // Mate genomic span [min, max] in reference coordinates. fragment
        // ref_positions are ascending and may include duplicates (insertions
        // map multiple read bases to one ref pos) — first/last are still
        // the bounding values.
        let mate_min = mate_ref_positions[0];
        let mate_max = *mate_ref_positions.last().unwrap();

        // Covered C reference position depends on which strand carries the
        // chemistry: CT → top-C at p; GA → bottom-C at p+1.
        // The CpG site "starts" at p (top-C) regardless of which strand
        // carried the call — that's how MethylDackel groups per-CpG calls.
        let table = methylation.table_for(haplotype_index);

        // Walk the ref CpGs in the mate's range. Use binary search to bound.
        let lo = ref_cpgs.partition_point(|&p| p < mate_min);
        let hi = ref_cpgs.partition_point(|&p| p <= mate_max);
        for &site_p in &ref_cpgs[lo..hi] {
            let covered_pos = if is_forward_fragment { site_p } else { site_p + 1 };
            // The strand-specific C must lie within the mate's genomic span.
            // (For CT fragments at site_p in [mate_min, mate_max] this is
            // always true; for GA fragments site_p+1 may fall just past
            // mate_max if the CpG sits at the right edge.)
            if covered_pos < mate_min || covered_pos > mate_max {
                continue;
            }
            // Skip CpGs that fall in deletion gaps on this haplotype: a
            // deleted reference position has no haplotype base, so the read
            // at this fragment doesn't actually carry that C. We can detect
            // this by looking up the position in the mate's ref_positions —
            // if it's missing, the C was deleted out.
            if mate_ref_positions.binary_search(&covered_pos).is_err() {
                continue;
            }
            let hap_pos = haplotype.hap_position_for(covered_pos);
            let is_meth = table.is_methylated(hap_pos, !is_forward_fragment);
            let entry = self.counts.entry((contig_idx, site_p)).or_insert((0, 0));
            if is_meth {
                entry.0 = entry.0.saturating_add(1);
            } else {
                entry.1 = entry.1.saturating_add(1);
            }
        }
    }

    /// Write the tally to `path` as a MethylDackel-format CpG bedGraph.
    ///
    /// Format (tab-separated): `chrom  start  end  rate  n_meth  n_unmeth`,
    /// preceded by a `track` header line. Sites are emitted in (contig
    /// index, reference position) order. `rate` is the integer percentage
    /// `round(100 * n_meth / (n_meth + n_unmeth))`. Sites with zero total
    /// coverage are skipped (they would not appear in MethylDackel's output
    /// either).
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be created, if a contig index
    /// in the tally is missing from `dict`, or if any I/O write fails.
    pub fn write_bedgraph(&self, dict: &SequenceDictionary, path: &Path) -> Result<()> {
        let file = File::create(path)
            .with_context(|| format!("creating CpG truth bedGraph: {}", path.display()))?;
        let mut w = BufWriter::new(file);
        writeln!(
            w,
            "track type=\"bedGraph\" description=\"holodeck CpG truth (per-haplotype methylation calls per simulated read)\""
        )?;
        let mut sorted: Vec<((usize, u32), (u32, u32))> =
            self.counts.iter().map(|(k, v)| (*k, *v)).collect();
        sorted.sort_by_key(|&(k, _)| k);
        for ((contig_idx, p), (n_meth, n_unmeth)) in sorted {
            let total = n_meth + n_unmeth;
            if total == 0 {
                continue;
            }
            let rate = ((f64::from(n_meth) / f64::from(total)) * 100.0).round();
            #[expect(clippy::cast_possible_truncation, reason = "rate is in [0, 100]")]
            #[expect(clippy::cast_sign_loss, reason = "rate is non-negative")]
            let rate = rate as u32;
            let contig_name = dict
                .get_by_index(contig_idx)
                .with_context(|| format!("missing contig index {contig_idx} in dictionary"))?
                .name();
            writeln!(w, "{contig_name}\t{p}\t{}\t{rate}\t{n_meth}\t{n_unmeth}", p + 1)?;
        }
        w.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand::rngs::SmallRng;

    use super::*;
    use crate::haplotype::build_haplotypes;
    use crate::meth::{ContigMethylation, MethylationTable};

    /// Build a single-haplotype `ContigMethylation` with a manually-set
    /// bitmap, for tests that need precise control over the methylation
    /// state without depending on `from_haplotypes`'s RNG.
    fn cm_with_top(top_meth_positions: &[u32], len: usize) -> ContigMethylation {
        let mut table = MethylationTable::empty(len);
        for &p in top_meth_positions {
            table.set_top(p as usize, true);
        }
        ContigMethylation::from_tables(vec![table])
    }

    fn cm_with_bottom(bottom_meth_positions: &[u32], len: usize) -> ContigMethylation {
        let mut table = MethylationTable::empty(len);
        for &p in bottom_meth_positions {
            table.set_bottom(p as usize, true);
        }
        ContigMethylation::from_tables(vec![table])
    }

    /// All-reference haplotype (no variants) so haplotype coords == reference coords.
    fn ref_haplotype() -> crate::haplotype::Haplotype {
        let haps = build_haplotypes(&[], 1, &mut SmallRng::seed_from_u64(0));
        haps.into_iter().next().unwrap()
    }

    #[test]
    fn test_record_mate_ct_fragment_methylated_top() {
        // Reference: ACGTACG (CpGs at 1, 5). Mate covers positions [0, 6].
        // CT fragment: tally top[hap_pos] at sites 1 and 5.
        // Bitmap: top[1]=true, top[5]=false.
        let cm = cm_with_top(&[1], 7);
        let hap = ref_haplotype();
        let cpgs = vec![1u32, 5];
        let mate_positions: Vec<u32> = (0u32..7).collect();

        let mut tally = CpgTruthTally::new();
        tally.record_mate(0, &mate_positions, &cpgs, &cm, &hap, 0, true);

        assert_eq!(tally.get(0, 1), Some((1, 0)), "site 1 → 1 meth");
        assert_eq!(tally.get(0, 5), Some((0, 1)), "site 5 → 1 unmeth");
        assert_eq!(tally.len(), 2);
    }

    #[test]
    fn test_record_mate_ga_fragment_methylated_bottom() {
        // GA fragment: queries bottom[site_p + 1].
        // Bitmap: bottom[2]=true (covers site_p=1), bottom[6]=false (covers site_p=5).
        let cm = cm_with_bottom(&[2], 7);
        let hap = ref_haplotype();
        let cpgs = vec![1u32, 5];
        let mate_positions: Vec<u32> = (0u32..7).collect();

        let mut tally = CpgTruthTally::new();
        tally.record_mate(0, &mate_positions, &cpgs, &cm, &hap, 0, false);

        assert_eq!(tally.get(0, 1), Some((1, 0)), "site 1 (GA) → 1 meth");
        assert_eq!(tally.get(0, 5), Some((0, 1)), "site 5 (GA) → 1 unmeth");
    }

    #[test]
    fn test_record_mate_only_counts_cpgs_in_span() {
        // Mate covers ref positions [3, 6] only. CpG at 1 is OUT of range,
        // CpG at 5 is IN range.
        let cm = cm_with_top(&[5], 10);
        let hap = ref_haplotype();
        let cpgs = vec![1u32, 5];
        let mate_positions: Vec<u32> = (3u32..7).collect();

        let mut tally = CpgTruthTally::new();
        tally.record_mate(0, &mate_positions, &cpgs, &cm, &hap, 0, true);

        assert_eq!(tally.get(0, 1), None, "out-of-span CpG must not be tallied");
        assert_eq!(tally.get(0, 5), Some((1, 0)));
    }

    #[test]
    fn test_record_mate_ga_excludes_site_when_bottom_c_falls_outside_span() {
        // GA fragment: site_p=5 has bottom-C at top-pos 6. If mate spans
        // [3, 5], the bottom-C (6) is NOT covered → must not tally.
        let cm = cm_with_bottom(&[6], 10);
        let hap = ref_haplotype();
        let cpgs = vec![5u32];
        let mate_positions: Vec<u32> = (3u32..6).collect(); // 3, 4, 5 — no 6.

        let mut tally = CpgTruthTally::new();
        tally.record_mate(0, &mate_positions, &cpgs, &cm, &hap, 0, false);

        assert_eq!(tally.get(0, 5), None);
    }

    #[test]
    fn test_record_mate_aggregates_multiple_calls_at_same_site() {
        let cm = cm_with_top(&[1], 7);
        let hap = ref_haplotype();
        let cpgs = vec![1u32];
        let mate_positions: Vec<u32> = (0u32..7).collect();

        let mut tally = CpgTruthTally::new();
        tally.record_mate(0, &mate_positions, &cpgs, &cm, &hap, 0, true);
        tally.record_mate(0, &mate_positions, &cpgs, &cm, &hap, 0, true);
        tally.record_mate(0, &mate_positions, &cpgs, &cm, &hap, 0, true);

        assert_eq!(tally.get(0, 1), Some((3, 0)), "three meth calls aggregate to 3/0");
    }

    #[test]
    fn test_write_bedgraph_format() {
        // Build a tally with two sites and write to a temp file. Verify the
        // exact MethylDackel-format output.
        let mut tally = CpgTruthTally::new();
        tally.counts.insert((0, 1), (3, 1));
        tally.counts.insert((0, 5), (0, 4));

        let dict =
            SequenceDictionary::from_entries(vec![crate::sequence_dict::SequenceMetadata::new(
                0,
                "chr1".to_string(),
                100,
            )]);

        let tmp =
            std::env::temp_dir().join(format!("holodeck_cpgtruth_{}.bedGraph", std::process::id()));
        tally.write_bedgraph(&dict, &tmp).unwrap();

        let written = std::fs::read_to_string(&tmp).unwrap();
        std::fs::remove_file(&tmp).ok();

        let lines: Vec<&str> = written.lines().collect();
        assert_eq!(lines.len(), 3, "header + 2 data rows");
        assert!(lines[0].starts_with("track "), "first line must be a bedGraph track header");
        // Site at p=1: rate = round(100 * 3/4) = 75.
        assert_eq!(lines[1], "chr1\t1\t2\t75\t3\t1");
        // Site at p=5: rate = round(100 * 0/4) = 0.
        assert_eq!(lines[2], "chr1\t5\t6\t0\t0\t4");
    }
}
