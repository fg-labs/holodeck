use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use anyhow::{Context, Result, bail};
use coitrees::{COITree, Interval, IntervalTree};
use rand::Rng;

use crate::sequence_dict::SequenceDictionary;

/// The gzip magic number (first two bytes of any gzip/bgzip file).
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// A collection of genomic target regions loaded from a BED file, indexed for
/// efficient overlap queries.
///
/// Internally stores one [`COITree`] per contig for O(log N) overlap detection.
/// BED coordinates are 0-based half-open `[start, end)`, but coitrees uses
/// end-inclusive intervals, so we store `[start, end-1]` internally.
///
/// `Debug` is implemented manually because `COITree` does not implement `Debug`.
pub struct TargetRegions {
    /// One interval tree per contig, indexed by the contig's position in the
    /// sequence dictionary. Empty trees for contigs with no targets.
    trees: Vec<COITree<(), u32>>,
    /// Total bases covered by all target regions.
    total_territory: u64,
    /// Per-contig target territory (bases), indexed by contig position.
    per_contig_territory: Vec<u64>,
    /// Sorted intervals per contig in 0-based half-open coordinates `[start, end)`,
    /// used to build padded sampling regions for fragment start position selection.
    sorted_intervals: Vec<Vec<(u32, u32)>>,
    /// Sequence dictionary for contig lookups.
    dict: SequenceDictionary,
}

impl std::fmt::Debug for TargetRegions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("TargetRegions")
            .field("num_contigs", &self.trees.len())
            .field("total_territory", &self.total_territory)
            .field("dict", &self.dict)
            .finish_non_exhaustive()
    }
}

impl TargetRegions {
    /// Load target regions from a BED file.
    ///
    /// The file may be plain text or gzipped. Coordinates are validated against
    /// `dict` — unknown contigs or out-of-range coordinates cause an error.
    ///
    /// # Errors
    /// Returns an error if the file cannot be read, contains invalid entries,
    /// or references contigs not present in the dictionary.
    pub fn from_path(path: &Path, dict: &SequenceDictionary) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open BED file: {}", path.display()))?;

        // Detect gzip by magic bytes, then re-open to read from the start.
        let mut magic = [0u8; 2];
        let is_gzipped = {
            let mut peek = BufReader::new(file);
            peek.read_exact(&mut magic).is_ok() && magic == GZIP_MAGIC
        };

        let file = File::open(path)?;
        let reader: Box<dyn BufRead> = if is_gzipped {
            Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        let mut intervals_by_contig: Vec<Vec<Interval<()>>> = vec![Vec::new(); dict.len()];
        let mut raw_intervals_by_contig: Vec<Vec<(u32, u32)>> = vec![Vec::new(); dict.len()];
        let mut total_territory: u64 = 0;
        let mut per_contig_territory: Vec<u64> = vec![0; dict.len()];

        for (line_num, line) in reader.lines().enumerate() {
            let line =
                line.with_context(|| format!("Failed to read line {} of BED file", line_num + 1))?;
            let line = line.trim();
            if line.is_empty()
                || line.starts_with('#')
                || line.starts_with("track ")
                || line.starts_with("browser ")
            {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 {
                bail!("BED line {} has fewer than 3 fields: {line}", line_num + 1);
            }

            let contig = fields[0];
            let start: u32 = fields[1].parse().with_context(|| {
                format!("Invalid start coordinate on BED line {}: {}", line_num + 1, fields[1])
            })?;
            let end: u32 = fields[2].parse().with_context(|| {
                format!("Invalid end coordinate on BED line {}: {}", line_num + 1, fields[2])
            })?;

            if start >= end {
                bail!("BED line {} has start >= end: {start} >= {end}", line_num + 1);
            }

            let meta = dict.get_by_name(contig).ok_or_else(|| {
                anyhow::anyhow!(
                    "BED line {} references unknown contig '{contig}'. \
                     Ensure the BED file matches the reference FASTA.",
                    line_num + 1
                )
            })?;

            #[expect(clippy::cast_possible_truncation, reason = "contig lengths fit in u32")]
            let contig_len = meta.length() as u32;
            if end > contig_len {
                bail!(
                    "BED line {} has end ({end}) > contig length ({contig_len}) for '{contig}'",
                    line_num + 1
                );
            }

            // coitrees uses end-inclusive: BED [start, end) → coitrees [start, end-1]
            #[expect(clippy::cast_possible_wrap, reason = "genomic coords < i32::MAX")]
            let iv = Interval::new(start as i32, (end - 1) as i32, ());
            intervals_by_contig[meta.index()].push(iv);
            raw_intervals_by_contig[meta.index()].push((start, end));
            let bases = u64::from(end - start);
            total_territory += bases;
            per_contig_territory[meta.index()] += bases;
        }

        let trees: Vec<COITree<(), u32>> = intervals_by_contig.iter().map(COITree::new).collect();

        // Sort raw intervals per contig for padded sampling region construction.
        let sorted_intervals: Vec<Vec<(u32, u32)>> = raw_intervals_by_contig
            .into_iter()
            .map(|mut ivs| {
                ivs.sort_unstable();
                ivs
            })
            .collect();

        Ok(Self {
            trees,
            total_territory,
            per_contig_territory,
            sorted_intervals,
            dict: dict.clone(),
        })
    }

    /// Return the total number of bases covered by all target regions.
    ///
    /// Note: if the BED file contains overlapping intervals, overlapping bases
    /// are counted multiple times. Callers that need exact territory should
    /// provide a non-overlapping BED.
    #[must_use]
    pub fn total_territory(&self) -> u64 {
        self.total_territory
    }

    /// Return the target territory (in bases) for a specific contig.
    #[must_use]
    pub fn contig_territory(&self, contig_index: usize) -> u64 {
        self.per_contig_territory.get(contig_index).copied().unwrap_or(0)
    }

    /// Check whether the interval `[start, end)` (0-based half-open) on the
    /// given contig overlaps any target region.
    #[must_use]
    #[expect(clippy::cast_possible_wrap, reason = "genomic coords < i32::MAX")]
    pub fn overlaps(&self, contig_index: usize, start: u32, end: u32) -> bool {
        self.trees
            .get(contig_index)
            .is_some_and(|tree| tree.query_count(start as i32, (end.saturating_sub(1)) as i32) > 0)
    }

    /// Return the sorted target intervals for a contig in 0-based half-open
    /// coordinates `[start, end)`.
    #[must_use]
    pub fn contig_intervals(&self, contig_index: usize) -> &[(u32, u32)] {
        self.sorted_intervals.get(contig_index).map_or(&[], Vec::as_slice)
    }

    /// Compute the effective territory for coverage calculation, accounting for
    /// the fact that fragments extend beyond target boundaries.
    ///
    /// For a target of width W, a fragment of length L placed uniformly to
    /// overlap the target has an on-target fraction of `W / (W + L - 1)`.
    /// The effective territory per target is therefore `W + L - 1` — the
    /// catchment zone of fragment start positions — and the total effective
    /// territory is the sum across all targets.  Using this as `effective_size`
    /// in the standard coverage formula `N = C * effective_size / bases_per_read`
    /// yields the correct number of reads for the desired mean target coverage.
    #[must_use]
    pub fn effective_territory(&self, fragment_mean: usize) -> u64 {
        let l_minus_1 = fragment_mean.saturating_sub(1) as u64;
        self.sorted_intervals
            .iter()
            .flat_map(|ivs| ivs.iter())
            .map(|&(start, end)| u64::from(end - start) + l_minus_1)
            .sum()
    }

    /// Compute the effective territory for a single contig.
    ///
    /// See [`effective_territory`](Self::effective_territory) for the rationale.
    #[must_use]
    pub fn contig_effective_territory(&self, contig_index: usize, fragment_mean: usize) -> u64 {
        let l_minus_1 = fragment_mean.saturating_sub(1) as u64;
        self.sorted_intervals.get(contig_index).map_or(0, |ivs| {
            ivs.iter().map(|&(start, end)| u64::from(end - start) + l_minus_1).sum()
        })
    }

    /// Return a reference to the underlying sequence dictionary.
    #[must_use]
    pub fn dict(&self) -> &SequenceDictionary {
        &self.dict
    }
}

/// A sampler that draws random fragment start positions from padded target regions.
///
/// Given a set of target intervals, each interval is padded on the left by a
/// specified amount (to include fragments that start before a target but extend
/// into it), overlapping padded intervals are merged, and start positions are
/// sampled uniformly across the merged regions.  The resulting fragments should
/// still be checked for overlap with the original unpadded targets — starts in
/// the pad zone whose drawn fragment length is too short to reach the target
/// will be rejected, but this is rare when the pad matches the expected max
/// fragment length.
pub struct PaddedIntervalSampler {
    /// Merged padded intervals in 0-based half-open coordinates.
    intervals: Vec<(u32, u32)>,
    /// Cumulative territory sums for binary-search sampling.
    /// `cumulative[j]` = total bases in intervals `0..=j`.
    cumulative: Vec<u64>,
    /// Total bases across all padded intervals.
    total: u64,
}

impl PaddedIntervalSampler {
    /// Build a sampler from sorted target intervals, padding each on the left.
    ///
    /// `pad` is typically the maximum expected fragment length (e.g.
    /// `fragment_mean + 4 * fragment_stddev`), ensuring that fragments starting
    /// before a target but overlapping it are represented in the sampling space.
    /// Padded intervals are merged where they overlap and clamped to
    /// `[0, contig_len)`.
    #[must_use]
    pub fn new(intervals: &[(u32, u32)], pad: u32, contig_len: u32) -> Self {
        if intervals.is_empty() {
            return Self { intervals: Vec::new(), cumulative: Vec::new(), total: 0 };
        }

        // Pad each interval on the left and clamp to [0, contig_len).
        let mut padded: Vec<(u32, u32)> = intervals
            .iter()
            .map(|&(start, end)| (start.saturating_sub(pad), end.min(contig_len)))
            .collect();
        padded.sort_unstable();

        // Merge overlapping or abutting intervals.
        let mut merged: Vec<(u32, u32)> = Vec::with_capacity(padded.len());
        for (start, end) in padded {
            if let Some(last) = merged.last_mut()
                && start <= last.1
            {
                last.1 = last.1.max(end);
                continue;
            }
            merged.push((start, end));
        }

        // Build cumulative territory sums.
        let mut cumulative = Vec::with_capacity(merged.len());
        let mut running = 0u64;
        for &(start, end) in &merged {
            running += u64::from(end - start);
            cumulative.push(running);
        }
        let total = running;

        Self { intervals: merged, cumulative, total }
    }

    /// Sample a random start position uniformly from the padded intervals.
    ///
    /// Returns `None` if there are no intervals.
    pub fn sample_start(&self, rng: &mut impl Rng) -> Option<u32> {
        if self.total == 0 {
            return None;
        }

        let r = rng.random_range(0..self.total);
        let idx = self.cumulative.partition_point(|&c| c <= r);
        let (start, _end) = self.intervals[idx];
        let base_before = if idx > 0 { self.cumulative[idx - 1] } else { 0 };
        let offset = r - base_before;

        #[expect(clippy::cast_possible_truncation, reason = "offset within interval fits u32")]
        Some(start + offset as u32)
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use rand::SeedableRng;
    use tempfile::NamedTempFile;

    use super::*;
    use crate::sequence_dict::SequenceMetadata;

    /// Build a dict for testing.
    fn test_dict() -> SequenceDictionary {
        // We need to construct a dict. Use the same test helper pattern.
        let sequences = vec![
            SequenceMetadata::new(0, "chr1".to_string(), 10000),
            SequenceMetadata::new(1, "chr2".to_string(), 5000),
        ];
        SequenceDictionary::from_entries(sequences)
    }

    /// Write BED content to a temp file and return the path.
    fn write_bed(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f.flush().unwrap();
        f
    }

    #[test]
    fn test_load_simple_bed() {
        let dict = test_dict();
        let bed = write_bed("chr1\t100\t200\nchr1\t300\t400\nchr2\t50\t150\n");
        let regions = TargetRegions::from_path(bed.path(), &dict).unwrap();

        assert_eq!(regions.total_territory(), 300); // 100 + 100 + 100
    }

    #[test]
    fn test_overlap_hit() {
        let dict = test_dict();
        let bed = write_bed("chr1\t100\t200\n");
        let regions = TargetRegions::from_path(bed.path(), &dict).unwrap();

        // Fragment fully within target
        assert!(regions.overlaps(0, 120, 180));
        // Fragment starts before, ends within
        assert!(regions.overlaps(0, 50, 150));
        // Fragment starts within, ends after
        assert!(regions.overlaps(0, 150, 250));
        // Fragment engulfs target
        assert!(regions.overlaps(0, 0, 300));
    }

    #[test]
    fn test_overlap_miss() {
        let dict = test_dict();
        let bed = write_bed("chr1\t100\t200\n");
        let regions = TargetRegions::from_path(bed.path(), &dict).unwrap();

        // Fragment entirely before target
        assert!(!regions.overlaps(0, 0, 100));
        // Fragment entirely after target
        assert!(!regions.overlaps(0, 200, 300));
        // Wrong contig
        assert!(!regions.overlaps(1, 100, 200));
    }

    #[test]
    fn test_overlap_single_base() {
        let dict = test_dict();
        let bed = write_bed("chr1\t100\t200\n");
        let regions = TargetRegions::from_path(bed.path(), &dict).unwrap();

        // Single-base overlap at target start
        assert!(regions.overlaps(0, 99, 101));
        // Single-base overlap at target end
        assert!(regions.overlaps(0, 199, 201));
        // Adjacent but not overlapping
        assert!(!regions.overlaps(0, 200, 201));
    }

    #[test]
    fn test_skips_comments_and_blank_lines() {
        let dict = test_dict();
        let bed = write_bed("# header\n\nchr1\t100\t200\n\n");
        let regions = TargetRegions::from_path(bed.path(), &dict).unwrap();
        assert_eq!(regions.total_territory(), 100);
    }

    #[test]
    fn test_error_unknown_contig() {
        let dict = test_dict();
        let bed = write_bed("chrZ\t100\t200\n");
        let result = TargetRegions::from_path(bed.path(), &dict);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("unknown contig"));
    }

    #[test]
    fn test_error_start_gte_end() {
        let dict = test_dict();
        let bed = write_bed("chr1\t200\t100\n");
        let result = TargetRegions::from_path(bed.path(), &dict);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("start >= end"));
    }

    #[test]
    fn test_error_end_exceeds_contig_length() {
        let dict = test_dict();
        let bed = write_bed("chr1\t9000\t20000\n");
        let result = TargetRegions::from_path(bed.path(), &dict);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("contig length"));
    }

    #[test]
    fn test_effective_territory_single_target() {
        let dict = test_dict();
        let bed = write_bed("chr1\t100\t200\n"); // 100bp target
        let regions = TargetRegions::from_path(bed.path(), &dict).unwrap();

        // effective = W + L - 1 = 100 + 375 - 1 = 474
        assert_eq!(regions.effective_territory(375), 474);
        // With fragment_mean = 1: effective = W + 0 = 100
        assert_eq!(regions.effective_territory(1), 100);
    }

    #[test]
    fn test_effective_territory_multiple_targets() {
        let dict = test_dict();
        // Two 100bp targets on chr1, one 50bp target on chr2.
        let bed = write_bed("chr1\t100\t200\nchr1\t500\t600\nchr2\t0\t50\n");
        let regions = TargetRegions::from_path(bed.path(), &dict).unwrap();

        // 3 targets: sum of (W_i + L - 1) = (100+374) + (100+374) + (50+374) = 1372
        assert_eq!(regions.effective_territory(375), 1372);
    }

    #[test]
    fn test_contig_effective_territory() {
        let dict = test_dict();
        let bed = write_bed("chr1\t100\t200\nchr1\t500\t600\nchr2\t0\t50\n");
        let regions = TargetRegions::from_path(bed.path(), &dict).unwrap();

        // chr1: (100+374) + (100+374) = 948
        assert_eq!(regions.contig_effective_territory(0, 375), 948);
        // chr2: (50+374) = 424
        assert_eq!(regions.contig_effective_territory(1, 375), 424);
    }

    #[test]
    fn test_contig_intervals_returns_sorted_intervals() {
        let dict = test_dict();
        // Intervals intentionally out of order in the BED file.
        let bed = write_bed("chr1\t300\t400\nchr1\t100\t200\nchr2\t50\t150\n");
        let regions = TargetRegions::from_path(bed.path(), &dict).unwrap();

        let chr1_ivs = regions.contig_intervals(0);
        assert_eq!(chr1_ivs, &[(100, 200), (300, 400)]);

        let chr2_ivs = regions.contig_intervals(1);
        assert_eq!(chr2_ivs, &[(50, 150)]);
    }

    #[test]
    fn test_contig_intervals_empty_contig() {
        let dict = test_dict();
        let bed = write_bed("chr1\t100\t200\n");
        let regions = TargetRegions::from_path(bed.path(), &dict).unwrap();
        assert!(regions.contig_intervals(1).is_empty());
    }

    // -- PaddedIntervalSampler tests --

    #[test]
    fn test_sampler_empty_intervals() {
        let sampler = PaddedIntervalSampler::new(&[], 100, 10000);
        let mut rng = rand::rngs::SmallRng::seed_from_u64(42);
        assert!(sampler.sample_start(&mut rng).is_none());
    }

    #[test]
    fn test_sampler_single_interval_no_pad() {
        let sampler = PaddedIntervalSampler::new(&[(100, 200)], 0, 10000);
        let mut rng = rand::rngs::SmallRng::seed_from_u64(42);

        for _ in 0..1000 {
            let pos = sampler.sample_start(&mut rng).unwrap();
            assert!((100..200).contains(&pos), "pos {pos} not in [100, 200)");
        }
    }

    #[test]
    fn test_sampler_padding_extends_left() {
        // Target at [500, 600), pad of 200 → sampling region [300, 600).
        let sampler = PaddedIntervalSampler::new(&[(500, 600)], 200, 10000);
        let mut rng = rand::rngs::SmallRng::seed_from_u64(42);

        let mut min_seen = u32::MAX;
        let mut max_seen = 0u32;
        for _ in 0..10_000 {
            let pos = sampler.sample_start(&mut rng).unwrap();
            assert!((300..600).contains(&pos), "pos {pos} not in [300, 600)");
            min_seen = min_seen.min(pos);
            max_seen = max_seen.max(pos);
        }

        // With 10k samples across 300 positions, we should see positions near
        // both edges.
        assert!(min_seen <= 310, "min_seen {min_seen} too high");
        assert!(max_seen >= 590, "max_seen {max_seen} too low");
    }

    #[test]
    fn test_sampler_padding_clamped_to_zero() {
        // Target near contig start: [50, 150), pad of 200 → [0, 150).
        let sampler = PaddedIntervalSampler::new(&[(50, 150)], 200, 10000);
        let mut rng = rand::rngs::SmallRng::seed_from_u64(42);

        for _ in 0..1000 {
            let pos = sampler.sample_start(&mut rng).unwrap();
            assert!(pos < 150, "pos {pos} not in [0, 150)");
        }
    }

    #[test]
    fn test_sampler_merges_overlapping_padded_intervals() {
        // Two intervals [200, 300) and [350, 450) with pad 100.
        // Padded: [100, 300) and [250, 450) → merged: [100, 450).
        let sampler = PaddedIntervalSampler::new(&[(200, 300), (350, 450)], 100, 10000);
        let mut rng = rand::rngs::SmallRng::seed_from_u64(42);

        for _ in 0..1000 {
            let pos = sampler.sample_start(&mut rng).unwrap();
            assert!((100..450).contains(&pos), "pos {pos} not in [100, 450)");
        }
    }

    #[test]
    fn test_sampler_keeps_disjoint_padded_intervals_separate() {
        // Two intervals [100, 150) and [1000, 1050) with pad 50.
        // Padded: [50, 150) and [950, 1050) — no overlap, should stay separate.
        let sampler = PaddedIntervalSampler::new(&[(100, 150), (1000, 1050)], 50, 10000);
        let mut rng = rand::rngs::SmallRng::seed_from_u64(42);

        for _ in 0..1000 {
            let pos = sampler.sample_start(&mut rng).unwrap();
            let in_first = (50..150).contains(&pos);
            let in_second = (950..1050).contains(&pos);
            assert!(in_first || in_second, "pos {pos} not in either padded interval");
        }
    }

    #[test]
    fn test_sampler_samples_proportional_to_interval_size() {
        // One large interval and one small interval.  Samples should be roughly
        // proportional to padded sizes.
        // [1000, 2000) pad 100 → [900, 2000) = 1100 bp
        // [5000, 5010) pad 100 → [4900, 5010) = 110 bp
        // Ratio should be ~10:1.
        let sampler = PaddedIntervalSampler::new(&[(1000, 2000), (5000, 5010)], 100, 10000);
        let mut rng = rand::rngs::SmallRng::seed_from_u64(42);

        let mut count_first = 0u32;
        let mut count_second = 0u32;
        for _ in 0..11_000 {
            let pos = sampler.sample_start(&mut rng).unwrap();
            if (900..2000).contains(&pos) {
                count_first += 1;
            } else {
                count_second += 1;
            }
        }

        let ratio = f64::from(count_first) / f64::from(count_second);
        assert!(
            (8.0..12.0).contains(&ratio),
            "ratio {ratio:.1} not near expected 10:1 (first={count_first}, second={count_second})"
        );
    }
}
