//! Methylation bitmaps and chemistry-conversion logic.
//!
//! Provides per-haplotype methylation state via [`MethylationTable`] and
//! [`ContigMethylation`], plus the [`apply_methylation_conversion`] free
//! function that simulates the per-base chemistry of either an em-seq /
//! bisulfite library (unmethylated C → T) or a TAPS library (methylated C → T).
//!
//! # Per-haplotype CpG detection and per-strand bitmaps
//!
//! Each haplotype gets its own pair of [`BitVec`]s indexed by haplotype
//! position (0..haplotype_length). For every `CG` dinucleotide on that
//! haplotype's materialized sequence (case-insensitive), two independent
//! Bernoulli draws decide whether each strand's C is methylated:
//!
//! - `top[h]` = "the top-strand C at haplotype position `h` is methylated."
//! - `bottom[h + 1]` = "the bottom-strand C at haplotype position `h + 1` is
//!   methylated."
//!
//! Indexing by haplotype position rather than reference position naturally
//! handles SNPs, insertions, and deletions that create or destroy CpG sites
//! on a particular haplotype: each haplotype's bitmap reflects the CpG
//! context that actually exists on that haplotype.
//!
//! Both bitmaps for a haplotype have length equal to the haplotype's
//! materialized length; positions that don't host a strand-specific C (or
//! that host a non-CpG cytosine) always read `false`. Hemimethylation is
//! allowed because the two strands' draws are independent. Allele-specific
//! methylation falls out naturally because each haplotype draws independently.
//!
//! Non-CpG cytosines are always treated as unmethylated.
//!
//! # Chemistry modes
//!
//! [`MethylationMode`] selects which class of cytosines is converted to
//! thymine during chemistry simulation:
//!
//! - [`MethylationMode::EmSeq`] -- unmethylated cytosines convert to thymine;
//!   methylated cytosines are preserved. Matches both classical bisulfite
//!   chemistry and enzymatic methyl-seq (em-seq, NEBNext) -- the conversion
//!   patterns are identical.
//! - [`MethylationMode::Taps`] -- methylated cytosines convert to thymine
//!   (TET oxidation + pyridine borane); unmethylated cytosines are preserved.
//!   The inverse of em-seq: a `C→T` event at a CpG signals methylation.

use rand::Rng;

use bitvec::vec::BitVec;

/// Per-haplotype methylation state, one bitmap per strand. Bitmaps are
/// indexed by **haplotype position** (which may differ from reference
/// position when the haplotype contains indels).
///
/// The strand state is stored in [`bitvec::vec::BitVec`] rather than
/// `Vec<bool>` because the bitmaps are sized to the full materialized
/// haplotype length — one bit per base, per strand, per haplotype. At
/// whole-chromosome scale (e.g. ~250 Mb × 2 strands × ploidy) `Vec<bool>`
/// would cost 8× the memory for the same information; `bitvec` packs it to
/// one bit each. That density justifies the extra direct dependency.
#[derive(Debug, Clone)]
pub struct MethylationTable {
    /// Top-strand methylation bitmap.
    /// `top[h]` is `true` iff the top-strand C at haplotype position `h` is
    /// methylated. Length equals the haplotype's materialized length;
    /// positions without a strand-specific C in CpG context hold `false`.
    top: BitVec,
    /// Bottom-strand methylation bitmap, same shape as `top` but for the
    /// reverse-complement strand.
    bottom: BitVec,
}

impl MethylationTable {
    /// Build an empty methylation table of the given length (no methylation
    /// anywhere). Test-only: production code always builds tables via
    /// [`Self::from_haplotype`] / [`ContigMethylation::from_haplotypes`].
    #[cfg(test)]
    #[must_use]
    pub(crate) fn empty(len: usize) -> Self {
        Self::with_len(len)
    }

    /// Build an empty methylation table with the given length (all bits
    /// `false`). Used by [`Self::from_haplotype`] and the VCF reader.
    #[must_use]
    pub(crate) fn with_len(len: usize) -> Self {
        Self { top: BitVec::repeat(false, len), bottom: BitVec::repeat(false, len) }
    }

    /// Build a methylation table for a single haplotype by materializing the
    /// haplotype's full sequence (via [`crate::haplotype::Haplotype::extract_fragment`]
    /// over the entire contig length) and scanning the result for CpG
    /// dinucleotides. The resulting bitmap length equals the materialized
    /// haplotype length, which may differ from the reference length when
    /// the haplotype carries indels.
    ///
    /// At every `CG` dinucleotide on the haplotype's top strand, the
    /// top-strand C and the bottom-strand C are independently drawn from a
    /// Bernoulli(`methylation_rate`) distribution.
    ///
    /// # Panics
    ///
    /// Panics if `methylation_rate` is not a finite value in `[0.0, 1.0]`.
    /// The `methylate` CLI validates this at startup, but this is a `pub`
    /// constructor reachable from library code; without the guard a `NaN`
    /// would silently behave like `0.0` (the `rng < rate` comparison is
    /// always false) and values `> 1.0` like `1.0`, quietly breaking the
    /// documented `Bernoulli(methylation_rate)` contract.
    pub fn from_haplotype(
        haplotype: &crate::haplotype::Haplotype,
        reference: &[u8],
        methylation_rate: f64,
        rng: &mut impl Rng,
    ) -> Self {
        assert!(
            methylation_rate.is_finite() && (0.0..=1.0).contains(&methylation_rate),
            "methylation_rate must be a finite value in [0.0, 1.0]; got {methylation_rate}"
        );
        // Materialize the entire haplotype as one large fragment. The cap
        // passed to `extract_fragment` is BOTH a pre-allocation hint AND a
        // truncation limit on the output base count, so it must be large
        // enough to hold the full materialized haplotype — including all
        // net-positive insertions. Using a tighter cap (e.g.
        // `reference.len() + 1`) silently drops haplotype suffix bases when
        // the haplotype contains insertions adding more than one base, which
        // can lose CpG sites entirely.
        //
        // `hap_position_for(reference.len())` returns the haplotype-coordinate
        // length of the materialized haplotype (sum of `(alt_len - ref_len)`
        // across all variants whose `var_end <= reference.len()`, plus
        // `reference.len()` itself), giving the exact required capacity for
        // sane VCFs whose variants do not extend past the chromosome end.
        #[expect(clippy::cast_possible_truncation, reason = "reference length fits in u32")]
        let cap = haplotype.hap_position_for(reference.len() as u32) as usize;
        let (hap_bases, _ref_positions, _hap_start) = haplotype.extract_fragment(reference, 0, cap);
        let len = hap_bases.len();
        let mut table = Self::with_len(len);
        if len < 2 {
            return table;
        }
        for i in 0..len - 1 {
            let c0 = hap_bases[i].to_ascii_uppercase();
            let c1 = hap_bases[i + 1].to_ascii_uppercase();
            if c0 == b'C' && c1 == b'G' {
                if rng.random::<f64>() < methylation_rate {
                    table.top.set(i, true);
                }
                if rng.random::<f64>() < methylation_rate {
                    table.bottom.set(i + 1, true);
                }
            }
        }
        table
    }

    /// Return whether the C at haplotype position `pos` is methylated on
    /// the strand the read came from. For a forward-strand read, queries
    /// the top-strand bitmap; for a negative-strand read, queries the
    /// bottom-strand bitmap. Returns `false` for out-of-range positions.
    #[must_use]
    pub fn is_methylated(&self, pos: u32, is_negative_strand: bool) -> bool {
        let bv = if is_negative_strand { &self.bottom } else { &self.top };
        bv.get(pos as usize).is_some_and(|b| *b)
    }

    /// Length of each per-strand bitmap (equals the haplotype's materialized
    /// length, or the length passed to [`Self::empty`]). Test-only.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn len(&self) -> usize {
        self.top.len()
    }

    /// Whether both bitmaps are zero length. Test-only; paired with
    /// [`Self::len`] to satisfy clippy's `len_without_is_empty`.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn is_empty(&self) -> bool {
        self.top.is_empty()
    }

    /// Set a top-strand methylation bit at position `pos`. Panics on
    /// out-of-range index. Used by the VCF reader and tests.
    pub(crate) fn set_top(&mut self, pos: usize, value: bool) {
        self.top.set(pos, value);
    }

    /// Set a bottom-strand methylation bit at position `pos`. Panics on
    /// out-of-range index. Used by the VCF reader and tests.
    pub(crate) fn set_bottom(&mut self, pos: usize, value: bool) {
        self.bottom.set(pos, value);
    }
}

/// Per-contig methylation state covering ALL haplotypes for one contig.
/// Indexed by [`crate::fragment::Fragment::haplotype_index`].
#[derive(Debug, Clone)]
pub struct ContigMethylation {
    /// One [`MethylationTable`] per haplotype, in haplotype-index order.
    per_haplotype: Vec<MethylationTable>,
}

impl ContigMethylation {
    /// Build per-haplotype methylation tables by scanning each haplotype's
    /// materialized sequence for CpG dinucleotides. Methylation draws are
    /// independent per (haplotype, strand, position), allowing both
    /// hemimethylation and allele-specific methylation.
    pub fn from_haplotypes(
        haplotypes: &[crate::haplotype::Haplotype],
        reference: &[u8],
        methylation_rate: f64,
        rng: &mut impl Rng,
    ) -> Self {
        let per_haplotype = haplotypes
            .iter()
            .map(|hap| MethylationTable::from_haplotype(hap, reference, methylation_rate, rng))
            .collect();
        Self { per_haplotype }
    }

    /// Construct from a pre-built per-haplotype table list. Used by the VCF
    /// reader and tests; production code uses [`Self::from_haplotypes`].
    #[must_use]
    pub(crate) fn from_tables(per_haplotype: Vec<MethylationTable>) -> Self {
        Self { per_haplotype }
    }

    /// Return the methylation table for the haplotype at the given index.
    /// Panics if the index is out of range -- `Fragment::haplotype_index`
    /// is always derived from a haplotype actually built for the contig,
    /// so an out-of-range index is a programming error.
    #[must_use]
    pub fn table_for(&self, haplotype_index: usize) -> &MethylationTable {
        &self.per_haplotype[haplotype_index]
    }

    /// Number of haplotypes covered by this `ContigMethylation`.
    #[must_use]
    pub fn len(&self) -> usize {
        self.per_haplotype.len()
    }

    /// Whether there are no haplotypes covered.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.per_haplotype.is_empty()
    }
}

/// Methylation chemistry. Selects which class of cytosines is converted
/// to thymine during chemistry simulation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum MethylationMode {
    /// Unmethylated cytosines convert to thymine; methylated cytosines
    /// are preserved. Matches classical bisulfite chemistry and enzymatic
    /// methyl-seq. Pass `--methylation-mode bisulfite` as an alias.
    #[value(alias = "bisulfite")]
    EmSeq,
    /// Methylated cytosines convert to thymine (TET oxidation + pyridine
    /// borane); unmethylated cytosines are preserved. The inverse of
    /// bisulfite/em-seq: a `C→T` event at a CpG signals methylation.
    Taps,
}

impl MethylationMode {
    /// Canonical seed-string form of this mode. Used in
    /// [`crate::commands::simulate::Simulate::compute_seed`] so changes to the
    /// CLI alias system don't accidentally shift seed determinism.
    #[must_use]
    pub fn as_seed_str(self) -> &'static str {
        match self {
            Self::EmSeq => "em-seq",
            Self::Taps => "taps",
        }
    }
}

/// Configuration bundle passed through the read-generation pipeline when
/// methylation chemistry simulation is enabled. Carries the per-contig,
/// per-haplotype methylation state plus the chemistry parameters.
#[derive(Debug, Clone, Copy)]
pub struct MethylationConfig<'a> {
    /// Per-contig methylation tables, one per haplotype.
    pub contig_methylation: &'a ContigMethylation,
    /// Which chemistry to apply (em-seq vs TAPS).
    pub mode: MethylationMode,
    /// Probability that a qualifying C is converted to T in a molecule that
    /// converted normally. Clamped to `[0.0, 1.0]` by the caller; not
    /// re-validated in the hot loop.
    pub conversion_rate: f64,
    /// Per-molecule probability that a fragment is a *conversion failure*.
    /// Real bisulfite/EM-seq conversion is effectively bimodal: most
    /// molecules convert near-completely, a small fraction escape conversion
    /// as a unit (fragments that fail to denature, or re-anneal too fast).
    /// A failed molecule converts its should-convert cytosines at
    /// `1.0 - conversion_rate`, so it coherently retains almost all of them
    /// as C. Drawn once per fragment so both mates agree. Clamped to
    /// `[0.0, 1.0]` by the caller.
    ///
    /// The "near-zero failed rate" intuition assumes `conversion_rate` is
    /// close to `1.0` (the realistic regime). At the degenerate setting
    /// `conversion_rate == 0.0` the relationship inverts — failed molecules
    /// convert at `1.0` (fully) while normal molecules don't convert at all —
    /// which is a deliberate consequence of pinning the failed rate to
    /// `1.0 - conversion_rate`, not a special case.
    pub failure_rate: f64,
}

/// Reference CpG positions: the 0-based position of the top-strand `C` in
/// each `CG` dinucleotide (case-insensitive), returned in ascending order.
///
/// Shared by the CpG-truth tally ([`crate::output::cpg_truth`]) and the
/// MT/MB classifier ([`crate::vcf::methylation`]); both need the identical
/// scan over an unmodified reference.
#[must_use]
pub(crate) fn find_reference_cpgs(reference: &[u8]) -> Vec<u32> {
    let mut out = Vec::new();
    if reference.len() < 2 {
        return out;
    }
    for i in 0..reference.len() - 1 {
        let c0 = reference[i].to_ascii_uppercase();
        let c1 = reference[i + 1].to_ascii_uppercase();
        if c0 == b'C' && c1 == b'G' {
            #[expect(clippy::cast_possible_truncation, reason = "ref position fits u32")]
            out.push(i as u32);
        }
    }
    out
}

/// Apply per-base methylation chemistry conversion to read bases in place.
///
/// `bases` are in read 5'->3' orientation (FASTQ orientation). `n_genomic`
/// is the number of leading bases that come from the haplotype (anything
/// past `n_genomic` is adapter sequence and is left untouched). When
/// `is_negative_strand` is `true`, `bases[i]` covers haplotype position
/// `hap_start + (n_genomic - 1 - i)`; otherwise `hap_start + i`.
///
/// `haplotype_index` selects which per-haplotype methylation table to use
/// from the [`ContigMethylation`] in `config`.
///
/// Per-molecule conversion failure is drawn once here, before the per-base
/// loop: with probability `config.failure_rate` the molecule is a failure
/// and converts its should-convert cytosines at `1.0 - config.conversion_rate`
/// (near-zero), otherwise at `config.conversion_rate`. Because this runs once
/// per fragment (via [`crate::read::apply_fragment_chemistry`]), both mates
/// derive from the same converted buffer and stay coherent. Returns whether
/// the molecule was drawn as a conversion failure so the caller can record it
/// as ground truth.
///
/// Called between `uppercase_in_place` and `apply_errors` in
/// [`crate::read::generate_read_pair`]'s `build_mate` helper. Chemistry
/// runs before sequencing errors are applied (mirroring the biological
/// order).
///
/// # Panics
///
/// Panics if `config.conversion_rate` or `config.failure_rate` is not a
/// finite value in `[0.0, 1.0]`. Mirrors the guard on
/// [`MethylationTable::from_haplotype`]: the CLI validates this, but the
/// function is `pub`, and an unguarded `NaN` / out-of-range rate would
/// silently distort every chemistry draw.
pub fn apply_methylation_conversion(
    bases: &mut [u8],
    n_genomic: usize,
    is_negative_strand: bool,
    hap_start: u32,
    haplotype_index: usize,
    config: &MethylationConfig<'_>,
    rng: &mut impl Rng,
) -> bool {
    // Validate at the public boundary (see `from_haplotype` for the same
    // guard). An out-of-range / NaN rate would silently corrupt the
    // `rng < rate` decision rather than fail loudly.
    assert!(
        config.conversion_rate.is_finite() && (0.0..=1.0).contains(&config.conversion_rate),
        "conversion_rate must be a finite value in [0.0, 1.0]; got {}",
        config.conversion_rate,
    );
    assert!(
        config.failure_rate.is_finite() && (0.0..=1.0).contains(&config.failure_rate),
        "failure_rate must be a finite value in [0.0, 1.0]; got {}",
        config.failure_rate,
    );
    // Draw the per-molecule camp once, before the per-base loop. A failed
    // molecule converts at `1 - conversion_rate` (near-zero), retaining
    // essentially all of its should-convert cytosines as a coherent unit.
    // Only consume an RNG draw when failures are enabled.
    let conversion_failed = config.failure_rate > 0.0 && rng.random::<f64>() < config.failure_rate;
    let rate =
        if conversion_failed { 1.0 - config.conversion_rate } else { config.conversion_rate };
    let table = config.contig_methylation.table_for(haplotype_index);
    for (i, base) in bases.iter_mut().enumerate().take(n_genomic) {
        let b = *base;
        if b != b'C' {
            continue;
        }
        // `uppercase_in_place` runs before this function in the simulator
        // pipeline, so lowercase 'c' is structurally impossible here.
        debug_assert!(
            !base.is_ascii_lowercase(),
            "uppercase_in_place runs before apply_methylation_conversion"
        );

        let pos_idx = if is_negative_strand { n_genomic - 1 - i } else { i };
        #[expect(clippy::cast_possible_truncation, reason = "pos_idx fits in u32")]
        let hap_pos = hap_start + pos_idx as u32;
        let is_meth = table.is_methylated(hap_pos, is_negative_strand);
        let should_convert = match config.mode {
            MethylationMode::EmSeq => !is_meth,
            MethylationMode::Taps => is_meth,
        };
        if should_convert && rng.random::<f64>() < rate {
            *base = b'T';
        }
    }
    conversion_failed
}

/// Which methylation conversion pattern a read displays when mapped to the
/// reference. Mirrors the Bismark `XG:Z` (genome-strand) tag convention --
/// it's a strand indicator and stays the same for both em-seq and TAPS
/// chemistries; only the biological meaning of an observed `C→T` differs
/// (em-seq: unmethylated at that site; TAPS: methylated at that site).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConversionType {
    /// Read derived from the top (forward) strand; shows `C->T` pattern when
    /// mapped to the reference.
    Ct,
    /// Read derived from the bottom (reverse) strand; shows `G->A` pattern
    /// when mapped to the reference.
    Ga,
}

impl ConversionType {
    /// Two-letter tag value (`"CT"` or `"GA"`) for the `XG:Z` BAM tag.
    #[must_use]
    pub fn as_tag_str(self) -> &'static str {
        match self {
            Self::Ct => "CT",
            Self::Ga => "GA",
        }
    }

    /// Top strand → `Ct`; bottom strand → `Ga`.
    #[must_use]
    pub fn from_strand(is_top: bool) -> Self {
        if is_top { Self::Ct } else { Self::Ga }
    }
}

/// Annotation captured during methylation read generation, used by
/// downstream writers (notably the golden BAM) to emit the full Bismark-
/// compatible methylation tag set: `XG:Z` (genome-strand indicator),
/// `XR:Z` (read-conversion direction; derived from the SAM flag at
/// emission time), `YS:Z` (pre-conversion bases), the holodeck `cf:i`
/// conversion-failure flag, plus the Bismark call tags `XM:Z` / `YM:Z` /
/// `NM:i` / `MD:Z` carried in `r1_call_tags` / `r2_call_tags` and computed
/// by [`crate::methylation_tags::populate_pair_call_tags`].
#[derive(Debug, Clone)]
pub struct MethylationAnnotation {
    /// Conversion direction for the source fragment (same for R1 and R2).
    pub conversion_type: ConversionType,
    /// Whether the source molecule was drawn as a conversion failure. A
    /// molecule property, so it is identical for R1 and R2; surfaced in the
    /// golden BAM as the `cf:i` ground-truth tag.
    pub conversion_failed: bool,
    /// R1 pre-conversion bases. `None` when `capture_pre_conversion` was
    /// false at simulation time.
    pub r1_pre_conversion_bases: Option<Vec<u8>>,
    /// R2 pre-conversion bases. `None` for SE reads or when
    /// `capture_pre_conversion` was false.
    pub r2_pre_conversion_bases: Option<Vec<u8>>,
    /// R1 Bismark-style methylation call tags (`XM`, `YM`, `NM`, `MD`).
    /// Populated by the simulator when the golden BAM is requested.
    pub r1_call_tags: Option<crate::methylation_tags::CallTags>,
    /// R2 Bismark-style methylation call tags. `None` for SE reads or
    /// when the golden BAM is not requested.
    pub r2_call_tags: Option<crate::methylation_tags::CallTags>,
}

impl MethylationAnnotation {
    /// Tags for the R1 BAM record. Returns `None` when the pre-conversion
    /// bases were not captured (e.g., the user requested bisulfite
    /// simulation but not a golden BAM).
    #[must_use]
    pub fn r1_tags(&self) -> Option<MethylationRecordTags<'_>> {
        self.r1_pre_conversion_bases.as_deref().map(|bases| MethylationRecordTags {
            conversion_type: self.conversion_type,
            conversion_failed: self.conversion_failed,
            pre_conversion_bases: bases,
            call_tags: self.r1_call_tags.as_ref(),
        })
    }

    /// Tags for the R2 BAM record. `None` for SE pairs OR when pre-conversion
    /// bases were not captured.
    #[must_use]
    pub fn r2_tags(&self) -> Option<MethylationRecordTags<'_>> {
        self.r2_pre_conversion_bases.as_deref().map(|bases| MethylationRecordTags {
            conversion_type: self.conversion_type,
            conversion_failed: self.conversion_failed,
            pre_conversion_bases: bases,
            call_tags: self.r2_call_tags.as_ref(),
        })
    }
}

/// Per-record methylation annotation passed into the golden BAM record
/// builder when the pair was generated with methylation chemistry enabled
/// (em-seq or TAPS). Carries the `XG:Z` strand indicator (shared across R1
/// and R2) and the read's own pre-conversion bases in read 5'->3' (FASTQ)
/// orientation. The record builder reverse-complements them when the record
/// is reverse-strand so `YS:Z` ends up in the same orientation as `SEQ`.
///
/// `YS:Z` is a holodeck-specific BAM tag — no bisulfite aligner produces or
/// consumes it. It exists so downstream evaluators can diff `SEQ` against
/// `YS` base-for-base to recover ground-truth chemistry events.
///
/// `call_tags`, when present, carries the precomputed `XM:Z`, `YM:Z`,
/// `NM:i`, and `MD:Z` values for the record.
#[derive(Debug, Clone, Copy)]
pub struct MethylationRecordTags<'a> {
    /// Strand indicator for the `XG:Z` tag (same for R1 and R2 of a pair).
    pub conversion_type: ConversionType,
    /// Whether the source molecule was a conversion failure; emitted as the
    /// `cf:i` ground-truth tag. Same for R1 and R2 of a pair.
    pub conversion_failed: bool,
    /// Pre-conversion bases (read 5'->3' orientation) for the `YS:Z` tag.
    pub pre_conversion_bases: &'a [u8],
    /// Precomputed Bismark methylation call tags. `None` when the
    /// simulator did not compute them (no golden BAM requested).
    pub call_tags: Option<&'a crate::methylation_tags::CallTags>,
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand::rngs::SmallRng;

    use super::*;
    use crate::haplotype::build_haplotypes;
    use crate::vcf::genotype::{Genotype, VariantRecord};

    #[test]
    fn test_find_reference_cpgs_basic() {
        // Reference "ACGTACG" → CpGs at positions 1 and 5.
        assert_eq!(find_reference_cpgs(b"ACGTACG"), vec![1, 5]);
    }

    #[test]
    fn test_find_reference_cpgs_case_insensitive() {
        assert_eq!(find_reference_cpgs(b"acgTaCg"), vec![1, 5]);
    }

    #[test]
    fn test_find_reference_cpgs_empty_and_short() {
        assert!(find_reference_cpgs(b"").is_empty());
        assert!(find_reference_cpgs(b"C").is_empty());
        assert!(find_reference_cpgs(b"AT").is_empty());
        assert_eq!(find_reference_cpgs(b"CG"), vec![0]);
    }

    /// Build a [`MethylationConfig`] wrapping a single-haplotype
    /// [`ContigMethylation`] for chemistry-only invocation tests. Returns
    /// the owning `ContigMethylation` so the borrow lives long enough.
    fn single_hap_cm(table: MethylationTable) -> ContigMethylation {
        ContigMethylation::from_tables(vec![table])
    }

    /// Build a single all-reference haplotype (no variants).
    fn ref_haplotype() -> crate::haplotype::Haplotype {
        let haps = build_haplotypes(&[], 1, &mut SmallRng::seed_from_u64(0));
        haps.into_iter().next().unwrap()
    }

    /// Build a SNP variant record helper for tests.
    fn snp_variant(pos: u32, ref_base: u8, alt_base: u8, gt: &str) -> VariantRecord {
        VariantRecord {
            position: pos,
            ref_allele: vec![ref_base],
            alt_alleles: vec![vec![alt_base]],
            genotype: Genotype::parse(gt).unwrap(),
        }
    }

    /// Build an indel variant record helper for tests.
    fn indel_variant(pos: u32, ref_allele: &[u8], alt_allele: &[u8], gt: &str) -> VariantRecord {
        VariantRecord {
            position: pos,
            ref_allele: ref_allele.to_vec(),
            alt_alleles: vec![alt_allele.to_vec()],
            genotype: Genotype::parse(gt).unwrap(),
        }
    }

    // --- from_haplotype tests ---

    #[test]
    fn test_from_haplotype_matches_reference_for_no_variants_haplotype() {
        // No variants → haplotype 0 is the reference. The bitmap shape and
        // the methylation marks should match a direct CpG scan of the ref.
        let reference = b"ACGTACGT";
        let hap = ref_haplotype();
        let mut rng = SmallRng::seed_from_u64(42);
        let table = MethylationTable::from_haplotype(hap_borrow(&hap), reference, 1.0, &mut rng);

        assert_eq!(table.len(), reference.len());
        for i in 0u32..8 {
            let want_top = i == 1 || i == 5;
            let want_bottom = i == 2 || i == 6;
            assert_eq!(table.is_methylated(i, false), want_top, "top[{i}] mismatch");
            assert_eq!(table.is_methylated(i, true), want_bottom, "bottom[{i}] mismatch");
        }
    }

    /// Tiny shim so the no-variants-haplotype helper can be reused without
    /// transferring ownership for each call site.
    fn hap_borrow(h: &crate::haplotype::Haplotype) -> &crate::haplotype::Haplotype {
        h
    }

    #[test]
    fn test_from_haplotype_snp_creates_cpg() {
        // Reference ATG → SNP T→C at pos 1 on the variant haplotype gives
        // ACG, which has a CpG at hap positions (1, 2). With rate 1.0 both
        // strand bits should be set for haplotype 1.
        let reference = b"ATG";
        let variants = vec![snp_variant(1, b'T', b'C', "0|1")];
        let haps = build_haplotypes(&variants, 2, &mut SmallRng::seed_from_u64(7));
        let mut rng = SmallRng::seed_from_u64(42);
        let var_hap_table = MethylationTable::from_haplotype(&haps[1], reference, 1.0, &mut rng);

        assert_eq!(var_hap_table.len(), 3);
        assert!(var_hap_table.is_methylated(1, false), "top-strand C at hap pos 1 must be set");
        assert!(var_hap_table.is_methylated(2, true), "bottom-strand C at hap pos 2 must be set");

        // Reference haplotype: no CpG at all.
        let mut rng2 = SmallRng::seed_from_u64(42);
        let ref_hap_table = MethylationTable::from_haplotype(&haps[0], reference, 1.0, &mut rng2);
        assert!(!ref_hap_table.is_methylated(1, false));
        assert!(!ref_hap_table.is_methylated(2, true));
    }

    #[test]
    fn test_from_haplotype_snp_destroys_cpg() {
        // Reference ACG (CpG at positions 1-2) → SNP C→T at pos 1 on the
        // variant haplotype gives ATG, no CpG. No methylation marks
        // anywhere on hap 1.
        let reference = b"ACG";
        let variants = vec![snp_variant(1, b'C', b'T', "0|1")];
        let haps = build_haplotypes(&variants, 2, &mut SmallRng::seed_from_u64(7));
        let mut rng = SmallRng::seed_from_u64(42);
        let table = MethylationTable::from_haplotype(&haps[1], reference, 1.0, &mut rng);

        assert!(!table.is_methylated(0, false));
        assert!(!table.is_methylated(1, false));
        assert!(!table.is_methylated(2, false));
        assert!(!table.is_methylated(0, true));
        assert!(!table.is_methylated(1, true));
        assert!(!table.is_methylated(2, true));
    }

    #[test]
    fn test_from_haplotype_deletion_bridges_cpg() {
        // Reference CTG: no CpG. Delete the T (CT → C, anchor at pos 0)
        // → haplotype CG with hap positions 0 (C) and 1 (G).
        // Assert top[0] and bottom[1] are both set.
        let reference = b"CTG";
        // VCF deletion: REF=CT (positions 0..2) → ALT=C; net -1 base.
        let variants = vec![indel_variant(0, b"CT", b"C", "0|1")];
        let haps = build_haplotypes(&variants, 2, &mut SmallRng::seed_from_u64(7));
        let mut rng = SmallRng::seed_from_u64(42);
        let table = MethylationTable::from_haplotype(&haps[1], reference, 1.0, &mut rng);

        // Materialized hap is "CG" (length 2).
        assert_eq!(table.len(), 2);
        assert!(table.is_methylated(0, false), "top-strand C at hap pos 0 must be set");
        assert!(table.is_methylated(1, true), "bottom-strand C at hap pos 1 must be set");
    }

    #[test]
    fn test_from_haplotype_inserted_cpg() {
        // Reference AG: no CpG. Insertion of C between A and G:
        //   VCF style: REF=A (pos 0) → ALT=AC. After insertion the
        //   haplotype reads "ACG" with hap positions 0 (A), 1 (C), 2 (G).
        // Assert top[1] and bottom[2] are set on the variant haplotype.
        let reference = b"AG";
        let variants = vec![indel_variant(0, b"A", b"AC", "0|1")];
        let haps = build_haplotypes(&variants, 2, &mut SmallRng::seed_from_u64(7));
        let mut rng = SmallRng::seed_from_u64(42);
        let table = MethylationTable::from_haplotype(&haps[1], reference, 1.0, &mut rng);

        assert_eq!(table.len(), 3, "haplotype should be 3 bases (ACG)");
        assert!(table.is_methylated(1, false), "inserted top-strand C must be methylated");
        assert!(table.is_methylated(2, true), "bottom-strand C at G must be methylated");
    }

    #[test]
    fn test_from_haplotype_multi_base_insertion_not_truncated() {
        // Reference AG (length 2). Insertion: A → ACGCG (alt_len=5, ref_len=1,
        // net delta = +4). Materialized haplotype is "ACGCGG" (length 6) with
        // CpGs at hap positions (1,2) and (3,4). Both pairs must register —
        // a previous over-tight cap truncated tail bases past
        // reference.len() + 1, losing the second CpG entirely.
        let reference = b"AG";
        let variants = vec![indel_variant(0, b"A", b"ACGCG", "0|1")];
        let haps = build_haplotypes(&variants, 2, &mut SmallRng::seed_from_u64(7));
        let mut rng = SmallRng::seed_from_u64(42);
        let table = MethylationTable::from_haplotype(&haps[1], reference, 1.0, &mut rng);

        assert_eq!(table.len(), 6, "haplotype must materialize all 6 bases (ACGCGG)");
        assert!(table.is_methylated(1, false), "top[1] must be set (first CpG)");
        assert!(table.is_methylated(2, true), "bottom[2] must be set (first CpG)");
        assert!(
            table.is_methylated(3, false),
            "top[3] must be set (second CpG, lost when truncated)"
        );
        assert!(
            table.is_methylated(4, true),
            "bottom[4] must be set (second CpG, lost when truncated)"
        );
    }

    #[test]
    #[should_panic(expected = "methylation_rate must be a finite value in [0.0, 1.0]")]
    fn test_from_haplotype_rejects_nan_rate() {
        let reference = b"ACGT";
        let hap = ref_haplotype();
        let mut rng = SmallRng::seed_from_u64(42);
        let _ = MethylationTable::from_haplotype(&hap, reference, f64::NAN, &mut rng);
    }

    #[test]
    #[should_panic(expected = "methylation_rate must be a finite value in [0.0, 1.0]")]
    fn test_from_haplotype_rejects_rate_above_one() {
        let reference = b"ACGT";
        let hap = ref_haplotype();
        let mut rng = SmallRng::seed_from_u64(42);
        let _ = MethylationTable::from_haplotype(&hap, reference, 1.5, &mut rng);
    }

    #[test]
    fn test_from_haplotype_zero_rate_no_methylation() {
        let reference = b"ACGTACGTACGT";
        let hap = ref_haplotype();
        let mut rng = SmallRng::seed_from_u64(42);
        let table = MethylationTable::from_haplotype(&hap, reference, 0.0, &mut rng);
        #[expect(clippy::cast_possible_truncation, reason = "test reference len fits in u32")]
        let len = reference.len() as u32;
        for i in 0..len {
            assert!(!table.is_methylated(i, false));
            assert!(!table.is_methylated(i, true));
        }
    }

    #[test]
    fn test_from_haplotype_case_insensitive() {
        // Lowercase "acgt" should still detect a CpG at (1, 2).
        let reference = b"acgt";
        let hap = ref_haplotype();
        let mut rng = SmallRng::seed_from_u64(42);
        let table = MethylationTable::from_haplotype(&hap, reference, 1.0, &mut rng);
        assert!(table.is_methylated(1, false), "lowercase 'cg' must register top-strand C");
        assert!(table.is_methylated(2, true), "lowercase 'cg' must register bottom-strand C");
        assert!(!table.is_methylated(0, false));
        assert!(!table.is_methylated(3, true));
    }

    #[test]
    fn test_from_haplotype_independent_strand_draws() {
        // 1000 CpG sites at rate 0.5 -- both strands should each show
        // ~50% methylation, drawn independently.
        //
        // Band derived empirically with `SmallRng::seed_from_u64(7)` on rand
        // 0.9. If `rand` updates `SmallRng`'s output stream, this band may
        // need widening or re-derivation.
        let mut reference = Vec::with_capacity(1000 * 4);
        for _ in 0..1000 {
            reference.extend_from_slice(b"ACGT");
        }
        let hap = ref_haplotype();
        let mut rng = SmallRng::seed_from_u64(7);
        let table = MethylationTable::from_haplotype(&hap, &reference, 0.5, &mut rng);

        let mut top_meth = 0usize;
        let mut bottom_meth = 0usize;
        let mut both_meth = 0usize;
        let mut top_only = 0usize;
        let mut bottom_only = 0usize;
        for site in 0u32..1000 {
            let top_pos = site * 4 + 1; // C of CpG
            let bottom_pos = site * 4 + 2; // G of CpG (bottom-strand C)
            let t = table.is_methylated(top_pos, false);
            let b = table.is_methylated(bottom_pos, true);
            if t {
                top_meth += 1;
            }
            if b {
                bottom_meth += 1;
            }
            if t && b {
                both_meth += 1;
            }
            if t && !b {
                top_only += 1;
            }
            if !t && b {
                bottom_only += 1;
            }
        }

        // Each strand independently ~50%.
        assert!((400..=600).contains(&top_meth), "top methylation count out of band: {top_meth}");
        assert!(
            (400..=600).contains(&bottom_meth),
            "bottom methylation count out of band: {bottom_meth}"
        );
        // If the two strands were perfectly correlated we'd see ~500 doubles;
        // independent at p=0.5 each gives ~250.
        assert!(
            (180..=320).contains(&both_meth),
            "both-strand methylation count {both_meth} suggests strands are not independent"
        );
        // Independence requires plenty of sites where one strand is methylated
        // and the other is not -- positive evidence the draws aren't tied.
        assert!(
            top_only > 100 && bottom_only > 100,
            "expected hemimethylated sites in both directions; got top_only={top_only} bottom_only={bottom_only}"
        );
    }

    #[test]
    fn test_from_haplotype_empty_and_short() {
        let hap = ref_haplotype();
        let mut rng = SmallRng::seed_from_u64(42);
        let table = MethylationTable::from_haplotype(&hap, b"", 1.0, &mut rng);
        assert!(table.is_empty());
        assert_eq!(table.len(), 0);

        let table = MethylationTable::from_haplotype(&hap, b"C", 1.0, &mut rng);
        assert_eq!(table.len(), 1);
        assert!(!table.is_methylated(0, false));
        assert!(!table.is_methylated(0, true));
    }

    #[test]
    fn test_is_methylated_strand_selection() {
        let mut table = MethylationTable::empty(10);
        table.set_top(3, true);
        table.set_bottom(7, true);
        assert!(table.is_methylated(3, false));
        assert!(!table.is_methylated(3, true));
        assert!(!table.is_methylated(7, false));
        assert!(table.is_methylated(7, true));
    }

    #[test]
    fn test_is_methylated_out_of_range() {
        let table = MethylationTable::empty(10);
        assert!(!table.is_methylated(99, false));
        assert!(!table.is_methylated(99, true));
        assert!(!table.is_methylated(u32::MAX, false));
    }

    #[test]
    fn test_empty_table_returns_false_everywhere() {
        let table = MethylationTable::empty(100);
        assert_eq!(table.len(), 100);
        assert!(!table.is_empty());
        for i in 0..100 {
            assert!(!table.is_methylated(i, false));
            assert!(!table.is_methylated(i, true));
        }
    }

    // --- ContigMethylation tests ---

    #[test]
    fn test_contig_methylation_per_haplotype_independent() {
        // Two haplotypes with different CpG content. ContigMethylation
        // should keep their bitmaps separate.
        let reference = b"ACG"; // hap 0 = ref ACG, hap 1 = AAG (SNP C→A at pos 1)
        let variants = vec![snp_variant(1, b'C', b'A', "0|1")];
        let haps = build_haplotypes(&variants, 2, &mut SmallRng::seed_from_u64(0));
        let mut rng = SmallRng::seed_from_u64(42);
        let cm = ContigMethylation::from_haplotypes(&haps, reference, 1.0, &mut rng);

        assert_eq!(cm.len(), 2);
        assert!(!cm.is_empty());

        // Haplotype 0 (ref): CpG at hap_pos (1, 2).
        let t0 = cm.table_for(0);
        assert!(t0.is_methylated(1, false));
        assert!(t0.is_methylated(2, true));

        // Haplotype 1 (AAG): no CpG.
        let t1 = cm.table_for(1);
        assert!(!t1.is_methylated(1, false));
        assert!(!t1.is_methylated(2, true));
    }

    // --- apply_methylation_conversion tests (em-seq mode) ---

    #[test]
    fn test_apply_methylation_conversion_em_seq_zero_meth_full_conversion_forward() {
        let cm = single_hap_cm(MethylationTable::empty(20));
        let mut bases = b"ACGTACGT".to_vec();
        let mut rng = SmallRng::seed_from_u64(42);

        let c = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::EmSeq,
            conversion_rate: 1.0,
            failure_rate: 0.0,
        };
        apply_methylation_conversion(&mut bases, 8, false, 10, 0, &c, &mut rng);

        // 0% methylated, 100% conversion rate -> every C becomes T.
        assert_eq!(&bases, b"ATGTATGT");
    }

    #[test]
    fn test_apply_methylation_conversion_em_seq_full_meth_no_conversion() {
        // C's in "ACGTACGT" at read positions 1 and 5 cover hap positions
        // 11 and 15 (hap_start = 10). Set the top-strand bits directly so
        // they're protected.
        let mut table = MethylationTable::empty(20);
        table.set_top(11, true);
        table.set_top(15, true);
        let cm = single_hap_cm(table);

        let mut bases = b"ACGTACGT".to_vec();
        let mut rng = SmallRng::seed_from_u64(42);

        let c = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::EmSeq,
            conversion_rate: 1.0,
            failure_rate: 0.0,
        };
        apply_methylation_conversion(&mut bases, 8, false, 10, 0, &c, &mut rng);

        // Both C's are methylated -> no conversion, even at full rate.
        assert_eq!(&bases, b"ACGTACGT");
    }

    #[test]
    fn test_apply_methylation_conversion_em_seq_negative_strand_uses_reversed_index() {
        // Methylate ONLY hap position 17 on the bottom strand.
        let mut table = MethylationTable::empty(20);
        table.set_bottom(17, true);
        let cm = single_hap_cm(table);

        // For a negative-strand read, bases[0] covers hap_start + (n-1).
        // hap_start = 13, n = 5 → bases[0] covers hap_pos 17.
        // The C at read position 0 should look up hap pos 17 and find
        // 100% methylation (on the bottom strand), so it must NOT convert.
        // The C at read position 4 (covering hap pos 13) is not methylated
        // → converts.
        let mut bases = b"CAAAC".to_vec();
        let mut rng = SmallRng::seed_from_u64(42);

        let c = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::EmSeq,
            conversion_rate: 1.0,
            failure_rate: 0.0,
        };
        apply_methylation_conversion(&mut bases, 5, true, 13, 0, &c, &mut rng);

        assert_eq!(&bases, b"CAAAT");
    }

    #[test]
    fn test_apply_methylation_conversion_skips_adapter_bases() {
        let cm = single_hap_cm(MethylationTable::empty(10));
        // 3 genomic bases + 5 adapter bases. The adapter contains C's that
        // must NOT be touched because they have no haplotype position.
        let mut bases = b"ACGCCNNN".to_vec();
        let mut rng = SmallRng::seed_from_u64(42);

        let c = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::EmSeq,
            conversion_rate: 1.0,
            failure_rate: 0.0,
        };
        apply_methylation_conversion(&mut bases, 3, false, 0, 0, &c, &mut rng);

        // bases[0..3] = ACG -> ATG (one C converted)
        // bases[3..] untouched -> still CCNNN
        assert_eq!(&bases, b"ATGCCNNN");
    }

    // --- per-molecule conversion-failure tests ---

    #[test]
    fn test_failed_molecule_retains_all_should_convert_cytosines() {
        // failure_rate = 1.0 forces the failed camp; with conversion_rate = 1.0
        // the failed camp converts at 1 - 1.0 = 0.0, so every should-convert C
        // is retained and the returned flag reports the failure.
        let cm = single_hap_cm(MethylationTable::empty(20));
        let mut bases = b"ACGTACGT".to_vec();
        let mut rng = SmallRng::seed_from_u64(42);

        let c = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::EmSeq,
            conversion_rate: 1.0,
            failure_rate: 1.0,
        };
        let failed = apply_methylation_conversion(&mut bases, 8, false, 10, 0, &c, &mut rng);

        assert!(failed, "molecule should be flagged as a conversion failure");
        assert_eq!(&bases, b"ACGTACGT", "failed molecule must retain every should-convert C");
    }

    #[test]
    fn test_failed_molecule_converts_at_one_minus_conversion_rate() {
        // conversion_rate = 0.0 means the failed camp converts at 1 - 0.0 = 1.0,
        // so a forced-failed molecule converts every should-convert C. Pins the
        // "failed rate = 1 - conversion_rate" relationship.
        let cm = single_hap_cm(MethylationTable::empty(20));
        let mut bases = b"ACGTACGT".to_vec();
        let mut rng = SmallRng::seed_from_u64(42);

        let c = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::EmSeq,
            conversion_rate: 0.0,
            failure_rate: 1.0,
        };
        let failed = apply_methylation_conversion(&mut bases, 8, false, 10, 0, &c, &mut rng);

        assert!(failed);
        assert_eq!(&bases, b"ATGTATGT", "failed camp at 1 - 0.0 = 1.0 converts every C");
    }

    #[test]
    fn test_failure_rate_zero_never_flags_failure() {
        let cm = single_hap_cm(MethylationTable::empty(20));
        let mut bases = b"ACGTACGT".to_vec();
        let mut rng = SmallRng::seed_from_u64(42);

        let c = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::EmSeq,
            conversion_rate: 1.0,
            failure_rate: 0.0,
        };
        let failed = apply_methylation_conversion(&mut bases, 8, false, 10, 0, &c, &mut rng);

        assert!(!failed, "failure_rate 0.0 must never flag a failure");
        assert_eq!(&bases, b"ATGTATGT", "non-failed molecule at rate 1.0 converts every C");
    }

    #[test]
    fn test_zero_genomic_bases_is_noop_but_still_draws_failure() {
        // A fully-adapter read (n_genomic = 0) must touch no bases. The
        // negative-strand index math `n_genomic - 1 - i` would underflow if
        // the per-base loop ran, so this guards that `.take(0)` keeps it out
        // of the loop. The per-molecule failure draw still happens (it is
        // independent of base count), so the flag is still reported.
        let cm = single_hap_cm(MethylationTable::empty(20));
        let mut bases = b"CCCCCCCC".to_vec();
        let mut rng = SmallRng::seed_from_u64(42);

        let c = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::EmSeq,
            conversion_rate: 1.0,
            failure_rate: 1.0,
        };
        let failed = apply_methylation_conversion(&mut bases, 0, true, 10, 0, &c, &mut rng);

        assert!(failed, "failure draw is independent of genomic base count");
        assert_eq!(&bases, b"CCCCCCCC", "zero genomic bases must leave the buffer untouched");
    }

    #[test]
    fn test_failure_rate_observed_fraction_matches() {
        // Each call models one molecule; ~50% should be flagged failed at
        // failure_rate = 0.5. Band derived empirically with
        // `SmallRng::seed_from_u64(42)` on rand 0.9; the [0.47, 0.53] window
        // is ~4 sigma from 0.5 at n = 5000, so it tolerates RNG-stream churn,
        // but may need re-derivation if `SmallRng`'s output stream changes.
        let cm = single_hap_cm(MethylationTable::empty(8));
        let mut rng = SmallRng::seed_from_u64(42);
        let c = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::EmSeq,
            conversion_rate: 0.999,
            failure_rate: 0.5,
        };

        let n = 5000;
        let mut failures = 0;
        for _ in 0..n {
            let mut bases = b"ACGTACGT".to_vec();
            if apply_methylation_conversion(&mut bases, 8, false, 0, 0, &c, &mut rng) {
                failures += 1;
            }
        }
        let frac = f64::from(failures) / f64::from(n);
        assert!((0.47..=0.53).contains(&frac), "observed failed fraction {frac} out of band");
    }

    #[test]
    fn test_apply_methylation_conversion_em_seq_partial_conversion_rate_empirical() {
        let cm = single_hap_cm(MethylationTable::empty(10_000));

        // 10_000 C's, conversion rate 0.5, no methylation -> expect ~50% T's.
        //
        // Band derived empirically with `SmallRng::seed_from_u64(42)` on
        // rand 0.9. If `rand` updates `SmallRng`'s output stream, this band
        // may need widening or re-derivation.
        let mut bases = vec![b'C'; 10_000];
        let mut rng = SmallRng::seed_from_u64(42);

        let c = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::EmSeq,
            conversion_rate: 0.5,
            failure_rate: 0.0,
        };
        apply_methylation_conversion(&mut bases, 10_000, false, 0, 0, &c, &mut rng);

        #[expect(clippy::naive_bytecount, reason = "test, no bytecount dep")]
        let t_count = bases.iter().filter(|&&b| b == b'T').count();
        let frac = t_count as f64 / 10_000.0;
        assert!((0.48..0.52).contains(&frac), "expected ~50% conversion, got {frac:.3}");
    }

    #[test]
    fn test_apply_methylation_conversion_em_seq_zero_conversion_rate_no_change() {
        let cm = single_hap_cm(MethylationTable::empty(10));
        let mut bases = b"CCCCCCCC".to_vec();
        let mut rng = SmallRng::seed_from_u64(42);

        let c = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::EmSeq,
            conversion_rate: 0.0,
            failure_rate: 0.0,
        };
        apply_methylation_conversion(&mut bases, 8, false, 0, 0, &c, &mut rng);

        assert_eq!(&bases, b"CCCCCCCC");
    }

    #[test]
    fn test_apply_methylation_conversion_only_converts_c_not_other_bases() {
        let cm = single_hap_cm(MethylationTable::empty(10));
        let mut bases = b"AGTNAGTN".to_vec();
        let mut rng = SmallRng::seed_from_u64(42);

        let c = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::EmSeq,
            conversion_rate: 1.0,
            failure_rate: 0.0,
        };
        apply_methylation_conversion(&mut bases, 8, false, 0, 0, &c, &mut rng);

        assert_eq!(&bases, b"AGTNAGTN");
    }

    #[test]
    fn test_apply_methylation_conversion_with_hap_start_offset() {
        // Verify the hap_start offset is applied correctly. Bitmap of length
        // 100 with top[50] set; apply conversion to a 10-base read at
        // hap_start = 45 with a C at read position 5 → maps to hap_pos 50,
        // which is methylated, so it must be preserved at em-seq mode.
        let mut table = MethylationTable::empty(100);
        table.set_top(50, true);
        let cm = single_hap_cm(table);

        let mut bases = b"AAAAACAAAA".to_vec();
        let mut rng = SmallRng::seed_from_u64(42);

        let c = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::EmSeq,
            conversion_rate: 1.0,
            failure_rate: 0.0,
        };
        apply_methylation_conversion(&mut bases, 10, false, 45, 0, &c, &mut rng);

        // The C at read pos 5 → hap pos 50 is methylated → preserved.
        assert_eq!(&bases, b"AAAAACAAAA");
    }

    // --- apply_methylation_conversion tests (TAPS mode) ---

    #[test]
    fn test_apply_methylation_conversion_taps_methylated_converts() {
        // TAPS: methylated cytosines convert. Set top[15]=true so the C at
        // read position 5 (covering hap pos 15) is the only one that
        // converts.
        let mut table = MethylationTable::empty(20);
        table.set_top(15, true);
        let cm = single_hap_cm(table);

        let mut bases = b"ACGTACGT".to_vec();
        let mut rng = SmallRng::seed_from_u64(42);

        let c = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::Taps,
            conversion_rate: 1.0,
            failure_rate: 0.0,
        };
        apply_methylation_conversion(&mut bases, 8, false, 10, 0, &c, &mut rng);

        // Read pos 1 -> hap 11: not methylated, preserved.
        // Read pos 5 -> hap 15: methylated under TAPS -> converts.
        assert_eq!(&bases, b"ACGTATGT");
    }

    #[test]
    fn test_apply_methylation_conversion_taps_unmethylated_preserved() {
        // TAPS with no methylation: nothing converts.
        let cm = single_hap_cm(MethylationTable::empty(20));
        let mut bases = b"ACGTACGT".to_vec();
        let mut rng = SmallRng::seed_from_u64(42);

        let c = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::Taps,
            conversion_rate: 1.0,
            failure_rate: 0.0,
        };
        apply_methylation_conversion(&mut bases, 8, false, 10, 0, &c, &mut rng);

        assert_eq!(&bases, b"ACGTACGT");
    }

    #[test]
    fn test_apply_methylation_conversion_taps_negative_strand() {
        // Methylate hap pos 17 on the bottom strand. With negative-strand
        // indexing and hap_start=13, n=5 → bases[0] covers hap pos 17.
        // TAPS mode: methylated -> converts. The C at read pos 4 covers
        // hap pos 13 (unmethylated under TAPS -> preserved).
        let mut table = MethylationTable::empty(20);
        table.set_bottom(17, true);
        let cm = single_hap_cm(table);

        let mut bases = b"CAAAC".to_vec();
        let mut rng = SmallRng::seed_from_u64(42);

        let c = MethylationConfig {
            contig_methylation: &cm,
            mode: MethylationMode::Taps,
            conversion_rate: 1.0,
            failure_rate: 0.0,
        };
        apply_methylation_conversion(&mut bases, 5, true, 13, 0, &c, &mut rng);

        // bases[0] (covers hap 17, methylated bottom) -> T
        // bases[4] (covers hap 13, unmethylated) -> stays C
        assert_eq!(&bases, b"TAAAC");
    }
}
