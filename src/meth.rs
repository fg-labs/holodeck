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

/// Default target methylation fraction for CpG-island-interior CpGs.
/// Islands are characteristically hypomethylated. See [`MethylationModel`].
pub const DEFAULT_ISLAND_RATE: f64 = 0.1;

/// Default target methylation fraction for CpG-island-shore CpGs
/// (intermediate). See [`MethylationModel`].
pub const DEFAULT_SHORE_RATE: f64 = 0.5;

/// Default target methylation fraction for open-sea CpGs (hypermethylated;
/// the bulk genomic default). See [`MethylationModel`].
pub const DEFAULT_OPEN_SEA_RATE: f64 = 0.85;

/// Default spatial correlation length (bp) for every context. Sets how far
/// methylation state persists between consecutive CpGs. See
/// [`ContextParams::correlation_length_bp`].
pub const DEFAULT_CORRELATION_LENGTH_BP: f64 = 1000.0;

/// Default sporadic hemimethylation probability. See
/// [`MethylationModel::hemi_rate`].
pub const DEFAULT_HEMI_RATE: f64 = 0.01;

// CpG-island detector thresholds (Gardiner-Garden & Frommer, 1987). Internal
// constants, not CLI flags: these are the canonical island-calling criteria,
// not something a simulation user normally retunes.

/// Minimum window length (bp) over which the island criteria are evaluated.
const ISLAND_MIN_WINDOW_BP: usize = 200;

/// Minimum GC fraction for a window to qualify as a CpG island.
const ISLAND_MIN_GC: f64 = 0.5;

/// Minimum observed/expected CpG ratio for a window to qualify as an island.
const ISLAND_MIN_OE_RATIO: f64 = 0.6;

/// Distance (bp) from an island within which a CpG is classified as "shore".
const SHORE_WIDTH_BP: u32 = 2000;

/// Genomic-context class of a CpG, which selects its methylation parameters.
///
/// `Island` / `Shore` / `OpenSea` are the standard CpG-island taxonomy
/// (Gardiner-Garden 1987; shores from Irizarry 2009). The canonical scheme
/// also has a "shelf" tier 2-4 kb out; this model folds shelf into open-sea.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CpgContext {
    /// Inside a detected CpG island — hypomethylated.
    Island,
    /// Within [`SHORE_WIDTH_BP`] of an island — intermediate.
    Shore,
    /// Everywhere else — hypermethylated.
    OpenSea,
}

/// Per-context methylation parameters: the stationary target rate and the
/// spatial correlation length.
#[derive(Debug, Clone, Copy)]
pub struct ContextParams {
    /// Stationary target methylation fraction in `[0.0, 1.0]`. Over any large
    /// region of this context the mean methylation converges to this value.
    pub rate: f64,
    /// Spatial correlation length L (bp). For two consecutive CpGs separated
    /// by `d` bp the second keeps the first's methylation state with
    /// probability `exp(-d / L)`, otherwise it is redrawn from
    /// `Bernoulli(rate)`. Larger L → longer runs of like-methylated CpGs.
    /// Must be finite and `> 0`.
    pub correlation_length_bp: f64,
}

/// Resolved per-context methylation model used by the `methylate` generator.
///
/// Built once from CLI flags and threaded read-only through
/// [`ContigMethylation::from_haplotypes`] →
/// [`MethylationTable::from_haplotype`]. Each CpG is classified into a
/// [`CpgContext`] from the haplotype sequence, then a two-state
/// (methylated/unmethylated) Markov chain walks the CpG list using that
/// context's [`ContextParams`]: the chain's stationary mean equals the
/// context rate, and spatial autocorrelation decays with genomic distance per
/// the correlation length. Methylation is symmetric (both strands) by
/// default; [`Self::hemi_rate`] introduces sporadic per-CpG hemimethylation.
#[derive(Debug, Clone, Copy)]
pub struct MethylationModel {
    /// Parameters for CpG-island-interior CpGs.
    pub island: ContextParams,
    /// Parameters for island-shore CpGs.
    pub shore: ContextParams,
    /// Parameters for open-sea CpGs.
    pub open_sea: ContextParams,
    /// Probability that a methylated CpG is made hemimethylated — exactly one
    /// randomly chosen strand is left unmethylated. In `[0.0, 1.0]`.
    pub hemi_rate: f64,
}

impl MethylationModel {
    /// Validate every context rate and the hemi rate are finite in
    /// `[0.0, 1.0]` and every correlation length is finite and `> 0`.
    ///
    /// Called at the CLI boundary; [`MethylationTable::from_haplotype`] also
    /// asserts a valid model as defense-in-depth on the `pub` boundary.
    ///
    /// # Errors
    ///
    /// Returns an error naming the offending flag if any bound is violated.
    pub fn validate(&self) -> anyhow::Result<()> {
        for (name, p) in
            [("island", &self.island), ("shore", &self.shore), ("open-sea", &self.open_sea)]
        {
            if !p.rate.is_finite() || !(0.0..=1.0).contains(&p.rate) {
                anyhow::bail!("--methylation-rate-{name} must be in [0.0, 1.0]");
            }
            if !p.correlation_length_bp.is_finite() || p.correlation_length_bp <= 0.0 {
                anyhow::bail!("--methylation-correlation-length-{name} must be a finite value > 0");
            }
        }
        if !self.hemi_rate.is_finite() || !(0.0..=1.0).contains(&self.hemi_rate) {
            anyhow::bail!("--hemimethylation-rate must be in [0.0, 1.0]");
        }
        Ok(())
    }

    /// The [`ContextParams`] for a given CpG context.
    fn params_for(&self, context: CpgContext) -> &ContextParams {
        match context {
            CpgContext::Island => &self.island,
            CpgContext::Shore => &self.shore,
            CpgContext::OpenSea => &self.open_sea,
        }
    }
}

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
    /// over the entire contig length) and assigning methylation with a
    /// context-aware Markov chain. The resulting bitmap length equals the
    /// materialized haplotype length, which may differ from the reference
    /// length when the haplotype carries indels.
    ///
    /// Each CpG is classified ([`classify_cpg_contexts`]) into island / shore
    /// / open-sea, then a two-state (methylated/unmethylated) chain walks the
    /// CpG list in order. For each CpG the [`ContextParams`] of its context
    /// give the stationary target `m` and correlation length `L`: the first
    /// CpG is drawn `Bernoulli(m)`; thereafter, with probability `exp(-d/L)`
    /// (where `d` is the bp distance to the previous CpG) the previous state
    /// is kept, otherwise it is redrawn `Bernoulli(m)`. This yields a
    /// stationary mean of `m` and autocorrelation that decays with genomic
    /// distance. Across a context boundary the marginal relaxes toward the
    /// new target over ~`L / spacing` CpGs — the realistic shore gradient —
    /// so per-context means hold in region interiors, not in transition bands.
    ///
    /// Methylation is **symmetric** by default (a methylated CpG sets both the
    /// top-strand C and the bottom-strand C); with probability
    /// [`MethylationModel::hemi_rate`] a methylated CpG is made hemimethylated
    /// by unsetting exactly one randomly chosen strand. The stationary mean `m`
    /// is the per-CpG *methylated-state* rate; when `hemi_rate > 0` the realized
    /// per-strand methylated fraction is `m * (1 - hemi_rate / 2)`.
    ///
    /// # Panics
    ///
    /// Panics if `model` fails [`MethylationModel::validate`] (NaN/out-of-range
    /// rate or non-positive correlation length). The `methylate` CLI validates
    /// the model at startup, but this is a `pub` constructor reachable from
    /// library code, so the invariant is asserted here as defense-in-depth.
    pub fn from_haplotype(
        haplotype: &crate::haplotype::Haplotype,
        reference: &[u8],
        model: &MethylationModel,
        rng: &mut impl Rng,
    ) -> Self {
        assert!(model.validate().is_ok(), "invalid MethylationModel: {model:?}");
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

        // CpG list (top-strand C positions) on this haplotype's sequence, then
        // per-CpG context from the same materialized bases so variant-created
        // / -destroyed CpGs and indel shifts are handled in haplotype coords.
        let cpg_top_c = find_reference_cpgs(&hap_bases);
        if cpg_top_c.is_empty() {
            return table;
        }
        let island = island_mask(&hap_bases);
        let contexts = classify_cpg_contexts(&cpg_top_c, &island);

        let mut prev_state: Option<bool> = None;
        let mut prev_pos: u32 = 0;
        for (k, &c) in cpg_top_c.iter().enumerate() {
            let params = model.params_for(contexts[k]);
            let state = match prev_state {
                None => rng.random::<f64>() < params.rate,
                Some(prev) => {
                    let d = f64::from(c - prev_pos);
                    let keep_prob = (-d / params.correlation_length_bp).exp();
                    if rng.random::<f64>() < keep_prob {
                        prev
                    } else {
                        rng.random::<f64>() < params.rate
                    }
                }
            };
            if state {
                let c = c as usize;
                table.top.set(c, true);
                table.bottom.set(c + 1, true);
                // Sporadic hemimethylation: drop exactly one strand.
                if rng.random::<f64>() < model.hemi_rate {
                    if rng.random::<bool>() {
                        table.top.set(c, false);
                    } else {
                        table.bottom.set(c + 1, false);
                    }
                }
            }
            prev_state = Some(state);
            prev_pos = c;
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
    /// Build per-haplotype methylation tables with the context-aware Markov
    /// model (see [`MethylationTable::from_haplotype`]). Each haplotype is
    /// walked independently with the shared `rng`, so allele-specific
    /// methylation arises naturally and the per-contig draw order is
    /// deterministic for a fixed seed.
    pub fn from_haplotypes(
        haplotypes: &[crate::haplotype::Haplotype],
        reference: &[u8],
        model: &MethylationModel,
        rng: &mut impl Rng,
    ) -> Self {
        let per_haplotype = haplotypes
            .iter()
            .map(|hap| MethylationTable::from_haplotype(hap, reference, model, rng))
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

/// Per-base CpG-island mask: `mask[p]` is `true` iff position `p` lies in a
/// window that satisfies the Gardiner-Garden island criteria
/// ([`ISLAND_MIN_WINDOW_BP`]-bp window with GC fraction > [`ISLAND_MIN_GC`]
/// and observed/expected CpG ratio > [`ISLAND_MIN_OE_RATIO`]).
///
/// `O(len)`: a single rolling window maintains C, G, and CpG counts, and
/// qualifying windows are painted into the mask with a monotone cursor so each
/// base is written at most once. Case-insensitive; non-`ACGT` bases count as
/// neither GC nor CpG. Sequences shorter than the window yield an all-`false`
/// mask (no island can be called).
#[expect(
    clippy::similar_names,
    reason = "is_c/is_g and n_c/n_g/n_cg mirror the C/G/CpG quantities they track"
)]
fn island_mask(seq: &[u8]) -> BitVec {
    let len = seq.len();
    let window = ISLAND_MIN_WINDOW_BP;
    let mut mask = BitVec::repeat(false, len);
    if len < window {
        return mask;
    }

    let is_c = |j: usize| seq[j].eq_ignore_ascii_case(&b'C');
    let is_g = |j: usize| seq[j].eq_ignore_ascii_case(&b'G');
    // A CpG occupies `j` and `j + 1`; both must exist.
    let is_cg = |j: usize| j + 1 < len && is_c(j) && is_g(j + 1);

    // Counts for the window starting at `a == 0`, covering [0, window).
    let mut n_c = (0..window).filter(|&j| is_c(j)).count();
    let mut n_g = (0..window).filter(|&j| is_g(j)).count();
    // CpG positions fully inside [a, a + window) are `j` in [a, a + window - 1).
    let mut n_cg = (0..window - 1).filter(|&j| is_cg(j)).count();

    let window_f = window as f64;
    let mut painted_end = 0usize;
    for a in 0..=(len - window) {
        let gc_ok = (n_c + n_g) as f64 > ISLAND_MIN_GC * window_f;
        // observed/expected = n_cg / (n_c * n_g / window) > threshold, written
        // without division so a zero C or G count simply fails the test.
        let oe_ok = n_c > 0
            && n_g > 0
            && (n_cg as f64) * window_f > ISLAND_MIN_OE_RATIO * (n_c as f64) * (n_g as f64);
        if gc_ok && oe_ok {
            let start = painted_end.max(a);
            for p in start..(a + window) {
                mask.set(p, true);
            }
            painted_end = a + window;
        }
        // Slide to the window starting at `a + 1`, covering [a+1, a+window+1).
        if a < len - window {
            let s = a + 1;
            n_c = n_c + usize::from(is_c(s + window - 1)) - usize::from(is_c(s - 1));
            n_g = n_g + usize::from(is_g(s + window - 1)) - usize::from(is_g(s - 1));
            // CpG range shifts from [a, a+window-1) to [s, s+window-1): drop
            // the pair leaving at `s-1`, add the pair entering at `s+window-2`.
            n_cg = n_cg + usize::from(is_cg(s + window - 2)) - usize::from(is_cg(s - 1));
        }
    }
    mask
}

/// Contiguous `[start, end)` runs of `true` bits in an island mask, ascending.
fn island_runs(mask: &BitVec) -> Vec<(u32, u32)> {
    let mut runs = Vec::new();
    let mut start: Option<usize> = None;
    for i in 0..mask.len() {
        if mask[i] {
            start.get_or_insert(i);
        } else if let Some(s) = start.take() {
            #[expect(clippy::cast_possible_truncation, reason = "positions fit u32")]
            runs.push((s as u32, i as u32));
        }
    }
    if let Some(s) = start {
        #[expect(clippy::cast_possible_truncation, reason = "positions fit u32")]
        runs.push((s as u32, mask.len() as u32));
    }
    runs
}

/// Whether `c` is within [`SHORE_WIDTH_BP`] (inclusive) of any island run, but
/// not inside one — callers test the mask for `Island` first. `runs` are
/// ascending and disjoint, so only the run immediately left and right of `c`
/// can be the nearest.
fn near_island(c: u32, runs: &[(u32, u32)]) -> bool {
    // First run whose start is strictly greater than `c` (the right neighbor).
    let idx = runs.partition_point(|&(s, _)| s <= c);
    if let Some(&(rs, _)) = runs.get(idx)
        && rs - c <= SHORE_WIDTH_BP
    {
        return true;
    }
    if idx > 0 {
        let (_, re) = runs[idx - 1];
        // `re` is exclusive; the last island base is `re - 1`. `c >= re` here
        // (otherwise `c` would be inside the run → classified Island already).
        if c >= re && c - (re - 1) <= SHORE_WIDTH_BP {
            return true;
        }
    }
    false
}

/// Classify each CpG (given its top-strand C position) as island / shore /
/// open-sea using the island `mask` produced by [`island_mask`].
fn classify_cpg_contexts(cpg_top_c: &[u32], mask: &BitVec) -> Vec<CpgContext> {
    let runs = island_runs(mask);
    cpg_top_c
        .iter()
        .map(|&c| {
            if mask.get(c as usize).is_some_and(|b| *b) {
                CpgContext::Island
            } else if near_island(c, &runs) {
                CpgContext::Shore
            } else {
                CpgContext::OpenSea
            }
        })
        .collect()
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

    /// A [`MethylationModel`] with all three contexts sharing `rate`,
    /// correlation length `l`, and hemimethylation `hemi`. Lets tests isolate
    /// the stationary-mean, autocorrelation, and hemi behaviours independently
    /// of context classification.
    fn model(rate: f64, l: f64, hemi: f64) -> MethylationModel {
        let ctx = ContextParams { rate, correlation_length_bp: l };
        MethylationModel { island: ctx, shore: ctx, open_sea: ctx, hemi_rate: hemi }
    }

    /// A [`MethylationModel`] with all three contexts at the same `rate`,
    /// default correlation length, and no hemimethylation. At `rate` 1.0/0.0
    /// the walk is fully deterministic (every / no CpG methylated, symmetric),
    /// so tests can pin exact bits regardless of context or correlation.
    fn uniform_model(rate: f64) -> MethylationModel {
        model(rate, DEFAULT_CORRELATION_LENGTH_BP, 0.0)
    }

    /// Build a reference of `units` copies of "ACGT" — one CpG every 4 bp
    /// (top-strand C at positions 1, 5, 9, …). Used by the statistical walk
    /// tests; "ACGT" is 50% GC so these CpGs classify as open-sea.
    fn acgt_repeat(units: usize) -> Vec<u8> {
        b"ACGT".repeat(units)
    }

    // --- from_haplotype tests ---

    #[test]
    fn test_from_haplotype_matches_reference_for_no_variants_haplotype() {
        // No variants → haplotype 0 is the reference. The bitmap shape and
        // the methylation marks should match a direct CpG scan of the ref.
        let reference = b"ACGTACGT";
        let hap = ref_haplotype();
        let mut rng = SmallRng::seed_from_u64(42);
        let table = MethylationTable::from_haplotype(
            hap_borrow(&hap),
            reference,
            &uniform_model(1.0),
            &mut rng,
        );

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
        let var_hap_table =
            MethylationTable::from_haplotype(&haps[1], reference, &uniform_model(1.0), &mut rng);

        assert_eq!(var_hap_table.len(), 3);
        assert!(var_hap_table.is_methylated(1, false), "top-strand C at hap pos 1 must be set");
        assert!(var_hap_table.is_methylated(2, true), "bottom-strand C at hap pos 2 must be set");

        // Reference haplotype: no CpG at all.
        let mut rng2 = SmallRng::seed_from_u64(42);
        let ref_hap_table =
            MethylationTable::from_haplotype(&haps[0], reference, &uniform_model(1.0), &mut rng2);
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
        let table =
            MethylationTable::from_haplotype(&haps[1], reference, &uniform_model(1.0), &mut rng);

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
        let table =
            MethylationTable::from_haplotype(&haps[1], reference, &uniform_model(1.0), &mut rng);

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
        let table =
            MethylationTable::from_haplotype(&haps[1], reference, &uniform_model(1.0), &mut rng);

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
        let table =
            MethylationTable::from_haplotype(&haps[1], reference, &uniform_model(1.0), &mut rng);

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
    #[should_panic(expected = "invalid MethylationModel")]
    fn test_from_haplotype_rejects_nan_rate() {
        let reference = b"ACGT";
        let hap = ref_haplotype();
        let mut rng = SmallRng::seed_from_u64(42);
        let _ =
            MethylationTable::from_haplotype(&hap, reference, &uniform_model(f64::NAN), &mut rng);
    }

    #[test]
    #[should_panic(expected = "invalid MethylationModel")]
    fn test_from_haplotype_rejects_rate_above_one() {
        let reference = b"ACGT";
        let hap = ref_haplotype();
        let mut rng = SmallRng::seed_from_u64(42);
        let _ = MethylationTable::from_haplotype(&hap, reference, &uniform_model(1.5), &mut rng);
    }

    #[test]
    fn test_methylation_model_validate_rejects_bad_fields() {
        // Each invalid field is named in the error message.
        let mut m = uniform_model(0.5);
        m.island.rate = 1.5;
        assert!(format!("{}", m.validate().unwrap_err()).contains("--methylation-rate-island"));

        let mut m = uniform_model(0.5);
        m.shore.correlation_length_bp = 0.0;
        assert!(
            format!("{}", m.validate().unwrap_err())
                .contains("--methylation-correlation-length-shore")
        );

        let mut m = uniform_model(0.5);
        m.hemi_rate = f64::NAN;
        assert!(format!("{}", m.validate().unwrap_err()).contains("--hemimethylation-rate"));
    }

    #[test]
    fn test_from_haplotype_zero_rate_no_methylation() {
        let reference = b"ACGTACGTACGT";
        let hap = ref_haplotype();
        let mut rng = SmallRng::seed_from_u64(42);
        let table =
            MethylationTable::from_haplotype(&hap, reference, &uniform_model(0.0), &mut rng);
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
        let table =
            MethylationTable::from_haplotype(&hap, reference, &uniform_model(1.0), &mut rng);
        assert!(table.is_methylated(1, false), "lowercase 'cg' must register top-strand C");
        assert!(table.is_methylated(2, true), "lowercase 'cg' must register bottom-strand C");
        assert!(!table.is_methylated(0, false));
        assert!(!table.is_methylated(3, true));
    }

    #[test]
    fn test_from_haplotype_symmetric_by_default() {
        // With hemi_rate 0 every methylated CpG sets BOTH strands and every
        // unmethylated CpG sets neither — no hemimethylation. Replaces the old
        // independent-strand-draw test: the model is now symmetric by default.
        let reference = acgt_repeat(200);
        let hap = ref_haplotype();
        let mut rng = SmallRng::seed_from_u64(7);
        let table =
            MethylationTable::from_haplotype(&hap, &reference, &model(0.5, 1000.0, 0.0), &mut rng);
        for site in 0u32..200 {
            let top = table.is_methylated(site * 4 + 1, false);
            let bot = table.is_methylated(site * 4 + 2, true);
            assert_eq!(top, bot, "site {site}: strands must agree with hemi_rate 0");
        }
    }

    #[test]
    fn test_walk_stationary_mean_matches_rate() {
        // Tiny correlation length → near-independent draws, so the mean over
        // many CpGs converges to the target rate. Decouples the stationary
        // distribution from the correlation knob. Empirical band, seed-pinned
        // (SmallRng on rand 0.9); widen if the RNG stream changes.
        let reference = acgt_repeat(5000); // 5000 CpGs, spaced 4 bp
        let hap = ref_haplotype();
        let mut rng = SmallRng::seed_from_u64(42);
        let table =
            MethylationTable::from_haplotype(&hap, &reference, &model(0.3, 1.0, 0.0), &mut rng);
        let meth = (0u32..5000).filter(|&s| table.is_methylated(s * 4 + 1, false)).count();
        let frac = meth as f64 / 5000.0;
        assert!((0.27..=0.33).contains(&frac), "stationary mean {frac} should be ~0.3");
    }

    #[test]
    fn test_walk_autocorrelation_present_with_long_l_absent_with_short_l() {
        let reference = acgt_repeat(4000);
        let hap = ref_haplotype();
        let agreement = |l: f64, seed: u64| {
            let mut rng = SmallRng::seed_from_u64(seed);
            let table =
                MethylationTable::from_haplotype(&hap, &reference, &model(0.5, l, 0.0), &mut rng);
            let states: Vec<bool> =
                (0u32..4000).map(|s| table.is_methylated(s * 4 + 1, false)).collect();
            let agree = states.windows(2).filter(|w| w[0] == w[1]).count();
            agree as f64 / (states.len() - 1) as f64
        };
        // Long L (>> 4 bp spacing) → neighbours almost always agree.
        assert!(agreement(10_000.0, 1) > 0.9, "long correlation length should give long runs");
        // Tiny L → near-independent → agreement ~0.5 at p=0.5.
        let indep = agreement(1.0, 2);
        assert!((0.42..=0.58).contains(&indep), "tiny L should decorrelate neighbours: {indep}");
    }

    #[test]
    fn test_walk_hemi_rate_produces_hemimethylation() {
        // rate 1.0 → every CpG methylated; hemi 0.3 → ~30% become hemi (one
        // strand). Both strand-drop directions should occur. Seed-pinned band.
        let reference = acgt_repeat(3000);
        let hap = ref_haplotype();
        let mut rng = SmallRng::seed_from_u64(42);
        let table =
            MethylationTable::from_haplotype(&hap, &reference, &model(1.0, 1.0, 0.3), &mut rng);
        let (mut hemi, mut top_only, mut bot_only) = (0usize, 0usize, 0usize);
        for s in 0u32..3000 {
            let top = table.is_methylated(s * 4 + 1, false);
            let bot = table.is_methylated(s * 4 + 2, true);
            if top != bot {
                hemi += 1;
                if top {
                    top_only += 1;
                } else {
                    bot_only += 1;
                }
            }
        }
        let frac = hemi as f64 / 3000.0;
        assert!((0.25..=0.35).contains(&frac), "hemi fraction {frac} should be ~0.3");
        assert!(top_only > 100 && bot_only > 100, "both hemi directions should occur");
    }

    #[test]
    fn test_walk_island_hypomethylated_relative_to_open_sea() {
        // Build a CpG island ("CG" repeats, 100% GC, CpG-dense) embedded in a
        // long AT-rich open-sea background (sparse CpGs, 20% GC). With the
        // realistic context defaults the island CpGs should be markedly less
        // methylated than the far open-sea CpGs. Small L so each region's mean
        // converges rather than locking into one long run.
        let os_unit = b"AATTCGAATT"; // CG at offset 4, 20% GC → open-sea
        let os_units = 2000usize; // 20 kb each flank
        let island_units = 1500usize; // 3 kb island of "CG"
        let mut reference = Vec::new();
        for _ in 0..os_units {
            reference.extend_from_slice(os_unit);
        }
        let island_start = reference.len();
        for _ in 0..island_units {
            reference.extend_from_slice(b"CG");
        }
        let island_end = reference.len();
        for _ in 0..os_units {
            reference.extend_from_slice(os_unit);
        }

        let hap = ref_haplotype();
        let mut rng = SmallRng::seed_from_u64(42);
        let m = MethylationModel {
            island: ContextParams { rate: 0.1, correlation_length_bp: 10.0 },
            shore: ContextParams { rate: 0.5, correlation_length_bp: 10.0 },
            open_sea: ContextParams { rate: 0.85, correlation_length_bp: 10.0 },
            hemi_rate: 0.0,
        };
        let table = MethylationTable::from_haplotype(&hap, &reference, &m, &mut rng);

        let cpgs = find_reference_cpgs(&reference);
        let (mut isl_m, mut isl_n, mut sea_m, mut sea_n) = (0usize, 0usize, 0usize, 0usize);
        for &c in &cpgs {
            let cu = c as usize;
            let methylated = table.is_methylated(c, false);
            // Deep island interior vs open-sea well beyond the 2 kb shore +
            // relaxation band — avoids transition-band ambiguity.
            if cu >= island_start + 200 && cu < island_end - 200 {
                isl_n += 1;
                isl_m += usize::from(methylated);
            } else if cu + 6000 < island_start || cu > island_end + 6000 {
                sea_n += 1;
                sea_m += usize::from(methylated);
            }
        }
        let isl_frac = isl_m as f64 / isl_n as f64;
        let sea_frac = sea_m as f64 / sea_n as f64;
        assert!(isl_frac < 0.35, "island interior should be hypomethylated, got {isl_frac}");
        assert!(sea_frac > 0.65, "open sea should be hypermethylated, got {sea_frac}");
        assert!(isl_frac < sea_frac, "island must be less methylated than open sea");
    }

    #[test]
    fn test_walk_determinism_fixed_seed() {
        let reference = acgt_repeat(500);
        let hap = ref_haplotype();
        let m = model(0.6, 500.0, 0.05);
        let run = || {
            let mut rng = SmallRng::seed_from_u64(123);
            let t = MethylationTable::from_haplotype(&hap, &reference, &m, &mut rng);
            (0u32..2000)
                .map(|p| (t.is_methylated(p, false), t.is_methylated(p, true)))
                .collect::<Vec<_>>()
        };
        assert_eq!(run(), run(), "same seed must produce identical methylation");
    }

    #[test]
    fn test_walk_no_cpg_and_single_cpg_do_not_panic() {
        let hap = ref_haplotype();
        let mut rng = SmallRng::seed_from_u64(1);
        // No CpG anywhere.
        let t = MethylationTable::from_haplotype(&hap, b"AAAATTTT", &uniform_model(1.0), &mut rng);
        for i in 0..8 {
            assert!(!t.is_methylated(i, false) && !t.is_methylated(i, true));
        }
        // Single CpG, full methylation → both strands set.
        let t = MethylationTable::from_haplotype(&hap, b"ACGT", &uniform_model(1.0), &mut rng);
        assert!(t.is_methylated(1, false) && t.is_methylated(2, true));
    }

    #[test]
    fn test_walk_shore_rate_is_intermediate_between_island_and_open_sea() {
        // Guards the Shore branch of `params_for`: shore CpGs (within 2 kb of
        // the island, but outside it) must take the shore rate, landing
        // between island (hypo) and open-sea (hyper). A bug routing Shore to
        // the wrong ContextParams would otherwise pass every other test, which
        // all exclude the shore band. Distinct rates + small L for tight means.
        let os_unit = b"AATTCGAATT"; // CG every 10 bp, 20% GC → open-sea
        let mut seq = b"CG".repeat(150); // 300 bp island at the start
        let island_end = seq.len();
        for _ in 0..1500 {
            seq.extend_from_slice(os_unit); // 15 kb of CpG-bearing open-sea
        }
        let hap = ref_haplotype();
        let mut rng = SmallRng::seed_from_u64(42);
        let m = MethylationModel {
            island: ContextParams { rate: 0.1, correlation_length_bp: 30.0 },
            shore: ContextParams { rate: 0.5, correlation_length_bp: 30.0 },
            open_sea: ContextParams { rate: 0.85, correlation_length_bp: 30.0 },
            hemi_rate: 0.0,
        };
        let table = MethylationTable::from_haplotype(&hap, &seq, &m, &mut rng);

        let cpgs = find_reference_cpgs(&seq);
        let mean_over = |lo: usize, hi: usize| {
            let v: Vec<bool> = cpgs
                .iter()
                .filter(|&&c| (c as usize) >= lo && (c as usize) < hi)
                .map(|&c| table.is_methylated(c, false))
                .collect();
            assert!(!v.is_empty(), "no CpGs in [{lo}, {hi})");
            v.iter().filter(|&&b| b).count() as f64 / v.len() as f64
        };
        // Island interior; shore band (within 2 kb, past the relaxation zone);
        // far open sea (> 2 kb past the island).
        let island = mean_over(50, island_end - 20);
        let shore = mean_over(island_end + 400, island_end + 1900);
        let open_sea = mean_over(island_end + 6000, seq.len());
        assert!(island < shore, "island {island} should be below shore {shore}");
        assert!(shore < open_sea, "shore {shore} should be below open sea {open_sea}");
        assert!((0.3..=0.7).contains(&shore), "shore mean {shore} should be ~0.5");
    }

    #[test]
    fn test_from_haplotypes_allele_specific_methylation() {
        // Two haplotypes with IDENTICAL CpG content (no variants → both are the
        // reference) must receive DIFFERENT methylation patterns, because each
        // haplotype is walked as an independent chain. Pins the advertised
        // allele-specific-methylation behaviour, which the rate-1.0/0.0 tests
        // (deterministic) cannot exercise.
        let reference = acgt_repeat(2000); // 2000 CpGs, identical on both haps
        let haps = build_haplotypes(&[], 2, &mut SmallRng::seed_from_u64(0));
        let mut rng = SmallRng::seed_from_u64(42);
        let cm =
            ContigMethylation::from_haplotypes(&haps, &reference, &model(0.5, 1.0, 0.0), &mut rng);
        let (h0, h1) = (cm.table_for(0), cm.table_for(1));
        let differ = (0u32..2000)
            .filter(|&s| h0.is_methylated(s * 4 + 1, false) != h1.is_methylated(s * 4 + 1, false))
            .count();
        // Independent Bernoulli(0.5) draws differ ~50% of the time; require a
        // large fraction to rule out shared/duplicated draws across haplotypes.
        assert!(
            differ > 600,
            "haplotypes should diverge at many CpGs; only {differ}/2000 differed"
        );
    }

    // --- CpG-island detector tests ---

    #[test]
    fn test_island_mask_detects_gc_cpg_dense_block() {
        // 300 bp "CG" island flanked by AT-rich sequence.
        let mut seq = b"AT".repeat(400); // 800 bp AT
        let island_start = seq.len();
        seq.extend_from_slice(&b"CG".repeat(150)); // 300 bp island
        let island_end = seq.len();
        seq.extend_from_slice(&b"AT".repeat(400));
        let mask = island_mask(&seq);
        // Interior of the CG block is island; AT flanks are not.
        assert!(mask[island_start + 150], "CG-dense interior should be island");
        assert!(!mask[10], "AT-rich flank should not be island");
        assert!(!mask[island_end + 400], "far AT flank should not be island");
    }

    #[test]
    fn test_island_mask_none_in_at_rich_or_short() {
        assert!(island_mask(&b"AT".repeat(500)).not_any(), "AT-rich → no island");
        assert!(island_mask(b"CGCGCG").not_any(), "sub-window sequence → no island");
    }

    #[test]
    fn test_classify_cpg_contexts_island_shore_open_sea() {
        // Island of CG at the start, then a long AT-rich tail with sparse CpGs.
        let mut seq = b"CG".repeat(150); // 300 bp island
        let island_end = seq.len();
        seq.extend_from_slice(&b"AATTCGAATT".repeat(1000)); // CpGs every 10 bp out to ~10 kb
        let cpgs = find_reference_cpgs(&seq);
        let mask = island_mask(&seq);
        let contexts = classify_cpg_contexts(&cpgs, &mask);
        // First CpG (inside island) is Island.
        assert_eq!(contexts[0], CpgContext::Island);
        // A CpG ~1 kb past the island is Shore; one ~5 kb past is OpenSea.
        let ctx_at = |target: usize| {
            let idx = cpgs.iter().position(|&c| c as usize >= target).unwrap();
            contexts[idx]
        };
        assert_eq!(ctx_at(island_end + 1000), CpgContext::Shore, "within 2 kb → shore");
        assert_eq!(ctx_at(island_end + 5000), CpgContext::OpenSea, "beyond 2 kb → open sea");
    }

    #[test]
    fn test_from_haplotype_empty_and_short() {
        let hap = ref_haplotype();
        let mut rng = SmallRng::seed_from_u64(42);
        let table = MethylationTable::from_haplotype(&hap, b"", &uniform_model(1.0), &mut rng);
        assert!(table.is_empty());
        assert_eq!(table.len(), 0);

        let table = MethylationTable::from_haplotype(&hap, b"C", &uniform_model(1.0), &mut rng);
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
        let cm =
            ContigMethylation::from_haplotypes(&haps, reference, &uniform_model(1.0), &mut rng);

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
