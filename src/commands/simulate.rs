use std::fs::File;
use std::io::BufWriter;

use anyhow::{Result, bail};
use clap::Parser;
use pooled_writer::PoolBuilder;
use pooled_writer::bgzf::BgzfCompressor;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, Normal};

use super::command::{Command, output_path};
use super::common::{BedOptions, OutputPrefixOptions, ReferenceOptions, SeedOptions, VcfOptions};
use crate::bed::{PaddedIntervalSampler, TargetRegions};
use crate::error_model::illumina::IlluminaErrorModel;
use crate::fasta::Fasta;
use crate::fragment::extract_fragment;
use crate::haplotype::build_haplotypes;
use crate::output::fastq::FastqWriter;
use crate::output::golden_bam::{GoldenBamMetadata, GoldenBamWriter};
use crate::read::generate_read_pair;
use crate::seed::{derive_seed, resolve_seed};
use crate::sequence_dict::SequenceDictionary;
use crate::version::VERSION;

/// Default sample name used in the golden BAM `@RG SM` field when the
/// simulation is not driven by a VCF sample.
const DEFAULT_SAMPLE_NAME: &str = "holodeck-simulation";

/// Default Illumina TruSeq adapter sequence for read 1.
const DEFAULT_ADAPTER_R1: &str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";

/// Default Illumina TruSeq adapter sequence for read 2.
const DEFAULT_ADAPTER_R2: &str = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

/// Default per-cytosine conversion rate for molecules in the converted camp
/// when `--methylation-mode` is set without `--methylation-conversion-rate`.
const DEFAULT_CONVERSION_RATE: f64 = 0.999;

/// Default fraction of whole-molecule conversion failures when
/// `--methylation-mode` is set without `--methylation-failure-rate`.
const DEFAULT_FAILURE_RATE: f64 = 0.01;

/// Simulate sequencing reads from a reference genome.
///
/// Generates paired-end or single-end FASTQ files with optional ground-truth
/// BAM and VCF outputs for benchmarking alignment and variant calling pipelines.
/// Variants are applied from an input VCF to construct haplotype sequences, and
/// reads are sampled with a position-dependent Illumina error model.
#[derive(Parser, Debug)]
#[command(after_long_help = "EXAMPLES:\n  \
    holodeck simulate -r ref.fa -o out --coverage 30\n  \
    holodeck simulate -r ref.fa -v vars.vcf -o out --coverage 30 --golden-bam\n  \
    holodeck simulate -r ref.fa -v vars.vcf -b targets.bed -o out --coverage 100")]
#[allow(clippy::struct_excessive_bools)] // CLI flags are naturally boolean
pub struct Simulate {
    #[command(flatten)]
    pub reference: ReferenceOptions,

    #[command(flatten)]
    pub vcf: VcfOptions,

    #[command(flatten)]
    pub bed: BedOptions,

    #[command(flatten)]
    pub output: OutputPrefixOptions,

    #[command(flatten)]
    pub seed: SeedOptions,

    /// Mean sequencing coverage depth.
    #[arg(short = 'c', long, default_value_t = 30.0, value_name = "FLOAT")]
    pub coverage: f64,

    /// Length of each read in bases.
    #[arg(short = 'l', long, default_value_t = 150, value_name = "INT")]
    pub read_length: usize,

    /// Mean insert size (outer distance between read pair ends).
    #[arg(short = 'd', long, default_value_t = 300, value_name = "INT")]
    pub fragment_mean: usize,

    /// Standard deviation of the insert size distribution.
    #[arg(short = 's', long, default_value_t = 50, value_name = "INT")]
    pub fragment_stddev: usize,

    /// Minimum fragment length in bases; sampled lengths below this are clamped
    /// up.  Fragments shorter than the read length are padded with adapter
    /// sequence.  Must be at least 1.
    #[arg(long, default_value_t = 20, value_name = "INT")]
    pub min_fragment_length: usize,

    /// Generate single-end reads instead of paired-end.
    #[arg(long)]
    pub single_end: bool,

    /// Adapter sequence appended to read 1 when the fragment is shorter than
    /// the read length.
    #[arg(long, default_value = DEFAULT_ADAPTER_R1, value_name = "SEQ")]
    pub adapter_r1: String,

    /// Adapter sequence appended to read 2 when the fragment is shorter than
    /// the read length.
    #[arg(long, default_value = DEFAULT_ADAPTER_R2, value_name = "SEQ")]
    pub adapter_r2: String,

    /// Minimum per-base error rate, applied at the start of reads.
    #[arg(long, default_value_t = 0.001, value_name = "FLOAT")]
    pub min_error_rate: f64,

    /// Maximum per-base error rate, applied at the end of reads.
    #[arg(long, default_value_t = 0.01, value_name = "FLOAT")]
    pub max_error_rate: f64,

    /// Maximum fraction of a read's bases that may come from ambiguous
    /// reference positions (IUPAC codes or `N`). Reference bases resolved
    /// from ambiguity codes are flagged at load time and counted per read;
    /// read pairs whose R1 or R2 exceeds this fraction are rejected and
    /// resampled. Setting this to `1.0` disables the filter; setting it to
    /// `0.0` requires every base to come from an unambiguous reference
    /// position. Value must be in `[0.0, 1.0]`.
    #[arg(long, default_value_t = 0.02, value_name = "FLOAT")]
    pub max_n_frac: f64,

    /// Probability in `[0.0, 1.0]` that a read's 5' end carries a terminal
    /// soft-clip artifact (end-repair fill-in, damaged or non-templated read
    /// ends). The clipped bases are replaced with random sequence so a
    /// downstream aligner soft-clips them, and the golden BAM records the clip
    /// in the CIGAR so ground truth stays exact. Real Twist EM-seq shows ~8%
    /// of reads 5'-clipped. Defaults to `0.0` (disabled). Protocol-agnostic:
    /// independent of `--methylation-mode`.
    #[arg(long, default_value_t = 0.0, value_name = "FLOAT")]
    pub clip_5p_rate: f64,

    /// Probability in `[0.0, 1.0]` that a read's 3' end carries a terminal
    /// soft-clip artifact, independent of (and in addition to) 3' adapter
    /// read-through. Real Twist EM-seq shows ~2% of reads 3'-clipped. Defaults
    /// to `0.0` (disabled).
    #[arg(long, default_value_t = 0.0, value_name = "FLOAT")]
    pub clip_3p_rate: f64,

    /// Mean length in bases of an injected terminal soft-clip. Lengths are
    /// drawn from a truncated geometric distribution centered here and capped
    /// by `--clip-length-max`. Only takes effect when a clip rate is non-zero.
    #[arg(long, default_value_t = 8, value_name = "INT")]
    pub clip_length_mean: usize,

    /// Maximum length in bases of an injected terminal soft-clip. Only takes
    /// effect when a clip rate is non-zero.
    #[arg(long, default_value_t = 20, value_name = "INT")]
    pub clip_length_max: usize,

    /// Enable methylation chemistry simulation. `em-seq` (or `bisulfite`)
    /// converts unmethylated cytosines to thymine and preserves methylated
    /// ones (matches both bisulfite and em-seq protocols). `taps` is the
    /// inverse: methylated cytosines convert to thymine while unmethylated
    /// ones are preserved. When omitted, no methylation chemistry is
    /// applied. Requires a methylation-annotated VCF (MT/MB FORMAT fields)
    /// produced by `holodeck methylate`.
    #[arg(long, value_enum, value_name = "MODE")]
    pub methylation_mode: Option<crate::meth::MethylationMode>,

    /// Probability that the converting class of cytosines (unmethylated
    /// for em-seq/bisulfite, methylated for TAPS) is converted to thymine
    /// in a molecule that converted normally. `1.0` is perfect chemistry;
    /// real protocols typically achieve `0.99`–`0.999`. Only takes effect
    /// when `--methylation-mode` is set; passing this flag without a mode is
    /// rejected. Defaults to `0.999` when `--methylation-mode` is set
    /// without this flag.
    #[arg(long, value_name = "FLOAT")]
    pub methylation_conversion_rate: Option<f64>,

    /// Fraction of molecules that are whole-molecule conversion failures.
    /// Real bisulfite/EM-seq conversion is bimodal: most molecules convert
    /// near-completely, a small fraction (fragments that fail to denature or
    /// re-anneal too fast) escape conversion as a coherent unit. A failed
    /// molecule converts its should-convert cytosines at
    /// `1 - methylation-conversion-rate` (near-zero), so it retains almost
    /// all of them as C across both mates. The golden BAM stamps each
    /// record with `cf:i:{0|1}` recording which camp the molecule was in.
    /// Only takes effect when `--methylation-mode` is set; passing this flag
    /// without a mode is rejected. Defaults to `0.01` (1%) when
    /// `--methylation-mode` is set without this flag; pass `0.0` to disable.
    #[arg(long, value_name = "FLOAT")]
    pub methylation_failure_rate: Option<f64>,

    /// Write a per-CpG ground-truth methylation tally in MethylDackel's
    /// `extract` CpG bedGraph format. Columns: `chrom  start  end  rate
    /// (0–100)  n_methylated  n_unmethylated`. For every reference CpG
    /// covered by at least one simulated read, counts how many reads called
    /// the site as methylated vs unmethylated according to the simulator's
    /// per-haplotype, per-strand methylation bitmap (NOT the post-conversion
    /// sequenced base). Lets downstream evaluators measure aligner-derived
    /// methylation calls against ground truth using the same code paths
    /// that consume MethylDackel output. Requires `--methylation-mode`.
    #[arg(long, value_name = "PATH")]
    pub cpg_truth_bedgraph: Option<std::path::PathBuf>,

    /// Write a ground-truth BAM file with correct alignments.
    #[arg(long)]
    pub golden_bam: bool,

    /// Write a ground-truth VCF annotated with simulated coverage.
    #[arg(long)]
    pub golden_vcf: bool,

    /// Use simple read names (`holodeck::N`) instead of encoding truth
    /// coordinates in the read name.
    #[arg(long)]
    pub simple_names: bool,

    /// BGZF compression level (0-12). Lower values are faster with larger
    /// output files; higher values produce smaller files at the cost of speed.
    /// Level 0 is no compression; 1 is fastest; 12 is maximum compression.
    #[arg(long, default_value_t = 1, value_name = "INT")]
    pub compression: u8,

    /// Number of threads for parallel BGZF output compression.
    #[arg(short = 't', long, default_value_t = 4, value_name = "INT")]
    pub threads: usize,
}

impl Command for Simulate {
    fn execute(&self) -> Result<()> {
        let resolved_sample = self.validate()?;
        self.run_simulation(resolved_sample.as_deref())
    }
}

impl Simulate {
    /// Validate command-line arguments before running and return the
    /// resolved VCF sample name (if a VCF was provided).  Resolving the
    /// sample here means the VCF header is read only once per run.
    fn validate(&self) -> Result<Option<String>> {
        if !self.coverage.is_finite() || self.coverage <= 0.0 {
            bail!("--coverage must be a finite positive number");
        }
        if self.read_length == 0 {
            bail!("--read-length must be > 0");
        }
        if self.min_fragment_length == 0 {
            bail!("--min-fragment-length must be at least 1");
        }
        if self.min_error_rate < 0.0 || self.max_error_rate < 0.0 {
            bail!("Error rates must be >= 0");
        }
        if self.min_error_rate > self.max_error_rate {
            bail!("--min-error-rate must be <= --max-error-rate");
        }
        if !self.max_n_frac.is_finite() || !(0.0..=1.0).contains(&self.max_n_frac) {
            bail!("--max-n-frac must be in [0.0, 1.0]");
        }
        if !self.clip_5p_rate.is_finite() || !(0.0..=1.0).contains(&self.clip_5p_rate) {
            bail!("--clip-5p-rate must be in [0.0, 1.0]");
        }
        if !self.clip_3p_rate.is_finite() || !(0.0..=1.0).contains(&self.clip_3p_rate) {
            bail!("--clip-3p-rate must be in [0.0, 1.0]");
        }
        // Validate clip-length parameters only when the model is enabled, so a
        // disabled model (both rates zero) never rejects on unused values.
        if self.clip_5p_rate > 0.0 || self.clip_3p_rate > 0.0 {
            if self.clip_length_max == 0 {
                bail!("--clip-length-max must be > 0 when a terminal clip rate is set");
            }
            if self.clip_length_mean > self.clip_length_max {
                bail!("--clip-length-mean must be <= --clip-length-max");
            }
        }
        if let Some(rate) = self.methylation_conversion_rate
            && (!rate.is_finite() || !(0.0..=1.0).contains(&rate))
        {
            bail!("--methylation-conversion-rate must be in [0.0, 1.0]");
        }
        // Reject orphaned methylation conversion-rate flag: passing a tuning
        // rate without enabling a chemistry mode is almost certainly a user
        // mistake, and silently dropping the request is worse than bailing.
        if self.methylation_mode.is_none() && self.methylation_conversion_rate.is_some() {
            bail!("--methylation-conversion-rate requires --methylation-mode");
        }
        if let Some(rate) = self.methylation_failure_rate
            && (!rate.is_finite() || !(0.0..=1.0).contains(&rate))
        {
            bail!("--methylation-failure-rate must be in [0.0, 1.0]");
        }
        if self.methylation_mode.is_none() && self.methylation_failure_rate.is_some() {
            bail!("--methylation-failure-rate requires --methylation-mode");
        }
        if self.cpg_truth_bedgraph.is_some() && self.methylation_mode.is_none() {
            bail!("--cpg-truth-bedgraph requires --methylation-mode");
        }

        // Validate `--cpg-truth-bedgraph` parent directory exists up front so
        // a typo'd path doesn't waste an entire simulation run before
        // `write_bedgraph` discovers it at the very end.
        if let Some(path) = &self.cpg_truth_bedgraph
            && let Some(parent) = path.parent()
            && !parent.as_os_str().is_empty()
            && !parent.exists()
        {
            bail!("CpG truth bedGraph directory does not exist: {}", parent.display());
        }

        // Validation matrix for methylation chemistry simulation:
        // - --methylation-mode requires the input VCF to actually carry MT/MB
        //   methylation records (run `holodeck methylate` first to produce one).
        // - A methylated VCF without --methylation-mode produces a warning and
        //   variants-only output (user almost certainly forgot the mode flag).
        //
        // Probe for real records, not just header declarations: a VCF that
        // declares MT/MB in the header but has no annotated records carries no
        // methylation truth, and accepting it here would defer the failure to
        // an internal-invariant panic deep in the per-contig loop.
        let vcf_has_methylation = if let Some(vcf_path) = &self.vcf.vcf {
            crate::vcf::methylation::vcf_has_mt_mb_records(vcf_path)
                .map_err(|e| anyhow::anyhow!("failed to inspect VCF for MT/MB records: {e}"))?
        } else {
            false
        };
        match (vcf_has_methylation, self.methylation_mode.is_some()) {
            (true, true) | (false, false) => {} // normal paths
            (true, false) => {
                log::warn!(
                    "VCF contains methylation truth (MT/MB) but --methylation-mode is not set; \
                     methylation chemistry will not be applied"
                );
            }
            (false, true) => {
                bail!(
                    "--methylation-mode requires a methylation-annotated VCF (MT/MB FORMAT \
                     fields); run `holodeck methylate` first"
                );
            }
        }
        if self.compression > 12 {
            bail!("--compression must be between 0 and 12");
        }

        // --sample without --vcf is nonsensical.
        if self.vcf.sample.is_some() && self.vcf.vcf.is_none() {
            bail!("--sample requires --vcf");
        }

        // Validate VCF sample configuration upfront so the user gets a clear
        // error before the simulation loop starts, and capture the resolved
        // sample name for use in downstream metadata (e.g. the golden BAM
        // `@RG` line).
        let resolved_sample = if let Some(vcf_path) = &self.vcf.vcf {
            Some(crate::vcf::validate_vcf_sample(vcf_path, self.vcf.sample.as_deref())?)
        } else {
            None
        };

        // Validate output parent directory exists.
        if let Some(parent) = self.output.output.parent()
            && !parent.as_os_str().is_empty()
            && !parent.exists()
        {
            bail!("Output directory does not exist: {}", parent.display());
        }

        Ok(resolved_sample)
    }

    /// Run the main simulation pipeline.
    ///
    /// `resolved_vcf_sample` is the sample name resolved from the VCF during
    /// validation (if any), used for the golden BAM `@RG SM` field.
    #[allow(clippy::too_many_lines)] // Top-level orchestrator.
    fn run_simulation(&self, resolved_vcf_sample: Option<&str>) -> Result<()> {
        let seed = self.compute_seed();
        let mut rng = SmallRng::seed_from_u64(seed);
        log::info!("Using random seed: {seed}");

        let mut fasta = Fasta::from_path(&self.reference.reference)?;
        let dict = fasta.dict().clone();
        log::info!(
            "Loaded reference with {} contigs, total {} bp",
            dict.len(),
            dict.total_length()
        );

        let targets = self.load_targets(&dict)?;
        let effective_size = targets
            .as_ref()
            .map_or(dict.total_length(), |t| t.effective_territory(self.fragment_mean));
        if effective_size == 0 {
            bail!("Effective genome size is 0; nothing to simulate");
        }

        let total_reads = self.compute_total_reads(effective_size);
        log::info!("Will generate {total_reads} read pairs for {:.1}x coverage", self.coverage);

        let error_model =
            IlluminaErrorModel::new(self.read_length, self.min_error_rate, self.max_error_rate);

        // Terminal soft-clip artifact model. Passed as `Some` only when at
        // least one end can clip; a disabled model is equivalent to `None`
        // (it draws no randomness) but keeping it `None` makes intent clear.
        let clip_config = (self.clip_5p_rate > 0.0 || self.clip_3p_rate > 0.0).then_some(
            crate::clip::TerminalClipConfig {
                rate_5p: self.clip_5p_rate,
                rate_3p: self.clip_3p_rate,
                length_mean: self.clip_length_mean,
                length_max: self.clip_length_max,
            },
        );
        let frag_dist = Normal::new(self.fragment_mean as f64, self.fragment_stddev as f64)
            .map_err(|e| anyhow::anyhow!("Invalid fragment distribution parameters: {e}"))?;

        // Normalize user-supplied adapter sequences to uppercase so they
        // don't register as ambiguity-resolved bases (lowercase marker) in
        // the per-read lowercase-fraction filter.
        let adapter_r1 = self.adapter_r1.to_ascii_uppercase();
        let adapter_r2 = self.adapter_r2.to_ascii_uppercase();

        let compression = self.compression;
        let use_pool = self.threads > 1;

        // When using multiple threads, create a shared compression pool that
        // handles BGZF block compression and writing across all output files.
        let mut pool_builder: Option<PoolBuilder<BufWriter<File>, BgzfCompressor>> = if use_pool {
            let pb = PoolBuilder::new()
                .threads(self.threads)
                .compression_level(compression)
                .map_err(|e| anyhow::anyhow!("failed to set compression level: {e}"))?;
            log::info!("Using {} threads for BGZF compression", self.threads);
            Some(pb)
        } else {
            None
        };

        let mut r1_writer = self.create_fastq_writer(".r1.fastq.gz", &mut pool_builder)?;
        let mut r2_writer = if self.single_end {
            None
        } else {
            Some(self.create_fastq_writer(".r2.fastq.gz", &mut pool_builder)?)
        };

        let mut golden_bam_writer = if self.golden_bam {
            let bam_path = output_path(&self.output.output, ".golden.bam");
            log::info!("Writing golden BAM to: {}", bam_path.display());
            let meta = Self::golden_bam_metadata(resolved_vcf_sample);
            if let Some(pb) = &mut pool_builder {
                let file = File::create(&bam_path)?;
                let pooled = pb.exchange(BufWriter::new(file));
                Some(GoldenBamWriter::from_writer(Box::new(pooled), &dict, &meta)?)
            } else {
                Some(GoldenBamWriter::new(&bam_path, &dict, compression, &meta)?)
            }
        } else {
            None
        };

        // Build the pool after all writers have been exchanged.
        let mut pool = pool_builder
            .map(PoolBuilder::build)
            .transpose()
            .map_err(|e| anyhow::anyhow!("failed to build compression pool: {e}"))?;

        if self.golden_vcf {
            log::warn!("--golden-vcf is not yet implemented; skipping");
        }

        // CpG truth tally is built across all contigs and written once at the
        // end. Validation in `validate()` ensures `--cpg-truth-bedgraph` is
        // only set when methylation chemistry is also enabled.
        let mut cpg_truth = self
            .cpg_truth_bedgraph
            .as_ref()
            .map(|_| crate::output::cpg_truth::CpgTruthTally::new());

        let mut read_num: u64 = 0;
        let contig_names: Vec<String> = dict.names().into_iter().map(String::from).collect();

        // Parse the methylation-annotated VCF (if any) once up front so the
        // per-contig loop can pluck out per-contig records without
        // re-decompressing or re-scanning the file every iteration.
        // Validation in `validate()` ensures that `--methylation-mode` is only
        // accepted when the VCF actually declares MT/MB.
        let methylation_records = if self.methylation_mode.is_some()
            && let Some(vcf_path) = &self.vcf.vcf
        {
            Some(crate::vcf::methylation::parse_methylation_vcf(vcf_path)?)
        } else {
            None
        };

        // Parse variants once, partitioning by contig. The per-contig loop
        // below looks up its slice via `HashMap::get` instead of re-scanning
        // the full VCF on every contig (the latter is O(contigs × records)
        // and dominates wall time on whole-genome inputs). The scan also
        // resolves the sample's ploidy from unfiltered GTs, so variant-free
        // contigs are sized the same as variant-bearing ones — a per-contig
        // fallback to `2` would silently mis-shape MT/MB on haploid /
        // triploid samples.
        let parsed_variants = if let Some(vcf_path) = &self.vcf.vcf {
            Some(crate::vcf::parse_variants_by_contig(vcf_path, resolved_vcf_sample, &dict)?)
        } else {
            None
        };
        let sample_ploidy = parsed_variants.as_ref().map_or(2, |p| p.sample_ploidy);

        for contig_name in &contig_names {
            let variants_for_contig: &[crate::vcf::genotype::VariantRecord] = parsed_variants
                .as_ref()
                .and_then(|p| p.by_contig.get(contig_name))
                .map_or(&[], Vec::as_slice);
            read_num += self.simulate_contig(
                contig_name,
                &dict,
                &mut fasta,
                targets.as_ref(),
                total_reads,
                &error_model,
                clip_config.as_ref(),
                &frag_dist,
                adapter_r1.as_bytes(),
                adapter_r2.as_bytes(),
                &mut r1_writer,
                &mut r2_writer,
                &mut golden_bam_writer,
                variants_for_contig,
                sample_ploidy,
                cpg_truth.as_mut(),
                methylation_records.as_ref(),
                read_num,
                seed,
                &mut rng,
            )?;
        }

        // Close writers first so pooled writers flush their buffers to the
        // pool, then stop the pool to wait for all compression/writing to
        // complete.
        r1_writer.close();
        if let Some(w) = r2_writer {
            w.close();
        }
        if let Some(w) = golden_bam_writer {
            w.close();
        }
        if let Some(ref mut p) = pool {
            p.stop_pool().map_err(|e| anyhow::anyhow!("failed to stop compression pool: {e}"))?;
        }

        if let (Some(path), Some(tally)) = (&self.cpg_truth_bedgraph, &cpg_truth) {
            log::info!("Writing CpG truth bedGraph to: {}", path.display());
            tally.write_bedgraph(&dict, path)?;
        }

        log::info!("Generated {read_num} total read pairs");
        Ok(())
    }

    /// Create a FASTQ writer, using the compression pool if available or
    /// single-threaded BGZF otherwise.
    fn create_fastq_writer(
        &self,
        suffix: &str,
        pool_builder: &mut Option<PoolBuilder<BufWriter<File>, BgzfCompressor>>,
    ) -> Result<FastqWriter> {
        let path = output_path(&self.output.output, suffix);
        if let Some(pb) = pool_builder {
            let file = File::create(&path)?;
            let pooled = pb.exchange(BufWriter::new(file));
            Ok(FastqWriter::from_writer(pooled))
        } else {
            FastqWriter::new(&path, self.compression)
        }
    }

    /// Build the `@PG`/`@RG` metadata for the golden BAM header.  The
    /// command line is captured verbatim from `std::env::args_os`, using
    /// lossy UTF-8 conversion so that non-Unicode arguments do not panic.
    /// `resolved_vcf_sample` should be the sample name returned by
    /// [`crate::vcf::validate_vcf_sample`] during validation; when absent,
    /// the sample defaults to [`DEFAULT_SAMPLE_NAME`].
    fn golden_bam_metadata(resolved_vcf_sample: Option<&str>) -> GoldenBamMetadata {
        let command_line = std::env::args_os()
            .map(|arg| arg.to_string_lossy().into_owned())
            .collect::<Vec<_>>()
            .join(" ");
        let sample =
            resolved_vcf_sample.map_or_else(|| DEFAULT_SAMPLE_NAME.to_string(), str::to_string);
        GoldenBamMetadata { command_line, version: VERSION.clone(), sample }
    }

    /// Compute the deterministic seed from simulation parameters.
    ///
    /// IMPORTANT: the legacy 7-field format
    /// `<reference>:<coverage>:<read_length>:<frag_mean>:<frag_stddev>:<min_err>:<max_err>`
    /// MUST NOT change. Existing users rely on the same params producing the
    /// same seed. New parameters (e.g. methylation chemistry) must only
    /// extend the string when explicitly enabled, and only as conditional
    /// suffixes.
    fn compute_seed(&self) -> u64 {
        let mut seed_desc = format!(
            "{}:{}:{}:{}:{}:{}:{}",
            self.reference.reference.display(),
            self.coverage,
            self.read_length,
            self.fragment_mean,
            self.fragment_stddev,
            self.min_error_rate,
            self.max_error_rate,
        );
        if let Some(mode) = self.methylation_mode {
            seed_desc.push(':');
            seed_desc.push_str(mode.as_seed_str());
            seed_desc.push(':');
            seed_desc.push_str(
                &self.methylation_conversion_rate.unwrap_or(DEFAULT_CONVERSION_RATE).to_string(),
            );
            seed_desc.push(':');
            seed_desc.push_str(
                &self.methylation_failure_rate.unwrap_or(DEFAULT_FAILURE_RATE).to_string(),
            );
        }
        resolve_seed(self.seed.seed, &seed_desc)
    }

    /// Load BED target regions if specified.
    fn load_targets(&self, dict: &SequenceDictionary) -> Result<Option<TargetRegions>> {
        match &self.bed.targets {
            Some(bed_path) => {
                let t = TargetRegions::from_path(bed_path, dict)?;
                log::info!("Loaded {} bp of target territory", t.total_territory());
                Ok(Some(t))
            }
            None => Ok(None),
        }
    }

    /// Compute total number of read pairs needed for the requested coverage.
    fn compute_total_reads(&self, effective_size: u64) -> u64 {
        let bases_per_read =
            if self.single_end { self.read_length as u64 } else { self.read_length as u64 * 2 };
        #[expect(clippy::cast_possible_truncation, reason = "read count fits u64")]
        #[expect(clippy::cast_sign_loss, reason = "coverage is positive")]
        let n = ((self.coverage * effective_size as f64) / bases_per_read as f64).round() as u64;
        n
    }

    /// Simulate reads for a single contig.
    ///
    /// `variants` is this contig's slice of the once-parsed
    /// `variants_by_contig` map (empty when no VCF is provided, or when the
    /// VCF has no variants on this contig).
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn simulate_contig(
        &self,
        contig_name: &str,
        dict: &SequenceDictionary,
        fasta: &mut Fasta,
        targets: Option<&TargetRegions>,
        total_reads: u64,
        error_model: &IlluminaErrorModel,
        clip_config: Option<&crate::clip::TerminalClipConfig>,
        frag_dist: &Normal<f64>,
        adapter_r1: &[u8],
        adapter_r2: &[u8],
        r1_writer: &mut FastqWriter,
        r2_writer: &mut Option<FastqWriter>,
        golden_bam: &mut Option<GoldenBamWriter>,
        variants: &[crate::vcf::genotype::VariantRecord],
        sample_ploidy: usize,
        cpg_truth: Option<&mut crate::output::cpg_truth::CpgTruthTally>,
        methylation_records: Option<&crate::vcf::methylation::MethylationVcfRecords>,
        start_read_num: u64,
        main_seed: u64,
        rng: &mut SmallRng,
    ) -> Result<u64> {
        let contig_meta = dict.get_by_name(contig_name).unwrap();
        let contig_len = contig_meta.length() as u64;
        let contig_idx = contig_meta.index();

        // Compute reads proportional to effective territory (if BED) or contig
        // size (whole genome).  For targeted mode, effective territory accounts
        // for the fact that fragments extend beyond targets — see
        // TargetRegions::effective_territory for the derivation.
        let contig_effective_size = targets
            .map_or(contig_len, |t| t.contig_effective_territory(contig_idx, self.fragment_mean));
        let effective_total =
            targets.map_or(dict.total_length(), |t| t.effective_territory(self.fragment_mean));

        if contig_effective_size == 0 || effective_total == 0 {
            return Ok(0);
        }

        #[expect(clippy::cast_possible_truncation, reason = "read count fits u64")]
        #[expect(clippy::cast_sign_loss, reason = "fraction is positive")]
        let contig_reads = (total_reads as f64 * contig_effective_size as f64
            / effective_total as f64)
            .round() as u64;

        if contig_reads == 0 {
            return Ok(0);
        }

        log::info!("Simulating {contig_reads} reads for contig {contig_name} ({contig_len} bp)");

        // Build a padded interval sampler when targets exist.  The pad covers
        // the catchment zone — fragment start positions outside a target whose
        // fragment still extends into the target.  Fragments whose drawn length
        // is too short to actually reach a target are caught by the overlap
        // check below (rare with this padding).
        #[expect(clippy::cast_possible_truncation, reason = "pad fits u32")]
        let sampler = targets.map(|tgt| {
            let pad = (self.fragment_mean + 4 * self.fragment_stddev) as u32;
            PaddedIntervalSampler::new(tgt.contig_intervals(contig_idx), pad, contig_len as u32)
        });

        // Seed a dedicated RNG for reference normalization so ambiguity-code
        // resolution is reproducible and independent of the main sampling RNG.
        let contig_seed = derive_seed(main_seed, contig_name);
        let mut ref_rng = SmallRng::seed_from_u64(contig_seed);
        let reference = fasta.load_contig(contig_name, &mut ref_rng)?;

        // Variants were partitioned by contig once in `run_simulation` —
        // here we just log how many landed on this contig.
        if !variants.is_empty() {
            log::info!("  Loaded {} variants for {contig_name}", variants.len());
        }

        // Size the haplotype set by `sample_ploidy` (resolved once over the
        // whole VCF by `run_simulation`), not by a per-contig `max + 2`
        // fallback. This keeps variant-free contigs aligned with the rest
        // of the sample on haploid / triploid genomes.
        let haplotypes = build_haplotypes(variants, sample_ploidy, rng);

        // Build per-haplotype methylation tables AFTER haplotypes exist so
        // that variant-driven CpG creations/destructions are reflected on
        // each haplotype's bitmap. Methylation truth must come from a VCF
        // produced by `holodeck methylate` (MT/MB FORMAT fields). The
        // validation matrix in `validate()` ensures that when
        // `--methylation-mode` is set, the VCF always declares MT/MB and
        // `methylation_records` is `Some(...)`; the `ok_or_else` below is
        // a defensive guard for that invariant.
        let methylation = if let Some(mode) = self.methylation_mode {
            let records = methylation_records.ok_or_else(|| {
                anyhow::anyhow!(
                    "internal invariant: --methylation-mode set but methylation \
                     records were not loaded (should be caught by validate())"
                )
            })?;
            let cm = crate::vcf::methylation::load_contig_methylation_from_records(
                records,
                contig_name,
                &reference,
                variants,
                sample_ploidy,
            )?
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "internal invariant: --methylation-mode set but VCF lacks MT/MB \
                     (should be caught by validate())"
                )
            })?;
            // The methylation reader derives its ploidy from `variants` the
            // same way we do above, so the two haplotype counts must agree.
            // Surface any drift loudly at contig-load time instead of as a
            // per-read `cm.table_for(hap_idx)` panic deep in the hot loop.
            anyhow::ensure!(
                cm.len() == haplotypes.len(),
                "internal invariant: methylation table count ({}) disagrees with \
                 haplotype count ({}) on contig {contig_name}",
                cm.len(),
                haplotypes.len(),
            );
            Some((cm, mode))
        } else {
            None
        };
        let methylation_config =
            methylation.as_ref().map(|(cm, mode)| crate::meth::MethylationConfig {
                contig_methylation: cm,
                mode: *mode,
                conversion_rate: self
                    .methylation_conversion_rate
                    .unwrap_or(DEFAULT_CONVERSION_RATE),
                failure_rate: self.methylation_failure_rate.unwrap_or(DEFAULT_FAILURE_RATE),
            });

        // Precompute reference CpG positions once per contig (only when the
        // CpG truth bedGraph is requested — otherwise we'd waste a scan).
        // The list is reused for every fragment's per-mate tally.
        let ref_cpgs: Vec<u32> = if cpg_truth.is_some() {
            crate::meth::find_reference_cpgs(&reference)
        } else {
            Vec::new()
        };
        let mut cpg_truth = cpg_truth;

        let mut generated: u64 = 0;
        let mut attempts: u64 = 0;
        let max_attempts = contig_reads * 100;

        while generated < contig_reads && attempts < max_attempts {
            attempts += 1;

            // Draw fragment length, clamped to [min_fragment_length, contig_length].
            // Fragments shorter than read_length are padded with adapter sequence
            // by the read extraction layer.
            #[expect(clippy::cast_possible_truncation, reason = "fragment length fits usize")]
            #[expect(clippy::cast_sign_loss, reason = "clamped to positive")]
            let frag_len = frag_dist
                .sample(rng)
                .round()
                .clamp(self.min_fragment_length as f64, contig_len as f64)
                as usize;

            if frag_len == 0 {
                continue;
            }

            // Pick a random start position.  When targets exist, sample from
            // the padded target regions so that nearly every draw overlaps a
            // target — vastly more efficient than rejection-sampling across
            // the whole contig.
            #[expect(clippy::cast_possible_truncation, reason = "position fits u32")]
            let ref_start = if let Some(samp) = &sampler {
                let s = samp.sample_start(rng).unwrap();
                // Ensure the fragment fits within the contig.
                s.min((contig_len - frag_len as u64) as u32)
            } else {
                let max_start = contig_len - frag_len as u64;
                if max_start > 0 { rng.random_range(0..=max_start) as u32 } else { 0 }
            };

            // Check BED target overlap — with padded sampling this rarely
            // rejects, but catches the occasional short fragment drawn from
            // the pad zone that doesn't reach the target.
            #[expect(clippy::cast_possible_truncation, reason = "frag end fits u32")]
            let frag_end = ref_start + frag_len as u32;
            if let Some(tgt) = targets
                && !tgt.overlaps(contig_idx, ref_start, frag_end)
            {
                continue;
            }

            let hap_idx = rng.random_range(0..haplotypes.len());
            let is_forward: bool = rng.random();
            let fragment =
                extract_fragment(&haplotypes[hap_idx], &reference, ref_start, frag_len, is_forward);

            let read_num = start_read_num + generated + 1;
            let Some(mut pair) = generate_read_pair(
                &fragment,
                contig_name,
                read_num,
                self.read_length,
                !self.single_end,
                adapter_r1,
                adapter_r2,
                self.max_n_frac,
                error_model,
                self.simple_names,
                methylation_config.as_ref(),
                // Pre-conversion bases are only consumed by the golden BAM's
                // YS:Z tag — skip the per-mate clone when no golden BAM is
                // requested.
                self.golden_bam,
                clip_config,
                rng,
            ) else {
                // Too many ambiguity-resolved bases in this read pair; resample.
                continue;
            };

            // When both methylation and a golden BAM are requested, compute
            // the per-record Bismark methylation call tags (XM/YM/NM/MD)
            // before writing. This is the single place where the reference,
            // haplotype, and methylation table are all in scope.
            if self.golden_bam
                && let (Some(annotation), Some((cm, mode))) =
                    (pair.methylation.as_mut(), methylation.as_ref())
            {
                crate::methylation_tags::populate_pair_call_tags(
                    &pair.read1,
                    pair.read2.as_ref(),
                    &pair.r1_truth,
                    pair.r2_truth.as_ref(),
                    &pair.r1_cigar,
                    pair.r2_cigar.as_ref(),
                    &reference,
                    cm.table_for(hap_idx),
                    &haplotypes[hap_idx],
                    *mode,
                    annotation,
                );
            }

            r1_writer.write_read(&pair.read1)?;
            if let Some(w) = r2_writer
                && let Some(r2) = &pair.read2
            {
                w.write_read(r2)?;
            }
            if let Some(bam_w) = golden_bam {
                bam_w.write_pair(&pair)?;
            }

            // Tally CpG truth from each mate's genomic ref-position slice.
            // Mirrors `read::build_mate`'s slicing: forward mates take the
            // left window (`[..genomic]`); negative-strand mates take the
            // right window (`[frag_len - genomic..frag_len]`).
            if let (Some(tally), Some((cm, _))) = (cpg_truth.as_deref_mut(), methylation.as_ref()) {
                let frag_len = fragment.bases.len();
                let genomic = frag_len.min(self.read_length);
                let right_start = frag_len.saturating_sub(genomic);
                let r1_neg = !fragment.is_forward;
                let r1_positions: &[u32] = if r1_neg {
                    &fragment.ref_positions[right_start..frag_len]
                } else {
                    &fragment.ref_positions[..genomic]
                };
                tally.record_mate(
                    contig_idx,
                    r1_positions,
                    &ref_cpgs,
                    cm,
                    &haplotypes[hap_idx],
                    hap_idx,
                    fragment.is_forward,
                );
                if pair.read2.is_some() {
                    let r2_neg = fragment.is_forward;
                    let r2_positions: &[u32] = if r2_neg {
                        &fragment.ref_positions[right_start..frag_len]
                    } else {
                        &fragment.ref_positions[..genomic]
                    };
                    tally.record_mate(
                        contig_idx,
                        r2_positions,
                        &ref_cpgs,
                        cm,
                        &haplotypes[hap_idx],
                        hap_idx,
                        fragment.is_forward,
                    );
                }
            }

            generated += 1;
        }

        if generated < contig_reads {
            log::warn!(
                "Only generated {generated}/{contig_reads} reads for {contig_name} \
                 after {max_attempts} attempts"
            );
        }

        Ok(generated)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::common::{
        BedOptions, OutputPrefixOptions, ReferenceOptions, SeedOptions, VcfOptions,
    };
    use std::path::PathBuf;

    /// Build a `Simulate` populated with the same defaults the CLI uses, so
    /// tests construct realistic instances without going through clap.
    fn make_default_simulate() -> Simulate {
        Simulate {
            reference: ReferenceOptions { reference: PathBuf::from("/tmp/ref.fa") },
            vcf: VcfOptions { vcf: None, sample: None },
            bed: BedOptions { targets: None },
            output: OutputPrefixOptions { output: PathBuf::from("out") },
            seed: SeedOptions { seed: None },
            coverage: 30.0,
            read_length: 150,
            fragment_mean: 300,
            fragment_stddev: 50,
            min_fragment_length: 20,
            single_end: false,
            adapter_r1: DEFAULT_ADAPTER_R1.to_string(),
            adapter_r2: DEFAULT_ADAPTER_R2.to_string(),
            min_error_rate: 0.001,
            max_error_rate: 0.01,
            max_n_frac: 0.02,
            clip_5p_rate: 0.0,
            clip_3p_rate: 0.0,
            clip_length_mean: 8,
            clip_length_max: 20,
            methylation_mode: None,
            methylation_conversion_rate: None,
            methylation_failure_rate: None,
            cpg_truth_bedgraph: None,
            golden_bam: false,
            golden_vcf: false,
            simple_names: false,
            compression: 1,
            threads: 4,
        }
    }

    /// Non-methylation invocations must produce the same seed they did
    /// before the methylation fields were added — i.e. the seed string is
    /// exactly
    /// `"<ref>:<coverage>:<read_length>:<fragment_mean>:<fragment_stddev>:<min_err>:<max_err>"`
    /// with no trailing fields. This pins the legacy 7-field format.
    #[test]
    fn test_compute_seed_no_methylation_matches_legacy() {
        let sim = make_default_simulate();
        let expected = resolve_seed(None, "/tmp/ref.fa:30:150:300:50:0.001:0.01");
        assert_eq!(sim.compute_seed(), expected);
    }

    /// Enabling methylation simulation extends the seed string with the
    /// three extra fields, which yields a different seed than the
    /// non-methylation case.
    #[test]
    fn test_compute_seed_methylation_differs_from_legacy() {
        let mut sim = make_default_simulate();
        let baseline = sim.compute_seed();
        sim.methylation_mode = Some(crate::meth::MethylationMode::EmSeq);
        assert_ne!(sim.compute_seed(), baseline);
    }

    /// An explicit `--seed` short-circuits the description-derived seed,
    /// so toggling methylation must not change the resolved value.
    #[test]
    fn test_compute_seed_explicit_seed_unaffected_by_methylation() {
        let mut sim = make_default_simulate();
        sim.seed.seed = Some(42);
        assert_eq!(sim.compute_seed(), 42);
        sim.methylation_mode = Some(crate::meth::MethylationMode::EmSeq);
        assert_eq!(sim.compute_seed(), 42);
    }

    /// Em-seq and TAPS modes feed different strings into the seed
    /// derivation, so they must produce different seeds with otherwise
    /// identical parameters.
    #[test]
    fn test_compute_seed_em_seq_and_taps_differ() {
        let mut sim_em = make_default_simulate();
        sim_em.methylation_mode = Some(crate::meth::MethylationMode::EmSeq);
        let mut sim_taps = make_default_simulate();
        sim_taps.methylation_mode = Some(crate::meth::MethylationMode::Taps);
        assert_ne!(sim_em.compute_seed(), sim_taps.compute_seed());
    }

    /// `--methylation-failure-rate` feeds the seed string, so two runs that
    /// differ only in failure rate must derive different seeds (otherwise
    /// reruns at a new failure rate would reuse the prior stream).
    #[test]
    fn test_compute_seed_differs_by_failure_rate() {
        let mut sim_a = make_default_simulate();
        sim_a.methylation_mode = Some(crate::meth::MethylationMode::EmSeq);
        sim_a.methylation_failure_rate = Some(0.0);
        let mut sim_b = make_default_simulate();
        sim_b.methylation_mode = Some(crate::meth::MethylationMode::EmSeq);
        sim_b.methylation_failure_rate = Some(0.5);
        assert_ne!(sim_a.compute_seed(), sim_b.compute_seed());
    }

    /// Non-default `--methylation-conversion-rate` without a mode is also
    /// rejected.
    #[test]
    fn test_methylation_conversion_rate_without_mode_rejected() {
        let mut sim = make_default_simulate();
        sim.methylation_mode = None;
        sim.methylation_conversion_rate = Some(0.5);
        let err =
            sim.validate().expect_err("validate must reject orphaned methylation conversion rate");
        let msg = format!("{err}");
        assert!(
            msg.contains("--methylation-conversion-rate requires --methylation-mode"),
            "error must mention required mode flag, got: {msg}"
        );
    }

    /// Default rates with no mode is the normal "no methylation" path and
    /// must not be rejected.
    #[test]
    fn test_default_rates_without_mode_accepted() {
        let sim = make_default_simulate();
        assert!(sim.validate().is_ok(), "defaults without mode must validate cleanly");
    }

    /// Explicit `--methylation-conversion-rate 1.0` without a mode is still a
    /// user error: the validator detects *flag presence*, not just numeric
    /// override.
    #[test]
    fn test_methylation_conversion_rate_explicit_one_without_mode_rejected() {
        let mut sim = make_default_simulate();
        sim.methylation_mode = None;
        sim.methylation_conversion_rate = Some(1.0);
        let err = sim.validate().expect_err(
            "validate must reject explicit --methylation-conversion-rate 1.0 without mode",
        );
        let msg = format!("{err}");
        assert!(
            msg.contains("--methylation-conversion-rate requires --methylation-mode"),
            "error must mention required mode flag, got: {msg}"
        );
    }

    /// `--methylation-failure-rate` without a mode is rejected, mirroring the
    /// conversion-rate orphan check.
    #[test]
    fn test_methylation_failure_rate_without_mode_rejected() {
        let mut sim = make_default_simulate();
        sim.methylation_mode = None;
        sim.methylation_failure_rate = Some(0.05);
        let err =
            sim.validate().expect_err("validate must reject orphaned methylation failure rate");
        let msg = format!("{err}");
        assert!(
            msg.contains("--methylation-failure-rate requires --methylation-mode"),
            "error must mention required mode flag, got: {msg}"
        );
    }

    /// `--methylation-failure-rate` outside `[0, 1]` is rejected even with a
    /// mode set.
    #[test]
    fn test_methylation_failure_rate_out_of_range_rejected() {
        let mut sim = make_default_simulate();
        sim.methylation_mode = Some(crate::meth::MethylationMode::EmSeq);
        sim.methylation_failure_rate = Some(1.5);
        let err = sim
            .validate()
            .expect_err("validate must reject --methylation-failure-rate outside [0, 1]");
        let msg = format!("{err}");
        assert!(
            msg.contains("--methylation-failure-rate must be in [0.0, 1.0]"),
            "error must mention the valid range, got: {msg}"
        );
    }

    /// `--methylation-conversion-rate` outside `[0, 1]` must be rejected with
    /// the mode set (mirrors the failure-rate guard, which is tested above).
    #[test]
    fn test_methylation_conversion_rate_out_of_range_rejected() {
        let mut sim = make_default_simulate();
        sim.methylation_mode = Some(crate::meth::MethylationMode::EmSeq);
        sim.methylation_conversion_rate = Some(1.5);
        let err = sim
            .validate()
            .expect_err("validate must reject --methylation-conversion-rate outside [0, 1]");
        let msg = format!("{err}");
        assert!(
            msg.contains("--methylation-conversion-rate must be in [0.0, 1.0]"),
            "error must mention the valid range, got: {msg}"
        );
    }

    /// A non-finite (NaN) `--methylation-conversion-rate` must also be rejected
    /// by the same range guard — `(0.0..=1.0).contains(&NaN)` is `false`.
    #[test]
    fn test_methylation_conversion_rate_nan_rejected() {
        let mut sim = make_default_simulate();
        sim.methylation_mode = Some(crate::meth::MethylationMode::EmSeq);
        sim.methylation_conversion_rate = Some(f64::NAN);
        let err =
            sim.validate().expect_err("validate must reject a NaN --methylation-conversion-rate");
        assert!(
            format!("{err}").contains("--methylation-conversion-rate must be in [0.0, 1.0]"),
            "error must mention the valid range"
        );
    }

    /// A `--cpg-truth-bedgraph` path whose parent directory does not exist
    /// must be rejected at validate time, not several minutes into the
    /// simulation when `write_bedgraph` tries to create the file. Mirrors
    /// the existing fail-fast check on `--output`.
    #[test]
    fn test_cpg_truth_bedgraph_missing_parent_rejected() {
        let mut sim = make_default_simulate();
        sim.methylation_mode = Some(crate::meth::MethylationMode::EmSeq);
        sim.cpg_truth_bedgraph =
            Some(PathBuf::from("/nonexistent-path-for-holodeck-test/sub/cpg.bg"));
        let err = sim
            .validate()
            .expect_err("validate must reject --cpg-truth-bedgraph with a non-existent parent");
        let msg = format!("{err}");
        assert!(
            msg.contains("CpG truth bedGraph directory does not exist"),
            "error must mention missing CpG truth bedGraph directory, got: {msg}"
        );
    }
}
