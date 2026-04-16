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
use crate::seed::resolve_seed;
use crate::sequence_dict::SequenceDictionary;
use crate::version::VERSION;

/// Default sample name used in the golden BAM `@RG SM` field when the
/// simulation is not driven by a VCF sample.
const DEFAULT_SAMPLE_NAME: &str = "holodeck-simulation";

/// Default Illumina TruSeq adapter sequence for read 1.
const DEFAULT_ADAPTER_R1: &str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";

/// Default Illumina TruSeq adapter sequence for read 2.
const DEFAULT_ADAPTER_R2: &str = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

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

    /// Write a ground-truth BAM file with correct alignments.
    #[arg(long)]
    pub golden_bam: bool,

    /// Write a ground-truth VCF annotated with simulated coverage.
    #[arg(long)]
    pub golden_vcf: bool,

    /// Use simple read names (`holodeck:N`) instead of encoding truth
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
        self.validate()?;
        self.run_simulation()
    }
}

impl Simulate {
    /// Validate command-line arguments before running.
    fn validate(&self) -> Result<()> {
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
        if self.compression > 12 {
            bail!("--compression must be between 0 and 12");
        }

        // --sample without --vcf is nonsensical.
        if self.vcf.sample.is_some() && self.vcf.vcf.is_none() {
            bail!("--sample requires --vcf");
        }

        // Validate VCF sample configuration upfront so the user gets a clear
        // error before the simulation loop starts.
        if let Some(vcf_path) = &self.vcf.vcf {
            crate::vcf::validate_vcf_sample(vcf_path, self.vcf.sample.as_deref())?;
        }

        // Validate output parent directory exists.
        if let Some(parent) = self.output.output.parent()
            && !parent.as_os_str().is_empty()
            && !parent.exists()
        {
            bail!("Output directory does not exist: {}", parent.display());
        }

        Ok(())
    }

    /// Run the main simulation pipeline.
    fn run_simulation(&self) -> Result<()> {
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
        let frag_dist = Normal::new(self.fragment_mean as f64, self.fragment_stddev as f64)
            .map_err(|e| anyhow::anyhow!("Invalid fragment distribution parameters: {e}"))?;

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
            let meta = self.golden_bam_metadata()?;
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

        let mut read_num: u64 = 0;
        let contig_names: Vec<String> = dict.names().into_iter().map(String::from).collect();

        for contig_name in &contig_names {
            read_num += self.simulate_contig(
                contig_name,
                &dict,
                &mut fasta,
                targets.as_ref(),
                total_reads,
                &error_model,
                &frag_dist,
                &mut r1_writer,
                &mut r2_writer,
                &mut golden_bam_writer,
                read_num,
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
    /// command line is captured verbatim from `std::env::args`.  The sample
    /// name is the resolved VCF sample when `--vcf` is given (either the
    /// value of `--sample`, or the sole sample in a single-sample VCF), and
    /// [`DEFAULT_SAMPLE_NAME`] otherwise.
    fn golden_bam_metadata(&self) -> Result<GoldenBamMetadata> {
        let command_line = std::env::args().collect::<Vec<_>>().join(" ");
        let sample = if let Some(vcf_path) = &self.vcf.vcf {
            crate::vcf::validate_vcf_sample(vcf_path, self.vcf.sample.as_deref())?
        } else {
            DEFAULT_SAMPLE_NAME.to_string()
        };
        Ok(GoldenBamMetadata { command_line, version: VERSION.clone(), sample })
    }

    /// Compute the deterministic seed from simulation parameters.
    fn compute_seed(&self) -> u64 {
        let seed_desc = format!(
            "{}:{}:{}:{}:{}:{}:{}",
            self.reference.reference.display(),
            self.coverage,
            self.read_length,
            self.fragment_mean,
            self.fragment_stddev,
            self.min_error_rate,
            self.max_error_rate,
        );
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
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn simulate_contig(
        &self,
        contig_name: &str,
        dict: &SequenceDictionary,
        fasta: &mut Fasta,
        targets: Option<&TargetRegions>,
        total_reads: u64,
        error_model: &IlluminaErrorModel,
        frag_dist: &Normal<f64>,
        r1_writer: &mut FastqWriter,
        r2_writer: &mut Option<FastqWriter>,
        golden_bam: &mut Option<GoldenBamWriter>,
        start_read_num: u64,
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

        let reference = fasta.load_contig(contig_name)?;

        // Load variants if VCF provided.
        let variants = if let Some(vcf_path) = &self.vcf.vcf {
            crate::vcf::load_variants_for_contig(
                vcf_path,
                contig_name,
                self.vcf.sample.as_deref(),
                dict,
            )?
        } else {
            Vec::new()
        };

        if !variants.is_empty() {
            log::info!("  Loaded {} variants for {contig_name}", variants.len());
        }

        let max_ploidy = variants.iter().map(|v| v.genotype.ploidy()).max().unwrap_or(2);
        let haplotypes = build_haplotypes(&variants, max_ploidy, rng);

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
            let pair = generate_read_pair(
                &fragment,
                contig_name,
                read_num,
                self.read_length,
                !self.single_end,
                self.adapter_r1.as_bytes(),
                self.adapter_r2.as_bytes(),
                error_model,
                self.simple_names,
                rng,
            );

            r1_writer.write_read(&pair.read1)?;
            if let Some(w) = r2_writer
                && let Some(r2) = &pair.read2
            {
                w.write_read(r2)?;
            }
            if let Some(bam_w) = golden_bam {
                bam_w.write_pair(&pair)?;
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
