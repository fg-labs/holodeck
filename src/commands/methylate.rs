//! `methylate` subcommand — given a reference (and optional variant VCF),
//! produce a methylation-annotated VCF that records per-haplotype
//! per-strand methylation state at every CpG.

use std::fs::File;
use std::io::{BufWriter, Write as _};
use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;
use rand::SeedableRng as _;

use super::command::Command;
use super::common::{ReferenceOptions, SeedOptions, VcfOptions};

/// Generate a methylation-annotated VCF.
///
/// Scans every CpG dinucleotide in the reference (optionally after applying
/// variants from a VCF) and records per-haplotype, per-strand methylation
/// state for each site. Methylation is assigned with a context-aware Markov
/// model: CpGs are classified into island / shore / open-sea (detected de
/// novo from the sequence), each with its own target rate and spatial
/// correlation length, so islands come out hypomethylated, open-sea
/// hypermethylated, with autocorrelated runs in between. Methylation is
/// symmetric by default with a low sporadic hemimethylation rate. Per-haplotype
/// draws give allele-specific methylation. All draws are deterministic from
/// `--seed`.
#[derive(Parser, Debug)]
#[command(after_long_help = "EXAMPLES:\n  \
    holodeck methylate -r ref.fa -o meth.vcf.gz\n  \
    holodeck methylate -r ref.fa -v variants.vcf -o meth.vcf.gz --seed 42\n  \
    holodeck methylate -r ref.fa -o meth.vcf.gz --methylation-rate-open-sea 0.9")]
pub struct Methylate {
    #[command(flatten)]
    pub reference: ReferenceOptions,

    #[command(flatten)]
    pub vcf: VcfOptions,

    #[command(flatten)]
    pub seed: SeedOptions,

    /// Target methylation fraction for CpG-island-interior CpGs. Islands are
    /// characteristically hypomethylated, so this defaults low. Must be in
    /// `[0.0, 1.0]`. To methylate uniformly (no island structure), set the
    /// three context rates equal.
    #[arg(long, default_value_t = crate::meth::DEFAULT_ISLAND_RATE, value_name = "FLOAT")]
    pub methylation_rate_island: f64,

    /// Target methylation fraction for island-shore CpGs (within 2 kb of an
    /// island; intermediate methylation). Must be in `[0.0, 1.0]`.
    #[arg(long, default_value_t = crate::meth::DEFAULT_SHORE_RATE, value_name = "FLOAT")]
    pub methylation_rate_shore: f64,

    /// Target methylation fraction for open-sea CpGs (the hypermethylated bulk
    /// of the genome, away from islands). Must be in `[0.0, 1.0]`.
    #[arg(long, default_value_t = crate::meth::DEFAULT_OPEN_SEA_RATE, value_name = "FLOAT")]
    pub methylation_rate_open_sea: f64,

    /// Spatial correlation length (bp) for island-interior CpGs: larger values
    /// give longer runs of like-methylated CpGs. Must be `> 0`.
    #[arg(long, default_value_t = crate::meth::DEFAULT_CORRELATION_LENGTH_BP, value_name = "BP")]
    pub methylation_correlation_length_island: f64,

    /// Spatial correlation length (bp) for shore CpGs. Must be `> 0`.
    #[arg(long, default_value_t = crate::meth::DEFAULT_CORRELATION_LENGTH_BP, value_name = "BP")]
    pub methylation_correlation_length_shore: f64,

    /// Spatial correlation length (bp) for open-sea CpGs. Must be `> 0`.
    #[arg(long, default_value_t = crate::meth::DEFAULT_CORRELATION_LENGTH_BP, value_name = "BP")]
    pub methylation_correlation_length_open_sea: f64,

    /// Probability that a methylated CpG is made hemimethylated — exactly one
    /// strand left unmethylated. Real hemimethylation is sporadic and rare, so
    /// this defaults low. Must be in `[0.0, 1.0]`.
    #[arg(long, default_value_t = crate::meth::DEFAULT_HEMI_RATE, value_name = "FLOAT")]
    pub hemimethylation_rate: f64,

    /// Output methylation-annotated VCF (BGZF-compressed).
    #[arg(long, short = 'o', value_name = "PATH")]
    pub output: PathBuf,

    /// Optional MethylDackel-format population-fraction bedGraph derived
    /// closed-form from the genome model. Independent of any read coverage.
    #[arg(long, value_name = "PATH")]
    pub bedgraph: Option<PathBuf>,
}

impl Methylate {
    /// Assemble the per-context methylation model from the CLI flags.
    fn methylation_model(&self) -> crate::meth::MethylationModel {
        use crate::meth::{ContextParams, MethylationModel};
        MethylationModel {
            island: ContextParams {
                rate: self.methylation_rate_island,
                correlation_length_bp: self.methylation_correlation_length_island,
            },
            shore: ContextParams {
                rate: self.methylation_rate_shore,
                correlation_length_bp: self.methylation_correlation_length_shore,
            },
            open_sea: ContextParams {
                rate: self.methylation_rate_open_sea,
                correlation_length_bp: self.methylation_correlation_length_open_sea,
            },
            hemi_rate: self.hemimethylation_rate,
        }
    }
}

impl Command for Methylate {
    fn execute(&self) -> Result<()> {
        use noodles_bgzf as bgzf;

        use crate::fasta::Fasta;
        use crate::haplotype::build_haplotypes;
        use crate::meth::ContigMethylation;
        use crate::seed::{derive_seed, resolve_seed};
        use crate::vcf::methylation::write_contig;
        use crate::vcf::methylation::write_vcf_header;
        use crate::version::VERSION;

        // 1. Build and validate the methylation model from the per-context flags.
        let model = self.methylation_model();
        model.validate()?;

        // 2. Resolve seed deterministically from args if not explicit. Fold
        //    every model parameter into the description so default-seed runs
        //    with different parameters get distinct streams.
        let seed_desc = format!(
            "{}:methylate:{}:{}:{}:{}:{}:{}:{}",
            self.reference.reference.display(),
            self.methylation_rate_island,
            self.methylation_rate_shore,
            self.methylation_rate_open_sea,
            self.methylation_correlation_length_island,
            self.methylation_correlation_length_shore,
            self.methylation_correlation_length_open_sea,
            self.hemimethylation_rate,
        );
        let seed = resolve_seed(self.seed.seed, &seed_desc);
        log::info!("Using random seed: {seed}");

        // 3. Open the reference and collect contig names.
        let mut fasta = Fasta::from_path(&self.reference.reference)?;
        let dict = fasta.dict().clone();
        let contig_names: Vec<String> = dict.names().into_iter().map(String::from).collect();

        // 4. Resolve the sample for any VCF, before per-contig variant loading.
        // Reject --sample without --vcf so a typo or missing flag fails fast
        // instead of running a reference-only job with the input silently
        // ignored. Mirrors the simulate guard.
        if self.vcf.sample.is_some() && self.vcf.vcf.is_none() {
            anyhow::bail!("--sample requires --vcf");
        }
        let resolved_sample = if let Some(vcf_path) = &self.vcf.vcf {
            Some(crate::vcf::validate_vcf_sample(vcf_path, self.vcf.sample.as_deref())?)
        } else {
            None
        };

        // 5. Open the output VCF, BGZF-compressed.
        let file = File::create(&self.output)?;
        let mut bgzf = bgzf::io::Writer::new(file);

        // 6. Write the VCF header.
        let cmd_line = capture_command_line();
        write_vcf_header(&mut bgzf, &dict, resolved_sample.as_deref(), &VERSION, &cmd_line)?;

        // 6b. Open the bedGraph output file and write the track header, if requested.
        let mut bedgraph_writer: Option<BufWriter<File>> =
            if let Some(bedgraph_path) = &self.bedgraph {
                let bg_file = File::create(bedgraph_path)?;
                let mut w = BufWriter::new(bg_file);
                crate::output::methylation_bedgraph::write_bedgraph_header(&mut w)?;
                Some(w)
            } else {
                None
            };

        // 6c. Parse the entire VCF once, partitioning variants by contig, so
        // the per-contig loop below performs O(1) lookups instead of
        // re-reading the full VCF per contig (O(contigs × records)). The scan
        // also resolves the sample's ploidy from unfiltered GTs and reuses it
        // everywhere — variant-bearing contigs and variant-free ones alike.
        // A per-contig `unwrap_or(2)` would give haploid / triploid samples
        // the wrong MT/MB shape on contigs with no ALT calls.
        let parsed_variants = if let Some(vcf_path) = &self.vcf.vcf {
            Some(crate::vcf::parse_variants_by_contig(vcf_path, resolved_sample.as_deref(), &dict)?)
        } else {
            None
        };
        let sample_ploidy = parsed_variants.as_ref().map_or(2, |p| p.sample_ploidy);

        // 7. For each contig, build haplotypes, compute methylation, emit rows.
        for contig_name in &contig_names {
            // Per-contig deterministic seed for reference normalization.
            // Use the bare contig name (no prefix) to match `simulate`'s
            // `derive_seed(main_seed, contig_name)` at simulate.rs line ~578,
            // so both commands resolve IUPAC codes identically for the same
            // `--seed` and reference.
            let ref_seed = derive_seed(seed, contig_name);
            let mut ref_rng = rand::rngs::SmallRng::seed_from_u64(ref_seed);
            let reference = fasta.load_contig(contig_name, &mut ref_rng)?;

            // Look this contig's variants up from the once-parsed map. An
            // absent key or `None` map (no VCF) both yield an empty slice.
            let variants: &[crate::vcf::genotype::VariantRecord] = parsed_variants
                .as_ref()
                .and_then(|p| p.by_contig.get(contig_name))
                .map_or(&[], Vec::as_slice);

            // Build haplotypes with their own deterministic sub-seed, sized
            // by the whole-VCF `sample_ploidy` so variant-free contigs land
            // the same number of haplotypes as variant-bearing ones.
            let hap_seed = derive_seed(seed, &format!("haps@{contig_name}"));
            let mut hap_rng = rand::rngs::SmallRng::seed_from_u64(hap_seed);
            let haplotypes = build_haplotypes(variants, sample_ploidy, &mut hap_rng);

            // Draw per-haplotype methylation with another deterministic sub-seed.
            let meth_seed = derive_seed(seed, &format!("meth@{contig_name}"));
            let mut meth_rng = rand::rngs::SmallRng::seed_from_u64(meth_seed);
            let methylation =
                ContigMethylation::from_haplotypes(&haplotypes, &reference, &model, &mut meth_rng);

            write_contig(
                &mut bgzf,
                contig_name,
                &reference,
                variants,
                &methylation,
                sample_ploidy,
            )?;

            // Append per-contig bedGraph records, if a bedGraph was requested.
            // Pass `&haplotypes` so the writer can map every reference CpG
            // through each haplotype's local coordinate system before reading
            // the methylation bitmap; variant-shifted and variant-destroyed
            // CpGs would otherwise be miscounted.
            if let Some(bg) = bedgraph_writer.as_mut() {
                crate::output::methylation_bedgraph::write_bedgraph_records(
                    bg,
                    contig_name,
                    &reference,
                    &haplotypes,
                    &methylation,
                )?;
            }
        }

        // 8. Flush and close the VCF.
        bgzf.flush()?;
        drop(bgzf);

        // 8b. Flush and close the bedGraph, if open.
        if let Some(mut bg) = bedgraph_writer {
            bg.flush()?;
            if let Some(path) = &self.bedgraph {
                log::info!("Wrote population-fraction bedGraph to {}", path.display());
            }
        }

        log::info!("Wrote methylation VCF to {}", self.output.display());
        Ok(())
    }
}

/// Capture the current process's command-line arguments as a space-joined
/// string, using lossy UTF-8 for non-Unicode arguments.
fn capture_command_line() -> String {
    std::env::args_os().map(|arg| arg.to_string_lossy().into_owned()).collect::<Vec<_>>().join(" ")
}
