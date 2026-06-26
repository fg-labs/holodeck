//! Alignment accuracy evaluation command.
//!
//! `holodeck eval` scores an aligner's BAM against holodeck's own truth.
//! Placement accuracy (true vs mapped position, by MAPQ bin) is always
//! reported from encoded read names. Optional truth inputs unlock further
//! metrics, each written to its own TSV alongside `<prefix>.eval.txt`:
//! - [`placement`] — placement accuracy (always; `<prefix>.eval.txt`).
//! - [`variants`] — variant-representation accuracy (`--variants` + `--truth`;
//!   `<prefix>.variants.tsv`).
//! - [`meth`] — methylation-level correlation (`--cpg-truth`;
//!   `<prefix>.meth.tsv`).

mod cigar;
mod golden;
mod meth;
mod placement;
mod variants;

use std::path::PathBuf;

use anyhow::{Result, bail};
use clap::Parser;

use super::command::Command;
use super::common::OutputPrefixOptions;

/// Evaluate alignment accuracy of simulated reads.
///
/// Always reports placement accuracy (true vs mapped position, by MAPQ bin)
/// from encoded read names. Given the truth VCF and golden BAM that `simulate`
/// emits, `--variants` additionally reports how faithfully aligned reads
/// represent the simulated substitutions, with `--meth` breaking the results
/// down by bisulfite substitution class.
#[derive(Parser, Debug)]
#[command(after_long_help = "EXAMPLES:\n  \
    holodeck eval --mapped aligned.bam -o eval_results\n  \
    holodeck eval --mapped aligned.bam --truth golden.bam \\\n    \
    --variants truth.vcf --meth -o eval_results")]
pub struct Eval {
    /// BAM file of mapped reads to evaluate.
    #[arg(short = 'm', long, value_name = "BAM")]
    pub mapped: PathBuf,

    /// Golden BAM (`simulate --golden-bam`) supplying each read's true span,
    /// haplotype, and MD/NM tags. Required by `--variants`.
    #[arg(long, value_name = "BAM")]
    pub truth: Option<PathBuf>,

    /// Truth VCF (`mutate` / `methylate`) of simulated variants. Enables
    /// variant-representation scoring; requires `--truth`.
    #[arg(long, value_name = "VCF")]
    pub variants: Option<PathBuf>,

    /// Sample name to resolve genotypes for in the truth VCF (defaults to the
    /// first sample).
    #[arg(long, value_name = "NAME")]
    pub sample: Option<String>,

    /// Break `--variants` results down by bisulfite substitution class
    /// (conversion, mirror, transversion, other). The bisulfite/EM-seq context
    /// comes from the golden BAM's true conversion strand, not from this flag.
    #[arg(long)]
    pub meth: bool,

    /// Per-CpG truth bedGraph (`simulate --cpg-truth-bedgraph`). Enables
    /// methylation-level correlation against the aligner's `XM` calls.
    #[arg(long, value_name = "BEDGRAPH")]
    pub cpg_truth: Option<PathBuf>,

    #[command(flatten)]
    pub output: OutputPrefixOptions,

    /// Maximum distance (in bases) between the true and mapped start positions
    /// of a read for it to be considered correctly mapped. Uses
    /// `|mapped_start - true_start| <= wiggle` on the same contig.
    #[arg(long, default_value_t = 5, value_name = "INT")]
    pub wiggle: u32,
}

impl Command for Eval {
    fn execute(&self) -> Result<()> {
        if self.variants.is_some() && self.truth.is_none() {
            bail!("--variants requires --truth (the golden BAM provides per-read truth spans)");
        }
        if self.meth && self.variants.is_none() {
            log::warn!("--meth has no effect without --variants");
        }

        placement::run(&self.mapped, &self.output.output, self.wiggle)?;

        if let Some(vcf) = &self.variants {
            // Safe: the guard above rejects --variants without --truth.
            let golden_path = self.truth.as_ref().expect("--variants requires --truth");
            let golden = golden::load(golden_path)?;
            let truth = variants::VariantTruth::from_vcf(vcf, self.sample.as_deref())?;
            variants::run(&self.mapped, &golden, &truth, self.meth, &self.output.output)?;
        } else if self.truth.is_some() {
            log::warn!("--truth is only used with --variants; placement uses encoded read names");
        }

        if let Some(cpg_truth) = &self.cpg_truth {
            meth::run(&self.mapped, cpg_truth, &self.output.output)?;
        }

        Ok(())
    }
}
