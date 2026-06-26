//! Alignment accuracy evaluation command.
//!
//! `holodeck eval` scores an aligner's BAM against holodeck's own truth.
//! Placement accuracy (true vs mapped position, by MAPQ bin) is always
//! reported from encoded read names. Optional truth inputs unlock further
//! metrics, each written to its own TSV alongside `<prefix>.eval.txt`:
//! - [`placement`] — placement accuracy (always).

mod placement;

use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;

use super::command::Command;
use super::common::OutputPrefixOptions;

/// Evaluate alignment accuracy of simulated reads.
///
/// Compares the true (simulated) positions of reads against their mapped
/// positions in a BAM file.  Reports mapping accuracy, mismapping rate, and
/// unmapped rate stratified by MAPQ bin.  Truth positions are parsed from
/// encoded read names (default holodeck format).
#[derive(Parser, Debug)]
#[command(after_long_help = "EXAMPLES:\n  \
    holodeck eval --mapped aligned.bam -o eval_results\n  \
    holodeck eval --mapped aligned.bam --truth golden.bam -o eval_results")]
pub struct Eval {
    /// BAM file of mapped reads to evaluate.
    #[arg(short = 'm', long, value_name = "BAM")]
    pub mapped: PathBuf,

    /// Optional golden BAM file with truth alignments. If omitted, truth
    /// positions are parsed from encoded read names.
    #[arg(long, value_name = "BAM")]
    pub truth: Option<PathBuf>,

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
        if self.truth.is_some() {
            log::warn!("--truth (golden BAM) is not yet implemented; using read names");
        }

        placement::run(&self.mapped, &self.output.output, self.wiggle)
    }
}
