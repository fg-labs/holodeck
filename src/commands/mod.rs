//! CLI subcommand implementations for the `holodeck` binary.
//!
//! Each subcommand lives in its own module and implements the [`command::Command`]
//! trait. [`common`] holds the shared CLI option groups (`--reference`, output,
//! VCF, BED, seed) flattened into every subcommand. `simulate` drives the
//! read-simulation pipeline; `mutate` produces random variant VCFs; `eval`
//! scores alignment accuracy; `methylate` builds the methylation-truth VCF
//! consumed by `simulate`'s bisulfite / TAPS chemistry.

pub mod command;
pub mod common;
pub mod eval;
/// `methylate` subcommand — generate per-haplotype CpG methylation truth
/// as an MT/MB-annotated VCF, optionally writing a population-fraction bedGraph.
pub mod methylate;
pub mod mutate;
pub mod simulate;
