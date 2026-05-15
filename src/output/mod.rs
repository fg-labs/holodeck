//! Output writers for simulated reads and ground truth files.
//!
//! Provides FASTQ output (always) and optional golden BAM output for
//! benchmarking alignment pipelines.  Supports single-threaded and
//! multi-threaded BGZF compression via pooled-writer.

/// Coverage-weighted per-CpG truth tally written as a MethylDackel-format
/// bedGraph from simulated reads.
pub mod cpg_truth;
pub mod fastq;
pub mod golden_bam;
/// Closed-form per-CpG population-fraction bedGraph writer derived from a
/// per-haplotype methylation bitmap; no simulated reads required.
pub mod methylation_bedgraph;
