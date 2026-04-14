//! Output writers for simulated reads and ground truth files.
//!
//! Provides FASTQ output (always) and optional golden BAM output for
//! benchmarking alignment pipelines.  Supports single-threaded and
//! multi-threaded BGZF compression via pooled-writer.

pub mod fastq;
pub mod golden_bam;
