//! Alignment accuracy evaluation command.

use std::collections::BTreeMap;
use std::io::Write;
use std::path::PathBuf;

use anyhow::{Context, Result};
use bstr::ByteSlice;
use clap::Parser;
use noodles::bam;

use super::command::{Command, output_path};
use super::common::OutputPrefixOptions;
use crate::read_naming::{parse_encoded_pe_name, parse_encoded_se_name};

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

/// Accuracy counts for a single MAPQ bin.
#[derive(Debug, Default, Clone)]
struct BinCounts {
    /// Reads mapped to the correct position (within wiggle).
    correct: u64,
    /// Reads mapped to a wrong position or wrong contig.
    mismapped: u64,
    /// Reads that are unmapped.
    unmapped: u64,
    /// Total reads in this bin.
    total: u64,
}

impl Command for Eval {
    fn execute(&self) -> Result<()> {
        if self.truth.is_some() {
            log::warn!("--truth (golden BAM) is not yet implemented; using read names");
        }

        let mut reader = bam::io::reader::Builder
            .build_from_path(&self.mapped)
            .with_context(|| format!("Failed to open BAM: {}", self.mapped.display()))?;

        let header = reader.read_header()?;

        // MAPQ bins: 0, 1-9, 10-19, 20-29, 30-39, 40-49, 50-59, 60+
        let mut bins: BTreeMap<u8, BinCounts> = BTreeMap::new();
        let mut total_reads: u64 = 0;
        let mut parse_failures: u64 = 0;

        for result in reader.records() {
            let record = result.with_context(|| "Failed to read BAM record")?;

            // Skip secondary and supplementary alignments before counting.
            let flags = record.flags();
            if flags.is_secondary() || flags.is_supplementary() {
                continue;
            }
            total_reads += 1;

            // Get read name.
            let name_bytes = record.name().map_or(&b""[..], |n| n.as_bytes());
            let name = name_bytes.to_str().unwrap_or("");

            // Parse truth from encoded read name.
            let truth = parse_encoded_pe_name(name)
                .map(|(_, r1, _)| r1)
                .or_else(|| parse_encoded_se_name(name).map(|(_, r1)| r1));

            let Some(truth) = truth else {
                parse_failures += 1;
                continue;
            };

            let mapq = record.mapping_quality().map_or(0, u8::from);
            let bin_key = mapq_bin(mapq);
            let counts = bins.entry(bin_key).or_default();
            counts.total += 1;

            if flags.is_unmapped() {
                counts.unmapped += 1;
                continue;
            }

            // Get mapped position.
            let mapped_contig_idx = record.reference_sequence_id().and_then(Result::ok);
            let mapped_pos = record.alignment_start().and_then(Result::ok).map(usize::from);

            let (Some(mapped_contig_idx), Some(mapped_pos_1based)) =
                (mapped_contig_idx, mapped_pos)
            else {
                counts.mismapped += 1;
                continue;
            };

            // Resolve mapped contig name.
            let Some((contig_name, _)) = header.reference_sequences().get_index(mapped_contig_idx)
            else {
                counts.mismapped += 1;
                continue;
            };
            let mapped_contig = String::from_utf8_lossy(contig_name.as_ref());

            // Compare truth vs mapped.
            #[expect(clippy::cast_possible_truncation, reason = "mapped positions fit u32")]
            let mapped_pos_u32 = mapped_pos_1based as u32;
            let is_correct = mapped_contig == truth.contig
                && mapped_pos_u32.abs_diff(truth.position) <= self.wiggle;

            if is_correct {
                counts.correct += 1;
            } else {
                counts.mismapped += 1;
            }
        }

        if parse_failures > 0 {
            log::warn!("{parse_failures} reads had unparseable names; skipped");
        }

        // Write results.
        let output_file = output_path(&self.output.output, ".eval.txt");
        let mut out = std::fs::File::create(&output_file)
            .with_context(|| format!("Failed to create {}", output_file.display()))?;

        writeln!(
            out,
            "mapq_bin\ttotal\tcorrect\tmismapped\tunmapped\tpct_correct\tpct_mismapped\tpct_unmapped"
        )?;

        let mut grand_total = BinCounts::default();
        for (&bin, counts) in &bins {
            write_bin_row(&mut out, &format_bin_label(bin), counts)?;
            grand_total.correct += counts.correct;
            grand_total.mismapped += counts.mismapped;
            grand_total.unmapped += counts.unmapped;
            grand_total.total += counts.total;
        }
        write_bin_row(&mut out, "ALL", &grand_total)?;

        log::info!(
            "Evaluated {total_reads} reads: {} correct, {} mismapped, {} unmapped",
            grand_total.correct,
            grand_total.mismapped,
            grand_total.unmapped
        );
        log::info!("Results written to: {}", output_file.display());

        Ok(())
    }
}

/// Map a MAPQ value to a bin key.
fn mapq_bin(mapq: u8) -> u8 {
    match mapq {
        0 => 0,
        1..=9 => 1,
        10..=19 => 10,
        20..=29 => 20,
        30..=39 => 30,
        40..=49 => 40,
        50..=59 => 50,
        _ => 60,
    }
}

/// Format a bin key as a label string.
fn format_bin_label(bin: u8) -> String {
    match bin {
        0 => "0".to_string(),
        1 => "1-9".to_string(),
        60 => "60+".to_string(),
        _ => format!("{}-{}", bin, bin + 9),
    }
}

/// Write one row of the evaluation results table.
fn write_bin_row(out: &mut impl Write, label: &str, counts: &BinCounts) -> Result<()> {
    let total = counts.total.max(1) as f64;
    writeln!(
        out,
        "{label}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}",
        counts.total,
        counts.correct,
        counts.mismapped,
        counts.unmapped,
        counts.correct as f64 / total * 100.0,
        counts.mismapped as f64 / total * 100.0,
        counts.unmapped as f64 / total * 100.0,
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mapq_bin() {
        assert_eq!(mapq_bin(0), 0);
        assert_eq!(mapq_bin(5), 1);
        assert_eq!(mapq_bin(10), 10);
        assert_eq!(mapq_bin(15), 10);
        assert_eq!(mapq_bin(30), 30);
        assert_eq!(mapq_bin(60), 60);
        assert_eq!(mapq_bin(255), 60);
    }

    #[test]
    fn test_format_bin_label() {
        assert_eq!(format_bin_label(0), "0");
        assert_eq!(format_bin_label(1), "1-9");
        assert_eq!(format_bin_label(10), "10-19");
        assert_eq!(format_bin_label(60), "60+");
    }

    #[test]
    fn test_write_bin_row() {
        let counts = BinCounts { correct: 90, mismapped: 8, unmapped: 2, total: 100 };
        let mut buf = Vec::new();
        write_bin_row(&mut buf, "30-39", &counts).unwrap();
        let line = String::from_utf8(buf).unwrap();
        assert!(line.starts_with("30-39\t100\t90\t8\t2\t"));
        assert!(line.contains("90.00"));
    }
}
