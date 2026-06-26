//! Placement accuracy: true vs mapped position, stratified by MAPQ bin.
//!
//! Truth positions are parsed from encoded holodeck read names. For each
//! primary mapped record, the mapped start is compared to the true start on
//! the same contig within a wiggle tolerance; reads are tallied as correct,
//! mismapped, or unmapped within their MAPQ bin.

use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use bstr::ByteSlice;
use noodles::bam;

use crate::commands::command::output_path;
use crate::read_naming::{parse_encoded_pe_name, parse_encoded_se_name};

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

/// Evaluate placement accuracy of `mapped` and write `<prefix>.eval.txt`.
///
/// # Errors
/// Returns an error if the BAM cannot be read or the output cannot be written.
pub fn run(mapped: &Path, output_prefix: &Path, wiggle: u32) -> Result<()> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(mapped)
        .with_context(|| format!("Failed to open BAM: {}", mapped.display()))?;

    let header = reader.read_header()?;

    // MAPQ bins: 0, 1-9, 10-19, 20-29, 30-39, 40-49, 50-59, 60+, NA (255)
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

        // Parse truth from encoded read name. For PE names, pick R1 or R2
        // based on the record's segment flag; mis-selecting here caused R2
        // alignments to be scored against the R1 truth position.
        let truth = if let Some((_, r1, r2)) = parse_encoded_pe_name(name) {
            if flags.is_last_segment() { Some(r2) } else { Some(r1) }
        } else {
            parse_encoded_se_name(name).map(|(_, truth)| truth)
        };

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

        let (Some(mapped_contig_idx), Some(mapped_pos_1based)) = (mapped_contig_idx, mapped_pos)
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
        let is_correct =
            mapped_contig == truth.contig && mapped_pos_u32.abs_diff(truth.position) <= wiggle;

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
    let output_file = output_path(output_prefix, ".eval.txt");
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
    log::info!("Placement results written to: {}", output_file.display());

    Ok(())
}

/// Map a MAPQ value to a bin key. MAPQ `255` is reserved by the SAM spec for
/// "mapping quality unavailable" and is binned separately so unknown-quality
/// reads are not counted with the highest-confidence `60+` alignments.
fn mapq_bin(mapq: u8) -> u8 {
    match mapq {
        0 => 0,
        1..=9 => 1,
        10..=19 => 10,
        20..=29 => 20,
        30..=39 => 30,
        40..=49 => 40,
        50..=59 => 50,
        255 => 255,
        _ => 60,
    }
}

/// Format a bin key as a label string. The `255` key is the SAM "unavailable"
/// MAPQ and is labelled `NA`.
fn format_bin_label(bin: u8) -> String {
    match bin {
        0 => "0".to_string(),
        1 => "1-9".to_string(),
        60 => "60+".to_string(),
        255 => "NA".to_string(),
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
        // 61..=254 are valid high-confidence MAPQs and stay in the 60+ bin.
        assert_eq!(mapq_bin(254), 60);
        // 255 is reserved by the SAM spec for "mapping quality unavailable" and
        // must not be counted as the highest-confidence bin.
        assert_eq!(mapq_bin(255), 255);
    }

    #[test]
    fn test_format_bin_label() {
        assert_eq!(format_bin_label(0), "0");
        assert_eq!(format_bin_label(1), "1-9");
        assert_eq!(format_bin_label(10), "10-19");
        assert_eq!(format_bin_label(60), "60+");
        assert_eq!(format_bin_label(255), "NA");
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
