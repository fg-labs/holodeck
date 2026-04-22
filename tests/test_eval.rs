#![allow(
    clippy::type_complexity,
    clippy::cast_possible_truncation,
    clippy::redundant_closure_for_method_calls
)]

mod helpers;

use std::path::PathBuf;

use helpers::{BamRecordSpec, TestEnv, non_repetitive_seq, run_eval, run_simulate, write_bam};
use noodles::sam::alignment::record::Flags;

/// Parse the eval output file and return (total, correct, mismapped, unmapped)
/// from the ALL row.
fn parse_eval_all_row(eval_path: &std::path::Path) -> (u64, u64, u64, u64) {
    let contents = std::fs::read_to_string(eval_path).unwrap();
    for line in contents.lines() {
        if line.starts_with("ALL\t") {
            let fields: Vec<&str> = line.split('\t').collect();
            let total: u64 = fields[1].parse().unwrap();
            let correct: u64 = fields[2].parse().unwrap();
            let mismapped: u64 = fields[3].parse().unwrap();
            let unmapped: u64 = fields[4].parse().unwrap();
            return (total, correct, mismapped, unmapped);
        }
    }
    panic!("ALL row not found in eval output: {}", eval_path.display());
}

/// Simulate single-end reads with a golden BAM, then run eval using the
/// golden BAM as the mapped BAM.  Since positions are truth, eval should
/// report 100% correct at MAPQ 60.
#[test]
fn test_eval_perfect_alignment() {
    let seq = non_repetitive_seq(2_000);
    let env = TestEnv::new(&[("chr1", &seq)]);
    let sim_out = env.output_prefix();

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        sim_out.to_str().unwrap(),
        "--coverage",
        "20",
        "--read-length",
        "50",
        "--fragment-mean",
        "150",
        "--fragment-stddev",
        "20",
        "--golden-bam",
        "--single-end",
        "--seed",
        "42",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", sim_out.display()));
    let eval_out = env.dir.path().join("eval");

    let (ok, _, stderr) = run_eval(&[
        "eval",
        "--mapped",
        bam_path.to_str().unwrap(),
        "-o",
        eval_out.to_str().unwrap(),
    ]);
    assert!(ok, "eval failed: {stderr}");

    let eval_file = PathBuf::from(format!("{}.eval.txt", eval_out.display()));
    assert!(eval_file.exists(), "Eval output should exist");

    let (total, correct, mismapped, unmapped) = parse_eval_all_row(&eval_file);
    assert!(total > 0, "Should have evaluated some reads");
    assert_eq!(mismapped, 0, "No reads should be mismapped with golden BAM");
    assert_eq!(unmapped, 0, "No reads should be unmapped with golden BAM");
    assert_eq!(correct, total, "All reads should be correct");
}

/// Paired-end variant of `test_eval_perfect_alignment`.
///
/// Regression test: previously eval compared every record against the R1
/// truth position, so R2 records in a PE golden BAM were reported as
/// mismapped even though the BAM was ground truth. Now eval picks the R1 or
/// R2 truth based on the record's `is_last_segment` flag.
#[test]
fn test_eval_perfect_alignment_paired_end() {
    let seq = non_repetitive_seq(4_000);
    let env = TestEnv::new(&[("chr1", &seq)]);
    let sim_out = env.output_prefix();

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        sim_out.to_str().unwrap(),
        "--coverage",
        "20",
        "--read-length",
        "50",
        "--fragment-mean",
        "150",
        "--fragment-stddev",
        "20",
        "--golden-bam",
        "--seed",
        "42",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", sim_out.display()));
    let eval_out = env.dir.path().join("eval");

    let (ok, _, stderr) = run_eval(&[
        "eval",
        "--mapped",
        bam_path.to_str().unwrap(),
        "-o",
        eval_out.to_str().unwrap(),
    ]);
    assert!(ok, "eval failed: {stderr}");

    let eval_file = PathBuf::from(format!("{}.eval.txt", eval_out.display()));
    let (total, correct, mismapped, unmapped) = parse_eval_all_row(&eval_file);
    assert!(total > 0, "Should have evaluated some reads");
    // The golden BAM contains both R1 and R2 records; every one must be scored
    // against the matching truth, not just R1's.
    assert_eq!(mismapped, 0, "No PE records should be mismapped against their own truth");
    assert_eq!(unmapped, 0, "No reads should be unmapped in a golden BAM");
    assert_eq!(correct, total, "All R1 and R2 records should be correct");
}

/// Create a BAM where all reads are unmapped.  Eval should report 100%
/// unmapped.
#[test]
fn test_eval_all_unmapped() {
    // First simulate to get valid read names, then create a BAM with those
    // names but unmapped flags.
    let seq = non_repetitive_seq(1_000);
    let env = TestEnv::new(&[("chr1", &seq)]);
    let sim_out = env.output_prefix();

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        sim_out.to_str().unwrap(),
        "--coverage",
        "10",
        "--read-length",
        "50",
        "--fragment-mean",
        "100",
        "--fragment-stddev",
        "20",
        "--golden-bam",
        "--single-end",
        "--seed",
        "42",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    // Read the golden BAM to extract read names, then write a new BAM where
    // every record is unmapped.
    let golden_path = PathBuf::from(format!("{}.golden.bam", sim_out.display()));
    let mut reader = noodles::bam::io::reader::Builder.build_from_path(&golden_path).unwrap();
    let header = reader.read_header().unwrap();

    let names: Vec<String> = reader
        .record_bufs(&header)
        .map(|r| r.unwrap())
        .map(|r| r.name().map(ToString::to_string).unwrap_or_default())
        .collect();

    assert!(!names.is_empty(), "Should have read names from golden BAM");

    // Write unmapped BAM — no contig index, no position, UNMAPPED flag.
    let unmapped_bam = env.dir.path().join("unmapped.bam");
    let bam_records: Vec<BamRecordSpec<'_>> =
        names.iter().map(|name| (name.as_str(), None, None, Flags::UNMAPPED, 0u8)).collect();
    write_bam(&unmapped_bam, &[("chr1", 1_000)], &bam_records);

    let eval_out = env.dir.path().join("eval");
    let (ok, _, stderr) = run_eval(&[
        "eval",
        "--mapped",
        unmapped_bam.to_str().unwrap(),
        "-o",
        eval_out.to_str().unwrap(),
    ]);
    assert!(ok, "eval failed: {stderr}");

    let eval_file = PathBuf::from(format!("{}.eval.txt", eval_out.display()));
    let (total, correct, mismapped, unmapped) = parse_eval_all_row(&eval_file);
    assert!(total > 0, "Should have evaluated some reads");
    assert_eq!(correct, 0, "No reads should be correct when all unmapped");
    assert_eq!(mismapped, 0, "No reads should be mismapped when all unmapped");
    assert_eq!(unmapped, total, "All reads should be unmapped");
}

/// Simulate SE reads, then create a BAM with positions shifted by a known
/// offset.  With wiggle=0 they should be mismapped; with wiggle >= offset they
/// should be correct.
#[test]
fn test_eval_wiggle_parameter() {
    let seq = non_repetitive_seq(2_000);
    let env = TestEnv::new(&[("chr1", &seq)]);
    let sim_out = env.output_prefix();

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        sim_out.to_str().unwrap(),
        "--coverage",
        "10",
        "--read-length",
        "50",
        "--fragment-mean",
        "150",
        "--fragment-stddev",
        "20",
        "--golden-bam",
        "--single-end",
        "--seed",
        "42",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    // Read golden BAM records and rewrite with positions shifted by +3.
    let golden_path = PathBuf::from(format!("{}.golden.bam", sim_out.display()));
    let mut reader = noodles::bam::io::reader::Builder.build_from_path(&golden_path).unwrap();
    let header = reader.read_header().unwrap();

    let shift: u32 = 3;
    let records: Vec<noodles::sam::alignment::RecordBuf> =
        reader.record_bufs(&header).map(|r| r.unwrap()).collect();

    // Build shifted BAM records.
    let shifted_bam = env.dir.path().join("shifted.bam");
    let shifted_records: Vec<BamRecordSpec<'_>> = records
        .iter()
        .filter_map(|r| {
            let name_str = r.name()?.to_string();
            let contig_idx = r.reference_sequence_id()?;
            let pos = r.alignment_start().map(usize::from)?;
            // Leak the name string into a &'static str to satisfy lifetime.
            let name: &str = Box::leak(name_str.into_boxed_str());
            Some((name, Some(contig_idx), Some(pos as u32 + shift), r.flags(), 60u8))
        })
        .collect();
    write_bam(&shifted_bam, &[("chr1", 2_000)], &shifted_records);

    // With wiggle=0, all should be mismapped (shift=3 > 0).
    let eval_out_0 = env.dir.path().join("eval_w0");
    let (ok, _, stderr) = run_eval(&[
        "eval",
        "--mapped",
        shifted_bam.to_str().unwrap(),
        "-o",
        eval_out_0.to_str().unwrap(),
        "--wiggle",
        "0",
    ]);
    assert!(ok, "eval (wiggle=0) failed: {stderr}");

    let eval_file_0 = PathBuf::from(format!("{}.eval.txt", eval_out_0.display()));
    let (total_0, correct_0, mismapped_0, _) = parse_eval_all_row(&eval_file_0);
    assert!(total_0 > 0, "Should have evaluated some reads");
    assert_eq!(correct_0, 0, "wiggle=0 with shift=3 should have 0 correct");
    assert_eq!(mismapped_0, total_0, "All should be mismapped with wiggle=0");

    // With wiggle=3, all should be correct (shift=3 <= 3).
    let eval_out_3 = env.dir.path().join("eval_w3");
    let (ok, _, stderr) = run_eval(&[
        "eval",
        "--mapped",
        shifted_bam.to_str().unwrap(),
        "-o",
        eval_out_3.to_str().unwrap(),
        "--wiggle",
        "3",
    ]);
    assert!(ok, "eval (wiggle=3) failed: {stderr}");

    let eval_file_3 = PathBuf::from(format!("{}.eval.txt", eval_out_3.display()));
    let (total_3, correct_3, mismapped_3, _) = parse_eval_all_row(&eval_file_3);
    assert_eq!(total_3, total_0, "Same reads evaluated in both runs");
    assert_eq!(correct_3, total_3, "wiggle=3 with shift=3 should be all correct");
    assert_eq!(mismapped_3, 0, "No mismapped with sufficient wiggle");
}
