#![allow(
    clippy::type_complexity,
    clippy::cast_possible_truncation,
    clippy::redundant_closure_for_method_calls
)]

mod helpers;

use std::path::PathBuf;

use helpers::{
    BamRecordSpec, TestEnv, VcfVariant, methylate_to_vcf, non_repetitive_seq, run_eval,
    run_simulate, write_bam,
};
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

/// Read the `all` (non-meth) class row of a `.variants.tsv` and return
/// `(n_expected, n_represented)`.
fn parse_variants_all_row(path: &std::path::Path) -> (u64, u64) {
    let contents = std::fs::read_to_string(path).unwrap();
    for line in contents.lines() {
        if let Some(rest) = line.strip_prefix("all\t") {
            let fields: Vec<&str> = rest.split('\t').collect();
            // class row: confounded, n_expected, n_represented, ...
            return (fields[1].parse().unwrap(), fields[2].parse().unwrap());
        }
    }
    panic!("`all` row not found in {}", path.display());
}

/// Read a `#key\tvalue` footer line from a TSV.
fn parse_footer<'a>(contents: &'a str, key: &str) -> &'a str {
    let needle = format!("#{key}\t");
    contents
        .lines()
        .find_map(|l| l.strip_prefix(&needle))
        .unwrap_or_else(|| panic!("footer #{key} not found"))
}

/// Simulate single-end reads carrying homozygous-alt SNVs with a golden BAM,
/// then run eval with the golden BAM as the mapped BAM. Every variant-bearing
/// read is perfectly placed, so all expected substitutions must be represented.
#[test]
fn test_eval_variant_representation_perfect() {
    let seq = non_repetitive_seq(2_000);
    let env = TestEnv::new(&[("chr1", &seq)]);

    // Hom-alt SNVs at known positions; ref base read from the sequence so the
    // VCF matches the reference, alt chosen to differ.
    let positions = [400usize, 800, 1200, 1600];
    let refs: Vec<String> = positions.iter().map(|&p| (seq[p] as char).to_string()).collect();
    let alts: Vec<String> =
        refs.iter().map(|r| if r == "A" { "C" } else { "A" }.to_string()).collect();
    let alt_arrays: Vec<[&str; 1]> = alts.iter().map(|a| [a.as_str()]).collect();
    let variants: Vec<VcfVariant<'_>> = positions
        .iter()
        .enumerate()
        .map(|(i, &p)| VcfVariant {
            chrom: "chr1",
            pos_1based: p as u32 + 1,
            ref_allele: refs[i].as_str(),
            alt_alleles: &alt_arrays[i],
            gt: "1|1",
        })
        .collect();
    let vcf = env.write_vcf("sample", &[("chr1", 2_000)], &variants);

    let sim_out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
        sim_out.to_str().unwrap(),
        "--coverage",
        "30",
        "--read-length",
        "50",
        "--fragment-mean",
        "150",
        "--fragment-stddev",
        "20",
        "--min-error-rate",
        "0",
        "--max-error-rate",
        "0",
        "--golden-bam",
        "--single-end",
        "--seed",
        "42",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let golden = PathBuf::from(format!("{}.golden.bam", sim_out.display()));
    let eval_out = env.dir.path().join("eval");
    let (ok, _, stderr) = run_eval(&[
        "eval",
        "--mapped",
        golden.to_str().unwrap(),
        "--truth",
        golden.to_str().unwrap(),
        "--variants",
        vcf.to_str().unwrap(),
        "-o",
        eval_out.to_str().unwrap(),
    ]);
    assert!(ok, "eval failed: {stderr}");

    let variants_tsv = PathBuf::from(format!("{}.variants.tsv", eval_out.display()));
    let (n_expected, n_represented) = parse_variants_all_row(&variants_tsv);
    assert!(n_expected > 0, "expected some variant-bearing reads");
    assert_eq!(n_represented, n_expected, "golden alignment must represent every variant");

    // A non-methylation golden BAM carries no MD/NM tags, so concordance has
    // nothing to compare against and must report NA (not 0%).
    let contents = std::fs::read_to_string(&variants_tsv).unwrap();
    assert_eq!(parse_footer(&contents, "md_concordant_pct"), "NA");
    assert_eq!(parse_footer(&contents, "nm_concordant_pct"), "NA");
}

/// With `--reference`, NM/MD concordance is the bisulfite-aware genomic edit
/// distance recomputed against the reference (not the raw tags). Grading the
/// golden BAM against itself must be perfect — every read's genomic edits match
/// its own — so both report 100%, not NA.
#[test]
fn test_eval_genomic_nm_md_concordance_with_reference() {
    let seq = non_repetitive_seq(2_000);
    let env = TestEnv::new(&[("chr1", &seq)]);

    let positions = [400usize, 800, 1200, 1600];
    let refs: Vec<String> = positions.iter().map(|&p| (seq[p] as char).to_string()).collect();
    let alts: Vec<String> =
        refs.iter().map(|r| if r == "A" { "C" } else { "A" }.to_string()).collect();
    let alt_arrays: Vec<[&str; 1]> = alts.iter().map(|a| [a.as_str()]).collect();
    let variants: Vec<VcfVariant<'_>> = positions
        .iter()
        .enumerate()
        .map(|(i, &p)| VcfVariant {
            chrom: "chr1",
            pos_1based: p as u32 + 1,
            ref_allele: refs[i].as_str(),
            alt_alleles: &alt_arrays[i],
            gt: "1|1",
        })
        .collect();
    let vcf = env.write_vcf("sample", &[("chr1", 2_000)], &variants);

    let sim_out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
        sim_out.to_str().unwrap(),
        "--coverage",
        "30",
        "--read-length",
        "50",
        "--fragment-mean",
        "150",
        "--fragment-stddev",
        "20",
        "--min-error-rate",
        "0",
        "--max-error-rate",
        "0",
        "--golden-bam",
        "--single-end",
        "--seed",
        "42",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let golden = PathBuf::from(format!("{}.golden.bam", sim_out.display()));
    let eval_out = env.dir.path().join("eval_ref");
    let (ok, _, stderr) = run_eval(&[
        "eval",
        "--mapped",
        golden.to_str().unwrap(),
        "--truth",
        golden.to_str().unwrap(),
        "--variants",
        vcf.to_str().unwrap(),
        "--reference",
        env.fasta_path.to_str().unwrap(),
        "-o",
        eval_out.to_str().unwrap(),
    ]);
    assert!(ok, "eval failed: {stderr}");

    let variants_tsv = PathBuf::from(format!("{}.variants.tsv", eval_out.display()));
    let contents = std::fs::read_to_string(&variants_tsv).unwrap();
    // Golden vs itself: genomic edits are identical → 100% concordance, and
    // crucially NOT "NA" (which is what the raw-tag path returned here).
    assert_eq!(parse_footer(&contents, "md_concordant_pct"), "100.00");
    assert_eq!(parse_footer(&contents, "nm_concordant_pct"), "100.00");
}

/// Simulate EM-seq reads with a methylation golden BAM and a cpg-truth
/// bedGraph, then correlate the golden BAM's own XM calls against that truth.
/// Because both derive from the same methylation draws, the correlation is
/// strong.
#[test]
fn test_eval_meth_correlation_on_golden() {
    let seq = non_repetitive_seq(4_000);
    let env = TestEnv::new(&[("chr1", &seq)]);
    // Mixed methylation (rate 0.5) so truth levels vary across CpGs and the
    // correlation is well-defined (a constant series would be undefined).
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 0.5, 7, "meth.vcf.gz");

    let sim_out = env.output_prefix();
    let bedgraph = env.dir.path().join("truth.bedGraph");
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
        sim_out.to_str().unwrap(),
        "--coverage",
        "30",
        "--read-length",
        "50",
        "--fragment-mean",
        "150",
        "--fragment-stddev",
        "20",
        "--min-error-rate",
        "0",
        "--max-error-rate",
        "0",
        "--methylation-mode",
        "em-seq",
        "--methylation-conversion-rate",
        "1.0",
        "--methylation-failure-rate",
        "0.0",
        "--cpg-truth-bedgraph",
        bedgraph.to_str().unwrap(),
        "--golden-bam",
        "--seed",
        "42",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let golden = PathBuf::from(format!("{}.golden.bam", sim_out.display()));
    let eval_out = env.dir.path().join("eval");
    let (ok, _, stderr) = run_eval(&[
        "eval",
        "--mapped",
        golden.to_str().unwrap(),
        "--cpg-truth",
        bedgraph.to_str().unwrap(),
        "-o",
        eval_out.to_str().unwrap(),
    ]);
    assert!(ok, "eval failed: {stderr}");

    let meth_tsv = std::fs::read_to_string(format!("{}.meth.tsv", eval_out.display())).unwrap();
    // Data row: n_cpg \t pearson_r \t rmse
    let row = meth_tsv.lines().nth(1).expect("meth.tsv data row");
    let fields: Vec<&str> = row.split('\t').collect();
    let n_cpg: u64 = fields[0].parse().unwrap();
    assert!(n_cpg > 0, "expected covered CpGs");
    assert_ne!(fields[1], "NA", "pearson_r should be defined");
    let r: f64 = fields[1].parse().unwrap();
    assert!(r > 0.8, "golden XM should track truth strongly, got r={r}");
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
