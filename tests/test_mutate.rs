#![allow(
    clippy::similar_names,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    clippy::manual_range_contains
)]

mod helpers;

use helpers::{TestEnv, non_repetitive_seq, read_vcf_records, run_mutate};

/// Run mutate on a small reference with default rates and verify the VCF is
/// valid: has variants, positions are in-bounds, and REF alleles match the
/// reference sequence.
#[test]
fn test_mutate_basic() {
    let seq = non_repetitive_seq(10_000);
    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.dir.path().join("output.vcf");

    let (ok, _, stderr) = run_mutate(&[
        "mutate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        vcf_path.to_str().unwrap(),
        "--seed",
        "42",
    ]);
    assert!(ok, "mutate failed: {stderr}");
    assert!(vcf_path.exists(), "VCF file should exist");

    let records = read_vcf_records(&vcf_path);
    assert!(!records.is_empty(), "Should have generated some variants");

    for rec in &records {
        assert_eq!(rec.chrom, "chr1", "All variants should be on chr1");
        assert!(rec.pos >= 1, "VCF positions are 1-based");
        assert!(rec.pos <= seq.len(), "Position should be within contig bounds");

        // REF allele should match the reference sequence.
        let ref_start = rec.pos - 1; // 0-based
        let ref_end = ref_start + rec.ref_allele.len();
        assert!(ref_end <= seq.len(), "REF should not extend beyond contig");
        let expected_ref: String = seq[ref_start..ref_end].iter().map(|&b| b as char).collect();
        assert_eq!(
            rec.ref_allele, expected_ref,
            "REF allele at pos {} should match reference",
            rec.pos
        );
    }
}

/// Verify that running mutate twice with the same seed produces identical VCF
/// output.
#[test]
fn test_mutate_seed_reproducibility() {
    let seq = non_repetitive_seq(5_000);
    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf1 = env.dir.path().join("run1.vcf");
    let vcf2 = env.dir.path().join("run2.vcf");
    let fasta = env.fasta_path.to_str().unwrap().to_string();

    for vcf_path in [&vcf1, &vcf2] {
        let vcf_str = vcf_path.to_str().unwrap();
        let (ok, _, stderr) = run_mutate(&["mutate", "-r", &fasta, "-o", vcf_str, "--seed", "99"]);
        assert!(ok, "mutate failed: {stderr}");
    }

    let content1 = std::fs::read_to_string(&vcf1).unwrap();
    let content2 = std::fs::read_to_string(&vcf2).unwrap();
    assert_eq!(content1, content2, "Same seed should produce identical VCF output");
}

/// Verify that all variants fall within BED target regions when a BED file is
/// provided.
#[test]
fn test_mutate_with_bed_targets() {
    let seq = non_repetitive_seq(10_000);
    let env = TestEnv::new(&[("chr1", &seq)]);
    let bed_path = env.write_bed(&[("chr1", 2000, 3000)]);
    let vcf_path = env.dir.path().join("output.vcf");

    let (ok, _, stderr) = run_mutate(&[
        "mutate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        vcf_path.to_str().unwrap(),
        "-b",
        bed_path.to_str().unwrap(),
        "--snp-rate",
        "0.01",
        "--seed",
        "7",
    ]);
    assert!(ok, "mutate failed: {stderr}");

    let records = read_vcf_records(&vcf_path);
    assert!(!records.is_empty(), "Should have generated variants in target");

    for rec in &records {
        let pos_0based = rec.pos - 1;
        assert!(
            pos_0based >= 2000 && pos_0based < 3000,
            "Variant at 0-based pos {pos_0based} should be within BED target [2000, 3000)"
        );
    }
}

/// Set indel and MNP rates to zero; verify all generated variants are
/// single-base SNPs (ref and alt are both single characters, and differ).
#[test]
fn test_mutate_snp_only() {
    let seq = non_repetitive_seq(10_000);
    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.dir.path().join("output.vcf");

    let (ok, _, stderr) = run_mutate(&[
        "mutate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        vcf_path.to_str().unwrap(),
        "--snp-rate",
        "0.005",
        "--indel-rate",
        "0",
        "--mnp-rate",
        "0",
        "--seed",
        "42",
    ]);
    assert!(ok, "mutate failed: {stderr}");

    let records = read_vcf_records(&vcf_path);
    assert!(!records.is_empty(), "Should have generated some SNPs");

    for rec in &records {
        assert_eq!(
            rec.ref_allele.len(),
            1,
            "SNP-only mode: REF should be single base, got '{}'",
            rec.ref_allele
        );
        assert_eq!(rec.alt_alleles.len(), 1, "Should have exactly one ALT allele");
        assert_eq!(
            rec.alt_alleles[0].len(),
            1,
            "SNP-only mode: ALT should be single base, got '{}'",
            rec.alt_alleles[0]
        );
        assert_ne!(rec.ref_allele, rec.alt_alleles[0], "REF and ALT should differ for a SNP");
    }
}

/// Run with a known high SNP rate and verify the observed variant count is
/// within a reasonable tolerance of the expected count.
#[test]
fn test_mutate_rates_respected() {
    let seq = non_repetitive_seq(10_000);
    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.dir.path().join("output.vcf");
    let snp_rate = 0.01;

    let (ok, _, stderr) = run_mutate(&[
        "mutate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        vcf_path.to_str().unwrap(),
        "--snp-rate",
        &snp_rate.to_string(),
        "--indel-rate",
        "0",
        "--mnp-rate",
        "0",
        "--seed",
        "42",
    ]);
    assert!(ok, "mutate failed: {stderr}");

    let records = read_vcf_records(&vcf_path);
    let count = records.len();

    // Expected: rate * bases = 0.01 * 10000 = 100 variants.
    // Allow generous tolerance for randomness.
    let expected = (snp_rate * seq.len() as f64) as usize;
    let lower = expected * 60 / 100; // -40%
    let upper = expected * 140 / 100; // +40%

    assert!(
        count >= lower && count <= upper,
        "Expected ~{expected} variants (tolerance [{lower}, {upper}]), got {count}"
    );
}

/// Verify that ploidy overrides work: chrX with ploidy=1 should produce
/// haploid genotypes, while chr1 with default ploidy=2 produces diploid.
#[test]
fn test_mutate_ploidy_override() {
    let seq1 = non_repetitive_seq(5_000);
    let seq_x = non_repetitive_seq(5_000);
    let env = TestEnv::new(&[("chr1", &seq1), ("chrX", &seq_x)]);
    let vcf_path = env.dir.path().join("output.vcf");

    let (ok, _, stderr) = run_mutate(&[
        "mutate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        vcf_path.to_str().unwrap(),
        "--snp-rate",
        "0.005",
        "--indel-rate",
        "0",
        "--mnp-rate",
        "0",
        "--ploidy",
        "2",
        "--ploidy-override",
        "chrX=1",
        "--seed",
        "42",
    ]);
    assert!(ok, "mutate failed: {stderr}");

    let records = read_vcf_records(&vcf_path);

    let chr1_records: Vec<_> = records.iter().filter(|r| r.chrom == "chr1").collect();
    let chrx_records: Vec<_> = records.iter().filter(|r| r.chrom == "chrX").collect();

    assert!(!chr1_records.is_empty(), "Should have chr1 variants");
    assert!(!chrx_records.is_empty(), "Should have chrX variants");

    // chr1 (diploid): genotypes should have a separator (/ or |).
    for rec in &chr1_records {
        assert!(
            rec.gt.contains('/') || rec.gt.contains('|'),
            "chr1 diploid GT '{}' should contain a separator",
            rec.gt
        );
    }

    // chrX (haploid): genotypes should be just "0" or "1" (no separator).
    for rec in &chrx_records {
        assert!(
            !rec.gt.contains('/') && !rec.gt.contains('|'),
            "chrX haploid GT '{}' should not contain a separator",
            rec.gt
        );
        // Since we filter hom-ref, the GT should be "1".
        assert_eq!(rec.gt, "1", "chrX haploid non-ref GT should be '1', got '{}'", rec.gt);
    }
}
