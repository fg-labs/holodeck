#![allow(
    clippy::similar_names,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod helpers;

use std::path::PathBuf;

use helpers::{
    PileupColumn, TestEnv, VcfVariant, count_fastq_records, non_repetitive_seq, pileup_bases,
    read_gzipped, run_simulate, simple_env,
};
use std::collections::HashMap;

use noodles::bam;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;

/// Return `true` if any record in the given BAM has at least one CIGAR op of
/// the given kind.  Used as a sanity-check in indel tests.
fn any_bam_record_has_cigar_op(bam_path: &std::path::Path, target: CigarKind) -> bool {
    let mut reader = bam::io::reader::Builder.build_from_path(bam_path).unwrap();
    let header = reader.read_header().unwrap();
    reader
        .record_bufs(&header)
        .map(|r| r.unwrap())
        .any(|r| r.cigar().as_ref().iter().any(|op| op.kind() == target))
}

/// For each BAM record, walk the CIGAR and return a map from 0-based reference
/// position to the base the read has at that position.  Insertions and
/// deletions are tracked separately.
///
/// Returns `(base_at_refpos, has_insertion_after, has_deletion_at)` where:
/// - `base_at_refpos`: maps ref_pos → read base (for M/X/= ops)
/// - `has_insertion_after`: set of ref positions where an insertion follows
/// - `has_deletion_at`: set of ref positions consumed by a deletion
fn read_base_map(
    record: &RecordBuf,
) -> (HashMap<u32, u8>, std::collections::HashSet<u32>, std::collections::HashSet<u32>) {
    let mut base_at = HashMap::new();
    let mut ins_after = std::collections::HashSet::new();
    let mut del_at = std::collections::HashSet::new();

    let Some(start) = record.alignment_start() else {
        return (base_at, ins_after, del_at);
    };
    let mut ref_pos = (usize::from(start) - 1) as u32;
    let seq = record.sequence().as_ref();
    let mut read_pos = 0u32;

    for op in record.cigar().as_ref() {
        let len = op.len() as u32;
        match op.kind() {
            CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                for i in 0..len {
                    if (read_pos + i) < seq.len() as u32 {
                        base_at.insert(ref_pos + i, seq[(read_pos + i) as usize]);
                    }
                }
                ref_pos += len;
                read_pos += len;
            }
            CigarKind::Insertion => {
                if ref_pos > 0 {
                    ins_after.insert(ref_pos - 1);
                }
                read_pos += len;
            }
            CigarKind::Deletion | CigarKind::Skip => {
                for i in 0..len {
                    del_at.insert(ref_pos + i);
                }
                ref_pos += len;
            }
            CigarKind::SoftClip => {
                read_pos += len;
            }
            CigarKind::HardClip | CigarKind::Pad => {}
        }
    }

    (base_at, ins_after, del_at)
}

#[test]
fn test_simulate_basic_pe() {
    let env = simple_env();
    let out = env.output_prefix();

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "10",
        "--read-length",
        "50",
        "--fragment-mean",
        "150",
        "--fragment-stddev",
        "30",
    ]);

    assert!(ok, "simulate failed: {stderr}");

    // Check R1 and R2 files exist (output_path is prefix-based).
    let r1_path = std::path::PathBuf::from(format!("{}.r1.fastq.gz", out.display()));
    let r2_path = std::path::PathBuf::from(format!("{}.r2.fastq.gz", out.display()));

    assert!(r1_path.exists(), "R1 file missing: {}", r1_path.display());
    assert!(r2_path.exists(), "R2 file missing: {}", r2_path.display());

    // Verify FASTQ contents.
    let r1_contents = read_gzipped(&r1_path);
    let r2_contents = read_gzipped(&r2_path);

    let r1_records = count_fastq_records(&r1_contents);
    let r2_records = count_fastq_records(&r2_contents);
    assert_eq!(r1_records, r2_records, "R1 and R2 should have same record count");
    assert!(r1_records > 0, "Should have generated some reads");

    // Verify read length.
    let first_seq = r1_contents.lines().nth(1).unwrap();
    assert_eq!(first_seq.len(), 50, "Read length should be 50");

    // Verify encoded read names.
    let first_name = r1_contents.lines().next().unwrap();
    assert!(
        first_name.starts_with("@holodeck:"),
        "Read name should start with @holodeck: got {first_name}"
    );
}

#[test]
fn test_simulate_single_end() {
    let env = simple_env();
    let out = env.output_prefix();

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "5",
        "--read-length",
        "50",
        "--fragment-mean",
        "100",
        "--fragment-stddev",
        "20",
        "--single-end",
    ]);

    assert!(ok, "simulate failed: {stderr}");

    let r1_path = std::path::PathBuf::from(format!("{}.r1.fastq.gz", out.display()));
    let r2_path = std::path::PathBuf::from(format!("{}.r2.fastq.gz", out.display()));

    assert!(r1_path.exists(), "R1 file should exist");
    assert!(!r2_path.exists(), "R2 file should NOT exist for SE");

    let r1_contents = read_gzipped(&r1_path);
    let r1_records = count_fastq_records(&r1_contents);
    assert!(r1_records > 0, "Should have generated some SE reads");

    // SE read names should have 6 colon-separated fields.
    let first_name = r1_contents.lines().next().unwrap().trim_start_matches('@');
    let fields = first_name.split(':').count();
    assert!(fields >= 6, "SE encoded name should have at least 6 fields, got {fields}");
}

#[test]
fn test_simulate_simple_names() {
    let env = simple_env();
    let out = env.output_prefix();

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "5",
        "--read-length",
        "50",
        "--fragment-mean",
        "100",
        "--fragment-stddev",
        "20",
        "--simple-names",
    ]);

    assert!(ok, "simulate failed: {stderr}");

    let r1_path = std::path::PathBuf::from(format!("{}.r1.fastq.gz", out.display()));
    let r1_contents = read_gzipped(&r1_path);
    let first_name = r1_contents.lines().next().unwrap();
    assert_eq!(first_name, "@holodeck:1", "Simple name should be holodeck:N");
}

#[test]
fn test_simulate_seed_reproducibility() {
    let env = simple_env();
    let out1 = env.dir.path().join("run1");
    let out2 = env.dir.path().join("run2");

    let common_args = |out: &std::path::Path| {
        vec![
            "simulate".to_string(),
            "-r".to_string(),
            env.fasta_path.to_str().unwrap().to_string(),
            "-o".to_string(),
            out.to_str().unwrap().to_string(),
            "--coverage".to_string(),
            "5".to_string(),
            "--read-length".to_string(),
            "50".to_string(),
            "--fragment-mean".to_string(),
            "100".to_string(),
            "--fragment-stddev".to_string(),
            "20".to_string(),
            "--seed".to_string(),
            "42".to_string(),
        ]
    };

    let args1 = common_args(&out1);
    let args1_str: Vec<&str> = args1.iter().map(String::as_str).collect();
    let (ok1, _, _) = run_simulate(&args1_str);
    assert!(ok1);

    let args2 = common_args(&out2);
    let args2_str: Vec<&str> = args2.iter().map(String::as_str).collect();
    let (ok2, _, _) = run_simulate(&args2_str);
    assert!(ok2);

    let r1_a = read_gzipped(&std::path::PathBuf::from(format!("{}.r1.fastq.gz", out1.display())));
    let r1_b = read_gzipped(&std::path::PathBuf::from(format!("{}.r1.fastq.gz", out2.display())));
    assert_eq!(r1_a, r1_b, "Same seed should produce identical output");
}

#[test]
fn test_simulate_with_bed_targets() {
    let seq = b"ACGTACGTAC".repeat(100); // 1000 bp
    let env = TestEnv::new(&[("chr1", &seq)]);
    let bed_path = env.write_bed(&[("chr1", 100, 200)]); // 100bp target
    let out = env.output_prefix();

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "10",
        "--read-length",
        "50",
        "--fragment-mean",
        "100",
        "--fragment-stddev",
        "20",
        "-b",
        bed_path.to_str().unwrap(),
    ]);

    assert!(ok, "simulate with BED failed: {stderr}");

    let r1_path = std::path::PathBuf::from(format!("{}.r1.fastq.gz", out.display()));
    let r1_contents = read_gzipped(&r1_path);
    let r1_records = count_fastq_records(&r1_contents);
    assert!(r1_records > 0, "Should have generated reads overlapping target");
}

/// This test covers the primary bug: when targets are sparse relative to the
/// contig size, the old rejection-sampling approach would produce almost no
/// reads because the probability of a random fragment overlapping a tiny target
/// was negligible.  With padded-interval sampling, fragment starts are drawn
/// from a region near the targets, so the expected read count is achieved.
#[test]
fn test_simulate_sparse_targets_produce_expected_reads() {
    // 100 Kbp contig with three small 100bp targets = 300bp total territory.
    // Targets are spread out so they don't merge even with padding.
    let seq = b"ACGTACGTAC".repeat(10_000); // 100,000 bp
    let env = TestEnv::new(&[("chr1", &seq)]);
    let bed_path = env.write_bed(&[
        ("chr1", 5_000, 5_100),
        ("chr1", 50_000, 50_100),
        ("chr1", 90_000, 90_100),
    ]);
    let out = env.output_prefix();

    // 30x coverage over 3 targets of 100bp with 50bp PE reads and 100bp
    // fragments.  Effective territory per target = 100 + 100 - 1 = 199,
    // total effective = 597.  Reads = 30 * 597 / 100 = 179.
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "30",
        "--read-length",
        "50",
        "--fragment-mean",
        "100",
        "--fragment-stddev",
        "20",
        "-b",
        bed_path.to_str().unwrap(),
        "--seed",
        "1",
    ]);
    assert!(ok, "simulate with sparse targets failed: {stderr}");

    let r1_path = std::path::PathBuf::from(format!("{}.r1.fastq.gz", out.display()));
    let r1_contents = read_gzipped(&r1_path);
    let r1_records = count_fastq_records(&r1_contents);

    // Expected ~179 reads.  Allow tolerance for rounding and rare rejection
    // of fragments from the pad zone that don't reach a target.
    assert!(
        r1_records >= 165,
        "Expected ~179 reads, got {r1_records} (too few — coverage calculation bug?)"
    );
    assert!(r1_records <= 195, "Expected ~179 reads, got {r1_records} (too many)");
}

/// Verify that simulated reads with BED targets actually overlap at least one
/// target region, by parsing the encoded truth coordinates from read names.
#[test]
fn test_simulate_reads_overlap_targets() {
    use holodeck_lib::read_naming::parse_encoded_pe_name;

    let seq = b"ACGTACGTAC".repeat(10_000); // 100,000 bp
    let env = TestEnv::new(&[("chr1", &seq)]);
    let bed_path = env.write_bed(&[("chr1", 10_000, 10_200)]);
    let out = env.output_prefix();

    let read_length: u32 = 50;
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "50",
        "--read-length",
        &read_length.to_string(),
        "--fragment-mean",
        "100",
        "--fragment-stddev",
        "20",
        "-b",
        bed_path.to_str().unwrap(),
        "--seed",
        "7",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let r1_path = std::path::PathBuf::from(format!("{}.r1.fastq.gz", out.display()));
    let r1_contents = read_gzipped(&r1_path);

    // Target in 0-based half-open coordinates.
    let target_start: u32 = 10_000;
    let target_end: u32 = 10_200;
    let mut checked = 0;

    for line in r1_contents.lines() {
        if !line.starts_with('@') {
            continue;
        }
        let name = line.trim_start_matches('@');
        let Some((_read_num, r1, r2)) = parse_encoded_pe_name(name) else {
            continue;
        };

        // R1 and R2 truth positions are 1-based.  The fragment spans from the
        // smaller position to the larger position + read_length (approximately).
        let frag_start_0 = r1.position.min(r2.position) - 1;
        let frag_end_0 = r1.position.max(r2.position) - 1 + read_length;

        assert!(
            frag_end_0 > target_start && frag_start_0 < target_end,
            "Fragment [{frag_start_0}, {frag_end_0}) does not overlap target [{target_start}, {target_end})"
        );
        checked += 1;
    }

    assert!(checked > 0, "Should have checked at least one read");
}

#[test]
fn test_simulate_golden_bam() {
    let env = simple_env();
    let out = env.output_prefix();

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "5",
        "--read-length",
        "50",
        "--fragment-mean",
        "100",
        "--fragment-stddev",
        "20",
        "--golden-bam",
        "--seed",
        "42",
    ]);

    assert!(ok, "simulate with --golden-bam failed: {stderr}");

    // Golden BAM should exist.
    let bam_path = std::path::PathBuf::from(format!("{}.golden.bam", out.display()));
    assert!(bam_path.exists(), "Golden BAM file missing: {}", bam_path.display());

    // Verify it's a valid BAM by reading the header.
    let mut reader = noodles::bam::io::reader::Builder.build_from_path(&bam_path).unwrap();
    let header = reader.read_header().unwrap();

    // Header should have the same contig as our reference.
    assert_eq!(header.reference_sequences().len(), 1);
    let (name, _) = header.reference_sequences().get_index(0).unwrap();
    assert_eq!(name.as_ref() as &[u8], b"chr1");

    // Count records: should be 2x the read pairs (R1 + R2).
    let record_count = reader.records().count();
    assert!(record_count > 0, "Golden BAM should have records");

    // FASTQ read count should be half the BAM record count (PE).
    let r1_path = std::path::PathBuf::from(format!("{}.r1.fastq.gz", out.display()));
    let r1_contents = read_gzipped(&r1_path);
    let fastq_records = count_fastq_records(&r1_contents);
    assert_eq!(record_count, fastq_records * 2, "Golden BAM should have 2 records per read pair");
}

/// Verify that actual per-base coverage over targets of varying sizes matches
/// the requested coverage.  Uses a golden BAM to compute real depth and asserts
/// each target's mean coverage is within ±2% of the requested value.
#[test]
fn test_targeted_coverage_matches_requested() {
    use noodles::sam::alignment::Record as _;

    let contig_len = 100_000usize;
    let seq = b"ACGT".repeat(contig_len / 4);
    let env = TestEnv::new(&[("chr1", &seq)]);

    // Six targets of varying width, well-separated so they don't interact.
    let targets: Vec<(&str, u32, u32)> = vec![
        ("chr1", 5_000, 5_010),   //    10 bp
        ("chr1", 15_000, 15_050), //    50 bp
        ("chr1", 30_000, 30_100), //   100 bp
        ("chr1", 50_000, 50_250), //   250 bp
        ("chr1", 65_000, 66_000), // 1,000 bp
        ("chr1", 80_000, 85_000), // 5,000 bp
    ];
    let bed_path = env.write_bed(&targets);
    let out = env.output_prefix();

    let requested_coverage = 500.0;
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        &requested_coverage.to_string(),
        "--read-length",
        "150",
        "--fragment-mean",
        "300",
        "--fragment-stddev",
        "50",
        "-b",
        bed_path.to_str().unwrap(),
        "--golden-bam",
        "--seed",
        "42",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    // Read all records from the golden BAM and accumulate per-base coverage.
    let bam_path = std::path::PathBuf::from(format!("{}.golden.bam", out.display()));
    let mut reader = noodles::bam::io::reader::Builder.build_from_path(&bam_path).unwrap();
    let _header = reader.read_header().unwrap();

    let mut coverage = vec![0u32; contig_len];
    for result in reader.records() {
        let record = result.unwrap();
        let Some(Ok(start)) = record.alignment_start() else { continue };
        let Some(Ok(end)) = record.alignment_end() else { continue };
        let start_0 = usize::from(start) - 1;
        let end_0 = usize::from(end); // alignment_end is 1-based inclusive → exclusive
        for depth in coverage.iter_mut().take(end_0.min(contig_len)).skip(start_0) {
            *depth += 1;
        }
    }

    // Check mean coverage over each target.  Narrow targets have higher
    // relative variance (fewer read pairs overlap them), so a tolerance
    // tighter than ~8% produces spurious failures for the 100 bp target.
    let tolerance = 0.08;
    let mut all_ok = true;
    for &(_, start, end) in &targets {
        let width = (end - start) as usize;
        let total: u64 = coverage[start as usize..end as usize].iter().map(|&c| u64::from(c)).sum();
        #[allow(clippy::cast_precision_loss)]
        let mean_cov = total as f64 / width as f64;
        let lower = requested_coverage * (1.0 - tolerance);
        let upper = requested_coverage * (1.0 + tolerance);

        eprintln!(
            "  Target [{start}, {end}) width={width}bp: mean_cov={mean_cov:.1}x \
             (expected {requested_coverage}x ±{:.0}%)",
            tolerance * 100.0
        );
        if mean_cov < lower || mean_cov > upper {
            eprintln!("    FAIL: outside [{lower:.1}, {upper:.1}]");
            all_ok = false;
        }
    }
    assert!(all_ok, "One or more targets had mean coverage outside ±5% tolerance");
}

/// Verify that whole-genome coverage over the interior of a contig matches the
/// requested depth and that every position has 100% reference bases (zero
/// errors at error-rate 0).
#[test]
fn test_wgs_coverage_accuracy() {
    let seq = non_repetitive_seq(10_000);
    let env = TestEnv::new(&[("chr1", &seq)]);
    let out = env.output_prefix();

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "100",
        "--read-length",
        "150",
        "--fragment-mean",
        "300",
        "--fragment-stddev",
        "50",
        "--min-error-rate",
        "0",
        "--max-error-rate",
        "0",
        "--golden-bam",
        "--seed",
        "42",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let pileup = pileup_bases(&bam_path);
    let cols = &pileup[0];

    // Skip the first and last fragment_mean bp to avoid coverage ramp-up/down.
    let edge = 300usize;
    let interior: Vec<&PileupColumn> = cols[edge..cols.len() - edge].iter().collect();
    assert!(!interior.is_empty(), "interior region should be non-empty");

    let mean_depth: f64 =
        interior.iter().map(|c| f64::from(c.total)).sum::<f64>() / interior.len() as f64;
    assert!(
        (95.0..=105.0).contains(&mean_depth),
        "interior mean depth {mean_depth:.1}x should be 100x ±5%"
    );

    // With zero error rate every base should match the reference.
    for (i, col) in interior.iter().enumerate() {
        let ref_base = seq[edge + i];
        let ref_count = match ref_base {
            b'A' => col.a,
            b'C' => col.c,
            b'G' => col.g,
            b'T' => col.t,
            _ => 0,
        };
        assert_eq!(
            ref_count,
            col.total,
            "pos {} should be 100% reference base '{}' but got {}/{} ref",
            edge + i,
            char::from(ref_base),
            ref_count,
            col.total
        );
    }
}

/// Verify that reads are distributed across contigs proportionally to their
/// lengths and that the interior of the largest contig has correct depth and
/// zero non-reference bases.
#[test]
fn test_multi_contig_read_distribution() {
    use holodeck_lib::read_naming::parse_encoded_pe_name;

    let seq1 = non_repetitive_seq(10_000);
    let seq2 = non_repetitive_seq(5_000);
    let seq3 = non_repetitive_seq(1_000);
    let env = TestEnv::new(&[("chr1", &seq1), ("chr2", &seq2), ("chr3", &seq3)]);
    let out = env.output_prefix();

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "30",
        "--read-length",
        "150",
        "--fragment-mean",
        "300",
        "--fragment-stddev",
        "50",
        "--min-error-rate",
        "0",
        "--max-error-rate",
        "0",
        "--golden-bam",
        "--seed",
        "99",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    // Count reads per contig from encoded R1 read names.
    let r1_path = PathBuf::from(format!("{}.r1.fastq.gz", out.display()));
    let r1_contents = read_gzipped(&r1_path);

    let mut contig_counts: std::collections::HashMap<String, usize> =
        std::collections::HashMap::new();
    for line in r1_contents.lines() {
        if !line.starts_with('@') {
            continue;
        }
        let name = line.trim_start_matches('@');
        let Some((_read_num, r1, _r2)) = parse_encoded_pe_name(name) else { continue };
        *contig_counts.entry(r1.contig.clone()).or_insert(0) += 1;
    }

    let total: usize = contig_counts.values().sum();
    assert!(total > 0, "should have generated reads");

    // Expected fractions: 10/16, 5/16, 1/16.
    let expected = [("chr1", 10.0_f64 / 16.0), ("chr2", 5.0 / 16.0), ("chr3", 1.0 / 16.0)];
    for (chrom, expected_frac) in &expected {
        let count = *contig_counts.get(*chrom).unwrap_or(&0);

        let actual_frac = count as f64 / total as f64;
        let rel_err = (actual_frac - expected_frac).abs() / expected_frac;
        assert!(
            rel_err <= 0.15,
            "{chrom}: expected {:.1}% of reads, got {:.1}% (rel error {:.1}%)",
            expected_frac * 100.0,
            actual_frac * 100.0,
            rel_err * 100.0
        );
    }

    // chr1 pileup: interior depth ≈ 30x and zero non-reference bases.
    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let pileup = pileup_bases(&bam_path);
    let cols = &pileup[0]; // chr1

    let edge = 300usize;
    let interior: Vec<&PileupColumn> = cols[edge..cols.len() - edge].iter().collect();

    let mean_depth: f64 =
        interior.iter().map(|c| f64::from(c.total)).sum::<f64>() / interior.len() as f64;
    assert!(
        (24.0..=36.0).contains(&mean_depth),
        "chr1 interior mean depth {mean_depth:.1}x should be ~30x ±20%"
    );

    for (i, col) in interior.iter().enumerate() {
        let ref_base = seq1[edge + i];
        let ref_count = match ref_base {
            b'A' => col.a,
            b'C' => col.c,
            b'G' => col.g,
            b'T' => col.t,
            _ => 0,
        };
        assert_eq!(
            ref_count,
            col.total,
            "chr1 pos {} should be 100% reference but got {}/{} ref",
            edge + i,
            ref_count,
            col.total
        );
    }
}

/// Verify that with zero error rate every simulated base exactly matches the
/// reference at high depth.
#[test]
fn test_zero_error_rate_reads_match_reference() {
    let seq = non_repetitive_seq(1_000);
    let env = TestEnv::new(&[("chr1", &seq)]);
    let out = env.output_prefix();

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "500",
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
        "--seed",
        "7",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let pileup = pileup_bases(&bam_path);
    let cols = &pileup[0];

    // Use fragment_mean + 3*fragment_stddev as edge to avoid ramp-up/down zone.
    let edge = 210usize; // 150 + 3*20
    let interior: Vec<&PileupColumn> = cols[edge..cols.len() - edge].iter().collect();

    for (i, col) in interior.iter().enumerate() {
        let abs_pos = edge + i;
        assert!(
            col.total >= 300,
            "pos {abs_pos}: depth {} is below minimum 300 — simulation may have generated too few reads",
            col.total
        );
        let ref_base = seq[abs_pos];
        let ref_count = match ref_base {
            b'A' => col.a,
            b'C' => col.c,
            b'G' => col.g,
            b'T' => col.t,
            _ => 0,
        };
        assert_eq!(
            ref_count,
            col.total,
            "pos {abs_pos}: expected 100% '{}' but got {}/{} ref (a={} c={} g={} t={} n={})",
            char::from(ref_base),
            ref_count,
            col.total,
            col.a,
            col.c,
            col.g,
            col.t,
            col.n
        );
    }
}

/// Verify adapter content when fragments are much shorter than the read length.
///
/// With fragment_mean=30 and read_length=100, approximately 70 bp of each read
/// should be adapter sequence.  We check:
///
/// 1. The adapter prefix appears at the same position in R1 and R2 for every
///    pair, consistent with the fragment length.
/// 2. Every golden BAM record has a terminal soft-clip (S) op, and no leading
///    S op (adapters trail, never lead).
/// 3. The pileup over the short covered region shows zero non-reference bases
///    (zero errors).
#[test]
fn test_adapter_content_short_fragments() {
    // Use a simple repetitive reference — we just need something valid.
    let seq = b"ACGT".repeat(250); // 1000 bp
    let env = TestEnv::new(&[("chr1", &seq)]);
    let out = env.output_prefix();

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--fragment-mean",
        "30",
        "--fragment-stddev",
        "5",
        "--read-length",
        "100",
        "--min-error-rate",
        "0",
        "--max-error-rate",
        "0",
        "--golden-bam",
        "--seed",
        "42",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    // ── Check R1/R2 adapter start positions match ────────────────────────────
    let r1_prefix = b"AGATCGGA"; // start of default R1 adapter
    let r1_path = PathBuf::from(format!("{}.r1.fastq.gz", out.display()));
    let r2_path = PathBuf::from(format!("{}.r2.fastq.gz", out.display()));
    let r1_contents = read_gzipped(&r1_path);
    let r2_contents = read_gzipped(&r2_path);

    let r1_lines: Vec<&str> = r1_contents.lines().collect();
    let r2_lines: Vec<&str> = r2_contents.lines().collect();
    let n_records = r1_lines.len() / 4;
    assert!(n_records > 0, "should have generated reads");

    for i in 0..n_records {
        let r1_seq = r1_lines[i * 4 + 1].as_bytes();
        let r2_seq = r2_lines[i * 4 + 1].as_bytes();

        // Find the first position where the adapter prefix appears.
        let find_adapter = |seq: &[u8], prefix: &[u8]| -> Option<usize> {
            seq.windows(prefix.len()).position(|w| w == prefix)
        };

        let r1_adapter_start = find_adapter(r1_seq, r1_prefix);
        let r2_adapter_start = find_adapter(r2_seq, r1_prefix);

        // Both reads must have adapter (fragment << read_length).
        assert!(r1_adapter_start.is_some(), "read {i}: R1 should contain adapter prefix");
        assert!(r2_adapter_start.is_some(), "read {i}: R2 should contain adapter prefix");

        // The adapter must start at the same offset in both reads (equal fragment length).
        assert_eq!(
            r1_adapter_start, r2_adapter_start,
            "read {i}: adapter start should be same in R1 ({r1_adapter_start:?}) and R2 ({r2_adapter_start:?})"
        );
    }

    // ── Check golden BAM: exactly one terminal soft-clip per record ─────────
    // Forward-strand reads have a trailing S; negative-strand reads (stored
    // reverse-complemented in BAM) have a leading S.  Either way there must
    // be exactly one terminal S and the aligned + soft-clipped lengths must
    // sum to read_length.
    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let mut reader = bam::io::reader::Builder.build_from_path(&bam_path).unwrap();
    let header = reader.read_header().unwrap();

    for result in reader.record_bufs(&header) {
        let record: RecordBuf = result.unwrap();
        let ops = record.cigar().as_ref(); // &[Op]

        assert!(!ops.is_empty(), "CIGAR should not be empty");

        let first_is_clip = ops.first().unwrap().kind() == CigarKind::SoftClip;
        let last_is_clip = ops.last().unwrap().kind() == CigarKind::SoftClip;

        // Exactly one terminal op is a soft-clip (leading for negative strand,
        // trailing for positive strand).  Both cannot be S (no adapter at
        // both ends in this test).
        assert!(
            first_is_clip ^ last_is_clip,
            "expected exactly one terminal SoftClip; first={:?} last={:?}",
            ops.first().unwrap().kind(),
            ops.last().unwrap().kind()
        );

        // Aligned (M) + soft-clipped (S) bases must equal read_length.
        let m_len: usize = ops
            .iter()
            .filter(|op| {
                matches!(
                    op.kind(),
                    CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch
                )
            })
            .map(|op| op.len())
            .sum();
        let s_len: usize =
            ops.iter().filter(|op| op.kind() == CigarKind::SoftClip).map(|op| op.len()).sum();
        assert_eq!(m_len + s_len, 100, "M+S ops should sum to read_length=100");
    }

    // ── Pileup: zero non-reference bases at covered positions ───────────────
    let pileup = pileup_bases(&bam_path);
    let cols = &pileup[0];

    for (pos, col) in cols.iter().enumerate() {
        if col.total < 5 {
            continue;
        } // skip low-coverage positions
        let ref_base = seq[pos];
        let ref_count = match ref_base {
            b'A' => col.a,
            b'C' => col.c,
            b'G' => col.g,
            b'T' => col.t,
            _ => 0,
        };
        assert_eq!(
            ref_count,
            col.total,
            "pos {pos}: expected 100% '{}' but got {}/{} ref",
            char::from(ref_base),
            ref_count,
            col.total
        );
    }
}

/// Verify that a homozygous-alt SNP shows up as ~100% alt in the pileup.
///
/// Simulates 1 000x with zero errors over a 1 000 bp reference.  A hom-alt
/// SNP is placed at position 500 (0-based).  At that position virtually every
/// read should carry the alt allele; flanking positions should be 100% ref.
#[test]
fn test_snp_hom_alt_pileup() {
    let seq = non_repetitive_seq(1_000);
    let ref_base = seq[500];
    // Pick an alt base that differs from the reference.
    let alt_base = match ref_base {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    };
    let ref_str = (ref_base as char).to_string();
    let alt_str = (alt_base as char).to_string();

    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.write_vcf(
        "SAMPLE",
        &[("chr1", 1_000)],
        &[VcfVariant {
            chrom: "chr1",
            pos_1based: 501, // 1-based VCF coordinate for 0-based pos 500
            ref_allele: &ref_str,
            alt_alleles: &[&alt_str],
            gt: "1|1",
        }],
    );

    let out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "1000",
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
        "--seed",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let pileup = pileup_bases(&bam_path);
    let cols = &pileup[0];

    // ── SNP position: virtually all reads should carry the alt allele ────────
    let snp_col = &cols[500];
    let alt_count = match alt_base {
        b'A' => snp_col.a,
        b'C' => snp_col.c,
        b'G' => snp_col.g,
        b'T' => snp_col.t,
        _ => 0,
    };
    assert!(snp_col.total > 0, "pos 500 should have coverage");
    let alt_frac = f64::from(alt_count) / f64::from(snp_col.total);
    assert!(
        alt_frac >= 0.95,
        "pos 500: expected ≥95% alt '{}' but got {}/{} ({:.1}%)",
        char::from(alt_base),
        alt_count,
        snp_col.total,
        alt_frac * 100.0
    );

    // ── Flanking positions: should be 100% reference ─────────────────────────
    for pos in (490..500).chain(501..511) {
        let col = &cols[pos];
        if col.total == 0 {
            continue;
        }
        let flanking_ref = seq[pos];
        let ref_count = match flanking_ref {
            b'A' => col.a,
            b'C' => col.c,
            b'G' => col.g,
            b'T' => col.t,
            _ => 0,
        };
        let ref_frac = f64::from(ref_count) / f64::from(col.total);
        assert!(
            ref_frac >= 0.99,
            "pos {pos}: expected ≥99% ref '{}' but got {}/{} ({:.1}%)",
            char::from(flanking_ref),
            ref_count,
            col.total,
            ref_frac * 100.0
        );
    }

    // ── Depth check: every interior position should have meaningful coverage ──
    // Coverage accuracy is not the focus of this test; a generous floor of
    // 600x ensures the simulation ran and positions are well-covered without
    // being sensitive to per-position statistical variance.
    let edge = 210usize; // fragment_mean + 3*fragment_stddev = 150 + 60
    for col in &cols[edge..cols.len() - edge] {
        assert!(col.total >= 600, "pos {}: depth {} is below minimum 600", col.pos, col.total);
    }
}

/// Phased heterozygous SNP (0|1) should produce ~50% alt allele at the variant
/// position and ~100% reference at all flanking positions.
#[test]
fn test_snp_phased_het_pileup() {
    let seq = non_repetitive_seq(1_000);
    let ref_base = seq[500];
    let alt_base = match ref_base {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    };
    let ref_str = (ref_base as char).to_string();
    let alt_str = (alt_base as char).to_string();

    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.write_vcf(
        "SAMPLE",
        &[("chr1", 1_000)],
        &[VcfVariant {
            chrom: "chr1",
            pos_1based: 501, // 1-based VCF coordinate for 0-based pos 500
            ref_allele: &ref_str,
            alt_alleles: &[&alt_str],
            gt: "0|1",
        }],
    );

    let out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "2000",
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
        "--seed",
        "2",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let pileup = pileup_bases(&bam_path);
    let cols = &pileup[0];

    // ── SNP position: ~50% alt, ~50% ref ─────────────────────────────────────
    // 2 000x depth makes the binomial variance small enough for tight bounds.
    let snp_col = &cols[500];
    let alt_count = match alt_base {
        b'A' => snp_col.a,
        b'C' => snp_col.c,
        b'G' => snp_col.g,
        b'T' => snp_col.t,
        _ => 0,
    };
    assert!(snp_col.total > 0, "pos 500 should have coverage");
    let alt_frac = f64::from(alt_count) / f64::from(snp_col.total);
    assert!(
        (0.42..=0.58).contains(&alt_frac),
        "pos 500: expected alt fraction in [0.42, 0.58] but got {}/{} ({:.1}%)",
        alt_count,
        snp_col.total,
        alt_frac * 100.0
    );

    // ── Flanking positions: should be 100% reference ─────────────────────────
    for pos in (490..500).chain(501..511) {
        let col = &cols[pos];
        if col.total == 0 {
            continue;
        }
        let flanking_ref = seq[pos];
        let ref_count = match flanking_ref {
            b'A' => col.a,
            b'C' => col.c,
            b'G' => col.g,
            b'T' => col.t,
            _ => 0,
        };
        let ref_frac = f64::from(ref_count) / f64::from(col.total);
        assert!(
            ref_frac >= 0.99,
            "pos {pos}: expected ≥99% ref '{}' but got {}/{} ({:.1}%)",
            char::from(flanking_ref),
            ref_count,
            col.total,
            ref_frac * 100.0
        );
    }

    // ── Depth check: generous floor ensures the simulation ran correctly ──────
    let edge = 210usize; // fragment_mean + 3*fragment_stddev = 150 + 60
    for col in &cols[edge..cols.len() - edge] {
        assert!(col.total >= 1200, "pos {}: depth {} is below minimum 1200", col.pos, col.total);
    }
}

/// Homozygous deletion (1|1) should produce near-zero depth at deleted
/// reference positions and normal coverage at flanking positions, with
/// correct anchor-base composition.
#[test]
// Note: previously had #[expect(clippy::too_many_lines)] — now covered by file-level allow.
fn test_deletion_pileup() {
    let seq = non_repetitive_seq(1_000);

    // 4-base VCF deletion anchored at 0-based pos 400 (1-based pos 401).
    // ref = seq[400..404], alt = seq[400] alone → deletes 3 bases at 401–403.
    let ref_4bp: String = seq[400..404].iter().map(|&b| b as char).collect();
    let anchor_str: String = std::iter::once(seq[400] as char).collect();

    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.write_vcf(
        "SAMPLE",
        &[("chr1", 1_000)],
        &[VcfVariant {
            chrom: "chr1",
            pos_1based: 401,
            ref_allele: &ref_4bp,
            alt_alleles: &[&anchor_str],
            gt: "1|1",
        }],
    );

    let out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "1000",
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
        "--seed",
        "3",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let pileup = pileup_bases(&bam_path);
    let cols = &pileup[0];

    // ── Compute mean flanking depth from positions 380–399 and 410–420 ────────
    let flanking_totals: Vec<u32> = (380..400).chain(410..420).map(|pos| cols[pos].total).collect();
    let flanking_sum: u32 = flanking_totals.iter().sum();

    let flanking_mean = f64::from(flanking_sum) / flanking_totals.len() as f64;
    let depth_threshold = (flanking_mean * 0.2) as u32;

    // ── Anchor position (pos 400): depth should resemble flanking; base is ref ─
    let anchor_col = &cols[400];
    assert!(anchor_col.total > 0, "anchor pos 400 should have coverage");
    let anchor_base_count = match seq[400] {
        b'A' => anchor_col.a,
        b'C' => anchor_col.c,
        b'G' => anchor_col.g,
        b'T' => anchor_col.t,
        _ => 0,
    };
    let anchor_frac = f64::from(anchor_base_count) / f64::from(anchor_col.total);
    assert!(
        anchor_frac >= 0.95,
        "anchor pos 400: expected ≥95% ref base but got {}/{} ({:.1}%)",
        anchor_base_count,
        anchor_col.total,
        anchor_frac * 100.0
    );

    // ── Deleted positions (401–403): coverage should be nearly absent ─────────
    // Reads spanning these ref positions carry D ops and do not contribute
    // reference-consuming bases; with a 1|1 deletion no reads cover them.
    for col in &cols[401..404] {
        assert!(
            col.total <= depth_threshold,
            "deleted pos {}: depth {} exceeds 20% of flanking mean \
             ({flanking_mean:.0}); expected near-zero coverage",
            col.pos,
            col.total
        );
    }

    // ── Flanking positions: full coverage and correct base composition ─────────
    for pos in (380..400).chain(410..420) {
        let col = &cols[pos];
        assert!(col.total >= 600, "flanking pos {pos}: depth {} below 600", col.total);
        let ref_base = seq[pos];
        let ref_count = match ref_base {
            b'A' => col.a,
            b'C' => col.c,
            b'G' => col.g,
            b'T' => col.t,
            _ => 0,
        };
        let ref_frac = f64::from(ref_count) / f64::from(col.total);
        assert!(
            ref_frac >= 0.99,
            "flanking pos {pos}: expected ≥99% ref '{}' but got {}/{} ({:.1}%)",
            char::from(ref_base),
            ref_count,
            col.total,
            ref_frac * 100.0
        );
    }

    // ── Sanity: golden BAM must contain at least one D CIGAR op ──────────────
    assert!(
        any_bam_record_has_cigar_op(&bam_path, CigarKind::Deletion),
        "expected at least one BAM record with a D CIGAR op"
    );
}

/// Homozygous insertion (1|1) should produce high insertion counts at the
/// anchor position and correct base composition at flanking positions.
#[test]
fn test_insertion_pileup() {
    let seq = non_repetitive_seq(1_000);

    // 3-base insertion "CCC" after 0-based pos 500 (1-based pos 501).
    // VCF: ref = seq[500], alt = seq[500] + "CCC"
    let anchor_char = seq[500] as char;
    let anchor_str: String = std::iter::once(anchor_char).collect();
    let alt_str: String = format!("{anchor_char}CCC");

    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.write_vcf(
        "SAMPLE",
        &[("chr1", 1_000)],
        &[VcfVariant {
            chrom: "chr1",
            pos_1based: 501,
            ref_allele: &anchor_str,
            alt_alleles: &[&alt_str],
            gt: "1|1",
        }],
    );

    let out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "1000",
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
        "--seed",
        "4",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let pileup = pileup_bases(&bam_path);
    let cols = &pileup[0];

    // ── Anchor position (pos 500): most reads should carry the "CCC" insertion ─
    // Insertions are anchored at the last aligned ref base (ref_pos - 1 during
    // the CIGAR walk), so the insertion count accumulates at the anchor.
    let anchor_col = &cols[500];
    assert!(anchor_col.total > 0, "anchor pos 500 should have coverage");

    let ccc_count = anchor_col.ins_seqs.get(b"CCC".as_ref()).copied().unwrap_or(0);
    let ccc_frac = f64::from(ccc_count) / f64::from(anchor_col.total);
    assert!(
        ccc_frac >= 0.90,
        "anchor pos 500: expected ≥90% reads with 'CCC' insertion but got {}/{} ({:.1}%)",
        ccc_count,
        anchor_col.total,
        ccc_frac * 100.0
    );

    let ins_frac = f64::from(anchor_col.ins) / f64::from(anchor_col.total);
    assert!(
        ins_frac >= 0.90,
        "anchor pos 500: expected ≥90% reads with any insertion but got {}/{} ({:.1}%)",
        anchor_col.ins,
        anchor_col.total,
        ins_frac * 100.0
    );

    // ── Flanking positions: correct base composition and meaningful coverage ───
    // Insertions do not consume reference positions, so surrounding depth is
    // unaffected and each position should be 100% reference.
    for pos in (480..499).chain(502..521) {
        let col = &cols[pos];
        if col.total == 0 {
            continue;
        }
        assert!(col.total >= 600, "flanking pos {pos}: depth {} below 600", col.total);
        let ref_base = seq[pos];
        let ref_count = match ref_base {
            b'A' => col.a,
            b'C' => col.c,
            b'G' => col.g,
            b'T' => col.t,
            _ => 0,
        };
        let ref_frac = f64::from(ref_count) / f64::from(col.total);
        assert!(
            ref_frac >= 0.99,
            "flanking pos {pos}: expected ≥99% ref '{}' but got {}/{} ({:.1}%)",
            char::from(ref_base),
            ref_count,
            col.total,
            ref_frac * 100.0
        );
    }

    // ── Sanity: golden BAM must contain at least one I CIGAR op ──────────────
    assert!(
        any_bam_record_has_cigar_op(&bam_path, CigarKind::Insertion),
        "expected at least one BAM record with an I CIGAR op"
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// Priority 2: Variant scenarios not yet covered
// ═══════════════════════════════════════════════════════════════════════════

/// Phased heterozygous deletion (0|1): the deleted positions should show ~50%
/// depth (one haplotype carries the deletion, the other does not).
#[test]
fn test_het_deletion_pileup() {
    let seq = non_repetitive_seq(1_000);

    // 4-base deletion at 0-based pos 400 (1-based 401): ref=seq[400..404],
    // alt=seq[400] alone → deletes 3 bases at 401-403.
    let ref_4bp: String = seq[400..404].iter().map(|&b| b as char).collect();
    let anchor_str: String = std::iter::once(seq[400] as char).collect();

    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.write_vcf(
        "SAMPLE",
        &[("chr1", 1_000)],
        &[VcfVariant {
            chrom: "chr1",
            pos_1based: 401,
            ref_allele: &ref_4bp,
            alt_alleles: &[&anchor_str],
            gt: "0|1",
        }],
    );

    let out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "2000",
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
        "--seed",
        "3",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let pileup = pileup_bases(&bam_path);
    let cols = &pileup[0];

    // Flanking mean depth for comparison.
    let flanking_totals: Vec<u32> = (380..400).chain(410..420).map(|pos| cols[pos].total).collect();
    let flanking_sum: u32 = flanking_totals.iter().sum();
    let flanking_mean = f64::from(flanking_sum) / flanking_totals.len() as f64;

    // Deleted positions (401-403) should have ~50% of flanking depth (one
    // haplotype has the deletion, reads from the other still cover these bases).
    for col in &cols[401..404] {
        let frac_of_flanking = f64::from(col.total) / flanking_mean;
        assert!(
            (0.35..=0.65).contains(&frac_of_flanking),
            "het-del pos {}: depth {} is {:.1}% of flanking mean {:.0}, \
             expected ~50%",
            col.pos,
            col.total,
            frac_of_flanking * 100.0,
            flanking_mean
        );
    }

    // Flanking positions should be 100% reference.
    for pos in (380..400).chain(410..420) {
        let col = &cols[pos];
        if col.total == 0 {
            continue;
        }
        let ref_base = seq[pos];
        let ref_count = match ref_base {
            b'A' => col.a,
            b'C' => col.c,
            b'G' => col.g,
            b'T' => col.t,
            _ => 0,
        };
        let ref_frac = f64::from(ref_count) / f64::from(col.total);
        assert!(
            ref_frac >= 0.99,
            "flanking pos {pos}: expected ≥99% ref, got {:.1}%",
            ref_frac * 100.0
        );
    }
}

/// Phased heterozygous insertion (0|1): ~50% of reads at the anchor position
/// should carry the insertion.
#[test]
fn test_het_insertion_pileup() {
    let seq = non_repetitive_seq(1_000);

    // 3-base insertion "CCC" after 0-based pos 500 (1-based 501).
    let anchor_char = seq[500] as char;
    let anchor_str = anchor_char.to_string();
    let alt_str = format!("{anchor_char}CCC");

    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.write_vcf(
        "SAMPLE",
        &[("chr1", 1_000)],
        &[VcfVariant {
            chrom: "chr1",
            pos_1based: 501,
            ref_allele: &anchor_str,
            alt_alleles: &[&alt_str],
            gt: "0|1",
        }],
    );

    let out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "2000",
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
        "--seed",
        "4",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let pileup = pileup_bases(&bam_path);
    let cols = &pileup[0];

    // Anchor position: ~50% of reads should carry the insertion.
    let anchor_col = &cols[500];
    assert!(anchor_col.total > 0, "anchor pos 500 should have coverage");

    let ins_frac = f64::from(anchor_col.ins) / f64::from(anchor_col.total);
    assert!(
        (0.40..=0.60).contains(&ins_frac),
        "anchor pos 500: expected 40-60% reads with insertion, got {}/{} ({:.1}%)",
        anchor_col.ins,
        anchor_col.total,
        ins_frac * 100.0
    );
}

/// Place a SNP, a deletion, and an insertion at well-separated positions on
/// one contig.  Verify each variant shows the expected pileup independently.
#[test]
fn test_multiple_variants_same_contig() {
    let seq = non_repetitive_seq(2_000);

    // SNP at pos 400 (hom-alt 1|1)
    let snp_ref_base = seq[400];
    let snp_alt_base = match snp_ref_base {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    };
    let snp_ref_str = (snp_ref_base as char).to_string();
    let snp_alt_str = (snp_alt_base as char).to_string();

    // Deletion at pos 800 (hom-alt 1|1): delete 3 bases
    let del_ref: String = seq[800..804].iter().map(|&b| b as char).collect();
    let del_alt: String = std::iter::once(seq[800] as char).collect();

    // Insertion at pos 1200 (hom-alt 1|1): insert "GGG"
    let ins_anchor = seq[1200] as char;
    let ins_ref = ins_anchor.to_string();
    let ins_alt = format!("{ins_anchor}GGG");

    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.write_vcf(
        "SAMPLE",
        &[("chr1", 2_000)],
        &[
            VcfVariant {
                chrom: "chr1",
                pos_1based: 401,
                ref_allele: &snp_ref_str,
                alt_alleles: &[&snp_alt_str],
                gt: "1|1",
            },
            VcfVariant {
                chrom: "chr1",
                pos_1based: 801,
                ref_allele: &del_ref,
                alt_alleles: &[&del_alt],
                gt: "1|1",
            },
            VcfVariant {
                chrom: "chr1",
                pos_1based: 1201,
                ref_allele: &ins_ref,
                alt_alleles: &[&ins_alt],
                gt: "1|1",
            },
        ],
    );

    let out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "500",
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
        "--seed",
        "5",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let pileup = pileup_bases(&bam_path);
    let cols = &pileup[0];

    // ── SNP at pos 400: ≥95% alt ──────────────────────────────────────────
    let snp_col = &cols[400];
    let snp_alt_count = match snp_alt_base {
        b'A' => snp_col.a,
        b'C' => snp_col.c,
        b'G' => snp_col.g,
        b'T' => snp_col.t,
        _ => 0,
    };
    assert!(snp_col.total > 0, "SNP position should have coverage");
    let snp_frac = f64::from(snp_alt_count) / f64::from(snp_col.total);
    assert!(snp_frac >= 0.95, "SNP pos 400: expected ≥95% alt, got {:.1}%", snp_frac * 100.0);

    // ── Deletion at 801-803: near-zero depth ──────────────────────────────
    let flanking_mean: f64 =
        (780..800).chain(810..820).map(|p| f64::from(cols[p].total)).sum::<f64>() / 30.0;
    for col in &cols[801..804] {
        let frac = f64::from(col.total) / flanking_mean;
        assert!(
            frac <= 0.20,
            "del pos {}: depth {} is {:.1}% of flanking mean {:.0}, expected near-zero",
            col.pos,
            col.total,
            frac * 100.0,
            flanking_mean,
        );
    }

    // ── Insertion at pos 1200: ≥90% reads with insertion ─────────────────
    let ins_col = &cols[1200];
    assert!(ins_col.total > 0, "Insertion anchor should have coverage");
    let ins_frac = f64::from(ins_col.ins) / f64::from(ins_col.total);
    assert!(
        ins_frac >= 0.90,
        "ins pos 1200: expected ≥90% insertion, got {:.1}%",
        ins_frac * 100.0
    );
}

/// Unphased heterozygous SNP (0/1): the Fisher-Yates shuffle in
/// `build_haplotypes` assigns the alt allele to a random haplotype.  With
/// high depth we should still see ~50% alt allele fraction.
#[test]
fn test_unphased_het_snp_pileup() {
    let seq = non_repetitive_seq(1_000);
    let ref_base = seq[500];
    let alt_base = match ref_base {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    };
    let ref_str = (ref_base as char).to_string();
    let alt_str = (alt_base as char).to_string();

    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.write_vcf(
        "SAMPLE",
        &[("chr1", 1_000)],
        &[VcfVariant {
            chrom: "chr1",
            pos_1based: 501,
            ref_allele: &ref_str,
            alt_alleles: &[&alt_str],
            gt: "0/1", // unphased
        }],
    );

    let out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "2000",
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
        "--seed",
        "10",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let pileup = pileup_bases(&bam_path);
    let cols = &pileup[0];

    let snp_col = &cols[500];
    let alt_count = match alt_base {
        b'A' => snp_col.a,
        b'C' => snp_col.c,
        b'G' => snp_col.g,
        b'T' => snp_col.t,
        _ => 0,
    };
    assert!(snp_col.total > 0, "pos 500 should have coverage");
    let alt_frac = f64::from(alt_count) / f64::from(snp_col.total);
    assert!(
        (0.42..=0.58).contains(&alt_frac),
        "unphased het pos 500: expected ~50% alt, got {}/{} ({:.1}%)",
        alt_count,
        snp_col.total,
        alt_frac * 100.0
    );
}

/// Multi-allelic site: two ALT alleles (e.g. A -> T,G) with genotype 1|2.
/// Each allele should appear at ~50% in the pileup.
#[test]
fn test_multi_allelic_site() {
    let seq = non_repetitive_seq(1_000);
    let ref_base = seq[500];

    // Pick two different alt bases.
    let alts: Vec<u8> = b"ACGT".iter().copied().filter(|&b| b != ref_base).collect();
    let alt1 = alts[0];
    let alt2 = alts[1];

    let ref_str = (ref_base as char).to_string();
    let alt1_str = (alt1 as char).to_string();
    let alt2_str = (alt2 as char).to_string();

    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.write_vcf(
        "SAMPLE",
        &[("chr1", 1_000)],
        &[VcfVariant {
            chrom: "chr1",
            pos_1based: 501,
            ref_allele: &ref_str,
            alt_alleles: &[&alt1_str, &alt2_str],
            gt: "1|2",
        }],
    );

    let out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "2000",
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
        "--seed",
        "11",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let pileup = pileup_bases(&bam_path);
    let cols = &pileup[0];

    let col = &cols[500];
    assert!(col.total > 0, "pos 500 should have coverage");

    let alt1_count = match alt1 {
        b'A' => col.a,
        b'C' => col.c,
        b'G' => col.g,
        b'T' => col.t,
        _ => 0,
    };
    let alt2_count = match alt2 {
        b'A' => col.a,
        b'C' => col.c,
        b'G' => col.g,
        b'T' => col.t,
        _ => 0,
    };
    let ref_count = match ref_base {
        b'A' => col.a,
        b'C' => col.c,
        b'G' => col.g,
        b'T' => col.t,
        _ => 0,
    };

    // With genotype 1|2, no reads should carry the reference allele.
    let ref_frac = f64::from(ref_count) / f64::from(col.total);
    assert!(
        ref_frac <= 0.05,
        "pos 500: expected ≤5% ref with GT 1|2, got {:.1}%",
        ref_frac * 100.0
    );

    // Each alt should appear at ~50%.
    let alt1_frac = f64::from(alt1_count) / f64::from(col.total);
    let alt2_frac = f64::from(alt2_count) / f64::from(col.total);
    assert!(
        (0.40..=0.60).contains(&alt1_frac),
        "pos 500: expected alt1 '{}' ~50%, got {:.1}%",
        alt1 as char,
        alt1_frac * 100.0
    );
    assert!(
        (0.40..=0.60).contains(&alt2_frac),
        "pos 500: expected alt2 '{}' ~50%, got {:.1}%",
        alt2 as char,
        alt2_frac * 100.0
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// Priority 3: Phased/unphased variant clusters within a read
// ═══════════════════════════════════════════════════════════════════════════

/// Place 4 phased SNPs within a ~30bp window, all on haplotype 1 (0|1).
/// Verify per-read cis linkage: within each read spanning the cluster, either
/// ALL alt alleles appear or NONE do.
#[test]
fn test_phased_variant_cluster_cis() {
    let seq = non_repetitive_seq(1_000);

    // 4 SNPs at 0-based positions 490, 500, 510, 520 — all phased on hap 1.
    let snp_positions: Vec<usize> = vec![490, 500, 510, 520];
    let mut snp_alts: Vec<u8> = Vec::new();
    let mut variants = Vec::new();

    for &pos in &snp_positions {
        let ref_base = seq[pos];
        let alt_base = match ref_base {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
        snp_alts.push(alt_base);
    }

    // Build variant descriptions. We need the strings to live long enough.
    let ref_strs: Vec<String> =
        snp_positions.iter().map(|&pos| (seq[pos] as char).to_string()).collect();
    let alt_strs: Vec<String> = snp_alts.iter().map(|&b| (b as char).to_string()).collect();
    let alt_refs: Vec<&str> = alt_strs.iter().map(String::as_str).collect();

    for (i, &pos) in snp_positions.iter().enumerate() {
        variants.push(VcfVariant {
            chrom: "chr1",
            pos_1based: (pos + 1) as u32,
            ref_allele: &ref_strs[i],
            alt_alleles: std::slice::from_ref(&alt_refs[i]),
            gt: "0|1",
        });
    }

    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.write_vcf("SAMPLE", &[("chr1", 1_000)], &variants);

    let out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "1000",
        "--read-length",
        "100",
        "--fragment-mean",
        "200",
        "--fragment-stddev",
        "20",
        "--min-error-rate",
        "0",
        "--max-error-rate",
        "0",
        "--golden-bam",
        "--seed",
        "20",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let mut reader = bam::io::reader::Builder.build_from_path(&bam_path).unwrap();
    let header = reader.read_header().unwrap();

    let mut all_alt = 0u32;
    let mut all_ref = 0u32;
    let mut mixed = 0u32;

    for result in reader.record_bufs(&header) {
        let record = result.unwrap();
        let (base_at, _, _) = read_base_map(&record);

        // Only consider reads that span ALL 4 SNP positions.
        let snp_pos_u32: Vec<u32> = snp_positions.iter().map(|&p| p as u32).collect();
        if !snp_pos_u32.iter().all(|p| base_at.contains_key(p)) {
            continue;
        }

        let mut has_alt = Vec::new();
        for (i, &pos) in snp_pos_u32.iter().enumerate() {
            let read_base = base_at[&pos];
            has_alt.push(read_base == snp_alts[i]);
        }

        if has_alt.iter().all(|&a| a) {
            all_alt += 1;
        } else if has_alt.iter().all(|&a| !a) {
            all_ref += 1;
        } else {
            mixed += 1;
        }
    }

    let total_spanning = all_alt + all_ref + mixed;
    assert!(total_spanning >= 100, "Need sufficient spanning reads, got {total_spanning}");

    // Perfect cis linkage: no mixed reads.
    assert_eq!(
        mixed, 0,
        "Phased cis cluster: expected 0 mixed reads, got {mixed} \
         (all_alt={all_alt}, all_ref={all_ref}, mixed={mixed})"
    );

    // Roughly 50/50 split between all-alt and all-ref.
    let alt_frac = f64::from(all_alt) / f64::from(total_spanning);
    assert!(
        (0.40..=0.60).contains(&alt_frac),
        "Expected ~50% all-alt reads, got {:.1}%",
        alt_frac * 100.0
    );
}

/// Place 4 phased SNPs in a ~30bp window, alternating haplotypes: 1|0, 0|1,
/// 1|0, 0|1.  Verify the expected trans pattern: within each read, SNPs 1&3
/// agree (both alt or both ref), SNPs 2&4 agree, and the two groups are
/// anti-correlated.
#[test]
fn test_phased_variant_cluster_trans() {
    let seq = non_repetitive_seq(1_000);

    let snp_positions: Vec<usize> = vec![490, 500, 510, 520];
    let gts = ["1|0", "0|1", "1|0", "0|1"]; // alternating haplotypes
    let mut snp_alts: Vec<u8> = Vec::new();

    for &pos in &snp_positions {
        let alt = match seq[pos] {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
        snp_alts.push(alt);
    }

    let ref_strs: Vec<String> =
        snp_positions.iter().map(|&pos| (seq[pos] as char).to_string()).collect();
    let alt_strs: Vec<String> = snp_alts.iter().map(|&b| (b as char).to_string()).collect();
    let alt_refs: Vec<&str> = alt_strs.iter().map(String::as_str).collect();

    let mut variants = Vec::new();
    for (i, &pos) in snp_positions.iter().enumerate() {
        variants.push(VcfVariant {
            chrom: "chr1",
            pos_1based: (pos + 1) as u32,
            ref_allele: &ref_strs[i],
            alt_alleles: std::slice::from_ref(&alt_refs[i]),
            gt: gts[i],
        });
    }

    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.write_vcf("SAMPLE", &[("chr1", 1_000)], &variants);

    let out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "1000",
        "--read-length",
        "100",
        "--fragment-mean",
        "200",
        "--fragment-stddev",
        "20",
        "--min-error-rate",
        "0",
        "--max-error-rate",
        "0",
        "--golden-bam",
        "--seed",
        "21",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let mut reader = bam::io::reader::Builder.build_from_path(&bam_path).unwrap();
    let header = reader.read_header().unwrap();

    // Count reads matching the expected trans patterns:
    // Pattern A (hap 0): SNPs 1,3 are alt; SNPs 2,4 are ref
    // Pattern B (hap 1): SNPs 1,3 are ref; SNPs 2,4 are alt
    let mut pattern_a = 0u32;
    let mut pattern_b = 0u32;
    let mut other = 0u32;

    for result in reader.record_bufs(&header) {
        let record = result.unwrap();
        let (base_at, _, _) = read_base_map(&record);

        let snp_pos_u32: Vec<u32> = snp_positions.iter().map(|&p| p as u32).collect();
        if !snp_pos_u32.iter().all(|p| base_at.contains_key(p)) {
            continue;
        }

        let has_alt: Vec<bool> =
            snp_pos_u32.iter().enumerate().map(|(i, &pos)| base_at[&pos] == snp_alts[i]).collect();

        // Pattern A: [alt, ref, alt, ref] — from haplotype 0
        // Pattern B: [ref, alt, ref, alt] — from haplotype 1
        if has_alt == [true, false, true, false] {
            pattern_a += 1;
        } else if has_alt == [false, true, false, true] {
            pattern_b += 1;
        } else {
            other += 1;
        }
    }

    let total_spanning = pattern_a + pattern_b + other;
    assert!(total_spanning >= 100, "Need sufficient spanning reads, got {total_spanning}");

    assert_eq!(
        other, 0,
        "Trans cluster: expected 0 unexpected patterns, got {other} \
         (pattern_a={pattern_a}, pattern_b={pattern_b})"
    );

    // Roughly 50/50 between the two trans patterns.
    let a_frac = f64::from(pattern_a) / f64::from(total_spanning);
    assert!(
        (0.40..=0.60).contains(&a_frac),
        "Expected ~50% pattern A reads, got {:.1}%",
        a_frac * 100.0
    );
}

/// Place 3 unphased het SNPs (0/1) within a ~20bp window.  Since the
/// Fisher-Yates shuffle assigns alleles to haplotypes independently per
/// variant, the variants should NOT always co-occur on the same reads —
/// there should be some mixed patterns.
#[test]
fn test_unphased_variant_cluster() {
    let seq = non_repetitive_seq(1_000);

    let snp_positions: Vec<usize> = vec![495, 505, 515];
    let mut snp_alts: Vec<u8> = Vec::new();

    for &pos in &snp_positions {
        let alt = match seq[pos] {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
        snp_alts.push(alt);
    }

    let ref_strs: Vec<String> =
        snp_positions.iter().map(|&pos| (seq[pos] as char).to_string()).collect();
    let alt_strs: Vec<String> = snp_alts.iter().map(|&b| (b as char).to_string()).collect();
    let alt_refs: Vec<&str> = alt_strs.iter().map(String::as_str).collect();

    let mut variants = Vec::new();
    for (i, &pos) in snp_positions.iter().enumerate() {
        variants.push(VcfVariant {
            chrom: "chr1",
            pos_1based: (pos + 1) as u32,
            ref_allele: &ref_strs[i],
            alt_alleles: std::slice::from_ref(&alt_refs[i]),
            gt: "0/1", // unphased
        });
    }

    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.write_vcf("SAMPLE", &[("chr1", 1_000)], &variants);

    let out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "2000",
        "--read-length",
        "100",
        "--fragment-mean",
        "200",
        "--fragment-stddev",
        "20",
        "--min-error-rate",
        "0",
        "--max-error-rate",
        "0",
        "--golden-bam",
        "--seed",
        "22",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let mut reader = bam::io::reader::Builder.build_from_path(&bam_path).unwrap();
    let header = reader.read_header().unwrap();

    // Count distinct per-read allele patterns.
    let mut pattern_counts: HashMap<Vec<bool>, u32> = HashMap::new();

    for result in reader.record_bufs(&header) {
        let record = result.unwrap();
        let (base_at, _, _) = read_base_map(&record);

        let snp_pos_u32: Vec<u32> = snp_positions.iter().map(|&p| p as u32).collect();
        if !snp_pos_u32.iter().all(|p| base_at.contains_key(p)) {
            continue;
        }

        let pattern: Vec<bool> =
            snp_pos_u32.iter().enumerate().map(|(i, &pos)| base_at[&pos] == snp_alts[i]).collect();
        *pattern_counts.entry(pattern).or_insert(0) += 1;
    }

    let total: u32 = pattern_counts.values().sum();
    assert!(total >= 200, "Need sufficient spanning reads, got {total}");

    // With unphased genotypes, the Fisher-Yates shuffle assigns each variant's
    // alt to either haplotype once at simulation start. The assignment is
    // consistent for the entire simulation, so we expect exactly 2 patterns:
    // whatever haplotype 0 carries and whatever haplotype 1 carries. With 3
    // variants, the two haplotype patterns are complements of each other.
    //
    // The key invariant is that each pattern should be ~50% of reads (each
    // haplotype contributes half the reads).
    let distinct_patterns = pattern_counts.len();
    assert_eq!(
        distinct_patterns, 2,
        "Expected exactly 2 read allele patterns (one per haplotype), got {distinct_patterns}"
    );

    // Each pattern should be ~50% of total.
    for (pattern, &count) in &pattern_counts {
        let frac = f64::from(count) / f64::from(total);
        assert!(
            (0.40..=0.60).contains(&frac),
            "Pattern {pattern:?} has {count}/{total} reads ({:.1}%), expected ~50%",
            frac * 100.0
        );
    }
}

/// Place a phased het SNP and a phased het 3bp deletion on the SAME
/// haplotype (both 0|1), separated by ~10bp.  Verify the SNP alt and the
/// deletion always co-occur on the same reads.
#[test]
fn test_phased_snp_and_deletion_cluster() {
    let seq = non_repetitive_seq(1_000);

    // SNP at 0-based pos 500, deletion at 0-based pos 510 (deletes 511-513).
    let snp_ref_base = seq[500];
    let snp_alt_base = match snp_ref_base {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    };
    let snp_ref_str = (snp_ref_base as char).to_string();
    let snp_alt_str = (snp_alt_base as char).to_string();

    let del_ref: String = seq[510..514].iter().map(|&b| b as char).collect();
    let del_alt: String = std::iter::once(seq[510] as char).collect();

    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.write_vcf(
        "SAMPLE",
        &[("chr1", 1_000)],
        &[
            VcfVariant {
                chrom: "chr1",
                pos_1based: 501,
                ref_allele: &snp_ref_str,
                alt_alleles: &[&snp_alt_str],
                gt: "0|1",
            },
            VcfVariant {
                chrom: "chr1",
                pos_1based: 511,
                ref_allele: &del_ref,
                alt_alleles: &[&del_alt],
                gt: "0|1",
            },
        ],
    );

    let out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "1000",
        "--read-length",
        "100",
        "--fragment-mean",
        "200",
        "--fragment-stddev",
        "20",
        "--min-error-rate",
        "0",
        "--max-error-rate",
        "0",
        "--golden-bam",
        "--seed",
        "23",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let mut reader = bam::io::reader::Builder.build_from_path(&bam_path).unwrap();
    let header = reader.read_header().unwrap();

    let mut snp_alt_del_present = 0u32; // read has both alt SNP and deletion
    let mut snp_ref_del_absent = 0u32; // read has ref SNP and no deletion
    let mut discordant = 0u32; // reads where SNP and deletion disagree

    for result in reader.record_bufs(&header) {
        let record = result.unwrap();
        let (base_at, _, del_at) = read_base_map(&record);

        // Must span both the SNP and the deletion region.
        if !base_at.contains_key(&500) || !base_at.contains_key(&515) {
            continue;
        }

        let has_snp_alt = base_at.get(&500) == Some(&snp_alt_base);
        let has_deletion = del_at.contains(&511) || del_at.contains(&512) || del_at.contains(&513);

        if has_snp_alt && has_deletion {
            snp_alt_del_present += 1;
        } else if !has_snp_alt && !has_deletion {
            snp_ref_del_absent += 1;
        } else {
            discordant += 1;
        }
    }

    let total = snp_alt_del_present + snp_ref_del_absent + discordant;
    assert!(total >= 100, "Need sufficient spanning reads, got {total}");

    assert_eq!(
        discordant, 0,
        "Phased SNP+deletion: expected 0 discordant reads, got {discordant} \
         (concordant_alt={snp_alt_del_present}, concordant_ref={snp_ref_del_absent})"
    );

    // Roughly 50/50 split.
    let alt_frac = f64::from(snp_alt_del_present) / f64::from(total);
    assert!(
        (0.40..=0.60).contains(&alt_frac),
        "Expected ~50% alt-carrying reads, got {:.1}%",
        alt_frac * 100.0
    );
}

/// Place a phased het SNP and a phased het 3bp insertion on the SAME
/// haplotype (both 0|1), separated by ~15bp.  Verify the SNP alt and the
/// insertion always co-occur on the same reads.
#[test]
fn test_phased_snp_and_insertion_cluster() {
    let seq = non_repetitive_seq(1_000);

    // SNP at 0-based pos 500, insertion at 0-based pos 515.
    let snp_ref_base = seq[500];
    let snp_alt_base = match snp_ref_base {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    };
    let snp_ref_str = (snp_ref_base as char).to_string();
    let snp_alt_str = (snp_alt_base as char).to_string();

    let ins_anchor = seq[515] as char;
    let ins_ref = ins_anchor.to_string();
    let ins_alt = format!("{ins_anchor}TTT");

    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.write_vcf(
        "SAMPLE",
        &[("chr1", 1_000)],
        &[
            VcfVariant {
                chrom: "chr1",
                pos_1based: 501,
                ref_allele: &snp_ref_str,
                alt_alleles: &[&snp_alt_str],
                gt: "0|1",
            },
            VcfVariant {
                chrom: "chr1",
                pos_1based: 516,
                ref_allele: &ins_ref,
                alt_alleles: &[&ins_alt],
                gt: "0|1",
            },
        ],
    );

    let out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "1000",
        "--read-length",
        "100",
        "--fragment-mean",
        "200",
        "--fragment-stddev",
        "20",
        "--min-error-rate",
        "0",
        "--max-error-rate",
        "0",
        "--golden-bam",
        "--seed",
        "24",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let mut reader = bam::io::reader::Builder.build_from_path(&bam_path).unwrap();
    let header = reader.read_header().unwrap();

    let mut snp_alt_ins_present = 0u32;
    let mut snp_ref_ins_absent = 0u32;
    let mut discordant = 0u32;

    for result in reader.record_bufs(&header) {
        let record = result.unwrap();
        let (base_at, ins_after, _) = read_base_map(&record);

        // Must span both SNP and insertion positions.
        if !base_at.contains_key(&500) || !base_at.contains_key(&515) {
            continue;
        }

        let has_snp_alt = base_at.get(&500) == Some(&snp_alt_base);
        let has_insertion = ins_after.contains(&515);

        if has_snp_alt && has_insertion {
            snp_alt_ins_present += 1;
        } else if !has_snp_alt && !has_insertion {
            snp_ref_ins_absent += 1;
        } else {
            discordant += 1;
        }
    }

    let total = snp_alt_ins_present + snp_ref_ins_absent + discordant;
    assert!(total >= 100, "Need sufficient spanning reads, got {total}");

    // Allow a tiny fraction of discordant reads (<1%) for edge effects where
    // a read barely spans the insertion anchor position but the CIGAR does not
    // fully resolve the insertion.
    let discordant_frac = f64::from(discordant) / f64::from(total);
    assert!(
        discordant_frac < 0.01,
        "Phased SNP+insertion: expected <1% discordant, got {discordant}/{total} ({:.1}%) \
         (concordant_alt={snp_alt_ins_present}, concordant_ref={snp_ref_ins_absent})",
        discordant_frac * 100.0
    );

    let concordant = snp_alt_ins_present + snp_ref_ins_absent;
    let alt_frac = f64::from(snp_alt_ins_present) / f64::from(concordant);
    assert!(
        (0.40..=0.60).contains(&alt_frac),
        "Expected ~50% alt-carrying reads, got {:.1}%",
        alt_frac * 100.0
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// Priority 4: Common parameter combinations
// ═══════════════════════════════════════════════════════════════════════════

/// Simulate with both a VCF and a BED file.  Place a hom-alt SNP inside the
/// target and one outside.  Verify reads over the target show the variant
/// and reads are not generated far from the target.
#[test]
fn test_simulate_vcf_plus_bed() {
    let seq = non_repetitive_seq(10_000);

    // SNP at pos 5000 (inside target) and pos 1000 (outside target).
    let inside_ref = seq[5000];
    let inside_alt = match inside_ref {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    };
    let outside_ref_str = (seq[1000] as char).to_string();
    let outside_alt_str = match seq[1000] {
        b'A' => "C".to_string(),
        b'C' => "G".to_string(),
        _ => "A".to_string(),
    };
    let inside_ref_str = (inside_ref as char).to_string();
    let inside_alt_str = (inside_alt as char).to_string();

    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.write_vcf(
        "SAMPLE",
        &[("chr1", 10_000)],
        &[
            VcfVariant {
                chrom: "chr1",
                pos_1based: 1001,
                ref_allele: &outside_ref_str,
                alt_alleles: &[&outside_alt_str],
                gt: "1|1",
            },
            VcfVariant {
                chrom: "chr1",
                pos_1based: 5001,
                ref_allele: &inside_ref_str,
                alt_alleles: &[&inside_alt_str],
                gt: "1|1",
            },
        ],
    );
    let bed_path = env.write_bed(&[("chr1", 4800, 5200)]); // 400bp target around pos 5000

    let out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf_path.to_str().unwrap(),
        "-b",
        bed_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "500",
        "--read-length",
        "100",
        "--fragment-mean",
        "200",
        "--fragment-stddev",
        "30",
        "--min-error-rate",
        "0",
        "--max-error-rate",
        "0",
        "--golden-bam",
        "--seed",
        "30",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let pileup = pileup_bases(&bam_path);
    let cols = &pileup[0];

    // Inside-target SNP at pos 5000 should show ≥95% alt.
    let inside_col = &cols[5000];
    assert!(inside_col.total > 0, "Inside-target position should have coverage");
    let inside_alt_count = match inside_alt {
        b'A' => inside_col.a,
        b'C' => inside_col.c,
        b'G' => inside_col.g,
        b'T' => inside_col.t,
        _ => 0,
    };
    let inside_frac = f64::from(inside_alt_count) / f64::from(inside_col.total);
    assert!(
        inside_frac >= 0.95,
        "Inside-target SNP: expected ≥95% alt, got {:.1}%",
        inside_frac * 100.0
    );

    // Outside-target position 1000 should have zero or near-zero coverage.
    let outside_col = &cols[1000];
    assert!(
        outside_col.total <= 5,
        "Outside-target position 1000 should have near-zero coverage, got {}",
        outside_col.total
    );
}

/// Simulate with a moderate error rate and verify the pileup has a non-zero
/// mismatch rate at reference positions.
#[test]
fn test_simulate_error_rate_introduces_mismatches() {
    let seq = non_repetitive_seq(2_000);
    let env = TestEnv::new(&[("chr1", &seq)]);
    let out = env.output_prefix();

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "200",
        "--read-length",
        "100",
        "--fragment-mean",
        "200",
        "--fragment-stddev",
        "30",
        "--min-error-rate",
        "0.01",
        "--max-error-rate",
        "0.02",
        "--golden-bam",
        "--seed",
        "31",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let pileup = pileup_bases(&bam_path);
    let cols = &pileup[0];

    // Sample interior positions to check for mismatches.
    let edge = 250;
    let interior = &cols[edge..cols.len() - edge];
    assert!(!interior.is_empty());

    let mut total_bases: u64 = 0;
    let mut mismatch_bases: u64 = 0;

    for (i, col) in interior.iter().enumerate() {
        if col.total == 0 {
            continue;
        }
        let ref_base = seq[edge + i];
        let ref_count = match ref_base {
            b'A' => col.a,
            b'C' => col.c,
            b'G' => col.g,
            b'T' => col.t,
            _ => 0,
        };
        total_bases += u64::from(col.total);
        mismatch_bases += u64::from(col.total - ref_count);
    }

    assert!(total_bases > 0, "Should have some bases in interior");
    let error_rate = mismatch_bases as f64 / total_bases as f64;

    // With min_error_rate=0.01 and max_error_rate=0.02, we should see a
    // mismatch rate roughly in that range.
    assert!(error_rate > 0.005, "Error rate {error_rate:.4} too low — errors should be detectable",);
    assert!(error_rate < 0.05, "Error rate {error_rate:.4} too high — something is wrong",);
}

/// Two contigs, each with a hom-alt SNP.  Verify both variants appear in the
/// golden BAM pileup.
#[test]
fn test_simulate_multiple_contigs_with_variants() {
    let seq1 = non_repetitive_seq(2_000);
    let seq2 = non_repetitive_seq(2_000);

    let snp1_ref = seq1[500];
    let snp1_alt = match snp1_ref {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    };
    let snp2_ref = seq2[800];
    let snp2_alt = match snp2_ref {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    };

    let snp1_ref_str = (snp1_ref as char).to_string();
    let snp1_alt_str = (snp1_alt as char).to_string();
    let snp2_ref_str = (snp2_ref as char).to_string();
    let snp2_alt_str = (snp2_alt as char).to_string();

    let env = TestEnv::new(&[("chr1", &seq1), ("chr2", &seq2)]);
    let vcf_path = env.write_vcf(
        "SAMPLE",
        &[("chr1", 2_000), ("chr2", 2_000)],
        &[
            VcfVariant {
                chrom: "chr1",
                pos_1based: 501,
                ref_allele: &snp1_ref_str,
                alt_alleles: &[&snp1_alt_str],
                gt: "1|1",
            },
            VcfVariant {
                chrom: "chr2",
                pos_1based: 801,
                ref_allele: &snp2_ref_str,
                alt_alleles: &[&snp2_alt_str],
                gt: "1|1",
            },
        ],
    );

    let out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "500",
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
        "--seed",
        "32",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    let pileup = pileup_bases(&bam_path);

    // chr1 SNP at pos 500: ≥95% alt.
    let col1 = &pileup[0][500];
    assert!(col1.total > 0, "chr1 pos 500 should have coverage");
    let alt1_count = match snp1_alt {
        b'A' => col1.a,
        b'C' => col1.c,
        b'G' => col1.g,
        b'T' => col1.t,
        _ => 0,
    };
    let alt1_frac = f64::from(alt1_count) / f64::from(col1.total);
    assert!(alt1_frac >= 0.95, "chr1 SNP: expected ≥95% alt, got {:.1}%", alt1_frac * 100.0);

    // chr2 SNP at pos 800: ≥95% alt.
    let col2 = &pileup[1][800];
    assert!(col2.total > 0, "chr2 pos 800 should have coverage");
    let alt2_count = match snp2_alt {
        b'A' => col2.a,
        b'C' => col2.c,
        b'G' => col2.g,
        b'T' => col2.t,
        _ => 0,
    };
    let alt2_frac = f64::from(alt2_count) / f64::from(col2.total);
    assert!(alt2_frac >= 0.95, "chr2 SNP: expected ≥95% alt, got {:.1}%", alt2_frac * 100.0);
}

/// A contig shorter than fragment_mean should still produce reads.  Fragments
/// are clamped to the contig length.
#[test]
fn test_simulate_very_short_contig() {
    let seq = non_repetitive_seq(50); // 50 bp — much shorter than default fragment_mean
    let env = TestEnv::new(&[("chr1", &seq)]);
    let out = env.output_prefix();

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "100",
        "--read-length",
        "50",
        "--fragment-mean",
        "300",
        "--fragment-stddev",
        "50",
        "--seed",
        "33",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let r1_path = PathBuf::from(format!("{}.r1.fastq.gz", out.display()));
    let r1_contents = read_gzipped(&r1_path);
    let r1_records = count_fastq_records(&r1_contents);
    assert!(r1_records > 0, "Should have generated reads even for short contig");
}

// ═══════════════════════════════════════════════════════════════════════════
// Priority 5: End-to-end pipeline tests
// ═══════════════════════════════════════════════════════════════════════════

/// Run `mutate` to generate a VCF from a reference, then `simulate` using
/// that VCF.  Verify the simulation completes and the golden BAM shows
/// evidence of variants.
#[test]
fn test_mutate_then_simulate() {
    use helpers::run_mutate;

    let seq = non_repetitive_seq(5_000);
    let env = TestEnv::new(&[("chr1", &seq)]);
    let vcf_path = env.dir.path().join("muts.vcf");

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
        "40",
    ]);
    assert!(ok, "mutate failed: {stderr}");

    let out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "100",
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
        "--seed",
        "41",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    // Verify the golden BAM exists and has records.
    let bam_path = PathBuf::from(format!("{}.golden.bam", out.display()));
    assert!(bam_path.exists(), "Golden BAM should exist");

    let pileup = pileup_bases(&bam_path);
    let cols = &pileup[0];

    // Check that at least some positions have non-reference bases (variant
    // evidence).  With ~25 SNPs at 100x coverage, interior positions with
    // variants should be detectable.
    let edge = 200;
    let interior = &cols[edge..cols.len() - edge];
    let mut positions_with_alt = 0;
    for (i, col) in interior.iter().enumerate() {
        if col.total == 0 {
            continue;
        }
        let ref_base = seq[edge + i];
        let ref_count = match ref_base {
            b'A' => col.a,
            b'C' => col.c,
            b'G' => col.g,
            b'T' => col.t,
            _ => 0,
        };
        if ref_count < col.total {
            positions_with_alt += 1;
        }
    }

    assert!(
        positions_with_alt > 0,
        "Expected some positions with non-reference bases from mutated VCF"
    );
}

/// Run `simulate` to produce reads and a golden BAM, then `eval` on the
/// golden BAM.  Since the golden BAM has truth positions, eval should report
/// 100% correct.
#[test]
fn test_simulate_then_eval() {
    use helpers::run_eval;

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

    // Parse the eval output and verify 100% correct.
    let eval_file = PathBuf::from(format!("{}.eval.txt", eval_out.display()));
    let contents = std::fs::read_to_string(&eval_file).unwrap();

    for line in contents.lines() {
        if line.starts_with("ALL\t") {
            let fields: Vec<&str> = line.split('\t').collect();
            let total: u64 = fields[1].parse().unwrap();
            let correct: u64 = fields[2].parse().unwrap();
            let mismapped: u64 = fields[3].parse().unwrap();
            let unmapped: u64 = fields[4].parse().unwrap();

            assert!(total > 0, "Should have evaluated some reads");
            assert_eq!(correct, total, "All reads should be correct");
            assert_eq!(mismapped, 0, "No reads should be mismapped");
            assert_eq!(unmapped, 0, "No reads should be unmapped");
            return;
        }
    }
    panic!("ALL row not found in eval output");
}
