//! End-to-end tests for the `methylate` subcommand.
#![allow(clippy::cast_precision_loss)]

mod helpers;

use helpers::TestEnv;

#[test]
fn methylate_command_runs_on_empty_reference() {
    // Build a tiny reference with a few CpGs (ACGT repeat gives C in CpG context).
    let seq = b"ACGTACGT".to_vec();
    let env = TestEnv::new(&[("chr1", &seq)]);
    let out_path = env.dir.path().join("methylated.vcf.gz");

    let (ok, _, stderr) = run_methylate(&[
        "methylate",
        "--reference",
        env.fasta_path.to_str().unwrap(),
        "--output",
        out_path.to_str().unwrap(),
        "--methylation-rate-island",
        "1.0",
        "--methylation-rate-shore",
        "1.0",
        "--methylation-rate-open-sea",
        "1.0",
        "--hemimethylation-rate",
        "0",
        "--seed",
        "42",
    ]);
    assert!(ok, "methylate exited with non-zero status: {stderr}");
    assert!(out_path.exists(), "output file was not created");
}

#[test]
fn methylate_produces_loadable_vcf_with_mt_mb() {
    // Two CpGs at top-C positions 1 and 5 on a single contig (ACGTACG has
    // CG at offsets 1–2 and 5–6).
    let seq = b"ACGTACG".to_vec();
    let env = TestEnv::new(&[("chr1", &seq)]);
    let out_path = env.dir.path().join("methylated.vcf.gz");

    let (ok, _, stderr) = run_methylate(&[
        "methylate",
        "--reference",
        env.fasta_path.to_str().unwrap(),
        "--output",
        out_path.to_str().unwrap(),
        "--methylation-rate-island",
        "1.0",
        "--methylation-rate-shore",
        "1.0",
        "--methylation-rate-open-sea",
        "1.0",
        "--hemimethylation-rate",
        "0",
        "--seed",
        "42",
    ]);
    assert!(ok, "methylate exited with non-zero status: {stderr}");

    // Read the BGZF-compressed VCF as raw text and check the header.
    let raw_text = helpers::read_gzipped(&out_path);
    assert!(raw_text.contains("##FORMAT=<ID=MT,"), "missing MT FORMAT header in:\n{raw_text}");
    assert!(raw_text.contains("##FORMAT=<ID=MB,"), "missing MB FORMAT header in:\n{raw_text}");

    // Count data lines (non-header lines).
    let data_lines: Vec<&str> = raw_text.lines().filter(|l| !l.starts_with('#')).collect();
    // ACGTACG has CpGs at ref positions 1 and 5 → expect 2 standalone records.
    assert_eq!(
        data_lines.len(),
        2,
        "expected one record per CpG; got {} records:\n{}",
        data_lines.len(),
        data_lines.join("\n")
    );
}

/// End-to-end pipeline test: `methylate` → `simulate` propagates methylation
/// truth into the golden BAM's `YM:Z` tag.
///
/// Reference: a 1 kb ACGT-repeat contig where every C is in a CpG context.
/// Pipeline:
///   1. `holodeck methylate` with all context rates 1.0 → `meth.vcf.gz`
///      (all CpGs unconditionally methylated on both haplotypes)
///   2. `holodeck simulate --vcf meth.vcf.gz --methylation-mode em-seq
///        --methylation-conversion-rate 0.0 --golden-bam`
///      (VCF has MT/MB → loader path; zero conversion so no C→T noise)
///   3. Parse the golden BAM. For every record, assert that every CpG C in
///      the `YM:Z` truth tag is `'Z'` (methylated) and no non-CpG C appears
///      as methylated.
#[test]
fn methylate_then_simulate_propagates_methylation_to_golden_bam() {
    use noodles::sam::alignment::record::data::field::Tag as DataTag;
    use noodles::sam::alignment::record_buf::data::field::Value as DataValue;

    // 1. Build a 1 kb reference where every C is in CpG context.
    let seq = b"ACGT".repeat(250); // 1000 bp
    let env = TestEnv::new(&[("chr1", &seq)]);
    let meth_vcf = env.dir.path().join("meth.vcf.gz");
    let out = env.dir.path().join("sim");

    // 2. Run `holodeck methylate` to produce a methylated VCF.
    let (ok, _, stderr) = run_methylate(&[
        "methylate",
        "--reference",
        env.fasta_path.to_str().unwrap(),
        "--output",
        meth_vcf.to_str().unwrap(),
        "--methylation-rate-island",
        "1.0",
        "--methylation-rate-shore",
        "1.0",
        "--methylation-rate-open-sea",
        "1.0",
        "--hemimethylation-rate",
        "0",
        "--seed",
        "42",
    ]);
    assert!(ok, "methylate failed: {stderr}");
    assert!(meth_vcf.exists(), "methylate did not produce output VCF");

    // 3. Run `holodeck simulate` with the methylated VCF, using the MT/MB
    //    loading path. Zero conversion rate keeps all CpG C's unconverted so
    //    the YM tag directly reflects the methylation truth. Zero error rate
    //    eliminates sequencing noise.
    let (ok, _, stderr) = helpers::run_simulate(&[
        "simulate",
        "--reference",
        env.fasta_path.to_str().unwrap(),
        "--vcf",
        meth_vcf.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--coverage",
        "5",
        "--read-length",
        "50",
        "--fragment-mean",
        "100",
        "--fragment-stddev",
        "10",
        "--min-error-rate",
        "0",
        "--max-error-rate",
        "0",
        "--methylation-mode",
        "em-seq",
        "--methylation-conversion-rate",
        "0.0",
        "--methylation-failure-rate",
        "0.0",
        "--golden-bam",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    // 4. Open the golden BAM and assert that every CpG in YM is 'Z'.
    let bam_path = out.with_extension("golden.bam");
    let records = helpers::read_bam_records(&bam_path);
    assert!(!records.is_empty(), "golden BAM must have records");

    let ym_tag = DataTag::new(b'Y', b'M');

    for rec in &records {
        let DataValue::String(ym) =
            rec.data().get(&ym_tag).expect("YM:Z tag must be present in golden BAM record")
        else {
            panic!("YM must be a String value");
        };
        let ym_bytes: &[u8] = ym.as_ref();
        let has_z = ym_bytes.contains(&b'Z');
        assert!(has_z, "expected at least one methylated CpG (Z) in YM tag");
        let has_lowercase_z = ym_bytes.contains(&b'z');
        assert!(
            !has_lowercase_z,
            "unexpected unmethylated CpG (z) in YM tag under 100% methylation-rate; YM={}",
            std::str::from_utf8(ym_bytes).unwrap_or("?")
        );
    }
}

#[test]
fn methylate_writes_bedgraph_when_requested() {
    // Reference ACGTACG has two CpGs: top-C at ref positions 1 and 5.
    // With all context rates 1.0 both are fully methylated on both haplotypes.
    let seq = b"ACGTACG".to_vec();
    let env = TestEnv::new(&[("chr1", &seq)]);
    let out_path = env.dir.path().join("meth.vcf.gz");
    let bg_path = env.dir.path().join("meth.bedgraph");

    let (ok, _, stderr) = run_methylate(&[
        "methylate",
        "--reference",
        env.fasta_path.to_str().unwrap(),
        "--output",
        out_path.to_str().unwrap(),
        "--bedgraph",
        bg_path.to_str().unwrap(),
        "--methylation-rate-island",
        "1.0",
        "--methylation-rate-shore",
        "1.0",
        "--methylation-rate-open-sea",
        "1.0",
        "--hemimethylation-rate",
        "0",
        "--seed",
        "42",
    ]);
    assert!(ok, "methylate exited with non-zero status: {stderr}");
    assert!(bg_path.exists(), "bedgraph file was not created");

    let contents = std::fs::read_to_string(&bg_path).unwrap();
    assert!(contents.starts_with("track type="), "missing track header: {contents}");
    // Two CpGs at top-C ref positions 1 and 5, both fully methylated (rate 100).
    assert!(contents.contains("chr1\t1\t2\t100"), "missing first CpG record: {contents}");
    assert!(contents.contains("chr1\t5\t6\t100"), "missing second CpG record: {contents}");
}

/// Parse a population-fraction bedGraph into `(start, rate)` pairs, skipping
/// the track header line.
fn parse_bedgraph(contents: &str) -> Vec<(usize, u32)> {
    contents
        .lines()
        .filter(|l| !l.starts_with("track"))
        .filter_map(|l| {
            let cols: Vec<&str> = l.split('\t').collect();
            Some((cols.get(1)?.parse().ok()?, cols.get(3)?.parse().ok()?))
        })
        .collect()
}

#[test]
fn methylate_default_is_no_longer_fully_methylated() {
    // The no-flags default is now context-aware (open-sea ~0.85), NOT 100%
    // methylated. Over a long open-sea reference the MEAN per-CpG rate should
    // sit near the open-sea target and clearly below 100 — a structural check,
    // not "at least one happens to be unmethylated". Guards against a
    // regression to the old flat-1.0 default. Seed-pinned (measured ~81.5 at
    // this seed; ~81–91 across seeds), so the band is tight enough to catch a
    // flat-default regression yet leaves headroom for the SmallRng-on-rand-0.9
    // caveat; widen if the RNG stream changes.
    let seq = b"ACGT".repeat(25_000); // 100 kb, CpG every 4 bp → all open-sea
    let env = TestEnv::new(&[("chr1", &seq)]);
    let out_path = env.dir.path().join("meth.vcf.gz");
    let bg_path = env.dir.path().join("meth.bedgraph");
    let (ok, _, stderr) = run_methylate(&[
        "methylate",
        "--reference",
        env.fasta_path.to_str().unwrap(),
        "--output",
        out_path.to_str().unwrap(),
        "--bedgraph",
        bg_path.to_str().unwrap(),
        "--seed",
        "7",
    ]);
    assert!(ok, "methylate failed: {stderr}");
    let rates = parse_bedgraph(&std::fs::read_to_string(&bg_path).unwrap());
    assert!(!rates.is_empty(), "expected bedgraph records");
    let mean = rates.iter().map(|&(_, r)| f64::from(r)).sum::<f64>() / rates.len() as f64;
    assert!(
        (70.0..=95.0).contains(&mean),
        "default open-sea mean rate {mean} should be near ~85 and clearly below the old flat 100"
    );
}

#[test]
fn methylate_context_flags_produce_island_hypomethylation() {
    // A CpG island ("CG" repeats) embedded in AT-rich open-sea. With island
    // rate low, open-sea rate high, and short correlation lengths (so each
    // region's mean converges), the island region's bedGraph rates should be
    // markedly lower than the far open-sea rates. Validates the full
    // CLI → context-classification → bedGraph wiring.
    let os_unit = b"AATTCGAATT"; // CG at offset 4, 20% GC → open-sea
    let os_units = 2000usize; // 20 kb flanks
    let island_units = 1500usize; // 3 kb island
    let mut seq = Vec::new();
    for _ in 0..os_units {
        seq.extend_from_slice(os_unit);
    }
    let island_start = seq.len();
    for _ in 0..island_units {
        seq.extend_from_slice(b"CG");
    }
    let island_end = seq.len();
    for _ in 0..os_units {
        seq.extend_from_slice(os_unit);
    }
    let env = TestEnv::new(&[("chr1", &seq)]);
    let out_path = env.dir.path().join("meth.vcf.gz");
    let bg_path = env.dir.path().join("meth.bedgraph");
    let (ok, _, stderr) = run_methylate(&[
        "methylate",
        "--reference",
        env.fasta_path.to_str().unwrap(),
        "--output",
        out_path.to_str().unwrap(),
        "--bedgraph",
        bg_path.to_str().unwrap(),
        "--methylation-rate-island",
        "0.05",
        "--methylation-rate-open-sea",
        "0.9",
        "--methylation-correlation-length-island",
        "20",
        "--methylation-correlation-length-shore",
        "20",
        "--methylation-correlation-length-open-sea",
        "20",
        "--hemimethylation-rate",
        "0",
        "--seed",
        "7",
    ]);
    assert!(ok, "methylate failed: {stderr}");
    let rates = parse_bedgraph(&std::fs::read_to_string(&bg_path).unwrap());

    let mean = |sel: &dyn Fn(usize) -> bool| {
        let vals: Vec<u32> = rates.iter().filter(|&&(s, _)| sel(s)).map(|&(_, r)| r).collect();
        assert!(!vals.is_empty(), "no CpGs selected");
        f64::from(vals.iter().sum::<u32>()) / vals.len() as f64
    };
    let island_mean = mean(&|s| s >= island_start + 200 && s < island_end - 200);
    let open_sea_mean = mean(&|s| s + 6000 < island_start || s > island_end + 6000);
    assert!(island_mean < 40.0, "island region should be hypomethylated, got {island_mean}");
    assert!(open_sea_mean > 60.0, "open-sea region should be hypermethylated, got {open_sea_mean}");
    assert!(island_mean < open_sea_mean, "island must be less methylated than open sea");
}

#[test]
fn methylate_rejects_methylation_rate_above_one() {
    let env = TestEnv::new(&[("chr1", b"ACGT")]);
    let out_path = env.dir.path().join("meth.vcf.gz");
    let (ok, _, stderr) = run_methylate(&[
        "methylate",
        "--reference",
        env.fasta_path.to_str().unwrap(),
        "--output",
        out_path.to_str().unwrap(),
        "--methylation-rate-island",
        "1.5",
    ]);
    assert!(!ok, "methylate should reject rate > 1.0");
    assert!(
        stderr.contains("--methylation-rate-island must be in [0.0, 1.0]"),
        "stderr did not mention rate range: {stderr}"
    );
}

#[test]
fn methylate_rejects_negative_methylation_rate() {
    let env = TestEnv::new(&[("chr1", b"ACGT")]);
    let out_path = env.dir.path().join("meth.vcf.gz");
    let (ok, _, stderr) = run_methylate(&[
        "methylate",
        "--reference",
        env.fasta_path.to_str().unwrap(),
        "--output",
        out_path.to_str().unwrap(),
        // Use `--flag=value` form so clap doesn't parse `-0.1` as a flag.
        "--methylation-rate-open-sea=-0.1",
    ]);
    assert!(!ok, "methylate should reject negative rate");
    assert!(
        stderr.contains("--methylation-rate-open-sea must be in [0.0, 1.0]"),
        "stderr did not mention rate range: {stderr}"
    );
}

#[test]
fn methylate_rejects_nonpositive_correlation_length() {
    let env = TestEnv::new(&[("chr1", b"ACGT")]);
    let out_path = env.dir.path().join("meth.vcf.gz");
    let (ok, _, stderr) = run_methylate(&[
        "methylate",
        "--reference",
        env.fasta_path.to_str().unwrap(),
        "--output",
        out_path.to_str().unwrap(),
        "--methylation-correlation-length-shore",
        "0",
    ]);
    assert!(!ok, "methylate should reject a non-positive correlation length");
    assert!(
        stderr.contains("--methylation-correlation-length-shore must be a finite value > 0"),
        "stderr did not mention correlation-length constraint: {stderr}"
    );
}

#[test]
fn methylate_rejects_sample_without_vcf() {
    let env = TestEnv::new(&[("chr1", b"ACGT")]);
    let out_path = env.dir.path().join("meth.vcf.gz");
    let (ok, _, stderr) = run_methylate(&[
        "methylate",
        "--reference",
        env.fasta_path.to_str().unwrap(),
        "--output",
        out_path.to_str().unwrap(),
        "--sample",
        "NA12878",
    ]);
    assert!(!ok, "methylate should reject --sample without --vcf");
    assert!(
        stderr.contains("--sample requires --vcf"),
        "stderr did not mention --sample requires --vcf: {stderr}"
    );
}

#[test]
fn methylate_rejects_nonfinite_methylation_rate() {
    let env = TestEnv::new(&[("chr1", b"ACGT")]);
    let out_path = env.dir.path().join("meth.vcf.gz");
    let (ok, _, stderr) = run_methylate(&[
        "methylate",
        "--reference",
        env.fasta_path.to_str().unwrap(),
        "--output",
        out_path.to_str().unwrap(),
        "--methylation-rate-island",
        "nan",
    ]);
    assert!(!ok, "methylate should reject NaN rate");
    assert!(
        stderr.contains("--methylation-rate-island must be in [0.0, 1.0]"),
        "stderr did not mention rate range: {stderr}"
    );
}

#[test]
fn methylate_seed_determinism() {
    // Same --seed should produce byte-identical output across runs.
    let seq = b"ACGTACGTACGTACGT".to_vec();
    let env = TestEnv::new(&[("chr1", &seq)]);
    let out_a = env.dir.path().join("a.vcf.gz");
    let out_b = env.dir.path().join("b.vcf.gz");

    for out in [&out_a, &out_b] {
        let (ok, _, stderr) = run_methylate(&[
            "methylate",
            "--reference",
            env.fasta_path.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--methylation-rate-island",
            "0.5",
            "--methylation-rate-shore",
            "0.5",
            "--methylation-rate-open-sea",
            "0.5",
            "--hemimethylation-rate",
            "0",
            "--seed",
            "42",
        ]);
        assert!(ok, "methylate failed: {stderr}");
    }

    // Compare decompressed contents so any BGZF block-boundary differences
    // (which are immaterial to VCF semantics) don't fail the test.
    let bytes_a = decompress_vcf(&out_a);
    let bytes_b = decompress_vcf(&out_b);
    // Strip the ##holodeckCommand line, which embeds the output path
    // (different across the two runs).
    let a = strip_command_line(&bytes_a);
    let b = strip_command_line(&bytes_b);
    assert_eq!(a, b, "same seed produced different VCF bodies");
}

#[test]
fn methylate_applies_variants_from_input_vcf() {
    // Input variant VCF has one SNP. methylate's output should preserve the
    // variant record (with GT) and emit standalone methylation records at
    // every reference CpG outside the variant span.
    let seq = b"ACGTACG".to_vec(); // CpGs at top-C ref pos 1 and 5
    let env = TestEnv::new(&[("chr1", &seq)]);
    let variants = vec![helpers::VcfVariant {
        chrom: "chr1",
        pos_1based: 3, // ref T -> alt A; no CpG impact
        ref_allele: "T",
        alt_alleles: &["A"],
        gt: "1|0",
    }];
    let vcf_in = env.write_vcf("SAMPLE", &[("chr1", seq.len())], &variants);
    let out_path = env.dir.path().join("meth.vcf.gz");

    let (ok, _, stderr) = run_methylate(&[
        "methylate",
        "--reference",
        env.fasta_path.to_str().unwrap(),
        "--vcf",
        vcf_in.to_str().unwrap(),
        "--output",
        out_path.to_str().unwrap(),
        "--methylation-rate-island",
        "1.0",
        "--methylation-rate-shore",
        "1.0",
        "--methylation-rate-open-sea",
        "1.0",
        "--hemimethylation-rate",
        "0",
        "--seed",
        "42",
    ]);
    assert!(ok, "methylate failed: {stderr}");

    let text = decompress_vcf(&out_path);
    // The variant record should appear with GT:MT:MB.
    assert!(
        text.lines().any(|l| l.starts_with("chr1\t3\t.\tT\tA") && l.contains("GT:MT:MB")),
        "variant record not preserved in methylate output:\n{text}"
    );
    // Standalone CpG records should appear at POS 2 (top-C pos 1, 0-based -> 2, 1-based)
    // and POS 6 (top-C pos 5 -> 6, 1-based), with MT:MB (no GT).
    assert!(
        text.lines().any(|l| l.starts_with("chr1\t2\t.\tC\t.") && l.contains("MT:MB")),
        "standalone CpG at POS 2 missing:\n{text}"
    );
    assert!(
        text.lines().any(|l| l.starts_with("chr1\t6\t.\tC\t.") && l.contains("MT:MB")),
        "standalone CpG at POS 6 missing:\n{text}"
    );
}

#[test]
fn methylate_output_is_identical_across_thread_counts() {
    // Several contigs of varied length so the per-contig parallel map spans
    // multiple work units; non-repetitive sequence carries CpGs on each. The
    // output must be byte-identical whether methylated on one thread or many:
    // each contig's RNGs are seeded purely from `(seed, contig)` with no state
    // carried between contigs, so processing order cannot affect the result.
    let c1 = helpers::non_repetitive_seq(8000);
    let c2 = helpers::non_repetitive_seq(5000);
    let c3 = helpers::non_repetitive_seq(3000);
    let c4 = helpers::non_repetitive_seq(1500);
    let env = TestEnv::new(&[("chr1", &c1), ("chr2", &c2), ("chr3", &c3), ("chr4", &c4)]);

    let run = |tag: &str, threads: usize| -> (String, Vec<u8>) {
        let vcf = env.dir.path().join(format!("{tag}.vcf.gz"));
        let bg = env.dir.path().join(format!("{tag}.bedgraph"));
        let (ok, _, stderr) = run_methylate_with_threads(
            threads,
            &[
                "methylate",
                "--reference",
                env.fasta_path.to_str().unwrap(),
                "--output",
                vcf.to_str().unwrap(),
                "--bedgraph",
                bg.to_str().unwrap(),
                "--seed",
                "42",
            ],
        );
        assert!(ok, "methylate (threads={threads}) failed: {stderr}");
        (strip_command_line(&decompress_vcf(&vcf)), std::fs::read(&bg).expect("read bedgraph"))
    };

    let (vcf_single, bg_single) = run("single", 1);
    let (vcf_multi, bg_multi) = run("multi", 8);

    assert_eq!(vcf_single, vcf_multi, "VCF content differs between 1 and 8 threads");
    assert_eq!(bg_single, bg_multi, "bedGraph differs between 1 and 8 threads");
    // Guard against a vacuous pass on two equal-but-empty outputs.
    assert!(
        vcf_single.lines().any(|l| l.starts_with("chr1\t") && l.contains("MT:MB")),
        "expected CpG methylation rows in the output"
    );
}

// ── Command runner ───────────────────────────────────────────────────────────

/// Run `holodeck methylate` with the given arguments and return
/// `(success, stdout, stderr)`.
fn run_methylate(args: &[&str]) -> (bool, String, String) {
    let output = std::process::Command::new(env!("CARGO_BIN_EXE_holodeck"))
        .args(args)
        .output()
        .expect("Failed to run holodeck");
    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();
    (output.status.success(), stdout, stderr)
}

/// Run `holodeck methylate` with `RAYON_NUM_THREADS` pinned, returning
/// `(success, stdout, stderr)`. Used to prove per-contig parallelism does not
/// change the output.
fn run_methylate_with_threads(threads: usize, args: &[&str]) -> (bool, String, String) {
    let output = std::process::Command::new(env!("CARGO_BIN_EXE_holodeck"))
        .env("RAYON_NUM_THREADS", threads.to_string())
        .args(args)
        .output()
        .expect("Failed to run holodeck");
    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();
    (output.status.success(), stdout, stderr)
}

/// Decompress a BGZF-compressed VCF into a String. The methylate command
/// always writes its output as BGZF; tests that need to inspect the text
/// use this helper.
fn decompress_vcf(path: &std::path::Path) -> String {
    use std::io::Read as _;
    let file = std::fs::File::open(path).expect("open VCF");
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(file).read_to_string(&mut text).expect("decompress VCF");
    text
}

/// Remove the `##holodeckCommand=...` line so that two runs that differ only
/// in their output path can be compared for content equivalence.
fn strip_command_line(text: &str) -> String {
    text.lines().filter(|l| !l.starts_with("##holodeckCommand=")).collect::<Vec<_>>().join("\n")
}
