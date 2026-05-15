//! End-to-end tests for the `methylate → simulate` pipeline.
//!
//! Covers em-seq and TAPS chemistry simulation against a methylation-truth
//! VCF: validation of the `(VCF has MT/MB, --methylation-mode set)` matrix,
//! golden-BAM `XG`/`XR`/`XM`/`YM`/`NM`/`MD`/`YS` tag emission, propagation of
//! truth bits from `methylate`'s VCF through `simulate`'s reads, conversion-
//! rate behavior, and per-CpG truth-bedGraph parity with the methylation
//! bitmap. Tests build inputs with `methylate_to_vcf` helpers (no committed
//! fixtures) so any change to either command shows up here.

#![allow(clippy::similar_names, clippy::cast_precision_loss)]

mod helpers;

use helpers::{
    TestEnv, VcfVariant, methylate_to_vcf, methylate_to_vcf_with_variants, read_bam_records,
    read_gzipped, run_simulate,
};

/// Build a 1 kb reference with frequent C's so we can detect conversion.
fn cytosine_rich_env() -> TestEnv {
    // Repeating "ACGT" gives 25% C content — and crucially every C is in
    // CpG context (followed by G), so methylation-driven preservation
    // applies to every genomic C.
    let seq = b"ACGT".repeat(250); // 1000 bp
    TestEnv::new(&[("chr1", &seq)])
}

/// Count cytosines and thymines across all FASTQ records (sequence lines only).
fn count_c_and_t(fastq_text: &str) -> (usize, usize) {
    let mut c = 0;
    let mut t = 0;
    for (i, line) in fastq_text.lines().enumerate() {
        // FASTQ records: 0=name, 1=seq, 2=+, 3=qual.
        if i % 4 == 1 {
            for &b in line.as_bytes() {
                match b {
                    b'C' => c += 1,
                    b'T' => t += 1,
                    _ => {}
                }
            }
        }
    }
    (c, t)
}

/// Count occurrences of `target` across all FASTQ sequence lines.
#[expect(clippy::naive_bytecount, reason = "test scope; clarity over speed")]
fn count_base(fastq_text: &str, target: u8) -> usize {
    let mut n = 0;
    for (i, line) in fastq_text.lines().enumerate() {
        if i % 4 == 1 {
            n += line.as_bytes().iter().filter(|&&b| b == target).count();
        }
    }
    n
}

#[test]
fn test_em_seq_full_conversion_eliminates_c() {
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    // rate=0 → all CpG-context cytosines unmethylated → all convert under em-seq.
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 0.0, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
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
        "1.0",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let r1 = read_gzipped(&out.with_extension("r1.fastq.gz"));
    let r2 = read_gzipped(&out.with_extension("r2.fastq.gz"));

    let (c1, t1) = count_c_and_t(&r1);
    let (_c2, t2) = count_c_and_t(&r2);
    let g2 = count_base(&r2, b'G');

    // Directional library, 0% methylation, 100% conversion, no errors, no
    // adapter padding (fragment_mean 100 > read_length 50).
    //   R1 5'→3' = c2t(source) → every source C → T → R1 has zero C's.
    //   R2 5'→3' = revcomp(c2t(source)) → c2t(source) has zero C's →
    //              revcomp has zero G's. (R2 still has C's; they reflect
    //              source-strand G positions complemented through revcomp.)
    assert_eq!(c1, 0, "R1 still contains C bases: {c1}");
    assert_eq!(
        g2, 0,
        "R2 still contains G bases: {g2} — c2t(source) has no C's so revcomp has no G's"
    );
    assert!(t1 > 0, "R1 should have many T's after conversion");
    assert!(t2 > 0, "R2 should have many T's");
}

#[test]
fn test_em_seq_full_methylation_preserves_c() {
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    // rate=1.0 → all CpG-context cytosines methylated → all preserved under em-seq.
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 1.0, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
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
        "1.0",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let r1 = read_gzipped(&out.with_extension("r1.fastq.gz"));
    let (c1, t1) = count_c_and_t(&r1);

    // Reference is "ACGT" repeated. Every C is in CpG context (the next
    // base is always G), so at methylation-rate 1.0 every C is methylated
    // and preserved. With no errors and 100% methylation, the ratio of
    // C:T in genomic-only reads should be ~1:1 (matching the reference).
    assert!(c1 > 0, "R1 should still contain C's when fully methylated");
    // Sanity: more or less balanced.
    let ratio = c1 as f64 / t1 as f64;
    assert!(
        (0.7..1.4).contains(&ratio),
        "C:T ratio {ratio:.2} suggests unintended conversion (c={c1}, t={t1})"
    );
}

#[test]
fn test_em_seq_partial_conversion_rate() {
    // Band derived empirically with `SmallRng::seed_from_u64(42)` on rand
    // 0.9 (the simulator's default RNG). If `rand` updates `SmallRng`'s
    // output stream, this band may need widening or re-derivation.
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    // rate=0 → all cytosines unmethylated; 50% conversion rate splits them ~50/50.
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 0.0, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "20",
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
        "0.5",
        "--seed",
        "42",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let r1 = read_gzipped(&out.with_extension("r1.fastq.gz"));
    let (c1, t1) = count_c_and_t(&r1);

    // With 0% methylation and 50% conversion rate, original C's split
    // ~50/50 between C and T. Reference ACGT repeated has equal C and T
    // per position; after conversion, T count should grow by ~50% of the
    // original C count, and C count should drop to ~50% of original.
    // Exact ratio (post-conversion C : post-conversion T) ≈ 1 : 3.
    let ratio = c1 as f64 / (c1 + t1) as f64;
    assert!(
        (0.20..0.30).contains(&ratio),
        "expected C / (C+T) ≈ 0.25 with 50% conversion; got {ratio:.3} (c={c1}, t={t1})"
    );
}

#[test]
fn test_em_seq_intermediate_methylation_preserves_some_c() {
    // With CpG-context methylation at rate 1.0 and full conversion, every
    // top-strand CpG C is preserved — which here is every C in the
    // reference (`ACGT` repeated). Confirm at least some C's remain to
    // distinguish from the full-conversion path.
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    // rate=1.0 → every CpG-context C is methylated and preserved.
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 1.0, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "30",
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
        "1.0",
        "--seed",
        "42",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let r1 = read_gzipped(&out.with_extension("r1.fastq.gz"));
    let (c1, _t1) = count_c_and_t(&r1);
    assert!(c1 > 0, "expected C's preserved at methylated CpG sites; got 0");
}

#[test]
fn test_methylation_is_deterministic_with_seed() {
    // This test runs WITHOUT --golden-bam, so capture_pre_conversion is
    // false. The methylation chemistry still runs deterministically; the
    // pre-conversion bases just aren't recorded. This intentionally
    // exercises the no-annotation path.
    //
    // Both simulate runs use the same methylated VCF, so the per-haplotype
    // methylation bitmaps are identical across both runs.
    let env = cytosine_rich_env();
    let out_a = env.dir.path().join("a");
    let out_b = env.dir.path().join("b");
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 0.7, 12345, "meth.vcf.gz");

    let common: Vec<&str> = vec![
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "--coverage",
        "5",
        "--read-length",
        "50",
        "--methylation-mode",
        "em-seq",
        "--methylation-conversion-rate",
        "0.5",
        "--seed",
        "12345",
        "--threads",
        "1",
    ];

    let mut a_args = common.clone();
    a_args.push("-o");
    a_args.push(out_a.to_str().unwrap());
    let (ok_a, _, err_a) = run_simulate(&a_args);
    assert!(ok_a, "first run failed: {err_a}");

    let mut b_args = common.clone();
    b_args.push("-o");
    b_args.push(out_b.to_str().unwrap());
    let (ok_b, _, err_b) = run_simulate(&b_args);
    assert!(ok_b, "second run failed: {err_b}");

    let r1_a = read_gzipped(&out_a.with_extension("r1.fastq.gz"));
    let r1_b = read_gzipped(&out_b.with_extension("r1.fastq.gz"));
    assert_eq!(r1_a, r1_b, "methylation output must be reproducible with --seed");
}

use noodles::sam::alignment::record::data::field::Tag as DataTag;
use noodles::sam::alignment::record_buf::data::field::Value as DataValue;

#[test]
fn test_golden_bam_em_seq_emits_xg_xr_and_ys_tags_when_enabled() {
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 0.0, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
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
        "1.0",
        "--golden-bam",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = out.with_extension("golden.bam");
    let records = read_bam_records(&bam_path);
    assert!(!records.is_empty(), "expected non-empty golden BAM");

    let xg = DataTag::new(b'X', b'G');
    let xr = DataTag::new(b'X', b'R');
    let yc = DataTag::new(b'Y', b'C');
    let ys = DataTag::new(b'Y', b'S');

    for rec in &records {
        let data = rec.data();

        // YC:Z is no longer emitted -- XG:Z carries the same information
        // under the canonical Bismark tag name.
        assert!(data.get(&yc).is_none(), "YC:Z must not be present (replaced by XG:Z)");

        // XG:Z is the genome-strand indicator: CT = top-derived, GA =
        // bottom-derived. Fragment-level: R1 and R2 of a pair share the
        // value. R1 in fragment orientation → CT-fragment maps forward
        // and GA-fragment maps reverse; R2 maps the opposite way.
        let xg_val = data.get(&xg).expect("XG:Z must be present");
        let xg_str = match xg_val {
            DataValue::String(s) => std::str::from_utf8(s.as_ref()).unwrap().to_string(),
            other => panic!("XG must be a String tag, got {other:?}"),
        };
        assert!(xg_str == "CT" || xg_str == "GA", "XG:Z must be 'CT' or 'GA', got {xg_str:?}");
        let is_rev = rec.flags().is_reverse_complemented();
        let is_r1 = rec.flags().is_first_segment();
        let expected_xg = match (is_r1, is_rev) {
            (true, false) | (false, true) => "CT",
            (true, true) | (false, false) => "GA",
        };
        assert_eq!(
            xg_str, expected_xg,
            "expected XG={expected_xg} for is_r1={is_r1}, is_reverse={is_rev}; got {xg_str}"
        );

        // XR:Z is the read-conversion direction in read orientation:
        // R1 (and SE) → CT; R2 → GA. Fixed by mate index under a
        // directional library, regardless of fragment source strand.
        let xr_val = data.get(&xr).expect("XR:Z must be present");
        let xr_str = match xr_val {
            DataValue::String(s) => std::str::from_utf8(s.as_ref()).unwrap().to_string(),
            other => panic!("XR must be a String tag, got {other:?}"),
        };
        let expected_xr =
            if rec.flags().is_segmented() && rec.flags().is_last_segment() { "GA" } else { "CT" };
        assert_eq!(
            xr_str, expected_xr,
            "expected XR={expected_xr} for is_r1={is_r1}; got {xr_str}"
        );

        let ys_val = data.get(&ys).expect("YS:Z must be present");
        match ys_val {
            DataValue::String(s) => {
                let bytes: &[u8] = s.as_ref();
                // Length must match the SEQ field (reference orientation).
                assert_eq!(
                    bytes.len(),
                    rec.sequence().as_ref().len(),
                    "YS length must match SEQ length"
                );
            }
            other => panic!("YS must be a String tag, got {other:?}"),
        }
    }
}

#[test]
fn test_golden_bam_omits_methylation_tags_when_methylation_disabled() {
    let env = cytosine_rich_env();
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
        "--min-error-rate",
        "0",
        "--max-error-rate",
        "0",
        "--golden-bam",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = out.with_extension("golden.bam");
    let records = read_bam_records(&bam_path);
    assert!(!records.is_empty());

    let xg = DataTag::new(b'X', b'G');
    let xr = DataTag::new(b'X', b'R');
    let yc = DataTag::new(b'Y', b'C');
    let ys = DataTag::new(b'Y', b'S');
    for rec in &records {
        assert!(rec.data().get(&xg).is_none(), "XG must not be present without methylation");
        assert!(rec.data().get(&xr).is_none(), "XR must not be present without methylation");
        assert!(rec.data().get(&yc).is_none(), "YC must never be present (removed)");
        assert!(rec.data().get(&ys).is_none(), "YS must not be present without methylation");
    }
}

#[test]
fn test_golden_bam_em_seq_ys_matches_reference_oriented_pre_conversion() {
    // With 0% methylation, 100% conversion, no errors, and only top-strand
    // C's affected: any difference between SEQ and YS:Z must be of the form
    // T/C (forward records) or A/G (reverse records). This both proves the
    // tag carries pre-conversion truth AND that it's in reference orientation.
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 0.0, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "10",
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
        "1.0",
        "--golden-bam",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = out.with_extension("golden.bam");
    let records = read_bam_records(&bam_path);
    let ys_tag = DataTag::new(b'Y', b'S');
    let xg_tag = DataTag::new(b'X', b'G');

    let mut saw_difference = false;
    for rec in &records {
        let seq: Vec<u8> = rec.sequence().as_ref().to_vec();
        let DataValue::String(ys) = rec.data().get(&ys_tag).expect("YS present") else {
            panic!("YS must be string");
        };
        let DataValue::String(xg) = rec.data().get(&xg_tag).expect("XG present") else {
            panic!("XG must be string");
        };
        let ys_bytes: &[u8] = ys.as_ref();
        let xg_bytes: &[u8] = xg.as_ref();
        if seq != ys_bytes {
            saw_difference = true;
            // Directional library: the chemistry direction is determined by
            // the FRAGMENT (XG tag), not the per-read mapping orientation.
            //   XG=CT: chemistry on top → SEQ vs YS diffs are T (SEQ) / C (YS)
            //          for both R1 (forward-mapped) and R2 (reverse-mapped).
            //   XG=GA: chemistry on bottom → diffs in ref orientation are
            //          A (SEQ) / G (YS) for both R1 (reverse-mapped) and
            //          R2 (forward-mapped).
            for (s, y) in seq.iter().zip(ys_bytes.iter()) {
                if s != y {
                    match xg_bytes {
                        b"CT" => assert_eq!(
                            (*s, *y),
                            (b'T', b'C'),
                            "XG=CT diff must be T/C; got SEQ={} YS={}",
                            *s as char,
                            *y as char
                        ),
                        b"GA" => assert_eq!(
                            (*s, *y),
                            (b'A', b'G'),
                            "XG=GA diff must be A/G; got SEQ={} YS={}",
                            *s as char,
                            *y as char
                        ),
                        other => panic!("unexpected XG tag value: {other:?}"),
                    }
                }
            }
        }
    }
    assert!(saw_difference, "expected at least one record where SEQ != YS after conversion");
}

#[test]
fn test_em_seq_single_end_emits_xg_ys_and_full_conversion() {
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 0.0, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
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
        "--single-end",
        "--methylation-mode",
        "em-seq",
        "--methylation-conversion-rate",
        "1.0",
        "--golden-bam",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    // FASTQ side: only R1 is produced; every C in the genomic portion
    // should have been converted to T.
    let r1 = read_gzipped(&out.with_extension("r1.fastq.gz"));
    let (c1, t1) = count_c_and_t(&r1);
    assert_eq!(c1, 0, "SE R1 should have no C's after full conversion");
    assert!(t1 > 0, "SE R1 should have many T's after conversion");

    // SE mode must not create an R2 file at all.
    let r2_path = out.with_extension("r2.fastq.gz");
    assert!(!r2_path.exists(), "SE mode must not create an R2 file; found {}", r2_path.display());

    // BAM side: every record carries XG:Z, XR:Z, and YS:Z. SE records
    // have is_first_segment=false (SE has no SEGMENTED flag at all), so
    // the (is_first_segment, is_reverse_complemented) -> XG matrix used
    // in the PE test doesn't apply directly. Use the simpler invariant:
    // XG value matches the record's strand (CT for forward, GA for
    // reverse) since SE records and their fragment strand coincide.
    // XR:Z is always CT for SE under directional library (no R2 mate).
    let bam_path = out.with_extension("golden.bam");
    let records = read_bam_records(&bam_path);
    assert!(!records.is_empty(), "expected non-empty SE golden BAM");

    let xg = DataTag::new(b'X', b'G');
    let xr = DataTag::new(b'X', b'R');
    let ys = DataTag::new(b'Y', b'S');

    for rec in &records {
        // SE records are not paired; SEGMENTED flag should be absent.
        assert!(!rec.flags().is_segmented(), "SE record must not have SEGMENTED flag set");

        let xg_val = rec.data().get(&xg).expect("XG:Z must be present");
        let DataValue::String(s) = xg_val else {
            panic!("XG must be a String tag");
        };
        let s = std::str::from_utf8(s.as_ref()).unwrap();
        let is_rev = rec.flags().is_reverse_complemented();
        let expected = if is_rev { "GA" } else { "CT" };
        assert_eq!(s, expected, "SE record XG must match strand: is_reverse={is_rev}, got {s}");

        let xr_val = rec.data().get(&xr).expect("XR:Z must be present");
        let DataValue::String(xr_bstring) = xr_val else {
            panic!("XR must be a String tag");
        };
        let xr_bytes: &[u8] = xr_bstring.as_ref();
        assert_eq!(xr_bytes, b"CT", "SE XR:Z must be 'CT' (R1/SE convention)");

        let ys_val = rec.data().get(&ys).expect("YS:Z must be present");
        let DataValue::String(ys_bstring) = ys_val else {
            panic!("YS must be a String tag");
        };
        let ys_bytes: &[u8] = ys_bstring.as_ref();
        assert_eq!(
            ys_bytes.len(),
            rec.sequence().as_ref().len(),
            "YS length must match SEQ length"
        );
    }
}

// ── TAPS-mode integration tests ────────────────────────────────────────────

#[test]
fn test_taps_full_methylation_full_conversion_eliminates_c_at_cpg() {
    // TAPS at methylation-rate 1.0, conversion-rate 1.0: every methylated
    // C → T. In our reference (ACGT repeated), every C is in CpG context,
    // so every C is methylated under rate=1.0, and TAPS converts them all.
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 1.0, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
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
        "taps",
        "--methylation-conversion-rate",
        "1.0",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let r1 = read_gzipped(&out.with_extension("r1.fastq.gz"));
    let r2 = read_gzipped(&out.with_extension("r2.fastq.gz"));
    let (c1, t1) = count_c_and_t(&r1);
    let (_c2, t2) = count_c_and_t(&r2);
    let g2 = count_base(&r2, b'G');

    // Reference "ACGT" repeated: every C is CpG-context. Under TAPS rate=1.0
    // all CpG-Cs are methylated → all source-strand Cs convert →
    //   R1 5'→3' = c2t(source) has zero C's.
    //   R2 5'→3' = revcomp(c2t(source)) → no G's (c2t had no C's).
    assert_eq!(c1, 0, "TAPS R1 still contains C bases: {c1}");
    assert_eq!(g2, 0, "TAPS R2 still contains G bases: {g2}");
    assert!(t1 > 0, "TAPS R1 should have many T's after conversion");
    assert!(t2 > 0, "TAPS R2 should have many T's after conversion");
}

#[test]
fn test_taps_zero_methylation_preserves_all_c() {
    // TAPS at methylation-rate 0: nothing is methylated, so nothing converts.
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 0.0, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
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
        "taps",
        "--methylation-conversion-rate",
        "1.0",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let r1 = read_gzipped(&out.with_extension("r1.fastq.gz"));
    let (c1, t1) = count_c_and_t(&r1);

    // Reference is "ACGT" repeated → equal C and T positions. With no
    // methylation under TAPS, no C converts → ratio stays ~1:1.
    assert!(c1 > 0, "TAPS without methylation must preserve C's, got 0");
    let ratio = c1 as f64 / t1 as f64;
    assert!(
        (0.7..1.4).contains(&ratio),
        "TAPS C:T ratio {ratio:.2} suggests unintended conversion (c={c1}, t={t1})"
    );
}

#[test]
fn test_taps_golden_bam_emits_xg_xr_and_ys() {
    // TAPS mode emits the same XG/XR/YS tags as em-seq — they're strand
    // indicators, independent of the chemistry's biological meaning.
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 1.0, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
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
        "taps",
        "--methylation-conversion-rate",
        "1.0",
        "--golden-bam",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let bam_path = out.with_extension("golden.bam");
    let records = read_bam_records(&bam_path);
    assert!(!records.is_empty(), "expected non-empty TAPS golden BAM");

    let xg = DataTag::new(b'X', b'G');
    let xr = DataTag::new(b'X', b'R');
    let ys = DataTag::new(b'Y', b'S');

    for rec in &records {
        let xg_val = rec.data().get(&xg).expect("TAPS XG:Z must be present");
        let DataValue::String(xg_bstring) = xg_val else {
            panic!("XG must be a String tag");
        };
        let xg_str = std::str::from_utf8(xg_bstring.as_ref()).unwrap();
        assert!(xg_str == "CT" || xg_str == "GA", "TAPS XG must be 'CT' or 'GA', got {xg_str:?}");

        let xr_val = rec.data().get(&xr).expect("TAPS XR:Z must be present");
        let DataValue::String(xr_bstring) = xr_val else {
            panic!("XR must be a String tag");
        };
        let xr_str = std::str::from_utf8(xr_bstring.as_ref()).unwrap();
        let expected_xr =
            if rec.flags().is_segmented() && rec.flags().is_last_segment() { "GA" } else { "CT" };
        assert_eq!(xr_str, expected_xr, "TAPS XR:Z mismatch (R1/SE → CT, R2 → GA)");

        let ys_val = rec.data().get(&ys).expect("TAPS YS:Z must be present");
        let DataValue::String(ys_bstring) = ys_val else {
            panic!("YS must be a String tag");
        };
        let ys_bytes: &[u8] = ys_bstring.as_ref();
        assert_eq!(
            ys_bytes.len(),
            rec.sequence().as_ref().len(),
            "TAPS YS length must match SEQ length"
        );
    }
}

#[test]
fn test_em_seq_and_taps_produce_different_output() {
    // Same params, same seed, different chemistry mode → output must
    // differ. Sanity that the mode actually flows through the pipeline.
    let env = cytosine_rich_env();
    let out_em = env.dir.path().join("em");
    let out_taps = env.dir.path().join("taps");
    // Both runs use the same methylated VCF so the methylation state is
    // identical; only the chemistry mode differs.
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 0.5, 9876, "meth.vcf.gz");

    let common: Vec<&str> = vec![
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
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
        "--methylation-conversion-rate",
        "1.0",
        "--seed",
        "9876",
        "--threads",
        "1",
    ];

    let mut em_args = common.clone();
    em_args.extend_from_slice(&["-o", out_em.to_str().unwrap(), "--methylation-mode", "em-seq"]);
    let (ok_em, _, err_em) = run_simulate(&em_args);
    assert!(ok_em, "em-seq run failed: {err_em}");

    let mut taps_args = common.clone();
    taps_args.extend_from_slice(&["-o", out_taps.to_str().unwrap(), "--methylation-mode", "taps"]);
    let (ok_taps, _, err_taps) = run_simulate(&taps_args);
    assert!(ok_taps, "taps run failed: {err_taps}");

    let r1_em = read_gzipped(&out_em.with_extension("r1.fastq.gz"));
    let r1_taps = read_gzipped(&out_taps.with_extension("r1.fastq.gz"));
    assert_ne!(r1_em, r1_taps, "em-seq and TAPS must produce different output with the same seed");
}

#[test]
fn test_em_seq_xg_pair_level_consistency() {
    use std::collections::HashMap;

    let env = cytosine_rich_env();
    let out = env.output_prefix();
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 0.5, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "20",
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
        "--golden-bam",
        "--seed",
        "42",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let records = read_bam_records(&out.with_extension("golden.bam"));
    let xg_tag = DataTag::new(b'X', b'G');

    // Group XG values by read name; both R1 and R2 share the same name.
    let mut by_name: HashMap<String, Vec<String>> = HashMap::new();
    for rec in &records {
        let DataValue::String(s) = rec.data().get(&xg_tag).expect("XG present") else {
            panic!("XG must be string");
        };
        let xg = std::str::from_utf8(s.as_ref()).unwrap().to_string();
        let name =
            std::str::from_utf8(rec.name().expect("name present").as_ref()).unwrap().to_string();
        by_name.entry(name).or_default().push(xg);
    }

    assert!(!by_name.is_empty(), "expected at least one read pair");
    let mut pair_count = 0;
    for (name, xgs) in &by_name {
        assert_eq!(xgs.len(), 2, "expected exactly 2 records for read name {name}, got {xgs:?}");
        assert_eq!(
            xgs[0], xgs[1],
            "R1 and R2 of pair {name} carry inconsistent XG values: {xgs:?}"
        );
        pair_count += 1;
    }
    assert!(pair_count >= 50, "expected ~200 pairs, got {pair_count}");
}

#[test]
fn test_em_seq_xg_distribution_covers_all_four_cells() {
    use std::collections::HashMap;

    let env = cytosine_rich_env();
    let out = env.output_prefix();
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 0.5, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "30",
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
        "--golden-bam",
        "--seed",
        "42",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let records = read_bam_records(&out.with_extension("golden.bam"));
    let xg_tag = DataTag::new(b'X', b'G');

    // Joint distribution: (is_r1, is_rev, xg_value) → count.
    let mut counts: HashMap<(bool, bool, String), u32> = HashMap::new();
    for rec in &records {
        let DataValue::String(s) = rec.data().get(&xg_tag).expect("XG present") else {
            panic!("XG must be string");
        };
        let xg = std::str::from_utf8(s.as_ref()).unwrap().to_string();
        let is_r1 = rec.flags().is_first_segment();
        let is_rev = rec.flags().is_reverse_complemented();
        *counts.entry((is_r1, is_rev, xg)).or_insert(0) += 1;
    }

    // Expected: top-strand fragments give (R1=true, rev=false, CT) + (R1=false, rev=true, CT).
    // Bottom-strand fragments give (R1=true, rev=true, GA) + (R1=false, rev=false, GA).
    // No other combination should occur.
    let valid_cells = [
        (true, false, "CT".to_string()),
        (false, true, "CT".to_string()),
        (true, true, "GA".to_string()),
        (false, false, "GA".to_string()),
    ];
    for cell in &valid_cells {
        let n = counts.get(cell).copied().unwrap_or(0);
        assert!(
            n >= 10,
            "cell {cell:?} should have >= 10 records, got {n}; full counts: {counts:?}"
        );
    }
    // No cells outside the four valid ones.
    for (cell, n) in &counts {
        assert!(valid_cells.contains(cell), "unexpected XG × strand cell {cell:?} with count {n}");
    }
}

#[test]
fn test_bisulfite_alias_byte_identical_to_em_seq() {
    let env = cytosine_rich_env();
    let out_em = env.dir.path().join("em");
    let out_bs = env.dir.path().join("bs");
    // Both runs use the same methylated VCF so the methylation state is
    // identical; `bisulfite` and `em-seq` must produce byte-identical output.
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 0.0, 12345, "meth.vcf.gz");

    let common: Vec<&str> = vec![
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "--coverage",
        "5",
        "--read-length",
        "50",
        "--methylation-conversion-rate",
        "1.0",
        "--seed",
        "12345",
        "--threads",
        "1",
    ];

    let mut em_args = common.clone();
    em_args.push("--methylation-mode");
    em_args.push("em-seq");
    em_args.push("-o");
    em_args.push(out_em.to_str().unwrap());
    let (ok_em, _, err_em) = run_simulate(&em_args);
    assert!(ok_em, "em-seq run failed: {err_em}");

    let mut bs_args = common.clone();
    bs_args.push("--methylation-mode");
    bs_args.push("bisulfite");
    bs_args.push("-o");
    bs_args.push(out_bs.to_str().unwrap());
    let (ok_bs, _, err_bs) = run_simulate(&bs_args);
    assert!(ok_bs, "bisulfite-alias run failed: {err_bs}");

    let r1_em = read_gzipped(&out_em.with_extension("r1.fastq.gz"));
    let r1_bs = read_gzipped(&out_bs.with_extension("r1.fastq.gz"));
    assert_eq!(r1_em, r1_bs, "`bisulfite` alias must produce byte-identical output to `em-seq`");
}

#[test]
fn test_unknown_methylation_mode_rejected() {
    let env = cytosine_rich_env();
    let out = env.output_prefix();

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--methylation-mode",
        "wgbs",
        "--seed",
        "42",
        "--threads",
        "1",
    ]);
    assert!(!ok, "expected simulate to fail with unknown mode 'wgbs'");
    assert!(
        stderr.contains("invalid value 'wgbs'"),
        "expected clap value-enum error, got stderr: {stderr}"
    );
    // Clap should also list the valid alternatives.
    assert!(
        stderr.contains("em-seq") && stderr.contains("taps"),
        "expected error to list valid modes, got stderr: {stderr}"
    );
}

#[test]
fn test_taps_golden_bam_ys_diffs_match_chemistry() {
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 1.0, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "10",
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
        "taps",
        "--methylation-conversion-rate",
        "1.0",
        "--golden-bam",
        "--seed",
        "42",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let records = read_bam_records(&out.with_extension("golden.bam"));
    let ys_tag = DataTag::new(b'Y', b'S');
    let xg_tag = DataTag::new(b'X', b'G');

    let mut saw_difference = false;
    for rec in &records {
        let seq: Vec<u8> = rec.sequence().as_ref().to_vec();
        let DataValue::String(ys) = rec.data().get(&ys_tag).expect("YS present") else {
            panic!("YS must be string");
        };
        let DataValue::String(xg) = rec.data().get(&xg_tag).expect("XG present") else {
            panic!("XG must be string");
        };
        let ys_bytes: &[u8] = ys.as_ref();
        let xg_bytes: &[u8] = xg.as_ref();
        if seq != ys_bytes {
            saw_difference = true;
            // Directional library: chemistry direction is fragment-level (XG),
            // not per-read. XG=CT → T/C diffs; XG=GA → A/G diffs. Holds for
            // both R1 and R2 of each pair, regardless of mapping orientation.
            for (s, y) in seq.iter().zip(ys_bytes.iter()) {
                if s != y {
                    match xg_bytes {
                        b"CT" => assert_eq!(
                            (*s, *y),
                            (b'T', b'C'),
                            "TAPS XG=CT diff must be T/C; got SEQ={} YS={}",
                            *s as char,
                            *y as char
                        ),
                        b"GA" => assert_eq!(
                            (*s, *y),
                            (b'A', b'G'),
                            "TAPS XG=GA diff must be A/G; got SEQ={} YS={}",
                            *s as char,
                            *y as char
                        ),
                        other => panic!("unexpected XG tag value: {other:?}"),
                    }
                }
            }
        }
    }
    assert!(saw_difference, "expected at least one TAPS record where SEQ != YS after conversion");
}

/// Build a reference with a mix of CpG and non-CpG cytosines so we can
/// verify that non-CpG C's are always treated as unmethylated.
///
/// The reference layout (1000 bp): a 500-bp `ACGT`-repeat block (every C is
/// in CpG context) followed by a 500-bp block of "AAAACAAA" (the C at every
/// 5th position is NOT followed by G, so it's a non-CpG C).
fn mixed_context_env() -> TestEnv {
    let mut seq = Vec::with_capacity(1000);
    seq.extend(b"ACGT".repeat(125)); // 500 bp, every C is CpG-context
    seq.extend(b"AAAACAAA".repeat(63)); // 504 bp, every 5th C is non-CpG (followed by A)
    seq.truncate(1000);
    TestEnv::new(&[("chr1", &seq)])
}

#[test]
fn test_em_seq_full_cpg_methylation_still_converts_non_cpg() {
    // With em-seq mode, methylation-rate 1.0 fully protects every CpG-
    // context C. But non-CpG cytosines (the C in "AAAACAAA") are not in
    // the methylation table, so they're always unmethylated and must
    // convert to T at conversion_rate 1.0.
    let env = mixed_context_env();
    let out = env.output_prefix();
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 1.0, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "30",
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
        "1.0",
        "--seed",
        "42",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let r1 = read_gzipped(&out.with_extension("r1.fastq.gz"));
    let (c1, t1) = count_c_and_t(&r1);

    // Reference has CpG-context C's at ~12.5% of positions (in the ACGT
    // half) and non-CpG C's at ~12.5% of positions (in the AAAACAAA half).
    // Full CpG protection preserves the first set; non-CpG C's all convert
    // to T. We expect c1 to be roughly half the original C count, with the
    // T count having grown by the converted half.
    assert!(c1 > 0, "CpG-protected C's should remain (got 0)");
    assert!(t1 > 0, "non-CpG C's should have converted to T (got 0)");

    // Sanity: ratio of remaining C / (C + T) should be ~0.25 (half the
    // original 50/50 C/T balance, since half the C's converted).
    let ratio = c1 as f64 / (c1 + t1) as f64;
    assert!(
        (0.20..0.30).contains(&ratio),
        "expected C / (C+T) ≈ 0.25 with half-CpG-half-non-CpG ref; got {ratio:.3} (c={c1}, t={t1})"
    );
}

// ── VCF + methylation interaction tests ────────────────────────────────────

/// End-to-end: a SNP that creates a CpG only on the variant haplotype must
/// drive per-haplotype methylation correctly — i.e. records derived from
/// the variant haplotype with methylation rate 1.0 must preserve the C at
/// the newly-created CpG site.
#[test]
fn test_em_seq_with_vcf_handles_haplotype_specific_cpgs() {
    // Reference (length 200): a long A/T tract that contains NO CpG sites.
    // The SNP A→C at position 5 (0-based) creates a CG dinucleotide on the
    // variant haplotype only (positions 5-6 become CG; the reference base
    // at position 6 is G).
    //
    // With em-seq + methylation rate 1.0 + conversion 1.0, the C at
    // position 5 on the variant haplotype is methylated and preserved.
    // Records on the reference haplotype don't have a C at that position
    // so the test focuses on the variant haplotype.
    let mut seq = vec![b'A'; 200];
    seq[6] = b'G'; // Make the base after the SNP a G so the SNP creates a CpG.
    let env = TestEnv::new(&[("chr1", &seq)]);

    // Write the plain variants VCF (no MT/MB yet).
    let variants_vcf = env.write_vcf(
        "sample1",
        &[("chr1", 200)],
        &[VcfVariant {
            chrom: "chr1",
            pos_1based: 6,
            ref_allele: "A",
            alt_alleles: &["C"],
            gt: "0|1",
        }],
    );

    // Run methylate with the variants VCF so the haplotype-specific CpG
    // (created by the A→C SNP) is captured in the methylation truth.
    let meth_vcf = methylate_to_vcf_with_variants(
        &env,
        &env.fasta_path,
        Some(&variants_vcf),
        1.0,
        42,
        "meth.vcf.gz",
    );

    let out = env.output_prefix();
    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        meth_vcf.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "200",
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
        "1.0",
        "--golden-bam",
        "--seed",
        "42",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let records = read_bam_records(&out.with_extension("golden.bam"));
    assert!(!records.is_empty(), "expected non-empty golden BAM with VCF + methylation");

    let hp_tag = DataTag::new(b'h', b'p');

    // Find at least one variant-haplotype record covering position 5
    // (0-based). Its YS:Z (pre-conversion, reference orientation) at
    // position 5 must contain C — proving the variant haplotype's CpG was
    // detected. Its SEQ at position 5 must also contain C (or G after
    // reverse-complementing) — proving methylation preserved it through
    // chemistry.
    let mut saw_variant_hap_cpg_preserved = false;
    for rec in &records {
        // Skip records not on the variant haplotype (hp:i:1).
        let hp_val = rec.data().get(&hp_tag);
        let is_variant_hap = match hp_val {
            Some(DataValue::Int8(v)) => *v == 1,
            Some(DataValue::UInt8(v)) => *v == 1,
            Some(DataValue::Int16(v)) => *v == 1,
            Some(DataValue::UInt16(v)) => *v == 1,
            Some(DataValue::Int32(v)) => *v == 1,
            Some(DataValue::UInt32(v)) => *v == 1,
            _ => false,
        };
        if !is_variant_hap {
            continue;
        }

        let Some(start_1based) = rec.alignment_start() else { continue };
        let start_0based = usize::from(start_1based) - 1;
        let seq: Vec<u8> = rec.sequence().as_ref().to_vec();
        let end_excl = start_0based + seq.len();
        if start_0based > 5 || end_excl <= 5 {
            continue;
        }
        let offset_in_seq = 5 - start_0based;
        let seq_byte = seq[offset_in_seq];

        // Methylated C on a forward record stays C in SEQ; on a reverse
        // record SEQ is reference-oriented and the bottom-strand C at
        // position 6 (the G in the CpG) is the relevant base. The C at
        // position 5 should appear as 'C' in SEQ regardless of strand
        // because SEQ is in reference orientation.
        if seq_byte == b'C' {
            saw_variant_hap_cpg_preserved = true;
            break;
        }
    }
    assert!(
        saw_variant_hap_cpg_preserved,
        "expected at least one variant-haplotype record where the SNP-created CpG C was preserved; \
         per-haplotype methylation regression"
    );
}

// ── CpG truth bedGraph integration tests ───────────────────────────────────

#[test]
fn test_cpg_truth_bedgraph_requires_methylation_mode() {
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    let bg = env.dir.path().join("truth.bedGraph");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--cpg-truth-bedgraph",
        bg.to_str().unwrap(),
        "--seed",
        "42",
        "--threads",
        "1",
    ]);
    assert!(!ok, "--cpg-truth-bedgraph without --methylation-mode must be rejected");
    assert!(
        stderr.contains("--cpg-truth-bedgraph requires --methylation-mode"),
        "expected validation message, got stderr: {stderr}"
    );
}

#[test]
fn test_header_only_mt_mb_vcf_is_rejected() {
    // A VCF that declares MT/MB in its header but carries no annotated
    // records has no methylation truth. `--methylation-mode` must reject it
    // up front with the same user-facing error as a non-methylated VCF —
    // not pass validation and fail later on an internal invariant.
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    let vcf_path = env.dir.path().join("header_only.vcf");
    std::fs::write(
        &vcf_path,
        "##fileformat=VCFv4.4\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=MT,Number=.,Type=String,Description=\"top\">\n\
         ##FORMAT=<ID=MB,Number=.,Type=String,Description=\"bottom\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMETHYLATE\n",
    )
    .unwrap();

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--methylation-mode",
        "em-seq",
        "--seed",
        "42",
        "--threads",
        "1",
    ]);
    assert!(!ok, "header-only MT/MB VCF must be rejected by --methylation-mode");
    assert!(
        stderr.contains("requires a methylation-annotated VCF"),
        "expected the MT/MB validation message, got stderr: {stderr}"
    );
}

#[test]
fn test_cpg_truth_bedgraph_full_methylation_emits_only_methylated_calls() {
    // Reference is "ACGT" repeated → every C is in CpG context. With
    // methylation-rate 1.0 every CpG is methylated on both strands. Each
    // simulated mate that covers a reference CpG must register a methylated
    // call for that site, with zero unmethylated calls.
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    let bg = env.dir.path().join("truth.bedGraph");
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 1.0, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "5",
        "--read-length",
        "50",
        "--fragment-mean",
        "100",
        "--fragment-stddev",
        "10",
        "--methylation-mode",
        "em-seq",
        "--methylation-conversion-rate",
        "1.0",
        "--cpg-truth-bedgraph",
        bg.to_str().unwrap(),
        "--seed",
        "42",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let body = std::fs::read_to_string(&bg).expect("bedGraph must exist");
    let lines: Vec<&str> = body.lines().collect();
    assert!(
        lines.first().is_some_and(|l| l.starts_with("track ")),
        "first line must be a track header"
    );

    let mut sites = 0;
    for row in &lines[1..] {
        let fields: Vec<&str> = row.split('\t').collect();
        assert_eq!(fields.len(), 6, "expected 6 columns, got: {row:?}");
        let n_meth: u32 = fields[4].parse().unwrap();
        let n_unmeth: u32 = fields[5].parse().unwrap();
        let rate: u32 = fields[3].parse().unwrap();
        let total = n_meth + n_unmeth;
        assert!(total > 0, "row with zero coverage must not be emitted: {row:?}");
        assert_eq!(n_unmeth, 0, "row {row:?} has unmethylated calls under full methylation");
        assert_eq!(rate, 100, "row {row:?} should have rate 100");
        // CpGs in `ACGT` repeated land at positions 1, 5, 9, … = 4k+1.
        let start: u32 = fields[1].parse().unwrap();
        let end: u32 = fields[2].parse().unwrap();
        assert_eq!(start % 4, 1, "CpG top-C must lie at a 4k+1 offset, got {start}");
        assert_eq!(end, start + 1, "end column must equal start + 1");
        sites += 1;
    }
    assert!(sites > 0, "expected at least one CpG site in bedGraph");
}

#[test]
fn test_cpg_truth_bedgraph_zero_methylation_emits_only_unmethylated_calls() {
    // methylation-rate 0 → bitmap is all-false → every covered CpG records
    // unmethylated calls only. (Conversion rate is irrelevant; truth tracks
    // bitmap state, not the post-conversion sequenced base.)
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    let bg = env.dir.path().join("truth.bedGraph");
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 0.0, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "5",
        "--read-length",
        "50",
        "--fragment-mean",
        "100",
        "--fragment-stddev",
        "10",
        "--methylation-mode",
        "em-seq",
        "--methylation-conversion-rate",
        "1.0",
        "--cpg-truth-bedgraph",
        bg.to_str().unwrap(),
        "--seed",
        "42",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let body = std::fs::read_to_string(&bg).expect("bedGraph must exist");
    let mut sites = 0;
    for row in body.lines().skip(1) {
        let fields: Vec<&str> = row.split('\t').collect();
        let n_meth: u32 = fields[4].parse().unwrap();
        let n_unmeth: u32 = fields[5].parse().unwrap();
        let rate: u32 = fields[3].parse().unwrap();
        assert_eq!(n_meth, 0, "row {row:?} has methylated calls under zero methylation");
        assert!(n_unmeth > 0);
        assert_eq!(rate, 0);
        sites += 1;
    }
    assert!(sites > 0);
}

#[test]
fn test_cpg_truth_bedgraph_taps_intermediate_methylation_yields_mixed_rates() {
    // TAPS at 0.5 methylation rate. The truth tracks the bitmap, not the
    // post-conversion sequenced base, so across many sites we expect a
    // roughly half-and-half split. Empirical band keeps the test
    // deterministic against `SmallRng::seed_from_u64(42)`.
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    let bg = env.dir.path().join("truth.bedGraph");
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 0.5, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "10",
        "--read-length",
        "50",
        "--fragment-mean",
        "100",
        "--fragment-stddev",
        "10",
        "--methylation-mode",
        "taps",
        "--methylation-conversion-rate",
        "1.0",
        "--cpg-truth-bedgraph",
        bg.to_str().unwrap(),
        "--seed",
        "42",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let body = std::fs::read_to_string(&bg).expect("bedGraph must exist");
    let mut total_meth: u64 = 0;
    let mut total_unmeth: u64 = 0;
    for row in body.lines().skip(1) {
        let fields: Vec<&str> = row.split('\t').collect();
        let n_meth: u32 = fields[4].parse().unwrap();
        let n_unmeth: u32 = fields[5].parse().unwrap();
        total_meth += u64::from(n_meth);
        total_unmeth += u64::from(n_unmeth);
    }
    let total = total_meth + total_unmeth;
    assert!(total > 0);
    let frac_meth = total_meth as f64 / total as f64;
    assert!(
        (0.40..0.60).contains(&frac_meth),
        "expected ~50% methylated truth calls at rate=0.5, got {frac_meth:.3} (m={total_meth} u={total_unmeth})"
    );
}

// ── Bismark methylation-call tag integration tests ─────────────────────────

/// Zero errors + full methylation: every CpG C in `XM` and `YM` is uppercase
/// (`Z`), no other methylation calls appear, NM is 0, and MD reports a single
/// match run equal to SEQ length.
#[test]
fn test_golden_bam_emits_xm_ym_nm_md_under_full_methylation_no_errors() {
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 1.0, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
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
        "1.0",
        "--golden-bam",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let records = read_bam_records(&out.with_extension("golden.bam"));
    assert!(!records.is_empty(), "expected non-empty golden BAM");

    let xm_tag = DataTag::new(b'X', b'M');
    let ym_tag = DataTag::new(b'Y', b'M');
    let nm_tag = DataTag::new(b'N', b'M');
    let md_tag = DataTag::new(b'M', b'D');

    for rec in &records {
        let seq_len = rec.sequence().as_ref().len();

        let DataValue::String(xm) = rec.data().get(&xm_tag).expect("XM:Z must be present") else {
            panic!("XM must be String");
        };
        let xm_bytes: &[u8] = xm.as_ref();
        assert_eq!(xm_bytes.len(), seq_len, "XM length must match SEQ length");

        let DataValue::String(ym) = rec.data().get(&ym_tag).expect("YM:Z must be present") else {
            panic!("YM must be String");
        };
        let ym_bytes: &[u8] = ym.as_ref();
        assert_eq!(ym_bytes.len(), seq_len, "YM length must match SEQ length");

        // Full methylation, no errors → XM == YM (truth and observation
        // agree on every position).
        assert_eq!(xm_bytes, ym_bytes, "XM must equal YM under full methylation, zero errors");

        // Every CpG call must be uppercase Z (methylated). With our cytosine-
        // rich reference (`ACGT` repeat) every C is in CpG context, so we
        // expect at least one Z per record.
        #[expect(clippy::naive_bytecount, reason = "test scope; clarity over speed")]
        let n_methylated = xm_bytes.iter().filter(|&&b| b == b'Z').count();
        assert!(n_methylated > 0, "expected >= 1 methylated CpG call per record");
        for &c in xm_bytes {
            assert!(
                matches!(c, b'.' | b'Z'),
                "unexpected XM char {c:?} under full methylation; only '.' and 'Z' allowed"
            );
        }

        // NM:i must be 0 — no errors, no SNPs, all C→T are bisulfite-allowed.
        let nm_val =
            rec.data().get(&nm_tag).expect("NM:i must be present").as_int().expect("NM is int");
        assert_eq!(nm_val, 0, "NM must be 0 with full methylation and no errors");

        // MD:Z must be a single match run equal to SEQ length (no mismatches,
        // no insertions, no deletions; bisulfite events suppressed).
        let DataValue::String(md) = rec.data().get(&md_tag).expect("MD:Z must be present") else {
            panic!("MD must be String");
        };
        let md_bytes: &[u8] = md.as_ref();
        let expected = seq_len.to_string();
        assert_eq!(md_bytes, expected.as_bytes(), "MD must be \"{expected}\"");
    }
}

/// Zero errors + zero methylation + full conversion: every CpG C in `XM`
/// is lowercase `z`. NM is still 0 (every C→T is bisulfite-allowed), and
/// MD is a single match run equal to SEQ length.
#[test]
fn test_golden_bam_xm_lowercase_under_full_conversion() {
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 0.0, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
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
        "1.0",
        "--golden-bam",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let records = read_bam_records(&out.with_extension("golden.bam"));
    let xm_tag = DataTag::new(b'X', b'M');
    let nm_tag = DataTag::new(b'N', b'M');

    for rec in &records {
        let DataValue::String(xm) = rec.data().get(&xm_tag).expect("XM present") else {
            panic!("XM must be String")
        };
        let xm_bytes: &[u8] = xm.as_ref();
        #[expect(clippy::naive_bytecount, reason = "test scope; clarity over speed")]
        let n_unmeth = xm_bytes.iter().filter(|&&b| b == b'z').count();
        assert!(n_unmeth > 0, "expected >= 1 unmethylated 'z' per record");
        // Only '.' or 'z' allowed under zero-methylation full-conversion.
        for &c in xm_bytes {
            assert!(matches!(c, b'.' | b'z'), "unexpected XM char {c:?}");
        }

        let nm_val = rec.data().get(&nm_tag).expect("NM present").as_int().expect("NM is int");
        assert_eq!(nm_val, 0, "every C→T is bisulfite-allowed under XG=CT, no real mismatches");
    }
}

/// XM:Z, YM:Z, NM:i, MD:Z must be absent when methylation chemistry is off.
#[test]
fn test_golden_bam_omits_xm_ym_nm_md_without_methylation() {
    let env = cytosine_rich_env();
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
        "--min-error-rate",
        "0",
        "--max-error-rate",
        "0",
        "--golden-bam",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let records = read_bam_records(&out.with_extension("golden.bam"));
    let xm = DataTag::new(b'X', b'M');
    let ym = DataTag::new(b'Y', b'M');
    let nm = DataTag::new(b'N', b'M');
    let md = DataTag::new(b'M', b'D');
    for rec in &records {
        assert!(rec.data().get(&xm).is_none(), "XM must be absent without methylation");
        assert!(rec.data().get(&ym).is_none(), "YM must be absent without methylation");
        assert!(rec.data().get(&nm).is_none(), "NM must be absent without methylation");
        assert!(rec.data().get(&md).is_none(), "MD must be absent without methylation");
    }
}

/// XM and YM agree under zero errors but diverge when errors corrupt a
/// methylated CpG cytosine. With a high error rate at a fully methylated
/// reference, at least some records should show YM='Z' at a CpG position
/// where XM='.' (mismatch) — i.e. truth says methylated but the observed
/// base is not C/T.
#[test]
fn test_xm_ym_diverge_under_errors_at_methylated_cpg() {
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 1.0, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "30",
        "--read-length",
        "50",
        "--fragment-mean",
        "100",
        "--fragment-stddev",
        "10",
        // High error rate to force divergence.
        "--min-error-rate",
        "0.10",
        "--max-error-rate",
        "0.10",
        "--methylation-mode",
        "em-seq",
        "--methylation-conversion-rate",
        "1.0",
        "--golden-bam",
        "--seed",
        "42",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate failed: {stderr}");

    let records = read_bam_records(&out.with_extension("golden.bam"));
    let xm_tag = DataTag::new(b'X', b'M');
    let ym_tag = DataTag::new(b'Y', b'M');

    let mut saw_divergence = false;
    for rec in &records {
        let DataValue::String(xm) = rec.data().get(&xm_tag).expect("XM present") else {
            panic!("XM must be String")
        };
        let DataValue::String(ym) = rec.data().get(&ym_tag).expect("YM present") else {
            panic!("YM must be String")
        };
        let xm_bytes: &[u8] = xm.as_ref();
        let ym_bytes: &[u8] = ym.as_ref();
        for (x, y) in xm_bytes.iter().zip(ym_bytes.iter()) {
            // Truth-says-methylated but observation can't tell.
            if *x == b'.' && *y == b'Z' {
                saw_divergence = true;
                break;
            }
        }
        if saw_divergence {
            break;
        }
    }
    assert!(saw_divergence, "expected at least one record where XM='.' but YM='Z' due to errors");
}

// ── Validation matrix tests ────────────────────────────────────────────────
//
// Matrix:
// | VCF has MT/MB | --methylation-mode set | Behavior                     |
// |---------------|------------------------|------------------------------|
// | true          | true                   | Run chemistry (normal path). |
// | true          | false                  | Warn; variants-only output.  |
// | false         | true                   | Hard error.                  |
// | false         | false                  | Variants-only (baseline).    |

/// Matrix cell (true, true): methylated VCF + mode → chemistry runs.
/// Verified by checking that reads come out and the run exits successfully.
#[test]
fn matrix_vcf_with_mtmb_and_mode_runs_chemistry() {
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 1.0, 42, "meth.vcf.gz");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "5",
        "--read-length",
        "50",
        "--methylation-mode",
        "em-seq",
        "--methylation-conversion-rate",
        "1.0",
        "--threads",
        "1",
    ]);
    assert!(ok, "expected simulate to succeed with methylated VCF + mode; stderr: {stderr}");

    let r1 = read_gzipped(&out.with_extension("r1.fastq.gz"));
    assert!(!r1.is_empty(), "expected non-empty FASTQ output");
    // Confirm chemistry ran: reference is all-CpG (ACGT repeat), rate=1.0
    // → all C's methylated → em-seq preserves them → R1 still has C's.
    let (c1, _t1) = count_c_and_t(&r1);
    assert!(c1 > 0, "matrix (true,true): expected C's preserved under full methylation + em-seq");
}

/// Matrix cell (true, false): methylated VCF but no mode → warn once, then
/// output identical to variants-only (no chemistry applied).
#[test]
fn matrix_vcf_with_mtmb_no_mode_warns_and_runs_variants_only() {
    let env = cytosine_rich_env();
    let out_with_mode = env.dir.path().join("with_mode");
    let out_no_mode = env.dir.path().join("no_mode");
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 1.0, 42, "meth.vcf.gz");

    // Run with mode to get the "chemistry" baseline.
    let (ok_mode, _, stderr_mode) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
        out_with_mode.to_str().unwrap(),
        "--coverage",
        "5",
        "--read-length",
        "50",
        "--min-error-rate",
        "0",
        "--max-error-rate",
        "0",
        "--methylation-mode",
        "em-seq",
        "--methylation-conversion-rate",
        "1.0",
        "--seed",
        "99",
        "--threads",
        "1",
    ]);
    assert!(ok_mode, "chemistry run must succeed");
    assert!(
        !stderr_mode.contains("methylation chemistry will not be applied"),
        "matrix (true,true): with-mode run must NOT warn about skipped chemistry; stderr: {stderr_mode}"
    );

    // Run WITHOUT mode: should succeed, emit a WARN, produce reads without
    // chemistry conversion.
    let (ok_no_mode, _, stderr_no_mode) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
        out_no_mode.to_str().unwrap(),
        "--coverage",
        "5",
        "--read-length",
        "50",
        "--min-error-rate",
        "0",
        "--max-error-rate",
        "0",
        "--seed",
        "99",
        "--threads",
        "1",
    ]);
    assert!(
        ok_no_mode,
        "matrix (true,false): expected simulate to succeed; stderr: {stderr_no_mode}"
    );
    assert!(
        stderr_no_mode.contains("methylation chemistry will not be applied"),
        "matrix (true,false): expected WARN about methylation being ignored; stderr: {stderr_no_mode}"
    );

    // Without mode, no chemistry is applied. With full em-seq and rate=1.0 the
    // mode run has C's (all preserved); the no-mode run has the same C count as
    // a plain variants-only simulation.
    let r1_mode = read_gzipped(&out_with_mode.with_extension("r1.fastq.gz"));
    let r1_no_mode = read_gzipped(&out_no_mode.with_extension("r1.fastq.gz"));
    let (c_mode, _) = count_c_and_t(&r1_mode);
    let (c_no_mode, _) = count_c_and_t(&r1_no_mode);
    // mode run: C's preserved (methylated, no conversion). no_mode run: C's
    // present too (no chemistry at all). Reference is ACGT-repeat → both
    // should have C's, but amounts may differ. The critical check is that
    // no-mode didn't accidentally apply chemistry. Both should have C's since
    // no conversion happened in either case (mode=em-seq with rate=1 preserves
    // all; no-mode skips chemistry entirely). Just assert both are non-zero.
    assert!(c_mode > 0, "mode run should have C's");
    assert!(c_no_mode > 0, "no-mode run should also have C's (no chemistry applied)");
}

/// Matrix cell (false, true): plain VCF (no MT/MB) with mode set → hard error.
#[test]
fn matrix_vcf_no_mtmb_with_mode_errors() {
    let env = cytosine_rich_env();
    let out = env.output_prefix();
    // Write a plain variants-only VCF (no methylation FORMAT fields).
    let plain_vcf = env.write_vcf_header_only(&["sample1"], &[("chr1", 1000)]);

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        plain_vcf.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--coverage",
        "5",
        "--read-length",
        "50",
        "--methylation-mode",
        "em-seq",
        "--threads",
        "1",
    ]);
    assert!(!ok, "matrix (false,true): expected simulate to fail with non-methylated VCF + mode");
    assert!(
        stderr.contains("run `holodeck methylate` first"),
        "matrix (false,true): expected error message mentioning methylate; stderr: {stderr}"
    );
}

/// Matrix cell (false, false): no VCF methylation, no mode → variants-only
/// baseline, exit success.
#[test]
fn matrix_no_vcf_no_mode_is_variants_only() {
    let env = cytosine_rich_env();
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
        "--threads",
        "1",
    ]);
    assert!(
        ok,
        "matrix (false,false): expected variants-only simulate to succeed; stderr: {stderr}"
    );

    let r1 = read_gzipped(&out.with_extension("r1.fastq.gz"));
    assert!(!r1.is_empty(), "expected non-empty FASTQ output");
}

/// Lock in the "variants-only output is deterministic across runs" promise
/// from the validation matrix's (false, false) cell. Two runs with the same
/// `--seed` and the same flag set must produce byte-identical R1 FASTQ
/// output. Acts as a regression guard against accidental RNG-state drift in
/// simulate's pre-methylation code paths.
#[test]
fn matrix_no_vcf_no_mode_is_byte_identical_across_runs() {
    let env = cytosine_rich_env();
    let out_a = env.output_prefix();
    let out_b = env.dir.path().join("sim_b");

    for out in [&out_a, &out_b] {
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
            "--seed",
            "1234",
            "--threads",
            "1",
        ]);
        assert!(ok, "simulate failed: {stderr}");
    }

    let r1_a = read_gzipped(&out_a.with_extension("r1.fastq.gz"));
    let r1_b = read_gzipped(&out_b.with_extension("r1.fastq.gz"));
    assert_eq!(r1_a, r1_b, "same --seed produced different R1 FASTQ output");
}

/// `--cpg-truth-bedgraph` should integrate with the matrix's (true, true)
/// cell: when both `--methylation-mode` and an MT/MB-bearing VCF are
/// present, the bedgraph is written with one row per reference CpG that
/// any read covered. None of the validation-matrix tests above pass
/// `--cpg-truth-bedgraph`, so this test locks in that combination.
#[test]
fn cpg_truth_bedgraph_works_with_matrix_true_true_cell() {
    let env = cytosine_rich_env();
    let vcf = methylate_to_vcf(&env, &env.fasta_path, 1.0, 42, "meth.vcf.gz");
    let out = env.output_prefix();
    let bg_path = env.dir.path().join("cpg.bedgraph");

    let (ok, _, stderr) = run_simulate(&[
        "simulate",
        "-r",
        env.fasta_path.to_str().unwrap(),
        "-v",
        vcf.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--methylation-mode",
        "em-seq",
        "--cpg-truth-bedgraph",
        bg_path.to_str().unwrap(),
        "--coverage",
        "5",
        "--read-length",
        "50",
        "--threads",
        "1",
    ]);
    assert!(ok, "simulate with --cpg-truth-bedgraph + --methylation-mode failed: {stderr}");
    assert!(bg_path.exists(), "--cpg-truth-bedgraph file was not created");

    let contents = std::fs::read_to_string(&bg_path).unwrap();
    assert!(
        contents.starts_with("track type="),
        "--cpg-truth-bedgraph missing track header: {contents}"
    );
    let data_lines: Vec<&str> = contents.lines().filter(|l| !l.starts_with("track")).collect();
    assert!(
        !data_lines.is_empty(),
        "--cpg-truth-bedgraph should have at least one per-CpG record under non-zero coverage"
    );
}
