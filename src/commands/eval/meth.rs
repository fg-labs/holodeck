//! Methylation-level correlation: per-CpG methylation fraction derived from
//! the aligner's Bismark `XM:Z` calls versus the simulated cpg-truth bedGraph.
//!
//! For each mapped read the `XM` string is walked alongside the CIGAR; every
//! `Z` (methylated CpG) or `z` (unmethylated CpG) call is tallied at its
//! reference position. The aligner methylation level at a site is
//! `n_methylated / coverage`, which is correlated against the truth level
//! (`rate / 100`) from the bedGraph that `simulate --cpg-truth-bedgraph`
//! writes. Reports Pearson r and RMSE over the shared, covered CpG sites in
//! `<prefix>.meth.tsv`.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};
use noodles::bam;
use noodles::sam::alignment::record_buf::Cigar;

use super::cigar;
use super::golden::{contig_name, string_tag};
use crate::commands::command::output_path;

/// Minimum aligner coverage at a CpG for it to enter the correlation.
const MIN_COVERAGE: u32 = 1;

/// Per-CpG methylated/unmethylated call counts, nested by contig then 0-based
/// position. Nesting keeps the contig name allocated once per read rather than
/// once per CpG call.
type AlignerTally = HashMap<String, HashMap<u32, (u32, u32)>>;

/// Per-CpG truth methylation level in `[0, 1]`, nested by contig then position.
type TruthLevels = HashMap<String, HashMap<u32, f64>>;

/// Correlation summary between aligner-called and truth methylation levels.
#[derive(Debug, Clone, Copy, PartialEq)]
struct MethCorr {
    n_cpg: usize,
    pearson_r: f64,
    rmse: f64,
}

/// Evaluate methylation-level correlation of `mapped` against `cpg_truth`.
///
/// # Errors
/// Returns an error if the BAM, the bedGraph, or the output cannot be read or
/// written.
pub fn run(mapped: &Path, cpg_truth: &Path, output_prefix: &Path) -> Result<()> {
    let tally = tally_aligner(mapped)?;
    let truth = parse_bedgraph(cpg_truth)?;
    let corr = correlate(&tally, &truth);

    let path = output_path(output_prefix, ".meth.tsv");
    let mut out =
        File::create(&path).with_context(|| format!("Failed to create {}", path.display()))?;
    writeln!(out, "n_cpg\tpearson_r\trmse")?;
    writeln!(out, "{}\t{}\t{}", corr.n_cpg, fmt_opt(corr.pearson_r), fmt_opt(corr.rmse))?;
    log::info!(
        "Methylation correlation over {} CpGs: r={}, rmse={}",
        corr.n_cpg,
        fmt_opt(corr.pearson_r),
        fmt_opt(corr.rmse)
    );
    log::info!("Methylation results written to: {}", path.display());
    Ok(())
}

/// Format a possibly-NaN statistic as a fixed-precision value or `NA`.
fn fmt_opt(value: f64) -> String {
    if value.is_nan() { "NA".to_string() } else { format!("{value:.4}") }
}

/// Tally per-CpG methylated/unmethylated calls from every mapped read's `XM`.
fn tally_aligner(mapped: &Path) -> Result<AlignerTally> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(mapped)
        .with_context(|| format!("Failed to open BAM: {}", mapped.display()))?;
    let header = reader.read_header()?;

    let mut tally: AlignerTally = HashMap::new();
    let mut saw_mapped_primary = false;
    let mut saw_xm = false;
    for result in reader.record_bufs(&header) {
        let record = result.context("Failed to read BAM record")?;
        let flags = record.flags();
        if flags.is_secondary() || flags.is_supplementary() || flags.is_unmapped() {
            continue;
        }
        saw_mapped_primary = true;
        let Some(xm) = string_tag(&record, b'X', b'M') else { continue };
        saw_xm = true;
        let Some(ref_id) = record.reference_sequence_id() else { continue };
        let Some(contig) = contig_name(&header, ref_id) else { continue };
        let Some(start) = record.alignment_start() else { continue };
        let start0 = u32::try_from(usize::from(start).saturating_sub(1)).unwrap_or(0);
        tally_read(&contig, record.cigar(), start0, xm.as_bytes(), &mut tally);
    }

    // A BAM with mapped reads but no XM tags carries no methylation calls — for
    // example a plain bisulfite aligner like bwameth, where calling is a
    // separate extractor (MethylDackel) step. Warn and return the empty tally so
    // the correlation is reported as NA rather than failing the whole eval;
    // placement and variant representation are still meaningful for such a BAM.
    if saw_mapped_primary && !saw_xm {
        log::warn!(
            "no XM methylation tags in mapped primaries of {}; \
             reporting NA methylation correlation",
            mapped.display()
        );
    }
    Ok(tally)
}

/// Tally one read's CpG calls into `tally`, mapping each `XM` symbol to its
/// reference position via the CIGAR. The contig key is allocated once here
/// rather than once per CpG call.
fn tally_read(contig: &str, cigar: &Cigar, start0: u32, xm: &[u8], tally: &mut AlignerTally) {
    let by_pos = tally.entry(contig.to_string()).or_default();
    cigar::for_each_aligned(cigar, start0, |read_off, ref_pos| match xm.get(read_off) {
        Some(b'Z') => by_pos.entry(ref_pos).or_default().0 += 1,
        Some(b'z') => by_pos.entry(ref_pos).or_default().1 += 1,
        _ => {}
    });
}

/// Parse a MethylDackel-format CpG bedGraph into per-site truth levels in
/// `[0, 1]`. Columns: `chrom start end rate(0-100) n_meth n_unmeth`.
///
/// Blank lines and `track` / `#` headers are skipped. A data row that is
/// truncated, has an unparseable start/rate, or a rate outside `0..=100` is a
/// hard error: silently dropping malformed rows would understate coverage and
/// distort the correlation rather than surfacing a bad input.
fn parse_bedgraph(path: &Path) -> Result<TruthLevels> {
    let file =
        File::open(path).with_context(|| format!("Failed to open bedGraph: {}", path.display()))?;
    let mut truth: TruthLevels = HashMap::new();
    for (idx, line) in BufReader::new(file).lines().enumerate() {
        let line = line.context("Failed to read bedGraph line")?;
        if line.is_empty() || line.starts_with("track") || line.starts_with('#') {
            continue;
        }
        let line_no = idx + 1;
        let mut fields = line.split_whitespace();
        let (Some(chrom), Some(start), Some(end), Some(rate)) =
            (fields.next(), fields.next(), fields.next(), fields.next())
        else {
            bail!("malformed bedGraph line {line_no} in {}: expected >=4 columns", path.display());
        };
        let start0 = start.parse::<u32>().with_context(|| {
            format!("invalid bedGraph start on line {line_no} in {}", path.display())
        })?;
        let end = end.parse::<u32>().with_context(|| {
            format!("invalid bedGraph end on line {line_no} in {}", path.display())
        })?;
        // Each CpG-truth row must span exactly one site; a merged or malformed
        // interval would otherwise be silently scored as a single CpG at start,
        // skewing coverage and the reported correlation.
        if end != start0.saturating_add(1) {
            bail!(
                "bedGraph interval must span exactly one CpG on line {line_no} in {}: {start0}-{end}",
                path.display()
            );
        }
        let rate = rate.parse::<f64>().with_context(|| {
            format!("invalid bedGraph rate on line {line_no} in {}", path.display())
        })?;
        if !(0.0..=100.0).contains(&rate) {
            bail!("bedGraph rate out of range on line {line_no} in {}: {rate}", path.display());
        }
        truth.entry(chrom.to_string()).or_default().insert(start0, rate / 100.0);
    }
    Ok(truth)
}

/// Correlate aligner methylation levels against truth over shared, covered
/// CpG sites.
fn correlate(tally: &AlignerTally, truth: &TruthLevels) -> MethCorr {
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for (contig, by_pos) in tally {
        let Some(truth_pos) = truth.get(contig) else { continue };
        for (pos, &(meth, unmeth)) in by_pos {
            let coverage = meth + unmeth;
            if coverage < MIN_COVERAGE {
                continue;
            }
            if let Some(&truth_level) = truth_pos.get(pos) {
                xs.push(f64::from(meth) / f64::from(coverage));
                ys.push(truth_level);
            }
        }
    }
    MethCorr { n_cpg: xs.len(), pearson_r: pearson(&xs, &ys), rmse: rmse(&xs, &ys) }
}

/// Pearson correlation coefficient; `NaN` when undefined (n < 2 or a constant
/// series).
fn pearson(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len();
    if n < 2 {
        return f64::NAN;
    }
    let nf = n as f64;
    let mean_x = xs.iter().sum::<f64>() / nf;
    let mean_y = ys.iter().sum::<f64>() / nf;
    let mut cov = 0.0;
    let mut var_x = 0.0;
    let mut var_y = 0.0;
    for (&x, &y) in xs.iter().zip(ys) {
        let (dx, dy) = (x - mean_x, y - mean_y);
        cov += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }
    let denom = (var_x * var_y).sqrt();
    if denom == 0.0 { f64::NAN } else { cov / denom }
}

/// Root-mean-square error between aligner and truth levels.
fn rmse(xs: &[f64], ys: &[f64]) -> f64 {
    if xs.is_empty() {
        return f64::NAN;
    }
    let sse: f64 = xs.iter().zip(ys).map(|(&x, &y)| (x - y).powi(2)).sum();
    (sse / xs.len() as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::record::cigar::op::{Kind, Op};

    fn cigar(ops: &[(Kind, usize)]) -> Cigar {
        Cigar::from(ops.iter().map(|&(k, n)| Op::new(k, n)).collect::<Vec<_>>())
    }

    #[test]
    fn tally_read_maps_cpg_calls_to_reference_positions() {
        // 5M at ref 100; XM "Zz.zZ": ref 100 -> Z, 101 -> z, 103 -> z, 104 -> Z.
        let mut tally = AlignerTally::new();
        tally_read("chr1", &cigar(&[(Kind::Match, 5)]), 100, b"Zz.zZ", &mut tally);
        let chr1 = &tally["chr1"];
        assert_eq!(chr1[&100], (1, 0));
        assert_eq!(chr1[&101], (0, 1));
        assert_eq!(chr1[&104], (1, 0));
        assert!(!chr1.contains_key(&102)); // '.' ignored
    }

    #[test]
    fn tally_read_handles_insertion_offset() {
        // 2M2I2M at ref 100: read offsets 0,1 -> ref 100,101; offsets 2,3 are
        // the insertion (no ref); offsets 4,5 -> ref 102,103.
        let mut tally = AlignerTally::new();
        tally_read(
            "chr1",
            &cigar(&[(Kind::Match, 2), (Kind::Insertion, 2), (Kind::Match, 2)]),
            100,
            b"zzZZzz",
            &mut tally,
        );
        let chr1 = &tally["chr1"];
        // Insertion calls (offsets 2,3 = "ZZ") must not land on a reference pos.
        assert_eq!(chr1[&102], (0, 1));
        assert_eq!(chr1[&103], (0, 1));
        assert_eq!(chr1.values().map(|&(m, _)| m).sum::<u32>(), 0);
    }

    #[test]
    fn pearson_perfect_positive() {
        let xs = [0.0, 0.5, 1.0];
        let ys = [0.0, 0.5, 1.0];
        assert!((pearson(&xs, &ys) - 1.0).abs() < 1e-9);
        assert!(rmse(&xs, &ys).abs() < 1e-9);
    }

    #[test]
    fn pearson_is_nan_for_constant_series() {
        let xs = [0.5, 0.5, 0.5];
        let ys = [0.1, 0.9, 0.5];
        assert!(pearson(&xs, &ys).is_nan());
    }

    #[test]
    fn correlate_joins_on_shared_covered_sites() {
        let mut tally = AlignerTally::new();
        tally.insert("chr1".to_string(), HashMap::from([(10, (3, 1)), (20, (0, 4)), (30, (4, 0))])); // levels 0.75, 0.0, and a site absent from truth
        let mut truth = TruthLevels::new();
        truth.insert("chr1".to_string(), HashMap::from([(10, 0.75), (20, 0.0), (99, 0.5)]));
        let corr = correlate(&tally, &truth);
        assert_eq!(corr.n_cpg, 2);
        assert!((corr.pearson_r - 1.0).abs() < 1e-9);
        assert!(corr.rmse.abs() < 1e-9);
    }

    #[test]
    fn fmt_opt_renders_na_for_nan() {
        assert_eq!(fmt_opt(f64::NAN), "NA");
        assert_eq!(fmt_opt(0.5), "0.5000");
    }

    fn write_temp(content: &str) -> tempfile::NamedTempFile {
        use std::io::Write as _;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f.flush().unwrap();
        f
    }

    #[test]
    fn parse_bedgraph_reads_valid_rows_and_skips_headers() {
        let f =
            write_temp("track type=\"bedGraph\"\nchr1\t10\t11\t75\t3\t1\nchr1\t20\t21\t0\t0\t4\n");
        let truth = parse_bedgraph(f.path()).unwrap();
        let chr1 = &truth["chr1"];
        assert_eq!(chr1.len(), 2);
        assert!((chr1[&10] - 0.75).abs() < 1e-9);
        assert!((chr1[&20] - 0.0).abs() < 1e-9);
    }

    #[test]
    fn parse_bedgraph_rejects_truncated_row() {
        let f = write_temp("chr1\t10\n");
        assert!(parse_bedgraph(f.path()).is_err());
    }

    #[test]
    fn parse_bedgraph_rejects_out_of_range_rate() {
        let f = write_temp("chr1\t10\t11\t150\t3\t1\n");
        assert!(parse_bedgraph(f.path()).is_err());
    }

    #[test]
    fn parse_bedgraph_rejects_multi_cpg_interval() {
        // A row spanning more than one base is a merged/malformed interval: it
        // must be rejected rather than silently scored as a single CpG at start.
        let f = write_temp("chr1\t10\t15\t75\t3\t1\n");
        assert!(parse_bedgraph(f.path()).is_err());
    }

    #[test]
    fn parse_bedgraph_rejects_unparseable_end() {
        let f = write_temp("chr1\t10\tnope\t75\t3\t1\n");
        assert!(parse_bedgraph(f.path()).is_err());
    }
}
