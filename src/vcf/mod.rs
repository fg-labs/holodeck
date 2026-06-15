//! VCF input/output for holodeck.
//!
//! Reads variants from a multi-sample VCF (resolving GT against a chosen
//! sample) for use by `simulate` and `methylate`, and houses the methylation
//! classifier and `MT`/`MB` FORMAT reader/writer that carries per-haplotype
//! per-strand CpG methylation truth between the two commands. Submodules:
//! - [`genotype`] — `Genotype` and `VariantRecord` types plus GT parsing.
//! - [`methylation`] — per-CpG classifier, `MT`/`MB` writer/reader, and
//!   header probe.

pub mod genotype;
/// Per-CpG methylation classifier and MT/MB FORMAT VCF reader/writer used by
/// the `methylate` subcommand and the methylation-aware `simulate` path.
pub mod methylation;

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Read};
use std::path::Path;

use anyhow::{Context, Result, bail};
use noodles::vcf;
use noodles::vcf::variant::record::AlternateBases;

use crate::sequence_dict::SequenceDictionary;
use genotype::{Genotype, VariantRecord};

/// Result of a single-pass VCF scan: per-contig ALT-bearing variants plus
/// the sample's ploidy resolved from *unfiltered* genotypes.
#[derive(Debug, Default)]
pub struct ParsedVariants {
    /// ALT-bearing variant records, partitioned by contig name, each sorted
    /// by position. Homozygous-reference / all-missing records are dropped
    /// from this map (downstream haplotype construction only needs ALT
    /// calls) but still counted toward [`Self::sample_ploidy`]. Contigs with
    /// no ALT-bearing records do not appear; treat a missing key as "no
    /// variants on that contig".
    pub by_contig: HashMap<String, Vec<VariantRecord>>,
    /// Max GT ploidy across every record with a parseable genotype for the
    /// selected sample — **including** hom-ref / all-missing records that are
    /// excluded from `by_contig`. Falls back to `2` only when the VCF has no
    /// genotyped records at all. This is the authoritative MT/MB entry count:
    /// deriving ploidy from `by_contig` alone would undercount on a sample
    /// whose only calls on some contig are hom-ref (e.g. a haploid sample
    /// with a hom-ref-only contig would otherwise fall back to diploid).
    pub sample_ploidy: usize,
}

/// Read every variant from a VCF in a single pass and partition by contig.
///
/// This is the one-pass replacement for the per-contig loader pattern that
/// re-scanned the entire VCF for every contig. Mirrors
/// [`crate::vcf::methylation::parse_methylation_vcf`]: read the file once
/// before the contig loop, then look up per-contig records by name. Trades
/// I/O for memory — every selected variant is held in the returned map until
/// the caller drops it — which is the same trade-off the methylation loader
/// already makes.
///
/// Ploidy is captured during this scan *before* the hom-ref/all-missing
/// filter, so [`ParsedVariants::sample_ploidy`] reflects the true sample
/// ploidy even on contigs (or whole samples) with no ALT calls.
///
/// # Arguments
/// * `path` — Path to the VCF file (plain text, gzip, or BGZF; detected from
///   the leading magic bytes).
/// * `sample_name` — Sample to resolve genotypes for, or `None` to use the
///   only sample in a single-sample VCF.
/// * `_dict` — Sequence dictionary (reserved for future cross-validation
///   against the reference).
///
/// # Errors
/// Returns an error if the VCF cannot be read, the sample is not found, or a
/// GT field is malformed.
pub fn parse_variants_by_contig(
    path: &Path,
    sample_name: Option<&str>,
    _dict: &SequenceDictionary,
) -> Result<ParsedVariants> {
    let (header, mut reader) = open_lenient_vcf_reader(path)?;
    let sample_index = resolve_sample_index(&header, sample_name)?;

    let mut by_contig: HashMap<String, Vec<VariantRecord>> = HashMap::new();
    // Track ploidy across ALL genotyped records, including hom-ref ones that
    // are filtered out of `by_contig` below.
    let mut max_ploidy: Option<usize> = None;

    for result in reader.record_bufs(&header) {
        let record = result.with_context(|| "Failed to read VCF record")?;

        // Parse position (VCF is 1-based, convert to 0-based).
        let Some(pos_1based) = record.variant_start() else {
            continue;
        };
        #[expect(clippy::cast_possible_truncation, reason = "genomic positions fit in u32")]
        let position = (usize::from(pos_1based) - 1) as u32;

        // Extract alleles.
        let ref_allele: Vec<u8> = record.reference_bases().as_bytes().to_vec();
        let alt_alleles: Vec<Vec<u8>> = record
            .alternate_bases()
            .iter()
            .filter_map(Result::ok)
            .map(|a| a.as_bytes().to_vec())
            .collect();

        // Parse genotype for the selected sample.
        let gt_str = extract_gt_string(&record, sample_index);
        let Some(gt_str) = gt_str else {
            continue;
        };
        let genotype = Genotype::parse(&gt_str)?;

        // Record ploidy from the unfiltered GT, BEFORE dropping hom-ref /
        // all-missing records — those still carry the sample's ploidy.
        let ploidy = genotype.ploidy();
        max_ploidy = Some(max_ploidy.map_or(ploidy, |m| m.max(ploidy)));

        if !genotype.has_alt() {
            continue;
        }

        let chrom = record.reference_sequence_name().to_string();
        by_contig.entry(chrom).or_default().push(VariantRecord {
            position,
            ref_allele,
            alt_alleles,
            genotype,
        });
    }

    // Sort each contig's variants by position. Valid VCFs are already sorted,
    // but a stray unsorted file would break downstream invariants
    // (`build_haplotypes`, overlap checks, etc.) silently if we trusted it.
    for v in by_contig.values_mut() {
        v.sort_by_key(|vr| vr.position);
    }

    Ok(ParsedVariants { by_contig, sample_ploidy: max_ploidy.unwrap_or(2) })
}

/// Validate that a VCF file has the expected sample configuration and
/// return the resolved sample name.
///
/// Opens the VCF, reads the header, and checks:
/// - If `sample_name` is `Some`, that the named sample exists.
/// - If `sample_name` is `None`, that the VCF has exactly one sample.
///
/// Call this during argument validation to surface sample errors before
/// the per-contig simulation loop begins.
///
/// # Errors
/// Returns an error if the VCF cannot be read or the sample configuration
/// is invalid.
pub(crate) fn validate_vcf_sample(path: &Path, sample_name: Option<&str>) -> Result<String> {
    let (header, _reader) = open_lenient_vcf_reader(path)?;
    let idx = resolve_sample_index(&header, sample_name)?;
    Ok(header
        .sample_names()
        .get_index(idx)
        .expect("resolve_sample_index returned a valid index")
        .clone())
}

/// Resolve the sample index from a VCF header.
///
/// If `sample_name` is `Some`, looks up that sample. If `None` and the VCF has
/// exactly one sample, uses index 0. Errors if the VCF has multiple samples
/// and no name was specified.
fn resolve_sample_index(header: &vcf::Header, sample_name: Option<&str>) -> Result<usize> {
    let sample_names = header.sample_names();

    match sample_name {
        Some(name) => {
            let idx = sample_names.iter().position(|s| s == name).ok_or_else(|| {
                if sample_names.is_empty() {
                    anyhow::anyhow!("Sample '{name}' not found in VCF (VCF has no sample columns)")
                } else {
                    let available: Vec<&str> =
                        sample_names.iter().map(String::as_str).take(10).collect();
                    anyhow::anyhow!(
                        "Sample '{name}' not found in VCF. Available samples: {}",
                        available.join(", ")
                    )
                }
            })?;
            Ok(idx)
        }
        None => {
            if sample_names.len() == 1 {
                Ok(0)
            } else if sample_names.is_empty() {
                bail!("VCF has no sample columns")
            } else {
                bail!(
                    "VCF has {} samples but no --sample was specified. \
                     Available samples: {}",
                    sample_names.len(),
                    sample_names.iter().take(10).cloned().collect::<Vec<_>>().join(", ")
                )
            }
        }
    }
}

/// Open a VCF and return its parsed header plus a record reader, tolerating
/// duplicate meta-information definitions that noodles would otherwise reject.
///
/// noodles' header parser errors on any redeclared `##INFO`/`##FORMAT`/
/// `##FILTER`/`##ALT`/`##contig` ID (e.g. `duplicate INFO ID: BREAKSIMLENGTH`).
/// Real VCFs from upstream tools (Manta, DELLY, hand-merged files) routinely
/// carry such duplicates, and bcftools reads them without complaint; holodeck
/// must not be stricter. This reads the raw header text first, drops every
/// repeated structured definition (keeping the first, matching bcftools), and
/// hands the cleaned header back to noodles for parsing. The record stream is
/// untouched — only the header is rewritten.
///
/// Compression is detected from the leading magic bytes (plain text, gzip, or
/// BGZF), so the reader is also indifferent to the file's extension.
///
/// # Errors
/// Returns an error if the file cannot be opened or the (cleaned) header fails
/// to parse.
fn open_lenient_vcf_reader(
    path: &Path,
) -> Result<(vcf::Header, vcf::io::Reader<Box<dyn BufRead>>)> {
    let mut inner = open_vcf_buf_reader(path)
        .with_context(|| format!("Failed to open VCF: {}", path.display()))?;

    // Read raw header lines up to and including the `#CHROM` column header,
    // dropping duplicate structured definitions. Everything after `#CHROM`
    // stays unread in `inner` and is streamed as records below.
    let mut header_text = String::new();
    let mut seen: HashSet<(String, String)> = HashSet::new();
    let mut line = String::new();
    loop {
        line.clear();
        let n = inner
            .read_line(&mut line)
            .with_context(|| format!("Failed to read VCF header: {}", path.display()))?;
        if n == 0 {
            break; // EOF before `#CHROM`; let noodles report the malformed header
        }
        let is_meta = line.starts_with("##");
        if is_meta
            && let Some(key) = structured_header_key(&line)
            && !seen.insert(key.clone())
        {
            log::warn!("dropping duplicate VCF header definition: ##{}=<ID={},...>", key.0, key.1);
            continue;
        }
        let is_column_header = !is_meta; // first non-`##` line is `#CHROM`
        header_text.push_str(&line);
        if is_column_header {
            break;
        }
    }

    // Replay the cleaned header in front of the still-unread record stream so
    // noodles parses the header and reads records from one continuous source.
    let chained: Box<dyn BufRead> =
        Box::new(std::io::Cursor::new(header_text.into_bytes()).chain(inner));
    let mut reader = vcf::io::Reader::new(chained);
    let header = reader
        .read_header()
        .with_context(|| format!("Failed to parse VCF header: {}", path.display()))?;
    Ok((header, reader))
}

/// Open a VCF for buffered reading, transparently decompressing gzip/BGZF.
///
/// Compression is detected from the leading two magic bytes (`1f 8b`) rather
/// than the file extension, and `MultiGzDecoder` handles both single-member
/// gzip and the multi-member BGZF blocks `bgzip`/`bcftools` produce.
///
/// # Errors
/// Returns an I/O error if the file cannot be opened or its first bytes read.
pub(crate) fn open_vcf_buf_reader(path: &Path) -> std::io::Result<Box<dyn BufRead>> {
    use std::io::BufReader;

    let mut peek = [0u8; 2];
    // A file shorter than two bytes cannot be gzip; only treat the input as
    // compressed when both magic bytes were actually read.
    let bytes_read = {
        let mut f = std::fs::File::open(path)?;
        f.read(&mut peek)?
    };
    let file = std::fs::File::open(path)?;
    if bytes_read == 2 && peek == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Extract the `(directive, ID)` identity of a structured meta-information
/// line, e.g. `##INFO=<ID=DP,...>` → `("INFO", "DP")`.
///
/// Returns `None` for lines without a leading `ID=` field — `##fileformat`,
/// free-form `##key=value`, or the rare case where `ID` is not the first
/// attribute. Such lines are never deduplicated, matching the categories
/// noodles enforces unique (`INFO`/`FORMAT`/`FILTER`/`ALT`/`contig`), whose
/// emitters all write `ID` first.
fn structured_header_key(line: &str) -> Option<(String, String)> {
    let rest = line.trim_start_matches('#').trim_end_matches(['\n', '\r']);
    let (directive, value) = rest.split_once('=')?;
    let id_field = value.strip_prefix('<')?.strip_prefix("ID=")?;
    let id = id_field.split([',', '>']).next()?;
    Some((directive.to_string(), id.to_string()))
}

/// Extract the GT field as a string for a specific sample from a VCF record.
///
/// Returns `None` if the GT field is absent for this sample.
fn extract_gt_string(record: &vcf::variant::RecordBuf, sample_index: usize) -> Option<String> {
    use noodles::vcf::variant::record::samples::keys::key;
    use noodles::vcf::variant::record_buf::samples::sample::Value;

    let samples = record.samples();
    let sample = samples.get_index(sample_index)?;
    let gt_value = sample.get(key::GENOTYPE)?;

    match gt_value {
        Some(Value::Genotype(gt)) => Some(format_genotype(gt)),
        _ => None,
    }
}

/// Format a noodles `Genotype` value back to a VCF GT string (e.g. "0/1" or
/// "1|0").
fn format_genotype(
    gt: &noodles::vcf::variant::record_buf::samples::sample::value::Genotype,
) -> String {
    use noodles::vcf::variant::record::samples::series::value::genotype::Phasing;

    let alleles = gt.as_ref();
    let mut parts = Vec::with_capacity(alleles.len());
    let mut is_phased = false;

    for (i, allele) in alleles.iter().enumerate() {
        if i > 0 && allele.phasing() == Phasing::Phased {
            is_phased = true;
        }
        match allele.position() {
            Some(idx) => parts.push(idx.to_string()),
            None => parts.push(".".to_string()),
        }
    }

    let sep = if is_phased { "|" } else { "/" };
    parts.join(sep)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    /// Write a temporary plain-text VCF (no BGZF) and return the
    /// `NamedTempFile` so the caller controls its lifetime.
    fn write_temp_vcf(body: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        f.write_all(body.as_bytes()).unwrap();
        f.flush().unwrap();
        f
    }

    /// Build a minimal `SequenceDictionary` covering the test contigs.
    fn dict_for(contigs: &[(&str, usize)]) -> SequenceDictionary {
        let entries: Vec<_> = contigs
            .iter()
            .enumerate()
            .map(|(i, &(name, len))| {
                crate::sequence_dict::SequenceMetadata::new(i, name.to_string(), len)
            })
            .collect();
        SequenceDictionary::from_entries(entries)
    }

    /// Variants on two contigs must be partitioned by contig name, sorted by
    /// position, and stripped of homozygous-reference / no-alt records. A
    /// caller that loops over contigs should be able to look up each
    /// contig's slice in `O(1)` after a single full-file pass.
    #[test]
    fn parse_variants_by_contig_partitions_and_sorts() {
        let vcf = "\
##fileformat=VCFv4.4
##contig=<ID=chr1,length=100>
##contig=<ID=chr2,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t30\t.\tA\tT\t.\t.\t.\tGT\t0|1
chr2\t10\t.\tC\tG\t.\t.\t.\tGT\t1|0
chr1\t10\t.\tG\tC\t.\t.\t.\tGT\t0|1
chr1\t50\t.\tT\tA\t.\t.\t.\tGT\t0|0
chr2\t20\t.\tA\tG\t.\t.\t.\tGT\t0|1
";
        let f = write_temp_vcf(vcf);
        let dict = dict_for(&[("chr1", 100), ("chr2", 100)]);
        let parsed = parse_variants_by_contig(f.path(), None, &dict).unwrap();
        let by_contig = &parsed.by_contig;

        // chr1: two records (POS 30 + POS 10), sorted by position. The
        // hom-ref `0|0` at POS 50 must be dropped.
        let chr1 = by_contig.get("chr1").expect("chr1 must be present");
        assert_eq!(chr1.len(), 2, "expected 2 chr1 variants, got {chr1:?}");
        assert_eq!(chr1[0].position, 9, "first chr1 variant should be POS 10 → 0-based 9");
        assert_eq!(chr1[1].position, 29, "second chr1 variant should be POS 30 → 0-based 29");

        // chr2: two records, sorted by position.
        let chr2 = by_contig.get("chr2").expect("chr2 must be present");
        assert_eq!(chr2.len(), 2, "expected 2 chr2 variants, got {chr2:?}");
        assert_eq!(chr2[0].position, 9);
        assert_eq!(chr2[1].position, 19);

        // Diploid sample → ploidy 2 (resolved from the GTs).
        assert_eq!(parsed.sample_ploidy, 2);
    }

    /// A VCF with no alt-bearing records on a given contig must produce no
    /// entry for that contig — callers treat a missing key as "no variants
    /// here", and an empty-vector value would be equally fine but wastes a
    /// `HashMap` allocation per absent contig.
    #[test]
    fn parse_variants_by_contig_skips_contigs_with_no_alt_records() {
        let vcf = "\
##fileformat=VCFv4.4
##contig=<ID=chr1,length=100>
##contig=<ID=chr2,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t10\t.\tA\tT\t.\t.\t.\tGT\t0|0
chr2\t10\t.\tC\tG\t.\t.\t.\tGT\t1|0
";
        let f = write_temp_vcf(vcf);
        let dict = dict_for(&[("chr1", 100), ("chr2", 100)]);
        let parsed = parse_variants_by_contig(f.path(), None, &dict).unwrap();
        assert!(!parsed.by_contig.contains_key("chr1"), "chr1 should be absent (only hom-ref)");
        assert_eq!(parsed.by_contig.get("chr2").map(Vec::len), Some(1));
    }

    /// Sample ploidy must be resolved from *unfiltered* genotypes. A haploid
    /// sample whose only calls are hom-ref (dropped from `by_contig`) must
    /// still report ploidy 1 — not the diploid fallback. Before this fix,
    /// deriving ploidy from the filtered map gave `2` and mis-shaped every
    /// downstream MT/MB string.
    #[test]
    fn parse_variants_by_contig_resolves_haploid_ploidy_from_homref_only_records() {
        let vcf = "\
##fileformat=VCFv4.4
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t10\t.\tA\tT\t.\t.\t.\tGT\t0
chr1\t20\t.\tC\tG\t.\t.\t.\tGT\t0
";
        let f = write_temp_vcf(vcf);
        let dict = dict_for(&[("chr1", 100)]);
        let parsed = parse_variants_by_contig(f.path(), None, &dict).unwrap();
        // All records are hom-ref haploid → no ALT-bearing variants…
        assert!(parsed.by_contig.is_empty(), "no ALT calls → empty map");
        // …but the sample is unambiguously haploid.
        assert_eq!(parsed.sample_ploidy, 1, "ploidy must come from unfiltered GTs");
    }

    /// A genotyped triploid ALT call resolves ploidy 3.
    #[test]
    fn parse_variants_by_contig_resolves_triploid_ploidy() {
        let vcf = "\
##fileformat=VCFv4.4
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t10\t.\tA\tT\t.\t.\t.\tGT\t0|0|1
";
        let f = write_temp_vcf(vcf);
        let dict = dict_for(&[("chr1", 100)]);
        let parsed = parse_variants_by_contig(f.path(), None, &dict).unwrap();
        assert_eq!(parsed.sample_ploidy, 3);
        assert_eq!(parsed.by_contig.get("chr1").map(Vec::len), Some(1));
    }

    /// A VCF that redeclares an `##INFO` ID must still parse. noodles' header
    /// parser rejects duplicate IDs outright (`duplicate INFO ID: ...`), but
    /// real-world VCFs from tools like Manta/DELLY carry repeated definitions
    /// that bcftools tolerates; holodeck must not be stricter than bcftools.
    #[test]
    fn parse_variants_by_contig_tolerates_duplicate_info_headers() {
        let vcf = "\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=100>
##INFO=<ID=BREAKSIMLENGTH,Number=1,Type=Integer,Description=\"first\">
##INFO=<ID=BREAKSIMLENGTH,Number=1,Type=Integer,Description=\"second\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t10\t.\tA\tT\t.\tPASS\tBREAKSIMLENGTH=5\tGT\t0/1
";
        let f = write_temp_vcf(vcf);
        let dict = dict_for(&[("chr1", 100)]);
        let parsed = parse_variants_by_contig(f.path(), None, &dict).unwrap();
        assert_eq!(parsed.by_contig.get("chr1").map(Vec::len), Some(1));
        assert_eq!(parsed.by_contig["chr1"][0].position, 9);
    }

    /// Duplicate `##FORMAT`/`##FILTER`/`##contig` declarations are tolerated
    /// alongside `##INFO`; all of noodles' duplicate-ID categories are covered
    /// by the same dedup pass.
    #[test]
    fn parse_variants_by_contig_tolerates_duplicate_format_filter_contig() {
        let vcf = "\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=100>
##contig=<ID=chr1,length=100>
##FILTER=<ID=q10,Description=\"first\">
##FILTER=<ID=q10,Description=\"second\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"first\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"second\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/1
";
        let f = write_temp_vcf(vcf);
        let dict = dict_for(&[("chr1", 100)]);
        let parsed = parse_variants_by_contig(f.path(), None, &dict).unwrap();
        assert_eq!(parsed.by_contig.get("chr1").map(Vec::len), Some(1));
    }

    /// `validate_vcf_sample` reads the header through the same lenient path, so
    /// it too must accept a VCF with duplicate definitions.
    #[test]
    fn validate_vcf_sample_tolerates_duplicate_info_headers() {
        let vcf = "\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=100>
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"first\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"second\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878
chr1\t10\t.\tA\tT\t.\tPASS\tDP=5\tGT\t0/1
";
        let f = write_temp_vcf(vcf);
        let sample = validate_vcf_sample(f.path(), None).unwrap();
        assert_eq!(sample, "NA12878");
    }

    #[test]
    fn structured_header_key_extracts_directive_and_id() {
        assert_eq!(
            structured_header_key("##INFO=<ID=DP,Number=1,Type=Integer>\n"),
            Some(("INFO".to_string(), "DP".to_string()))
        );
        assert_eq!(
            structured_header_key("##contig=<ID=chr1,length=100>"),
            Some(("contig".to_string(), "chr1".to_string()))
        );
        // No ID field -> not a dedup candidate.
        assert_eq!(structured_header_key("##fileformat=VCFv4.3"), None);
        assert_eq!(structured_header_key("##source=holodeck-mutate"), None);
        // ID not the first field -> conservatively not deduplicated.
        assert_eq!(structured_header_key("##INFO=<Number=1,ID=DP>"), None);
    }
}
