pub mod genotype;

use std::path::Path;

use anyhow::{Context, Result, bail};
use noodles::vcf;
use noodles::vcf::variant::record::AlternateBases;

use crate::sequence_dict::SequenceDictionary;
use genotype::{Genotype, VariantRecord};

/// Load variant records from a VCF for a given sample and contig.
///
/// Reads all records for `contig_name`, parses the GT field for the selected
/// sample, and returns a sorted list of [`VariantRecord`]s.  Records where the
/// sample's genotype is homozygous reference or entirely missing are skipped.
///
/// # Arguments
/// * `path` — Path to the VCF file.
/// * `contig_name` — Name of the contig to load variants for.
/// * `sample_name` — Sample name, or `None` to use the only sample in the VCF.
/// * `_dict` — Sequence dictionary (reserved for future validation).
///
/// # Errors
/// Returns an error if the VCF cannot be read, the sample is not found, or GT
/// fields are missing or malformed.
pub fn load_variants_for_contig(
    path: &Path,
    contig_name: &str,
    sample_name: Option<&str>,
    _dict: &SequenceDictionary,
) -> Result<Vec<VariantRecord>> {
    let mut reader = vcf::io::reader::Builder::default()
        .build_from_path(path)
        .with_context(|| format!("Failed to open VCF: {}", path.display()))?;

    let header = reader.read_header()?;
    let sample_index = resolve_sample_index(&header, sample_name)?;

    let mut variants = Vec::new();

    // TODO: Use tabix-indexed random access to read only the target contig.
    // Currently reads every record and filters by contig name, which is
    // O(total_variants) per contig call.
    for result in reader.record_bufs(&header) {
        let record = result.with_context(|| "Failed to read VCF record")?;

        // Filter to the target contig.
        let chrom = record.reference_sequence_name();
        if chrom != contig_name {
            continue;
        }

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

        // Skip genotypes with no alt alleles (hom-ref, all-missing, etc.).
        if !genotype.has_alt() {
            continue;
        }

        variants.push(VariantRecord { position, ref_allele, alt_alleles, genotype });
    }

    // Sort by position (should already be sorted in a valid VCF, but ensure).
    variants.sort_by_key(|v| v.position);

    Ok(variants)
}

/// Validate that a VCF file has the expected sample configuration.
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
pub(crate) fn validate_vcf_sample(path: &Path, sample_name: Option<&str>) -> Result<()> {
    let mut reader = vcf::io::reader::Builder::default()
        .build_from_path(path)
        .with_context(|| format!("Failed to open VCF: {}", path.display()))?;
    let header = reader.read_header()?;
    resolve_sample_index(&header, sample_name)?;
    Ok(())
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
