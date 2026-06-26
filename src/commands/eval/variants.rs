//! Variant-representation accuracy: do aligned reads carry the simulated
//! variants they should, and how confidently?
//!
//! Truth comes entirely from holodeck's own outputs: the truth VCF enumerates
//! the simulated single-base substitutions and the golden BAM gives each read's
//! true span and the actual base it carries at each site (so whether a read is
//! expected to show the alt is read from the golden sequence, independent of
//! whether the VCF is phased). For every expected substitution this pass walks
//! the *mapped* read's CIGAR to the variant's reference position and checks
//! whether the observed base matches the alternate allele, accumulating the
//! represented fraction together with the read's `MAPQ` and `AS` per
//! substitution class.
//!
//! ## Methylation framing
//!
//! Under bisulfite/EM-seq chemistry the `C->T` substitution (on a `CT`-strand
//! read) is indistinguishable from an unconverted/converted cytosine: it is
//! *intrinsically confounded* with the conversion and is reported as such
//! rather than treated as a true accuracy signal. The discriminating classes
//! are the mirror (`T->C`) and the transversions, where a methylation-aware
//! scoring mode should neither over- nor under-penalize relative to the
//! genomic truth. Classes are assigned from the read's TRUE conversion
//! direction, taken from the golden BAM (its `XG`) rather than the mapped
//! record's tags, so an aligner that omits or rewrites `XG`/`XR` cannot shift a
//! variant into the wrong class.

use std::collections::{BTreeMap, HashMap};
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use noodles::bam;

use super::cigar;
use super::edits::{self, RefCache};
use super::golden::{ConvDir, GoldenInfo, ReadKey, contig_name, int_tag};
use crate::commands::command::output_path;
use crate::sequence_dict::SequenceDictionary;
use crate::vcf::{ParsedVariants, parse_variants_by_contig};

/// Classification of a single-base substitution under a conversion direction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubClass {
    /// The freed cell, intrinsically confounded with bisulfite conversion.
    Conversion,
    /// The mirror of the freed cell (its reverse direction).
    Mirror,
    /// A transversion (purine<->pyrimidine).
    Transversion,
    /// A transition that is neither the freed nor the mirror cell.
    Other,
}

impl SubClass {
    /// Stable lowercase label used as the TSV class column.
    fn label(self) -> &'static str {
        match self {
            SubClass::Conversion => "conversion",
            SubClass::Mirror => "mirror",
            SubClass::Transversion => "transversion",
            SubClass::Other => "other",
        }
    }

    /// Whether this class is confounded with bisulfite conversion and so
    /// carries no independent accuracy signal.
    fn is_confounded(self) -> bool {
        matches!(self, SubClass::Conversion)
    }
}

/// Classify a single-base substitution under a bisulfite conversion direction.
///
/// Bases are compared case-insensitively. The freed cell (`C->T` for
/// [`ConvDir::Ct`], `G->A` for [`ConvDir::Ga`]) is [`SubClass::Conversion`];
/// its reverse is [`SubClass::Mirror`]; other transitions are
/// [`SubClass::Other`]; everything else is [`SubClass::Transversion`].
#[must_use]
pub fn classify_substitution(ref_base: u8, alt_base: u8, conv: ConvDir) -> SubClass {
    let (r, a) = (ref_base.to_ascii_uppercase(), alt_base.to_ascii_uppercase());
    let freed = match conv {
        ConvDir::Ct => (b'C', b'T'),
        ConvDir::Ga => (b'G', b'A'),
    };
    if (r, a) == freed {
        return SubClass::Conversion;
    }
    if (r, a) == (freed.1, freed.0) {
        return SubClass::Mirror;
    }
    let transition = matches!((r, a), (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C'));
    if transition { SubClass::Other } else { SubClass::Transversion }
}

/// A single-base substitution a read actually carries, with the alternate base
/// read from the golden truth sequence (not inferred from VCF phasing).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExpectedSnv {
    /// 0-based reference position.
    pub pos0: u32,
    /// Uppercased reference base.
    pub ref_base: u8,
    /// Uppercased alternate base the golden read carries at this site.
    pub alt_base: u8,
}

/// One truth substitution site. The position and reference base, plus the
/// single-base ALT base(s) the sample's genotype declares here (uppercased).
/// Which copy carries the alt — and thus which read shows it — is read per-read
/// from the golden BAM, so the VCF need not be phased; the ALT set is retained
/// only so a bisulfite conversion or sequencing error that yields some *other*
/// non-reference base is not mistaken for the variant.
#[derive(Debug, Clone)]
struct SnvSite {
    pos0: u32,
    ref_base: u8,
    alts: Vec<u8>,
}

/// Truth SNVs indexed by contig for per-read span queries.
#[derive(Debug, Default)]
pub struct VariantTruth {
    by_contig: BTreeMap<String, Vec<SnvSite>>,
}

impl VariantTruth {
    /// Load SNV truth from a VCF, resolving genotypes for `sample`.
    ///
    /// # Errors
    /// Returns an error if the VCF cannot be read or parsed.
    pub fn from_vcf(path: &Path, sample: Option<&str>) -> Result<Self> {
        // `parse_variants_by_contig` ignores the sequence dictionary argument
        // (it is reserved for future contig validation); eval has no reference
        // FASTA, so pass an empty dictionary.
        let dict = SequenceDictionary::from_entries(Vec::new());
        let parsed = parse_variants_by_contig(path, sample, &dict)
            .with_context(|| format!("Failed to read truth VCF: {}", path.display()))?;
        Ok(Self::from_parsed(&parsed))
    }

    /// Build SNV truth from already-parsed variants, keeping each site that is
    /// a single-base substitution the sample carries (single-base reference and
    /// at least one single-base ALT in the genotype). The realized allele per
    /// read is resolved later from the golden sequence, so no phasing is stored.
    fn from_parsed(parsed: &ParsedVariants) -> Self {
        let mut by_contig: BTreeMap<String, Vec<SnvSite>> = BTreeMap::new();
        for (contig, records) in &parsed.by_contig {
            let mut sites = Vec::new();
            for record in records {
                if record.ref_allele.len() != 1 {
                    continue; // SNV requires a single reference base.
                }
                // The single-base ALT base(s) the genotype carries here. Keep
                // the site only if it has at least one (i.e. the sample carries
                // a substitution); skip ref-only and non-SNV-ALT records.
                let mut alts: Vec<u8> = record
                    .genotype
                    .alleles()
                    .iter()
                    .filter_map(|allele| match allele {
                        Some(idx) if *idx > 0 => record
                            .allele_bases(*idx)
                            .filter(|b| b.len() == 1)
                            .map(|b| b[0].to_ascii_uppercase()),
                        _ => None,
                    })
                    .collect();
                alts.sort_unstable();
                alts.dedup();
                if !alts.is_empty() {
                    sites.push(SnvSite {
                        pos0: record.position,
                        ref_base: record.ref_allele[0].to_ascii_uppercase(),
                        alts,
                    });
                }
            }
            sites.sort_by_key(|s| s.pos0);
            if !sites.is_empty() {
                by_contig.insert(contig.clone(), sites);
            }
        }
        Self { by_contig }
    }

    /// Substitutions a read actually carries, read from the golden truth.
    ///
    /// For every truth SNV site within the read's true span, the golden read's
    /// own base at that site is the per-read oracle: when it matches one of the
    /// site's declared ALT bases, this read carries that variant (the value the
    /// simulator placed on the copy it was sequenced from). A read sequenced
    /// from the reference copy shows the reference base and yields nothing — so
    /// the result is correct whether or not the truth VCF is phased. Requiring
    /// an ALT match (rather than merely "non-reference") keeps a bisulfite
    /// conversion or sequencing error at the site from posing as the variant.
    #[must_use]
    pub fn expected_for_read(&self, golden: &GoldenInfo) -> Vec<ExpectedSnv> {
        let Some(sites) = self.by_contig.get(&golden.contig) else {
            return Vec::new();
        };
        let end0 = golden.end0();
        let lo = sites.partition_point(|s| s.pos0 < golden.start0);
        let mut out = Vec::new();
        for site in &sites[lo..] {
            if site.pos0 >= end0 {
                break;
            }
            if let Some(base) = golden.base_at(site.pos0)
                && site.alts.contains(&base)
            {
                out.push(ExpectedSnv { pos0: site.pos0, ref_base: site.ref_base, alt_base: base });
            }
        }
        out
    }
}

/// Per-class accumulator of expected vs represented substitutions.
#[derive(Debug, Default, Clone)]
struct ClassAcc {
    n_expected: u64,
    n_represented: u64,
    sum_mapq: u64,
    sum_as: i64,
    n_as: u64,
    confounded: bool,
}

/// Aggregated variant-representation results across all reads.
///
/// MD/NM concordance is tracked over *comparable* reads — variant-bearing
/// reads whose golden record actually carries the tag — since the golden BAM
/// only stamps MD/NM on methylation runs. Reporting concordance over all
/// variant-bearing reads would understate it to 0% whenever truth tags are
/// absent.
#[derive(Debug, Default)]
struct VariantReport {
    by_class: BTreeMap<&'static str, ClassAcc>,
    variant_bearing_reads: u64,
    md_comparable_reads: u64,
    md_concordant_reads: u64,
    nm_comparable_reads: u64,
    nm_concordant_reads: u64,
}

impl VariantReport {
    /// Record one expected substitution observed (or not) in a mapped read.
    ///
    /// In meth mode `class` is `None` when the read carried no `XG`/`XR`
    /// conversion direction, so the substitution cannot be classified and is
    /// counted under `unclassified` rather than silently mislabelled.
    fn record(
        &mut self,
        class: Option<SubClass>,
        meth: bool,
        represented: bool,
        mapq: u8,
        as_score: Option<i64>,
    ) {
        // Non-meth runs collapse to a single class; only meth distinguishes
        // conversion-confounded cells from the discriminating ones.
        let (label, confounded) = match (meth, class) {
            (false, _) => ("all", false),
            (true, Some(c)) => (c.label(), c.is_confounded()),
            (true, None) => ("unclassified", false),
        };
        let acc = self.by_class.entry(label).or_default();
        acc.confounded = confounded;
        acc.n_expected += 1;
        acc.n_represented += u64::from(represented);
        acc.sum_mapq += u64::from(mapq);
        if let Some(a) = as_score {
            acc.sum_as += a;
            acc.n_as += 1;
        }
    }

    /// Write the `<prefix>.variants.tsv` table.
    fn write_tsv(&self, output_prefix: &Path) -> Result<()> {
        let path = output_path(output_prefix, ".variants.tsv");
        let mut out = std::fs::File::create(&path)
            .with_context(|| format!("Failed to create {}", path.display()))?;
        writeln!(
            out,
            "class\tconfounded\tn_expected\tn_represented\trepresented_pct\tmean_mapq\tmean_as"
        )?;
        for (label, acc) in &self.by_class {
            let denom = acc.n_expected.max(1) as f64;
            let mean_as = if acc.n_as > 0 {
                format!("{:.2}", acc.sum_as as f64 / acc.n_as as f64)
            } else {
                "NA".to_string()
            };
            writeln!(
                out,
                "{label}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{mean_as}",
                acc.confounded,
                acc.n_expected,
                acc.n_represented,
                acc.n_represented as f64 / denom * 100.0,
                acc.sum_mapq as f64 / denom,
            )?;
        }

        // Per-read MD/NM concordance against golden, over the reads whose
        // golden record carries the tag (NA when none do).
        writeln!(out, "#variant_bearing_reads\t{}", self.variant_bearing_reads)?;
        writeln!(
            out,
            "#md_concordant_pct\t{}",
            pct_or_na(self.md_concordant_reads, self.md_comparable_reads)
        )?;
        writeln!(
            out,
            "#nm_concordant_pct\t{}",
            pct_or_na(self.nm_concordant_reads, self.nm_comparable_reads)
        )?;
        log::info!("Variant results written to: {}", path.display());
        Ok(())
    }
}

/// Format `numerator / denominator` as a percentage, or `NA` when there is
/// nothing to compare.
fn pct_or_na(numerator: u64, denominator: u64) -> String {
    if denominator == 0 {
        "NA".to_string()
    } else {
        format!("{:.2}", numerator as f64 / denominator as f64 * 100.0)
    }
}

/// Evaluate variant representation of `mapped` against `golden` + `truth`.
///
/// # Errors
/// Returns an error if the mapped BAM cannot be read or the output written.
pub fn run(
    mapped: &Path,
    golden: &HashMap<ReadKey, GoldenInfo>,
    truth: &VariantTruth,
    meth: bool,
    mut reference: Option<&mut RefCache>,
    output_prefix: &Path,
) -> Result<()> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(mapped)
        .with_context(|| format!("Failed to open BAM: {}", mapped.display()))?;
    let header = reader.read_header()?;
    let mut report = VariantReport::default();

    for result in reader.record_bufs(&header) {
        let record = result.context("Failed to read BAM record")?;
        let flags = record.flags();
        if flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        let Some(name) = record.name() else { continue };
        let key: ReadKey = (name.to_vec(), flags.is_last_segment());
        let Some(truth_aln) = golden.get(&key) else { continue };

        let expected = truth.expected_for_read(truth_aln);
        if expected.is_empty() {
            continue;
        }
        report.variant_bearing_reads += 1;

        let mapq = record.mapping_quality().map_or(0, u8::from);
        let as_score = int_tag(&record, b'A', b'S');
        let conv = if meth { truth_aln.conv_dir } else { None };

        // The mapped record represents a variant only if it is aligned to the
        // variant's contig; otherwise (unmapped / mismapped) it cannot.
        let mapped_contig = record.reference_sequence_id().and_then(|id| contig_name(&header, id));
        let mapped_start0 = record
            .alignment_start()
            .map(|p| u32::try_from(usize::from(p).saturating_sub(1)).unwrap_or(0));

        for snv in &expected {
            let represented = match (mapped_contig.as_deref(), mapped_start0) {
                (Some(contig), Some(start0)) if contig == truth_aln.contig => {
                    cigar::ref_pos_to_read_offset(record.cigar(), start0, snv.pos0)
                        .and_then(|off| record.sequence().as_ref().get(off).copied())
                        .is_some_and(|base| base.to_ascii_uppercase() == snv.alt_base)
                }
                _ => false,
            };
            // In meth mode the class needs the read's TRUE conversion direction
            // (from the golden truth, so an aligner that omits or rewrites XG/XR
            // cannot move a variant into the wrong class); a non-meth read yields
            // None and is counted as unclassified.
            let class = conv.map(|dir| classify_substitution(snv.ref_base, snv.alt_base, dir));
            report.record(class, meth, represented, mapq, as_score);
        }

        // Bisulfite-aware genomic NM/MD concordance against the golden read,
        // computed from the reference (not the aligner's tags) using the read's
        // TRUE strand, so it is comparable across aligners regardless of their
        // NM/MD convention. Requires --reference; NA without it.
        if let (Some(ref_cache), Some(m_contig), Some(m_start0)) =
            (reference.as_deref_mut(), mapped_contig.as_deref(), mapped_start0)
        {
            // Borrow each contig in turn: both calls return owned GenomicEdits,
            // so the &[u8] borrows do not overlap.
            let golden_edits = ref_cache.contig(&truth_aln.contig).map(|r| {
                edits::genomic_edits(
                    &truth_aln.sequence,
                    &truth_aln.cigar,
                    truth_aln.start0,
                    r,
                    truth_aln.conv_dir,
                )
            });
            let aligned_edits = ref_cache.contig(m_contig).map(|r| {
                edits::genomic_edits(
                    record.sequence().as_ref(),
                    record.cigar(),
                    m_start0,
                    r,
                    truth_aln.conv_dir,
                )
            });
            if let (Some(g), Some(a)) = (golden_edits, aligned_edits) {
                report.nm_comparable_reads += 1;
                report.nm_concordant_reads += u64::from(a.nm == g.nm);
                report.md_comparable_reads += 1;
                report.md_concordant_reads += u64::from(a.positions == g.positions);
            }
        }
    }

    report.write_tsv(output_prefix)
}

#[cfg(test)]
mod tests {
    use noodles::sam::alignment::record::cigar::op::{Kind, Op};
    use noodles::sam::alignment::record_buf::Cigar;

    use super::*;
    use crate::vcf::genotype::{Genotype, VariantRecord};

    #[test]
    fn classify_conversion_ct_and_ga() {
        assert_eq!(classify_substitution(b'C', b'T', ConvDir::Ct), SubClass::Conversion);
        assert_eq!(classify_substitution(b'G', b'A', ConvDir::Ga), SubClass::Conversion);
    }

    #[test]
    fn classify_mirror_is_the_reverse_of_the_freed_cell() {
        assert_eq!(classify_substitution(b'T', b'C', ConvDir::Ct), SubClass::Mirror);
        assert_eq!(classify_substitution(b'A', b'G', ConvDir::Ga), SubClass::Mirror);
    }

    #[test]
    fn classify_transversion_and_other_transition() {
        assert_eq!(classify_substitution(b'C', b'A', ConvDir::Ct), SubClass::Transversion);
        // A->G is a transition but neither the CT freed cell nor its mirror.
        assert_eq!(classify_substitution(b'A', b'G', ConvDir::Ct), SubClass::Other);
    }

    #[test]
    fn classify_is_case_insensitive() {
        assert_eq!(classify_substitution(b'c', b't', ConvDir::Ct), SubClass::Conversion);
    }

    fn snv_record(pos0: u32, ref_b: &str, alt_b: &str, gt: &str) -> VariantRecord {
        VariantRecord {
            position: pos0,
            ref_allele: ref_b.as_bytes().to_vec(),
            alt_alleles: vec![alt_b.as_bytes().to_vec()],
            genotype: Genotype::parse(gt).unwrap(),
        }
    }

    fn truth_from(records: Vec<VariantRecord>) -> VariantTruth {
        let mut parsed = ParsedVariants::default();
        parsed.by_contig.insert("chr1".to_string(), records);
        VariantTruth::from_parsed(&parsed)
    }

    fn golden_read(start0: u32, seq: &[u8], ops: &[(Kind, usize)]) -> GoldenInfo {
        let cigar = Cigar::from(ops.iter().map(|&(k, n)| Op::new(k, n)).collect::<Vec<_>>());
        GoldenInfo {
            contig: "chr1".to_string(),
            start0,
            ref_len: cigar::reference_len(&cigar),
            sequence: seq.to_vec(),
            cigar,
            conv_dir: None,
        }
    }

    #[test]
    fn expected_reads_alt_from_golden_even_when_vcf_is_unphased() {
        // SNV at ref 100 (ref C). The VCF is UNPHASED (0/1) — eval must not guess
        // a haplotype; it reads the read's actual base from the golden sequence.
        let truth = truth_from(vec![snv_record(100, "C", "T", "0/1")]);
        // Golden read covering [90, 200) as 110M; base at ref 100 = offset 10 = T.
        let mut seq = vec![b'A'; 110];
        seq[10] = b'T';
        let exp = truth.expected_for_read(&golden_read(90, &seq, &[(Kind::Match, 110)]));
        assert_eq!(exp.len(), 1);
        assert_eq!((exp[0].pos0, exp[0].ref_base, exp[0].alt_base), (100, b'C', b'T'));
    }

    #[test]
    fn expected_skips_reads_carrying_the_reference_allele() {
        // Same unphased site, but this golden read shows the REFERENCE base C at
        // 100 — sequenced from the reference copy, so nothing is expected.
        let truth = truth_from(vec![snv_record(100, "C", "T", "0/1")]);
        let mut seq = vec![b'A'; 110];
        seq[10] = b'C';
        assert!(truth.expected_for_read(&golden_read(90, &seq, &[(Kind::Match, 110)])).is_empty());
    }

    #[test]
    fn expected_honors_span_and_skips_indels() {
        let truth = truth_from(vec![
            snv_record(50, "A", "G", "1/1"),
            snv_record(100, "AT", "A", "1/1"), // deletion: not an SNV site
            snv_record(150, "C", "A", "0/1"),
        ]);
        // Golden read [60, 200) as 140M: pos 50 is before the span, the indel at
        // 100 is dropped, only site 150 (offset 90 = alt A) remains.
        let mut seq = vec![b'C'; 140];
        seq[90] = b'A';
        let exp = truth.expected_for_read(&golden_read(60, &seq, &[(Kind::Match, 140)]));
        assert_eq!(exp.len(), 1);
        assert_eq!((exp[0].pos0, exp[0].alt_base), (150, b'A'));
    }

    #[test]
    fn expected_for_read_on_unknown_contig_is_empty() {
        let truth = truth_from(vec![snv_record(10, "C", "T", "0/1")]);
        let g = GoldenInfo {
            contig: "chrX".to_string(),
            start0: 0,
            ref_len: 100,
            sequence: vec![b'T'; 100],
            cigar: Cigar::from(vec![Op::new(Kind::Match, 100)]),
            conv_dir: None,
        };
        assert!(truth.expected_for_read(&g).is_empty());
    }

    #[test]
    fn report_record_collapses_to_all_when_not_meth() {
        let mut r = VariantReport::default();
        r.record(Some(SubClass::Mirror), false, true, 60, Some(70));
        r.record(Some(SubClass::Transversion), false, false, 0, Some(40));
        assert_eq!(r.by_class.len(), 1);
        let acc = &r.by_class["all"];
        assert_eq!((acc.n_expected, acc.n_represented), (2, 1));
        assert_eq!(acc.sum_as, 110);
    }

    #[test]
    fn pct_or_na_handles_empty_denominator() {
        assert_eq!(pct_or_na(0, 0), "NA");
        assert_eq!(pct_or_na(3, 4), "75.00");
    }

    #[test]
    fn report_record_splits_classes_under_meth() {
        let mut r = VariantReport::default();
        r.record(Some(SubClass::Conversion), true, true, 10, None);
        r.record(Some(SubClass::Mirror), true, true, 60, Some(60));
        assert_eq!(r.by_class.len(), 2);
        assert!(r.by_class["conversion"].confounded);
        assert!(!r.by_class["mirror"].confounded);
    }

    #[test]
    fn report_record_buckets_meth_without_conv_dir_as_unclassified() {
        let mut r = VariantReport::default();
        r.record(None, true, true, 30, Some(50));
        assert_eq!(r.by_class.len(), 1);
        let acc = &r.by_class["unclassified"];
        assert_eq!(acc.n_expected, 1);
        assert!(!acc.confounded);
    }
}
