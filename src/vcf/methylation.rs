//! VCF MT/MB FORMAT-field schema, per-CpG ownership classifier, and
//! reader/writer for methylation-annotated VCFs.
//!
//! Each CpG on each haplotype is routed to exactly one VCF record by
//! [`classify_cpgs`]: methylation-only records for reference-coordinate
//! CpGs outside any variant span, and variant-record annotations for CpGs
//! inside or straddling a variant's alt allele.

use rand::SeedableRng;

use crate::haplotype::build_haplotypes;
use crate::vcf::genotype::VariantRecord;

/// Error returned by [`classify_cpgs`] when the input violates the
/// classifier's preconditions: phased genotypes only, and no two variants
/// whose REF spans overlap while both alt alleles land on a shared
/// haplotype.
#[derive(Debug, thiserror::Error)]
pub enum ClassifyError {
    /// A variant record has an unphased GT. Methylation truth is
    /// haplotype-specific and meaningless without phasing.
    #[error("methylation requires phased genotypes; variant at position {pos} has unphased GT")]
    UnphasedGenotype {
        /// 0-based reference position of the offending variant.
        pos: u32,
    },
    /// Two variant records have overlapping REF spans and at least one
    /// haplotype carries the alt allele at both sites, which would
    /// materialize incompatible alt spans on that haplotype. Phased
    /// overlaps on disjoint haplotypes (e.g. `1|0` + `0|1`) are accepted.
    #[error(
        "variants at positions {a_pos} and {b_pos} have overlapping REF spans on a shared haplotype"
    )]
    OverlappingVariants {
        /// 0-based reference position of the first (upstream) variant.
        a_pos: u32,
        /// 0-based reference position of the second (downstream) variant.
        b_pos: u32,
    },
}

/// Where a single CpG (haplotype-coord pair) is recorded in the
/// methylated VCF.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CpgPlacement {
    /// Stored as a standalone methylation-only record (REF=C, ALT='.')
    /// at the given reference coordinate of the top-strand C.
    Standalone { ref_pos: u32 },
    /// Stored on the variant record at the given index in the variants
    /// slice, at the given offset within that haplotype's alt allele.
    OnVariant { variant_index: usize, hap_offset: u32 },
}

/// Classify every CpG on every haplotype into [`CpgPlacement`]s.
///
/// Walks every haplotype's materialized sequence, attributes each CpG
/// dinucleotide to either a standalone methylation record (both bases come
/// from reference on at least one haplotype) or to a variant record (at
/// least one base of the pair comes from an alt allele). For CpGs where
/// the two bases come from different variants, the upstream variant (lower
/// reference position) owns the placement.
///
/// `reference` is the contig's reference sequence (uppercase). `variants`
/// are the sorted variant records for the contig.
///
/// # Errors
///
/// Returns [`ClassifyError::UnphasedGenotype`] if any variant record has an
/// unphased GT field. Returns [`ClassifyError::OverlappingVariants`] if any
/// two adjacent variant records have overlapping REF spans and at least one
/// haplotype carries the alt allele at both sites; phased overlaps on
/// disjoint haplotypes are accepted.
pub fn classify_cpgs(
    reference: &[u8],
    variants: &[VariantRecord],
) -> Result<Vec<CpgPlacement>, ClassifyError> {
    // Single-contig public API: derive ploidy from the supplied variants.
    // Multi-contig callers should use `resolve_sample_ploidy` over the full
    // VCF and pass the result to `write_contig` / `read_contig_methylation`
    // so variant-free contigs don't fall back to a default ploidy that
    // disagrees with the rest of the sample.
    let sample_ploidy = variants.iter().map(|v| v.genotype.ploidy()).max().unwrap_or(2);
    let haplotypes = build_methylation_haplotypes(variants, sample_ploidy);
    classify_cpgs_with_haplotypes(reference, variants, &haplotypes)
}

/// Build the canonical haplotype layout used everywhere methylation truth is
/// (de)serialized — `seed = 0`, max ploidy derived from `variants`. Writer
/// (`write_contig`) and reader (`read_contig_methylation`) must agree on
/// per-haplotype ordering or MT/MB bits land on the wrong indices, and
/// since `classify_cpgs` rejects unphased genotypes upfront, every layout
/// reaching this helper is deterministic regardless of seed. Centralizing
/// the construction here lets one contig produce a single
/// `Vec<Haplotype>` that all downstream methylation helpers borrow, instead
/// of each call site rebuilding it from scratch with its own
/// `SmallRng::seed_from_u64(0)`.
pub(crate) fn build_methylation_haplotypes(
    variants: &[VariantRecord],
    sample_ploidy: usize,
) -> Vec<crate::haplotype::Haplotype> {
    let mut rng = rand::rngs::SmallRng::seed_from_u64(0);
    build_haplotypes(variants, sample_ploidy, &mut rng)
}

/// Inner classifier used by [`write_contig`] when haplotypes are already
/// materialized for the same contig. Skips the redundant `build_haplotypes`
/// the public [`classify_cpgs`] entry point does for callers that don't
/// have a haplotype slice in hand. Behavior is otherwise identical and
/// honors the same `(seed=0, phased-only)` contract.
fn classify_cpgs_with_haplotypes(
    reference: &[u8],
    variants: &[VariantRecord],
    haplotypes: &[crate::haplotype::Haplotype],
) -> Result<Vec<CpgPlacement>, ClassifyError> {
    // Reject unphased GTs. Methylation truth is haplotype-specific and
    // requires a definite allele assignment for each haplotype.
    for v in variants {
        if !v.genotype.is_phased() {
            return Err(ClassifyError::UnphasedGenotype { pos: v.position });
        }
    }

    check_overlaps_on_shared_haplotypes(variants)?;

    // Fast path: no variants → every haplotype is the reference; scan once.
    // This is O(L) with zero allocations; the haplotype-aware path below
    // allocates per-variant range tables and materializes each haplotype.
    if variants.is_empty() {
        return Ok(crate::meth::find_reference_cpgs(reference)
            .into_iter()
            .map(|ref_pos| CpgPlacement::Standalone { ref_pos })
            .collect());
    }

    let mut placements: Vec<CpgPlacement> = Vec::new();

    for haplotype in haplotypes {
        // Materialize the full haplotype sequence.
        #[expect(clippy::cast_possible_truncation, reason = "reference length fits in u32")]
        let cap = haplotype.hap_position_for(reference.len() as u32) as usize;
        let (hap_bases, ref_positions, _hap_start) = haplotype.extract_fragment(reference, 0, cap);
        let len = hap_bases.len();
        if len < 2 {
            continue;
        }

        // Build per-variant hap-coord ranges for variants this haplotype carries.
        // Each entry: (variant_index_in_variants_slice, hap_start, hap_end_exclusive)
        // where hap_end = hap_start + alt_bases.len().
        //
        // A haplotype carries a variant when its allele index for that variant is
        // non-zero (alt). We use `haplotype.hap_position_for(v.position)` to map
        // the variant's reference position to haplotype coordinates.
        let hap_allele_index = haplotype.allele_index();
        // Precondition: no two variants on the same haplotype share a reference
        // position (overlapping variants). Enforced upfront in classify_cpgs's
        // validation pass (returns ClassifyError::OverlappingVariants).
        let mut var_hap_ranges: Vec<(usize, u32, u32)> = Vec::with_capacity(variants.len());
        for (vi, v) in variants.iter().enumerate() {
            // Get this haplotype's allele at this variant site. Alleles are
            // ordered by allele_index within the genotype. For a phased diploid
            // "1|0", allele_index 0 has allele Some(1), allele_index 1 has Some(0).
            let allele_num = v.genotype.alleles().get(hap_allele_index).copied().flatten();
            let Some(allele_num) = allele_num else { continue };
            if allele_num == 0 {
                // This haplotype carries the reference allele at this site.
                continue;
            }
            let Some(alt_bases) = v.allele_bases(allele_num) else { continue };
            let hap_start = haplotype.hap_position_for(v.position);
            #[expect(clippy::cast_possible_truncation, reason = "alt allele len fits u32")]
            let hap_end = hap_start + alt_bases.len() as u32;
            var_hap_ranges.push((vi, hap_start, hap_end));
        }

        // Scan materialized haplotype for CpG dinucleotides.
        for h in 0..len - 1 {
            let c0 = hap_bases[h].to_ascii_uppercase();
            let c1 = hap_bases[h + 1].to_ascii_uppercase();
            if c0 != b'C' || c1 != b'G' {
                continue;
            }

            // Cast the loop index to u32 once. Haplotype lengths are bounded
            // by reference length + net indel size, both of which fit in u32.
            #[expect(clippy::cast_possible_truncation, reason = "haplotype length fits in u32")]
            let hpos = h as u32;

            // Determine source of each base: variant or reference.
            let src_h = base_source(hpos, &var_hap_ranges);
            let src_h1 = base_source(hpos + 1, &var_hap_ranges);

            let placement = match (src_h, src_h1) {
                // Both bases come from reference → standalone at the ref
                // position of the top-strand C.
                (BaseSource::Reference, BaseSource::Reference) => {
                    CpgPlacement::Standalone { ref_pos: ref_positions[h] }
                }

                // C is from a variant, G is from reference → the variant owns.
                (BaseSource::Variant { variant_index, hap_start }, BaseSource::Reference) => {
                    let hap_offset = hpos - hap_start;
                    CpgPlacement::OnVariant { variant_index, hap_offset }
                }

                // C is from reference, G is from a variant → the variant owns.
                (BaseSource::Reference, BaseSource::Variant { variant_index, hap_start }) => {
                    let hap_offset = (hpos + 1) - hap_start;
                    CpgPlacement::OnVariant { variant_index, hap_offset }
                }

                // Both bases from the same variant → that variant owns.
                (
                    BaseSource::Variant { variant_index: vi0, hap_start: hs0 },
                    BaseSource::Variant { variant_index: vi1, hap_start: _ },
                ) if vi0 == vi1 => {
                    let hap_offset = hpos - hs0;
                    CpgPlacement::OnVariant { variant_index: vi0, hap_offset }
                }

                // Both bases from different variants → upstream wins (lower
                // reference position), which is the variant with the smaller
                // index because variants are sorted by position.
                (
                    BaseSource::Variant { variant_index: vi0, hap_start: hs0 },
                    BaseSource::Variant { variant_index: vi1, hap_start: hs1 },
                ) => {
                    // variants slice is sorted by position; smaller index ⟹ upstream.
                    // vi0 == vi1 is handled by the same-variant arm above.
                    if vi0 < vi1 {
                        let hap_offset = hpos - hs0;
                        CpgPlacement::OnVariant { variant_index: vi0, hap_offset }
                    } else {
                        let hap_offset = (hpos + 1) - hs1;
                        CpgPlacement::OnVariant { variant_index: vi1, hap_offset }
                    }
                }
            };

            placements.push(placement);
        }
    }

    // Deduplicate and sort. Sort key: (effective_ref_pos, discriminant, hap_offset).
    // Standalone sorts by ref_pos; OnVariant sorts by variants[vi].position
    // then hap_offset, so placements interleave in reference order.
    // Discriminant 0 (Standalone) < 1 (OnVariant) breaks ties at the same position.
    placements.sort_unstable_by_key(|p| sort_key(p, variants));
    placements.dedup();
    Ok(placements)
}

/// Reject overlapping REF spans only when the two variants can coexist on
/// the same haplotype. Phased records like `1|0` + `0|1` can overlap in
/// reference coordinates without ever materializing both alt spans on one
/// haplotype, and the haplotype builder handles them correctly. The
/// genuine conflict is one haplotype carrying both alt alleles.
///
/// Variants are assumed sorted by position (current API contract — see
/// [`crate::vcf::parse_variants_by_contig`]). For each haplotype, this
/// tracks the furthest ALT-span end seen so far; a new variant whose
/// haplotype is ALT and whose start falls before that end is a genuine
/// overlap. Tracking per-haplotype catches non-adjacent overlaps where an
/// intervening variant on a different haplotype hides a long upstream
/// REF span from a simple adjacent-pair scan.
fn check_overlaps_on_shared_haplotypes(variants: &[VariantRecord]) -> Result<(), ClassifyError> {
    let max_ploidy = variants.iter().map(|v| v.genotype.ploidy()).max().unwrap_or(0);
    // Per haplotype: (end_exclusive, start_pos) of the most recent ALT span.
    let mut last_alt_span: Vec<Option<(u32, u32)>> = vec![None; max_ploidy];

    for v in variants {
        #[expect(clippy::cast_possible_truncation, reason = "ref allele len fits u32")]
        let v_end = v.position + v.ref_allele.len() as u32;
        for (hi, slot) in last_alt_span.iter_mut().enumerate() {
            let carries_alt =
                v.genotype.alleles().get(hi).copied().flatten().is_some_and(|x| x != 0);
            if !carries_alt {
                continue;
            }
            if let Some((prev_end, prev_pos)) = *slot
                && prev_end > v.position
            {
                return Err(ClassifyError::OverlappingVariants {
                    a_pos: prev_pos,
                    b_pos: v.position,
                });
            }
            *slot = Some((v_end, v.position));
        }
    }
    Ok(())
}

/// The source of a single base at haplotype coordinate `h`.
#[derive(Debug, Clone, Copy)]
enum BaseSource {
    /// Base comes from the reference. The reference position is carried by
    /// [`ref_positions`] in the caller, not stored in this variant.
    Reference,
    /// Base comes from variant `variant_index`'s alt allele. `hap_start` is
    /// the haplotype-coordinate start of that variant's alt span.
    Variant { variant_index: usize, hap_start: u32 },
}

/// Determine whether haplotype coordinate `h` falls inside any variant's
/// alt-allele span. Returns `BaseSource::Variant` for the first matching
/// range, or `BaseSource::Reference` if no variant spans `h`.
fn base_source(h: u32, var_hap_ranges: &[(usize, u32, u32)]) -> BaseSource {
    for &(vi, hap_start, hap_end) in var_hap_ranges {
        if h >= hap_start && h < hap_end {
            return BaseSource::Variant { variant_index: vi, hap_start };
        }
    }
    BaseSource::Reference
}

/// Compute a sort key for a [`CpgPlacement`] so that placements are ordered
/// by their effective reference position. Standalone records use their own
/// `ref_pos`; `OnVariant` records use the variant's reference position as the
/// primary key and `hap_offset` as the tertiary key.
///
/// The discriminant (second tuple element) is 0 for `Standalone` and 1 for
/// `OnVariant`, so when both appear at the same reference position the
/// reference-coordinate methylation sorts before the variant-allele methylation.
fn sort_key(placement: &CpgPlacement, variants: &[VariantRecord]) -> (u32, u8, u32) {
    match placement {
        CpgPlacement::Standalone { ref_pos } => (*ref_pos, 0, 0),
        CpgPlacement::OnVariant { variant_index, hap_offset } => {
            (variants[*variant_index].position, 1, *hap_offset)
        }
    }
}

/// Default sample column name used by `methylate` when no VCF sample is provided.
pub const DEFAULT_METHYLATE_SAMPLE: &str = "METHYLATE";

/// Write a VCF header for a methylation-annotated VCF to `writer`.
///
/// The header includes:
/// - `##fileformat=VCFv4.4`
/// - `##holodeckVersion=<version>`
/// - `##holodeckCommand=<command_line>`
/// - `##FORMAT` lines for GT, MT, and MB
/// - One `##contig` line per entry in `dict`
/// - The `#CHROM` column header with one sample column
///
/// `sample` is the sample name written in the column header.  If `None`,
/// [`DEFAULT_METHYLATE_SAMPLE`] is used.
///
/// # Errors
///
/// Returns an error if any write to `writer` fails.
pub fn write_vcf_header<W: std::io::Write>(
    writer: &mut W,
    dict: &crate::sequence_dict::SequenceDictionary,
    sample: Option<&str>,
    version: &str,
    command_line: &str,
) -> anyhow::Result<()> {
    let sample_name = sample.unwrap_or(DEFAULT_METHYLATE_SAMPLE);
    writeln!(writer, "##fileformat=VCFv4.4")?;
    writeln!(writer, "##holodeckVersion={version}")?;
    writeln!(writer, "##holodeckCommand={command_line}")?;
    writeln!(writer, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
    writeln!(
        writer,
        "##FORMAT=<ID=MT,Number=.,Type=String,\
         Description=\"Methylation state, top strand. \
         Per-haplotype pipe-separated bitstring: 1=methylated, \
         0=unmethylated, .=haplotype carries REF or no owned CpG.\">"
    )?;
    writeln!(
        writer,
        "##FORMAT=<ID=MB,Number=.,Type=String,\
         Description=\"Methylation state, bottom strand. \
         Per-haplotype pipe-separated bitstring: 1=methylated, \
         0=unmethylated, .=haplotype carries REF or no owned CpG.\">"
    )?;
    for seq in dict.iter() {
        writeln!(writer, "##contig=<ID={},length={}>", seq.name(), seq.length())?;
    }
    writeln!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}")?;
    Ok(())
}

/// Discriminated union for `write_contig`'s internal emit list.
///
/// Carries either the 0-based reference position of a standalone CpG record
/// or the index into the variants slice for a variant record. Used as the
/// payload of the `(pos_1based, discriminant, RecordKind)` tuples that
/// `write_contig` sorts before emitting.
enum WriteKind {
    /// Case-1 methylation-only record. Payload is the 0-based reference
    /// position of the top-strand C.
    Standalone(u32),
    /// Case-2 variant record. Payload is the variant's index in the variants
    /// slice.
    Variant(usize),
}

/// Write all methylation-only records for one contig as VCF text rows.
///
/// Each row has the form:
/// - Standalone (Case 1, REF=C ALT=.):
///   `<chrom>\t<POS>\t.\tC\t.\t.\t.\t.\tMT:MB\t<mt>:<mb>`
/// - Variant (Case 2, with GT):
///   `<chrom>\t<POS>\t.\t<REF>\t<ALT>\t.\t.\t.\tGT:MT:MB\t<gt>:<mt>:<mb>`
///
/// `POS` is 1-based. `<mt>`/`<mb>` are `|`-separated per-haplotype bitstrings:
/// `"1"` = methylated CpG, `"0"` = unmethylated CpG, `"."` = haplotype
/// carries REF or no CpGs in the alt allele. For standalone records the
/// per-haplotype value is always a single bit (`"0"` or `"1"`).
///
/// Standalone and variant records are interleaved in ascending reference
/// position order. Header writing is the caller's responsibility.
///
/// # Errors
///
/// Returns an error if any variant has an unphased genotype, two variants
/// have overlapping REF spans (see [`classify_cpgs`]), or any write to
/// `writer` fails.
pub fn write_contig<W: std::io::Write>(
    writer: &mut W,
    chrom: &str,
    reference: &[u8],
    variants: &[VariantRecord],
    methylation: &crate::meth::ContigMethylation,
    sample_ploidy: usize,
) -> anyhow::Result<()> {
    // One canonical haplotype build per contig — see
    // [`build_methylation_haplotypes`]. The classifier, the per-variant
    // CpG-offset table, and the per-haplotype bit-string emitters below
    // all borrow this same slice instead of rebuilding it from scratch
    // with their own `SmallRng::seed_from_u64(0)`. `sample_ploidy` (passed
    // from the caller's whole-VCF resolution) sizes the haplotype slice
    // even on variant-free contigs, so haploid/triploid samples don't
    // silently get diploid MT/MB shapes here.
    let haplotypes = build_methylation_haplotypes(variants, sample_ploidy);

    let placements = classify_cpgs_with_haplotypes(reference, variants, &haplotypes)?;

    // Pre-compute per-variant, per-haplotype absolute CpG hap-coords so we can
    // look them up when writing variant records.
    let var_hap_coords =
        per_variant_per_hap_cpg_offsets_with_haplotypes(reference, variants, &haplotypes);

    // Step 1: Build sorted emit list.
    //
    // Build the full interleaved emit order: variant rows and standalone CpG
    // rows, sorted by their 1-based VCF POS. Ties are broken by discriminant:
    // standalone (0) before variant (1) — matching classify_cpgs's sort_key.
    //
    // All variants get a row in the output (even those with no owned CpGs),
    // so that simulate can read the GT from every variant's record.
    let mut records: Vec<(u32, u8, WriteKind)> = Vec::new();

    // Collect standalone CpG positions (deduplicated by classify_cpgs already)
    // and the first occurrence of each variant index from OnVariant placements.
    let mut seen_variant: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for placement in &placements {
        match placement {
            CpgPlacement::Standalone { ref_pos } => {
                records.push((*ref_pos + 1, 0, WriteKind::Standalone(*ref_pos)));
            }
            CpgPlacement::OnVariant { variant_index, .. } => {
                if seen_variant.insert(*variant_index) {
                    let pos_1based = variants[*variant_index].position + 1;
                    records.push((pos_1based, 1, WriteKind::Variant(*variant_index)));
                }
            }
        }
    }
    // Add any variants that have no owned CpGs (zero OnVariant placements).
    // These still need a row in the output for simulate to read the GT.
    for (vi, v) in variants.iter().enumerate() {
        if !seen_variant.contains(&vi) {
            let pos_1based = v.position + 1;
            records.push((pos_1based, 1, WriteKind::Variant(vi)));
        }
    }
    // Sort by (1-based POS, discriminant): standalone (0) before variant (1).
    records.sort_unstable_by_key(|(pos, disc, _)| (*pos, *disc));

    // Step 2: Emit each record.
    for (pos_1based, _, kind) in &records {
        match kind {
            WriteKind::Standalone(ref_pos) => {
                let mt = per_hap_top_bit_string(methylation, &haplotypes, *ref_pos);
                let mb = per_hap_bot_bit_string(methylation, &haplotypes, *ref_pos);
                writeln!(writer, "{chrom}\t{pos_1based}\t.\tC\t.\t.\t.\t.\tMT:MB\t{mt}:{mb}")?;
            }
            WriteKind::Variant(vi) => {
                let v = &variants[*vi];
                let ref_str = std::str::from_utf8(&v.ref_allele).unwrap_or("N");
                let alt_str = v
                    .alt_alleles
                    .iter()
                    .map(|a| std::str::from_utf8(a).unwrap_or("N"))
                    .collect::<Vec<_>>()
                    .join(",");
                let gt_str = format_genotype(&v.genotype);
                // Per-haplotype absolute hap-coords for this variant
                // (hap_coords[variant_index][hap_idx] → sorted absolute hap coords).
                let hap_coords_for_var = &var_hap_coords[*vi];
                let top_meth = format_variant_meth_field(v, hap_coords_for_var, methylation, false);
                let bot_meth = format_variant_meth_field(v, hap_coords_for_var, methylation, true);
                writeln!(
                    writer,
                    "{chrom}\t{pos_1based}\t.\t{ref_str}\t{alt_str}\t.\t.\t.\tGT:MT:MB\t{gt_str}:{top_meth}:{bot_meth}"
                )?;
            }
        }
    }
    Ok(())
}

/// Composite identity of one variant row in the methylation VCF — the
/// quadruple `(POS_1based, REF, ALT-list, GT)`. Used as the key in the
/// writer-side index that [`read_contig_methylation`] consults when decoding
/// MT/MB rows: keying by POS alone collapses multi-allelic / decomposed sites
/// at the same position and routes downstream rows through the wrong
/// `VariantRecord`. Strings here mirror exactly what [`write_contig`] emits
/// in cols 4 (REF), 5 (ALT) and the GT subfield of col 10, so writer and
/// reader hash to the same key without any normalization step.
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
struct VariantKey {
    /// 1-based VCF POS.
    pos_1based: u32,
    /// REF allele as written in column 4.
    ref_allele: String,
    /// ALT allele list as written in column 5 (comma-joined for multi-allelic).
    alt_alleles: String,
    /// Phased GT string as written by [`format_genotype`] (e.g. `1|0`).
    gt: String,
}

/// Build a [`VariantKey`] for a [`VariantRecord`] using the exact same
/// string conversions [`write_contig`] applies when emitting the row.
/// Centralizing this keeps the writer's index construction and the writer's
/// row formatting in lock-step.
fn variant_key_from_record(v: &VariantRecord) -> VariantKey {
    let ref_allele = std::str::from_utf8(&v.ref_allele).unwrap_or("N").to_string();
    let alt_alleles = v
        .alt_alleles
        .iter()
        .map(|a| std::str::from_utf8(a).unwrap_or("N"))
        .collect::<Vec<_>>()
        .join(",");
    VariantKey {
        pos_1based: v.position + 1,
        ref_allele,
        alt_alleles,
        gt: format_genotype(&v.genotype),
    }
}

/// Format a [`Genotype`] as a phased `allele1|allele2|...` string.
/// Missing alleles (`.`) are represented as `"."`. All GTs in the methylation
/// VCF are phased (enforced by [`classify_cpgs`]), so `|` is always used.
fn format_genotype(gt: &crate::vcf::genotype::Genotype) -> String {
    gt.alleles()
        .iter()
        .map(|a| match a {
            Some(idx) => idx.to_string(),
            None => ".".to_string(),
        })
        .collect::<Vec<_>>()
        .join("|")
}

/// For variant record `v` at index `vi`, return the per-haplotype `MT` or
/// `MB` field value as a `|`-separated string. Each haplotype's entry is:
/// - `"."` if it carries REF (allele 0) or has a missing allele
/// - A bitstring of `'0'`/`'1'` characters (one per owned CpG, 5'→3' order)
///   if it carries ALT, or `"."` if the alt has no owned CpGs
///
/// `hap_coords_for_var` is indexed by haplotype; each entry is the sorted
/// list of **absolute** hap-coords of owned top-C positions for this variant
/// on that haplotype (see [`per_variant_per_hap_cpg_offsets`]).
///
/// `bottom_strand` selects which bitmap to read: `false` = top (MT),
/// `true` = bottom (MB). Query top-strand at `hap_coord`, bottom-strand at
/// `hap_coord + 1`.
fn format_variant_meth_field(
    v: &VariantRecord,
    hap_coords_for_var: &[Vec<u32>],
    methylation: &crate::meth::ContigMethylation,
    bottom_strand: bool,
) -> String {
    (0..methylation.len())
        .map(|hi| {
            // Determine this haplotype's allele at this variant.
            let allele = v.genotype.alleles().get(hi).copied().flatten();
            match allele {
                None | Some(0) => ".".to_string(), // REF or missing → no methylation info
                Some(_) => {
                    // ALT haplotype — read methylation bits at owned CpG hap-coords.
                    let hap_coords = if hi < hap_coords_for_var.len() {
                        hap_coords_for_var[hi].as_slice()
                    } else {
                        &[]
                    };
                    if hap_coords.is_empty() {
                        // No CpGs owned by this variant on this haplotype.
                        ".".to_string()
                    } else {
                        let table = methylation.table_for(hi);
                        hap_coords
                            .iter()
                            .map(|&hap_coord| {
                                // hap_coord is the absolute hap-coord of the top-C;
                                // query top-strand at hap_coord, bottom-strand at hap_coord + 1.
                                let query_coord =
                                    if bottom_strand { hap_coord + 1 } else { hap_coord };
                                if table.is_methylated(query_coord, bottom_strand) {
                                    '1'
                                } else {
                                    '0'
                                }
                            })
                            .collect()
                    }
                }
            }
        })
        .collect::<Vec<_>>()
        .join("|")
}

/// For each (variant_index, haplotype_index) pair, return the sorted list of
/// **absolute** haplotype-coordinate positions where a CpG's top-C lives in
/// that haplotype's alt allele at that variant.
///
/// **Encoding note:** these are absolute hap-coords, NOT offsets relative to
/// the variant's `hap_start`. This differs from [`classify_cpgs`]'s
/// `OnVariant { hap_offset }` payload, which uses relative offsets. The
/// absolute encoding is used here so the writer can look up methylation
/// bits via `MethylationTable::is_methylated(hap_coord, ...)` directly,
/// without re-computing `hap_position_for` per CpG.
///
/// Returns `hap_coords[variant_index][hap_index]` = sorted `Vec<u32>` of
/// absolute hap-coord top-C positions owned by that variant on that haplotype
/// (i.e., CpGs where the upstream-wins rule attributes them to `variant_index`).
/// Haplotypes carrying REF (allele 0) or missing alleles have empty lists.
///
/// Mirrors the inner loop of [`classify_cpgs`] but tracks per-haplotype
/// attribution rather than producing a deduplicated placement list. Takes a
/// pre-built `&[Haplotype]` slice from [`build_methylation_haplotypes`] —
/// callers share that slice with [`classify_cpgs_with_haplotypes`] and the
/// writer's per-haplotype bit-string helpers so the `(seed=0)` haplotype
/// layout is materialized exactly once per contig.
fn per_variant_per_hap_cpg_offsets_with_haplotypes(
    reference: &[u8],
    variants: &[VariantRecord],
    haplotypes: &[crate::haplotype::Haplotype],
) -> Vec<Vec<Vec<u32>>> {
    if variants.is_empty() {
        return Vec::new();
    }

    // Result: outer indexed by variant_index, inner by hap_allele_index.
    let n_vars = variants.len();
    let n_haps = haplotypes.len();
    // hap_coords[vi][hi] = sorted Vec<u32> of absolute hap coords of owned top-C
    let mut hap_coords: Vec<Vec<Vec<u32>>> = vec![vec![Vec::new(); n_haps]; n_vars];

    for haplotype in haplotypes {
        let hi = haplotype.allele_index(); // haplotype's index in the allele order

        #[expect(clippy::cast_possible_truncation, reason = "reference length fits in u32")]
        let cap = haplotype.hap_position_for(reference.len() as u32) as usize;
        let (hap_bases, _ref_positions, _hap_start) = haplotype.extract_fragment(reference, 0, cap);
        let len = hap_bases.len();
        if len < 2 {
            continue;
        }

        // Build per-variant hap-coord ranges for variants this haplotype carries.
        let mut var_hap_ranges: Vec<(usize, u32, u32)> = Vec::with_capacity(variants.len());
        for (vi, v) in variants.iter().enumerate() {
            let allele_num = v.genotype.alleles().get(hi).copied().flatten();
            let Some(allele_num) = allele_num else { continue };
            if allele_num == 0 {
                continue;
            }
            let Some(alt_bases) = v.allele_bases(allele_num) else { continue };
            let hap_start = haplotype.hap_position_for(v.position);
            #[expect(clippy::cast_possible_truncation, reason = "alt allele len fits u32")]
            let hap_end = hap_start + alt_bases.len() as u32;
            var_hap_ranges.push((vi, hap_start, hap_end));
        }

        // Scan for CpG dinucleotides and assign them to the owning variant
        // using the same upstream-wins rule as classify_cpgs.
        for h in 0..len - 1 {
            let c0 = hap_bases[h].to_ascii_uppercase();
            let c1 = hap_bases[h + 1].to_ascii_uppercase();
            if c0 != b'C' || c1 != b'G' {
                continue;
            }

            #[expect(clippy::cast_possible_truncation, reason = "haplotype length fits in u32")]
            let hpos = h as u32;

            let src_h = base_source(hpos, &var_hap_ranges);
            let src_h1 = base_source(hpos + 1, &var_hap_ranges);

            // Determine the owning variant index using the same upstream-wins
            // rule as classify_cpgs. None means both bases are from reference
            // (standalone CpG — skip here, handled by write_contig directly).
            let owning_vi: Option<usize> = match (src_h, src_h1) {
                (BaseSource::Reference, BaseSource::Reference) => None,
                (BaseSource::Variant { variant_index, .. }, BaseSource::Reference)
                | (BaseSource::Reference, BaseSource::Variant { variant_index, .. }) => {
                    Some(variant_index)
                }
                (
                    BaseSource::Variant { variant_index: vi0, .. },
                    BaseSource::Variant { variant_index: vi1, .. },
                ) if vi0 == vi1 => Some(vi0),
                (
                    BaseSource::Variant { variant_index: vi0, .. },
                    BaseSource::Variant { variant_index: vi1, .. },
                ) => Some(vi0.min(vi1)), // upstream wins = smaller index
            };

            if let Some(vi) = owning_vi
                && vi < n_vars
                && hi < n_haps
            {
                hap_coords[vi][hi].push(hpos); // absolute hap coord of top-C
            }
        }
    }

    // Sort each inner vec (they're appended in scan order which is already
    // ascending, but sort for correctness guarantee).
    for var_hap_coords in &mut hap_coords {
        for hc in var_hap_coords.iter_mut() {
            hc.sort_unstable();
        }
    }

    hap_coords
}

/// Error returned by [`read_contig_methylation`] when the input VCF stream
/// violates the writer's encoding contract.
#[derive(Debug, thiserror::Error)]
pub enum ReadError {
    /// A record's `MT`/`MB` bitstring length disagrees with the expected
    /// number of CpGs for that haplotype's alt allele.
    #[error("variant at pos {pos} hap {hap}: MT/MB length {actual} != expected {expected}")]
    MtMbLengthMismatch {
        /// 1-based VCF POS of the offending variant.
        pos: u32,
        /// Zero-based haplotype index.
        hap: usize,
        /// Length of the parsed bitstring.
        actual: usize,
        /// Expected length (number of owned CpGs for this haplotype).
        expected: usize,
    },
    /// A record's `MT`/`MB` contains a non-`0`/`1`/`.` character.
    #[error("variant at pos {pos} hap {hap}: invalid MT/MB character {ch:?}")]
    InvalidMtMbChar {
        /// 1-based VCF POS of the offending variant.
        pos: u32,
        /// Zero-based haplotype index.
        hap: usize,
        /// The offending character.
        ch: char,
    },
    /// A record's `MT` or `MB` field had a `|`-separated entry count that
    /// did not match the per-sample ploidy. Silent truncation or padding
    /// could corrupt per-haplotype methylation, so the reader rejects the
    /// record up front.
    #[error("variant at pos {pos}: {field} entry count {actual} != expected ploidy {expected}")]
    PloidyEntryCountMismatch {
        /// 1-based VCF POS of the offending record.
        pos: u32,
        /// Which FORMAT field carried the bad count (`"MT"` or `"MB"`).
        field: &'static str,
        /// Number of `|`-separated entries actually present.
        actual: usize,
        /// Expected number of entries (per-sample ploidy).
        expected: usize,
    },
    /// Underlying I/O or parse error.
    #[error("malformed VCF record at line {line}: {message}")]
    MalformedRecord {
        /// 1-based line number in the byte slice.
        line: usize,
        /// Human-readable description of what was malformed.
        message: String,
    },
}

/// Apply a single `'0'`/`'1'`/`'.'`-encoded bit from a standalone MT/MB
/// record to the appropriate table position.
///
/// `bit_str` is the per-haplotype value from the `|`-separated field (`"0"`,
/// `"1"`, or `"."`). `hap_top_c_pos` is the 0-based **haplotype**-coordinate
/// position of the top-strand C of the CpG (already translated from the
/// VCF's reference POS via `Haplotype::hap_position_for`). `bottom_strand`
/// selects which bitmap and offset to use: `false` → top at `hap_top_c_pos`;
/// `true` → bottom at `hap_top_c_pos + 1`. Returns an error if `bit_str`
/// contains any character other than `'0'`, `'1'`, or `'.'`.
fn apply_standalone_bit(
    tables: &mut [crate::meth::MethylationTable],
    hi: usize,
    hap_top_c_pos: usize,
    bit_str: &str,
    bottom_strand: bool,
    pos_1based: u32,
) -> Result<(), ReadError> {
    match bit_str {
        "1" => {
            if bottom_strand {
                tables[hi].set_bottom(hap_top_c_pos + 1, true);
            } else {
                tables[hi].set_top(hap_top_c_pos, true);
            }
        }
        "0" | "." => {} // unmethylated or missing — no-op
        other => {
            return Err(ReadError::InvalidMtMbChar {
                pos: pos_1based,
                hap: hi,
                ch: other.chars().next().unwrap_or('?'),
            });
        }
    }
    Ok(())
}

/// Apply one haplotype's per-variant bitstring to `tables`, validating length
/// and character validity. `expected_coords` holds the absolute hap-coords of
/// each owned CpG's top-strand C. `bottom_strand` selects the bitmap: `false`
/// → top (set at `hap_coord`); `true` → bottom (set at `hap_coord + 1`).
/// The strand offset is applied internally — callers pass the top-C coordinate
/// regardless of strand.
///
/// Returns `Ok(())` if `bits` is `"."` (skip) or a valid `'0'`/`'1'`
/// bitstring of the right length. Errors on length mismatch or bad character.
fn apply_variant_bits(
    tables: &mut [crate::meth::MethylationTable],
    hi: usize,
    bits: &str,
    expected_coords: &[u32],
    bottom_strand: bool,
    pos_1based: u32,
) -> Result<(), ReadError> {
    if bits == "." {
        // `.` is only valid when this haplotype owns zero CpGs at this
        // variant (REF allele, missing allele, or an ALT whose owned-CpG
        // list is empty). When the variants slice says the haplotype does
        // own CpGs, treat `.` as a length-mismatch parse error instead of
        // silently dropping the methylation truth.
        return if expected_coords.is_empty() {
            Ok(())
        } else {
            Err(ReadError::MtMbLengthMismatch {
                pos: pos_1based,
                hap: hi,
                actual: 0,
                expected: expected_coords.len(),
            })
        };
    }
    // Validate characters.
    for ch in bits.chars() {
        if ch != '0' && ch != '1' {
            return Err(ReadError::InvalidMtMbChar { pos: pos_1based, hap: hi, ch });
        }
    }
    // Validate length.
    if bits.len() != expected_coords.len() {
        return Err(ReadError::MtMbLengthMismatch {
            pos: pos_1based,
            hap: hi,
            actual: bits.len(),
            expected: expected_coords.len(),
        });
    }
    // Set bits.
    for (bit_idx, ch) in bits.chars().enumerate() {
        if ch == '1' {
            let base_coord = expected_coords[bit_idx] as usize;
            if bottom_strand {
                tables[hi].set_bottom(base_coord + 1, true);
            } else {
                tables[hi].set_top(base_coord, true);
            }
        }
    }
    Ok(())
}

/// Parse MT/MB records from a per-contig VCF byte slice back into a
/// [`crate::meth::ContigMethylation`].
///
/// The byte slice should contain only data lines for `chrom` (no header
/// lines), one record per line, formatted by [`write_contig`]. Blank lines
/// are skipped.
///
/// Supports both record kinds produced by the writer:
///
/// - **Standalone** (`ALT='.'`, `FORMAT=MT:MB`): a single reference-coordinate
///   CpG row. Each per-haplotype field in `MT`/`MB` is a single `0` or `1`.
/// - **Variant** (`FORMAT=GT:MT:MB`): a variant row. Per-haplotype entries
///   in `MT`/`MB` are either `"."` (haplotype carries REF or has no owned
///   CpGs) or a multi-character bitstring of `'0'`/`'1'` chars.
///
/// `reference` is the contig's full reference sequence (uppercase).
/// `variants` is the sorted variant list for the contig (used to reconstruct
/// per-haplotype CpG layouts).
///
/// Ploidy is derived from `variants` (the max across all records) to match
/// what [`write_contig`] uses. For a variant-free contig the function falls
/// back to ploidy `2`, again matching the writer's `unwrap_or(2)`. Keeping
/// the derivation internal removes a parameter that callers were always
/// computing the same way and that — when supplied — was silently overridden
/// for variant-bearing contigs anyway.
///
/// # Errors
///
/// Returns [`ReadError::MtMbLengthMismatch`] if a variant record's bitstring
/// length disagrees with the expected number of owned CpGs. Returns
/// [`ReadError::InvalidMtMbChar`] if a bitstring contains a character other
/// than `'0'` or `'1'`. Returns [`ReadError::PloidyEntryCountMismatch`] if a
/// record's `MT` or `MB` field does not contain exactly one `|`-separated
/// entry per haplotype. Returns [`ReadError::MalformedRecord`] if a line
/// cannot be parsed at all.
pub fn read_contig_methylation(
    vcf_bytes: &[u8],
    chrom: &str,
    reference: &[u8],
    variants: &[VariantRecord],
    sample_ploidy: usize,
) -> Result<crate::meth::ContigMethylation, ReadError> {
    use crate::meth::{ContigMethylation, MethylationTable};

    // One canonical haplotype build per contig (same `(seed=0)` contract
    // `write_contig` uses), shared with the offset table below.
    // `sample_ploidy` is the caller-supplied whole-VCF ploidy; on
    // variant-free contigs it sizes the haplotype slice so we don't fall
    // back to `2` and silently get the wrong MT/MB shape for haploid /
    // triploid samples.
    let haplotypes = build_methylation_haplotypes(variants, sample_ploidy);

    // Build empty per-haplotype tables sized to each haplotype's materialized
    // length. When there are no variants, every haplotype is the reference.
    let mut tables: Vec<MethylationTable> = if haplotypes.is_empty() {
        (0..sample_ploidy).map(|_| MethylationTable::with_len(reference.len())).collect()
    } else {
        #[expect(clippy::cast_possible_truncation, reason = "reference length fits in u32")]
        let ref_len_u32 = reference.len() as u32;
        haplotypes
            .iter()
            .map(|hap| MethylationTable::with_len(hap.hap_position_for(ref_len_u32) as usize))
            .collect()
    };

    // Pre-compute per-variant, per-haplotype absolute CpG hap-coords using
    // the haplotype slice we just built — same `(seed=0)` contract.
    let var_hap_coords =
        per_variant_per_hap_cpg_offsets_with_haplotypes(reference, variants, &haplotypes);

    // Pre-build a composite `(POS, REF, ALT, GT)` → variant index map for
    // O(1) lookups in `apply_variant_record`. Keying by POS alone collapses
    // multi-allelic / decomposed sites that share a position (e.g. one row
    // per ALT on different haplotypes), routing every record at that POS
    // through the wrong `VariantRecord` and corrupting per-haplotype state.
    // The composite key matches the exact `(REF, ALT, GT)` strings emitted by
    // [`write_contig`], so writer / reader round-trips disambiguate even when
    // multiple records share POS.
    let pos_to_vi: std::collections::HashMap<VariantKey, usize> =
        variants.iter().enumerate().map(|(i, v)| (variant_key_from_record(v), i)).collect();

    // Fail fast on invalid UTF-8. Silently coercing to an empty body would
    // return an all-zeros methylation table, indistinguishable from a
    // legitimately unmethylated contig — exactly the kind of "silent truth
    // erasure" we work hard to surface everywhere else.
    let text = std::str::from_utf8(vcf_bytes).map_err(|e| ReadError::MalformedRecord {
        line: 0,
        message: format!("VCF body is not valid UTF-8: {e}"),
    })?;
    for (line_idx, line) in text.lines().enumerate() {
        let line_no = line_idx + 1;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        // Minimal tab-split: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE
        let cols: Vec<&str> = line.splitn(10, '\t').collect();
        if cols.len() < 10 {
            return Err(ReadError::MalformedRecord {
                line: line_no,
                message: format!("expected 10 tab-separated columns, found {}", cols.len()),
            });
        }
        if cols[0] != chrom {
            continue; // defensive: skip records for other contigs
        }
        let pos_1based: u32 = cols[1].parse().map_err(|_| ReadError::MalformedRecord {
            line: line_no,
            message: format!("invalid POS field {:?}", cols[1]),
        })?;
        // VCF POS is 1-based; standalone-record decoding does `pos_1based - 1`
        // and the per-haplotype lookup expects a non-negative reference offset.
        // Reject `0` here so a malformed input becomes a clean parse error
        // rather than a `u32` wraparound and a wildly out-of-range hap-coord.
        if pos_1based == 0 {
            return Err(ReadError::MalformedRecord {
                line: line_no,
                message: "POS must be 1-based; got 0".to_string(),
            });
        }

        let ref_str = cols[3];
        let alt = cols[4];
        let format_col = cols[8];
        let sample = cols[9];

        if alt == "." {
            apply_standalone_record(
                &mut tables,
                &haplotypes,
                pos_1based,
                format_col,
                sample,
                line_no,
            )?;
        } else {
            apply_variant_record(
                &mut tables,
                pos_1based,
                ref_str,
                alt,
                format_col,
                sample,
                &pos_to_vi,
                &var_hap_coords,
                line_no,
            )?;
        }
    }

    Ok(ContigMethylation::from_tables(tables))
}

/// Process one standalone (`ALT='.'`) methylation record, updating `tables`
/// in place. `pos_1based` is the 1-based VCF POS. `haplotypes` carries the
/// per-haplotype variant layout so the reference-coordinate POS can be
/// translated to each haplotype's local coordinate before the bit is set —
/// otherwise upstream indels would shift every downstream standalone CpG and
/// the bits would land on the wrong positions.
///
/// Returns a [`ReadError`] on malformed input.
fn apply_standalone_record(
    tables: &mut [crate::meth::MethylationTable],
    haplotypes: &[crate::haplotype::Haplotype],
    pos_1based: u32,
    format_col: &str,
    sample: &str,
    line_no: usize,
) -> Result<(), ReadError> {
    if format_col != "MT:MB" {
        return Err(ReadError::MalformedRecord {
            line: line_no,
            message: format!("standalone record expected FORMAT=MT:MB, got {format_col:?}"),
        });
    }
    let mut parts = sample.splitn(2, ':');
    let top_field = parts.next().unwrap_or("");
    let bot_field = parts.next().unwrap_or("");

    let top_parts: Vec<&str> = top_field.split('|').collect();
    let bot_parts: Vec<&str> = bot_field.split('|').collect();
    if top_parts.len() != tables.len() {
        return Err(ReadError::PloidyEntryCountMismatch {
            pos: pos_1based,
            field: "MT",
            actual: top_parts.len(),
            expected: tables.len(),
        });
    }
    if bot_parts.len() != tables.len() {
        return Err(ReadError::PloidyEntryCountMismatch {
            pos: pos_1based,
            field: "MB",
            actual: bot_parts.len(),
            expected: tables.len(),
        });
    }

    // ref_pos is 0-based; hap_top_c_pos translates it through each haplotype.
    let ref_pos = pos_1based - 1;

    for (hi, &bit_str) in top_parts.iter().enumerate() {
        let hap_top_c_pos =
            haplotypes.get(hi).map_or(ref_pos, |h| h.hap_position_for(ref_pos)) as usize;
        apply_standalone_bit(tables, hi, hap_top_c_pos, bit_str, false, pos_1based)?;
    }
    for (hi, &bit_str) in bot_parts.iter().enumerate() {
        let hap_top_c_pos =
            haplotypes.get(hi).map_or(ref_pos, |h| h.hap_position_for(ref_pos)) as usize;
        apply_standalone_bit(tables, hi, hap_top_c_pos, bit_str, true, pos_1based)?;
    }
    Ok(())
}

/// Process one variant (`ALT != '.'`) methylation record, updating `tables`
/// in place. Looks the record up via `pos_to_vi` (a pre-built
/// `(POS, REF, ALT, GT)` → variant index map; see [`VariantKey`]) and
/// applies the per-haplotype MT/MB bitstrings using `var_hap_coords`. Keying
/// on the full record identity rather than POS alone is what lets multiple
/// VCF rows at the same position resolve to their correct `VariantRecord`.
#[allow(clippy::too_many_arguments)] // tight internal helper; one call site
fn apply_variant_record(
    tables: &mut [crate::meth::MethylationTable],
    pos_1based: u32,
    ref_str: &str,
    alt_str: &str,
    format_col: &str,
    sample: &str,
    pos_to_vi: &std::collections::HashMap<VariantKey, usize>,
    var_hap_coords: &[Vec<Vec<u32>>],
    line_no: usize,
) -> Result<(), ReadError> {
    if format_col != "GT:MT:MB" {
        return Err(ReadError::MalformedRecord {
            line: line_no,
            message: format!("variant record expected FORMAT=GT:MT:MB, got {format_col:?}"),
        });
    }

    // Split GT:MT:MB; capture GT so we can build the composite lookup key.
    let mut parts = sample.splitn(3, ':');
    let gt_str = parts.next().unwrap_or("");
    let top_field = parts.next().unwrap_or("");
    let bot_field = parts.next().unwrap_or("");

    let key = VariantKey {
        pos_1based,
        ref_allele: ref_str.to_string(),
        alt_alleles: alt_str.to_string(),
        gt: gt_str.to_string(),
    };
    let vi = pos_to_vi.get(&key).copied().ok_or_else(|| ReadError::MalformedRecord {
        line: line_no,
        message: format!(
            "variant at POS {pos_1based} (REF={ref_str}, ALT={alt_str}, GT={gt_str}) \
             not found in provided variants slice"
        ),
    })?;

    let top_parts: Vec<&str> = top_field.split('|').collect();
    let bot_parts: Vec<&str> = bot_field.split('|').collect();
    if top_parts.len() != tables.len() {
        return Err(ReadError::PloidyEntryCountMismatch {
            pos: pos_1based,
            field: "MT",
            actual: top_parts.len(),
            expected: tables.len(),
        });
    }
    if bot_parts.len() != tables.len() {
        return Err(ReadError::PloidyEntryCountMismatch {
            pos: pos_1based,
            field: "MB",
            actual: bot_parts.len(),
            expected: tables.len(),
        });
    }

    let hap_coords_for_var =
        if vi < var_hap_coords.len() { var_hap_coords[vi].as_slice() } else { &[] };

    for (hi, &bits) in top_parts.iter().enumerate() {
        let expected: &[u32] =
            if hi < hap_coords_for_var.len() { &hap_coords_for_var[hi] } else { &[] };
        apply_variant_bits(tables, hi, bits, expected, false, pos_1based)?;
    }
    for (hi, &bits) in bot_parts.iter().enumerate() {
        let expected: &[u32] =
            if hi < hap_coords_for_var.len() { &hap_coords_for_var[hi] } else { &[] };
        apply_variant_bits(tables, hi, bits, expected, true, pos_1based)?;
    }
    Ok(())
}

/// Build the `|`-separated per-haplotype top-strand bit string for a
/// standalone CpG whose top-C lives at reference position `ref_pos`.
/// `"1"` for methylated, `"0"` for not.
///
/// Per-haplotype methylation tables are indexed in *haplotype* coordinates
/// (built from the materialized haplotype sequence in
/// [`crate::meth::MethylationTable::from_haplotype`]), so this maps
/// `ref_pos` through each haplotype's `hap_position_for` before querying.
/// Without that mapping, every standalone CpG downstream of an indel would
/// read the wrong bit. When `haplotypes` is empty (no variants → no
/// indels), the function falls back to `ref_pos` directly since
/// `hap_position_for(r) == r` in that case.
fn per_hap_top_bit_string(
    cm: &crate::meth::ContigMethylation,
    haplotypes: &[crate::haplotype::Haplotype],
    ref_pos: u32,
) -> String {
    (0..cm.len())
        .map(|i| {
            let hap_pos = haplotypes.get(i).map_or(ref_pos, |h| h.hap_position_for(ref_pos));
            if cm.table_for(i).is_methylated(hap_pos, false) { "1" } else { "0" }
        })
        .collect::<Vec<_>>()
        .join("|")
}

/// Build the `|`-separated per-haplotype bottom-strand bit string for the
/// CpG whose top-C lives at reference position `ref_pos`. The
/// bottom-strand C lives at the G's position — `ref_pos + 1` in reference
/// coordinates, `hap_pos + 1` in each haplotype's coordinates (a standalone
/// CpG is undisturbed by variants, so the G is adjacent to the C on every
/// haplotype). See [`per_hap_top_bit_string`] for the rationale behind the
/// `hap_position_for` mapping.
fn per_hap_bot_bit_string(
    cm: &crate::meth::ContigMethylation,
    haplotypes: &[crate::haplotype::Haplotype],
    ref_pos: u32,
) -> String {
    (0..cm.len())
        .map(|i| {
            let hap_pos = haplotypes.get(i).map_or(ref_pos, |h| h.hap_position_for(ref_pos));
            if cm.table_for(i).is_methylated(hap_pos + 1, true) { "1" } else { "0" }
        })
        .collect::<Vec<_>>()
        .join("|")
}

/// Open a VCF file as a buffered line reader, transparently decompressing
/// BGZF-/gzip-compressed inputs detected by the gzip magic bytes
/// `0x1f 0x8b` at the start of the file. Returns a boxed `BufRead` so
/// header-only probes and full-body readers can share the open/peek logic.
fn open_vcf_buf_reader(path: &std::path::Path) -> std::io::Result<Box<dyn std::io::BufRead>> {
    use std::io::{BufReader, Read as _};
    let mut peek_buf = [0u8; 2];
    {
        let mut f = std::fs::File::open(path)?;
        f.read_exact(&mut peek_buf)?;
    }
    let file = std::fs::File::open(path)?;
    if peek_buf == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Check whether a VCF actually carries methylation truth: the header must
/// declare both `MT` and `MB` FORMAT fields **and** at least one data record
/// must list `MT` (or `MB`) in its FORMAT column.
///
/// Header declarations alone are not enough — a VCF can declare `MT`/`MB`
/// in the header yet contain zero annotated records (e.g. an aborted
/// `holodeck methylate` run, or a hand-edited header). Such a file has no
/// methylation truth to apply; treating it as methylated lets `simulate`
/// past `validate()` only to hit an internal-invariant failure deep in the
/// per-contig loop. Requiring a real record turns that into a clean,
/// up-front user-facing error instead.
///
/// Streams line-by-line (no full-body buffering) and early-exits: it returns
/// `false` immediately at the column-header line if the header lacked
/// `MT`/`MB`, and returns `true` on the first record whose FORMAT names
/// `MT`/`MB`. Worst case (header declares the fields but no record uses them)
/// scans the whole body, which is unavoidable to prove a negative. Supports
/// plain-text and BGZF-compressed inputs.
///
/// Used by `Simulate::validate` to decide upfront whether methylation
/// chemistry can be applied.
///
/// # Errors
///
/// Returns an I/O error if the file cannot be opened or read.
pub fn vcf_has_mt_mb_records(path: &std::path::Path) -> std::io::Result<bool> {
    use std::io::BufRead as _;
    let mut reader = open_vcf_buf_reader(path)?;
    let mut saw_top_strand_field = false;
    let mut saw_bot_strand_field = false;
    let mut line = String::new();
    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break; // EOF — header declared MT/MB but no record used it.
        }
        let trimmed = line.trim_end_matches(['\n', '\r']);
        if trimmed.starts_with("##") {
            if trimmed.contains("ID=MT,") {
                saw_top_strand_field = true;
            }
            if trimmed.contains("ID=MB,") {
                saw_bot_strand_field = true;
            }
            continue;
        }
        // Column-header (`#CHROM…`) or a data line: if the header never
        // declared both fields, no record can carry methylation truth.
        if !(saw_top_strand_field && saw_bot_strand_field) {
            return Ok(false);
        }
        if trimmed.starts_with('#') {
            continue; // the `#CHROM` line itself carries no FORMAT
        }
        // Data line: FORMAT is the 9th tab-separated column (index 8), a
        // colon-separated list of keys. A real methylation record names
        // `MT`/`MB` there.
        if let Some(format_col) = trimmed.split('\t').nth(8)
            && format_col.split(':').any(|k| k == "MT" || k == "MB")
        {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Pre-parsed methylation VCF body: header MT/MB detection plus per-contig
/// record bytes ready for [`read_contig_methylation`].
///
/// Built once per VCF by [`parse_methylation_vcf`]. Reusing it across the
/// per-contig loop avoids re-decompressing and re-scanning the full VCF on
/// every contig — the cost previously paid by the now-removed
/// `load_contig_methylation_if_present`.
#[derive(Debug, Default)]
pub struct MethylationVcfRecords {
    /// True iff the VCF header declared both `MT` and `MB` FORMAT fields.
    pub has_mt_mb: bool,
    /// Per-contig record lines (header excluded), each `\n`-terminated and
    /// concatenated, ready to pass as the byte slice to
    /// [`read_contig_methylation`]. Contigs absent from the file map to an
    /// empty `Vec`.
    pub per_contig: std::collections::HashMap<String, Vec<u8>>,
}

/// Read and split a methylation-annotated VCF into [`MethylationVcfRecords`]
/// in a single pass: decompresses (BGZF or plain) once, detects whether the
/// header declares MT/MB, and partitions data lines by contig.
///
/// Callers that need per-contig methylation truth (e.g. `simulate`) should
/// call this once before iterating contigs and then call
/// [`load_contig_methylation_from_records`] inside the loop, rather than
/// re-reading the VCF per contig.
///
/// # Errors
///
/// Returns an I/O error if the file cannot be opened or read.
pub fn parse_methylation_vcf(path: &std::path::Path) -> std::io::Result<MethylationVcfRecords> {
    use std::io::BufRead as _;

    // Stream the VCF line-by-line into the per-contig record map. The
    // previous implementation called `read_to_string` first, doubling peak
    // memory: the whole-file `String` and the per-contig `Vec<u8>`s held
    // the body bytes twice. Streaming halves that to roughly one body
    // worth (the `HashMap` values) plus a single `read_line` buffer.
    let mut reader = open_vcf_buf_reader(path)?;

    let mut saw_top_strand_field = false;
    let mut saw_bot_strand_field = false;
    let mut per_contig: std::collections::HashMap<String, Vec<u8>> =
        std::collections::HashMap::new();
    let mut past_header = false;
    let mut line = String::new();
    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break; // EOF
        }
        let trimmed = line.trim_end_matches(['\n', '\r']);
        if trimmed.starts_with("##") {
            // Two independent checks (not `else if`) so a hand-edited header
            // that crams both fields onto one `##FORMAT` line is still
            // recognized. The writer emits them on separate lines, but
            // external tools and manual edits don't always.
            if trimmed.contains("ID=MT,") {
                saw_top_strand_field = true;
            }
            if trimmed.contains("ID=MB,") {
                saw_bot_strand_field = true;
            }
        } else if trimmed.starts_with('#') {
            past_header = true;
            if !(saw_top_strand_field && saw_bot_strand_field) {
                // No MT/MB → no records to collect; return early with the flag.
                return Ok(MethylationVcfRecords::default());
            }
        } else if past_header {
            let Some(chrom) = trimmed.split('\t').next() else {
                log::debug!("skipping malformed VCF line: {trimmed}");
                continue;
            };
            let entry = per_contig.entry(chrom.to_string()).or_default();
            entry.extend_from_slice(trimmed.as_bytes());
            entry.push(b'\n');
        }
    }

    let has_mt_mb = saw_top_strand_field && saw_bot_strand_field;
    Ok(MethylationVcfRecords {
        has_mt_mb,
        per_contig: if has_mt_mb { per_contig } else { std::collections::HashMap::new() },
    })
}

/// Look up the parsed methylation truth for a single contig from the
/// pre-parsed [`MethylationVcfRecords`].
///
/// Returns `Ok(None)` if the source VCF lacked MT/MB FORMAT (signalling that
/// the file carries no methylation truth). Contigs absent from the file
/// produce an empty-but-valid [`crate::meth::ContigMethylation`].
///
/// The MT/MB encoding is the round-trip pair to [`write_contig`]; the
/// reader expects the writer's format exactly.
///
/// `reference` is the contig's full reference sequence (uppercase).
/// `variants` is the sorted variant list for the contig (same contract as
/// [`read_contig_methylation`]). `sample_ploidy` is the whole-VCF ploidy
/// resolved via [`resolve_sample_ploidy`]; threading it down ensures
/// variant-free contigs are sized the same as variant-bearing ones on
/// haploid / triploid / etc. samples.
///
/// # Errors
///
/// Returns an error if [`read_contig_methylation`] rejects a record.
pub fn load_contig_methylation_from_records(
    records: &MethylationVcfRecords,
    contig_name: &str,
    reference: &[u8],
    variants: &[crate::vcf::genotype::VariantRecord],
    sample_ploidy: usize,
) -> anyhow::Result<Option<crate::meth::ContigMethylation>> {
    if !records.has_mt_mb {
        return Ok(None);
    }
    let empty: Vec<u8> = Vec::new();
    let bytes = records.per_contig.get(contig_name).unwrap_or(&empty);
    let cm = read_contig_methylation(bytes, contig_name, reference, variants, sample_ploidy)
        .map_err(|e| anyhow::anyhow!("failed to parse MT/MB for {contig_name}: {e}"))?;
    log::debug!("Loaded methylation truth for {contig_name} from VCF MT/MB");
    Ok(Some(cm))
}

#[cfg(test)]
mod roundtrip_tests {
    use super::*;
    use crate::meth::{ContigMethylation, MethylationTable};

    #[test]
    fn standalone_record_round_trip() {
        // Reference with one CpG at position 1 (0-based). No variants.
        // Hap 0 has the top-strand C at ref pos 1 methylated; hap 1 is
        // fully unmethylated. Neither haplotype has any bottom-strand
        // methylation.
        let reference = b"ACGT";
        let mut h0 = MethylationTable::empty(4);
        h0.set_top(1, true);
        let h1 = MethylationTable::empty(4);
        let cm_in = ContigMethylation::from_tables(vec![h0, h1]);

        let mut buf = Vec::new();
        write_contig(&mut std::io::Cursor::new(&mut buf), "chr1", reference, &[], &cm_in, 2)
            .unwrap();

        let cm_out = read_contig_methylation(&buf, "chr1", reference, &[], 2).unwrap();
        assert_eq!(cm_out.len(), 2);
        // Hap 0: top-strand bit at ref pos 1 must survive the round trip.
        assert!(cm_out.table_for(0).is_methylated(1, false), "hap0 top should be methylated");
        // Hap 1: all bits remain false.
        assert!(!cm_out.table_for(1).is_methylated(1, false), "hap1 top should not be methylated");
        // Bottom-strand bits were not set by the writer, so both should be false.
        assert!(
            !cm_out.table_for(0).is_methylated(2, true),
            "hap0 bottom should not be methylated"
        );
        assert!(
            !cm_out.table_for(1).is_methylated(2, true),
            "hap1 bottom should not be methylated"
        );
    }
}

#[cfg(test)]
mod fuzz_tests {
    //! Closed-loop fuzz round-trip tests for the methylation VCF format.
    //!
    //! Confirms that [`write_contig`] + [`read_contig_methylation`] form a
    //! faithful round-trip pair across random inputs, covering both the
    //! no-variant (standalone-record-only) path and the variant-bearing path.

    use rand::Rng as _;
    use rand::SeedableRng;
    use rand::rngs::SmallRng;

    use super::{read_contig_methylation, write_contig};
    use crate::haplotype::build_haplotypes;
    use crate::meth::{ContigMethylation, MethylationTable};
    use crate::vcf::genotype::{Genotype, VariantRecord};

    /// Assert per-position, per-strand, per-haplotype equality between two
    /// [`ContigMethylation`] values of equal ploidy.
    ///
    /// `hap_lengths[hi]` is the number of positions in haplotype `hi`'s bitmap.
    fn assert_methylation_eq(
        cm_in: &ContigMethylation,
        cm_out: &ContigMethylation,
        hap_lengths: &[usize],
    ) {
        assert_eq!(cm_in.len(), cm_out.len(), "haplotype count mismatch");
        for (hap, &len) in hap_lengths.iter().enumerate().take(cm_in.len()) {
            #[expect(clippy::cast_possible_truncation, reason = "hap length fits u32")]
            for pos in 0..len as u32 {
                assert_eq!(
                    cm_in.table_for(hap).is_methylated(pos, false),
                    cm_out.table_for(hap).is_methylated(pos, false),
                    "top mismatch at hap {hap} pos {pos}",
                );
                assert_eq!(
                    cm_in.table_for(hap).is_methylated(pos, true),
                    cm_out.table_for(hap).is_methylated(pos, true),
                    "bottom mismatch at hap {hap} pos {pos}",
                );
            }
        }
    }

    #[test]
    fn roundtrip_random_bitmap_no_variants() {
        let mut rng = SmallRng::seed_from_u64(42);
        // 10 kb of random ACGT, no variants → every haplotype is the reference.
        let reference: Vec<u8> =
            (0..10_000).map(|_| b"ACGT"[rng.random_range(0..4usize)]).collect();

        // Set random methylation on two diploid haplotypes at every CpG.
        let mut h0 = MethylationTable::with_len(reference.len());
        let mut h1 = MethylationTable::with_len(reference.len());
        for i in 0..reference.len() - 1 {
            if reference[i] == b'C' && reference[i + 1] == b'G' {
                h0.set_top(i, rng.random_bool(0.7));
                h0.set_bottom(i + 1, rng.random_bool(0.7));
                h1.set_top(i, rng.random_bool(0.7));
                h1.set_bottom(i + 1, rng.random_bool(0.7));
            }
        }
        let cm_in = ContigMethylation::from_tables(vec![h0, h1]);

        let mut buf = Vec::new();
        write_contig(&mut std::io::Cursor::new(&mut buf), "chr1", &reference, &[], &cm_in, 2)
            .unwrap();
        let cm_out = read_contig_methylation(&buf, "chr1", &reference, &[], 2).unwrap();

        assert_methylation_eq(&cm_in, &cm_out, &[reference.len(), reference.len()]);
    }

    #[test]
    fn roundtrip_random_bitmap_with_phased_variants() {
        // 500 bp reference with several non-overlapping phased SNPs. For each
        // haplotype, materialise the full sequence, scan for CpGs, and set
        // random top/bottom methylation bits. Write the result to a VCF byte
        // buffer via write_contig, read it back via read_contig_methylation,
        // and assert exact bit-for-bit equality at every haplotype position.
        //
        // The writer and reader both call build_haplotypes(variants, 2,
        // seed_from_u64(0)) internally, so the haplotype layout is identical
        // on both sides.

        // --- Build reference ---
        let reference: Vec<u8> = b"ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT"
            .to_vec();

        // Hand-crafted non-overlapping phased SNPs, all well away from each
        // other and away from the reference ends.
        // GT "1|0" → hap 0 carries ALT, hap 1 carries REF.
        // GT "0|1" → hap 0 carries REF, hap 1 carries ALT.
        // None introduce CpGs that straddle variant boundaries in a way that
        // changes haplotype length (these are all single-base SNPs).
        let variants: Vec<VariantRecord> = vec![
            VariantRecord {
                position: 10,
                ref_allele: b"A".to_vec(),
                alt_alleles: vec![b"T".to_vec()],
                genotype: Genotype::parse("1|0").unwrap(),
            },
            VariantRecord {
                position: 50,
                ref_allele: b"T".to_vec(),
                alt_alleles: vec![b"A".to_vec()],
                genotype: Genotype::parse("0|1").unwrap(),
            },
            VariantRecord {
                position: 100,
                ref_allele: b"G".to_vec(),
                alt_alleles: vec![b"C".to_vec()],
                genotype: Genotype::parse("1|0").unwrap(),
            },
            VariantRecord {
                position: 200,
                ref_allele: b"C".to_vec(),
                alt_alleles: vec![b"A".to_vec()],
                genotype: Genotype::parse("0|1").unwrap(),
            },
            VariantRecord {
                position: 350,
                ref_allele: b"A".to_vec(),
                alt_alleles: vec![b"G".to_vec()],
                genotype: Genotype::parse("1|0").unwrap(),
            },
        ];

        // Materialise haplotypes with the same seed the writer/reader use.
        let haplotypes = build_haplotypes(&variants, 2, &mut SmallRng::seed_from_u64(0));

        // Determine per-haplotype bitmap lengths.
        #[expect(clippy::cast_possible_truncation, reason = "reference length fits u32")]
        let ref_len_u32 = reference.len() as u32;
        let hap_lengths: Vec<usize> =
            haplotypes.iter().map(|h| h.hap_position_for(ref_len_u32) as usize).collect();

        // Build per-haplotype methylation tables with random CpG bits.
        let mut rng = SmallRng::seed_from_u64(7);
        let tables: Vec<MethylationTable> = haplotypes
            .iter()
            .zip(hap_lengths.iter())
            .map(|(hap, &len)| {
                let cap = len;
                let (hap_bases, _ref_positions, _hap_start) =
                    hap.extract_fragment(&reference, 0, cap);
                let mut table = MethylationTable::with_len(len);
                for i in 0..hap_bases.len().saturating_sub(1) {
                    let c0 = hap_bases[i].to_ascii_uppercase();
                    let c1 = hap_bases[i + 1].to_ascii_uppercase();
                    if c0 == b'C' && c1 == b'G' {
                        table.set_top(i, rng.random_bool(0.6));
                        table.set_bottom(i + 1, rng.random_bool(0.6));
                    }
                }
                table
            })
            .collect();
        let cm_in = ContigMethylation::from_tables(tables);

        // Round-trip through write_contig → read_contig_methylation.
        let mut buf = Vec::new();
        write_contig(&mut std::io::Cursor::new(&mut buf), "chr1", &reference, &variants, &cm_in, 2)
            .unwrap();
        let cm_out = read_contig_methylation(&buf, "chr1", &reference, &variants, 2).unwrap();

        assert_methylation_eq(&cm_in, &cm_out, &hap_lengths);
    }

    #[test]
    fn roundtrip_triploid_with_phased_variants() {
        // Ploidy 3 exercises the per-haplotype index path beyond the diploid
        // cases above: three distinct bitmaps, phased GTs that put ALT on
        // different subsets of the three haplotypes, and an insertion that
        // shifts downstream coordinates on only the haplotypes carrying it.
        let reference: Vec<u8> = b"ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT\
              ACGTCGATCGATCGCGATCGACGT"
            .to_vec();

        // "1|0|1" → haps 0 and 2 carry ALT; "0|1|0" → only hap 1; the
        // insertion "0|0|1" shifts downstream coordinates on hap 2 alone.
        let variants: Vec<VariantRecord> = vec![
            VariantRecord {
                position: 12,
                ref_allele: b"A".to_vec(),
                alt_alleles: vec![b"T".to_vec()],
                genotype: Genotype::parse("1|0|1").unwrap(),
            },
            VariantRecord {
                position: 40,
                ref_allele: b"T".to_vec(),
                alt_alleles: vec![b"A".to_vec()],
                genotype: Genotype::parse("0|1|0").unwrap(),
            },
            VariantRecord {
                position: 70,
                ref_allele: b"A".to_vec(),
                alt_alleles: vec![b"AGGG".to_vec()],
                genotype: Genotype::parse("0|0|1").unwrap(),
            },
        ];

        let haplotypes = build_haplotypes(&variants, 3, &mut SmallRng::seed_from_u64(0));
        #[expect(clippy::cast_possible_truncation, reason = "reference length fits u32")]
        let ref_len_u32 = reference.len() as u32;
        let hap_lengths: Vec<usize> =
            haplotypes.iter().map(|h| h.hap_position_for(ref_len_u32) as usize).collect();

        let mut rng = SmallRng::seed_from_u64(11);
        let tables: Vec<MethylationTable> = haplotypes
            .iter()
            .zip(hap_lengths.iter())
            .map(|(hap, &len)| {
                let (hap_bases, _ref_positions, _hap_start) =
                    hap.extract_fragment(&reference, 0, len);
                let mut table = MethylationTable::with_len(len);
                for i in 0..hap_bases.len().saturating_sub(1) {
                    if hap_bases[i].eq_ignore_ascii_case(&b'C')
                        && hap_bases[i + 1].eq_ignore_ascii_case(&b'G')
                    {
                        table.set_top(i, rng.random_bool(0.6));
                        table.set_bottom(i + 1, rng.random_bool(0.6));
                    }
                }
                table
            })
            .collect();
        let cm_in = ContigMethylation::from_tables(tables);
        assert_eq!(cm_in.len(), 3, "fixture must be triploid");

        let mut buf = Vec::new();
        write_contig(&mut std::io::Cursor::new(&mut buf), "chr1", &reference, &variants, &cm_in, 3)
            .unwrap();
        let cm_out = read_contig_methylation(&buf, "chr1", &reference, &variants, 3).unwrap();

        assert_methylation_eq(&cm_in, &cm_out, &hap_lengths);
    }

    #[test]
    fn standalone_cpg_downstream_of_indel_round_trips_per_haplotype() {
        // Regression: an earlier writer/reader treated `ref_pos` as a
        // haplotype coordinate when emitting/parsing standalone CpG records.
        // On a haplotype with an upstream insertion, every downstream
        // standalone CpG's bit was written and read at the wrong bitmap
        // index. This test pins the per-haplotype mapping in place.
        //
        // Reference layout (0-based):
        //   pos  0    5    10   15   20   25
        //   ref  ATATA T    CGATA CGATA CGATA CG  (CpGs at 6, 11, 16, 21, 24)
        //
        // Variant: 3 bp insertion at position 3 on hap 0 only ("0|1" → hap 1
        // carries ALT, hap 0 carries REF). On hap 1, every downstream
        // standalone CpG's haplotype coordinate is shifted by +3.
        let reference: Vec<u8> = b"ATATATCGATACGATACGATACGCG".to_vec();
        let variants: Vec<VariantRecord> = vec![VariantRecord {
            position: 3,
            ref_allele: b"T".to_vec(),
            alt_alleles: vec![b"TGGG".to_vec()],
            genotype: Genotype::parse("0|1").unwrap(),
        }];

        // Build haplotypes with the same seed the writer/reader use.
        let haplotypes = build_haplotypes(&variants, 2, &mut SmallRng::seed_from_u64(0));
        #[expect(clippy::cast_possible_truncation, reason = "reference length fits u32")]
        let ref_len_u32 = reference.len() as u32;
        let hap_lengths: Vec<usize> =
            haplotypes.iter().map(|h| h.hap_position_for(ref_len_u32) as usize).collect();

        // Set distinctive per-haplotype methylation: hap 0 has top-strand
        // bits set at every CpG-context C in its (insertion-free) sequence;
        // hap 1 has bottom-strand bits set at every G of every CpG. Both
        // patterns rely on per-haplotype coordinates being correct.
        let mut tables: Vec<MethylationTable> = haplotypes
            .iter()
            .zip(hap_lengths.iter())
            .map(|(hap, &len)| {
                let (hap_bases, _ref_positions, _hap_start) =
                    hap.extract_fragment(&reference, 0, len);
                let mut table = MethylationTable::with_len(len);
                for i in 0..hap_bases.len().saturating_sub(1) {
                    if hap_bases[i].eq_ignore_ascii_case(&b'C')
                        && hap_bases[i + 1].eq_ignore_ascii_case(&b'G')
                    {
                        table.set_top(i, true);
                        table.set_bottom(i + 1, true);
                    }
                }
                table
            })
            .collect();
        // Pin a hap-specific asymmetry so the round-trip would catch any
        // off-by-three between writer and reader: clear hap 1's first CpG's
        // top-strand bit.
        let hap1_first_cpg_hap_pos = haplotypes[1].hap_position_for(6) as usize;
        tables[1].set_top(hap1_first_cpg_hap_pos, false);
        let cm_in = ContigMethylation::from_tables(tables);

        let mut buf = Vec::new();
        write_contig(&mut std::io::Cursor::new(&mut buf), "chr1", &reference, &variants, &cm_in, 2)
            .unwrap();
        let cm_out = read_contig_methylation(&buf, "chr1", &reference, &variants, 2).unwrap();

        assert_methylation_eq(&cm_in, &cm_out, &hap_lengths);

        // Independent positive check: the wrong-path code would write hap
        // 1's first standalone CpG bit at hap-coord 6 (treating ref_pos as
        // hap_pos). Confirm the bit actually landed at the shifted
        // hap-coordinate, not at the reference position.
        let hap1_table = cm_out.table_for(1);
        let cpg_ref_pos: u32 = 11; // second reference CpG
        let cpg_hap_pos = haplotypes[1].hap_position_for(cpg_ref_pos);
        assert_ne!(cpg_hap_pos, cpg_ref_pos, "test setup assumes indel shifts hap 1");
        assert!(
            hap1_table.is_methylated(cpg_hap_pos, false),
            "hap 1 CpG at hap-coord {cpg_hap_pos} (ref pos {cpg_ref_pos}) should be methylated",
        );
        assert!(
            !hap1_table.is_methylated(cpg_ref_pos, false),
            "hap 1 should NOT have a top-strand bit at ref pos {cpg_ref_pos} \
             (that would indicate ref_pos was used as a hap coord)",
        );
    }
}

#[cfg(test)]
mod writer_tests {
    use super::*;
    use crate::meth::{ContigMethylation, MethylationTable};
    use std::io::Cursor;

    #[test]
    fn writes_methylation_only_record_for_reference_cpg() {
        // Reference with one CpG at position 1 (0-based), no variants.
        // Hap 0 is fully methylated on the top strand; hap 1 is unmethylated
        // everywhere. Both haplotypes are unmethylated on the bottom strand.
        let reference = b"ACGT";
        let mut h0 = MethylationTable::empty(4);
        h0.set_top(1, true);
        let h1 = MethylationTable::empty(4);
        let cm = ContigMethylation::from_tables(vec![h0, h1]);
        let mut buf = Vec::new();
        write_contig(&mut Cursor::new(&mut buf), "chr1", reference, &[], &cm, 2).unwrap();
        let s = String::from_utf8(buf).unwrap();
        // VCF POS is 1-based; CpG top-C is at ref pos 1 (0-based) → POS 2.
        assert!(s.contains("chr1\t2\t.\tC\t.\t.\t.\t.\tMT:MB\t1|0:0|0"), "got: {s}");
    }

    #[test]
    fn writes_variant_record_with_mt_mb_for_alt_cpg() {
        use crate::vcf::genotype::Genotype;
        // ref: AAATTAA; SNPs at 3 (T->C) and 4 (T->G) on hap0 -> CpG in alt on hap0.
        let reference = b"AAATTAA";
        let variants = vec![
            VariantRecord {
                position: 3,
                ref_allele: b"T".to_vec(),
                alt_alleles: vec![b"C".to_vec()],
                genotype: Genotype::parse("1|0").unwrap(),
            },
            VariantRecord {
                position: 4,
                ref_allele: b"T".to_vec(),
                alt_alleles: vec![b"G".to_vec()],
                genotype: Genotype::parse("1|0").unwrap(),
            },
        ];
        // Build a per-haplotype methylation table for hap0 with the alt-CpG
        // methylated on both strands; hap1 has nothing.
        //
        // Haplotype 0 (hap_allele_index=0) carries allele Some(1) at both
        // variants (GT="1|0", so allele_index 0 maps to allele 1). When
        // build_haplotypes runs, hap0 carries T->C at position 3 and T->G at
        // position 4, materializing as AAACGAA (CG at hap-coords 3,4).
        // So top-strand C is at hap-coord 3, bottom-strand C is at hap-coord 4.
        let mut h0 = MethylationTable::empty(7);
        h0.set_top(3, true);
        h0.set_bottom(4, true);
        let h1 = MethylationTable::empty(7);
        let cm = ContigMethylation::from_tables(vec![h0, h1]);

        let mut buf = Vec::new();
        write_contig(&mut Cursor::new(&mut buf), "chr1", reference, &variants, &cm, 2).unwrap();
        let s = String::from_utf8(buf).unwrap();
        // Variant 0 (position 3, 1-based POS=4) owns the CpG via upstream-wins.
        // GT=1|0; hap0 carries the alt 'C', which is the top-C of the CpG.
        // MT for hap0 = top-strand bit at hap-coord 3 = 1. MB for hap0 = bottom
        // bit at hap-coord 4 = 1. Hap1 carries REF → MT/MB = '.'.
        assert!(
            s.contains("chr1\t4\t.\tT\tC\t.\t.\t.\tGT:MT:MB\t1|0:1|.:1|."),
            "expected variant 0 row, got:\n{s}"
        );
        // Variant 1 (position 4, 1-based POS=5) owns no CpGs (variant 0 won).
        // MT/MB = '.|.' for both haplotypes.
        assert!(
            s.contains("chr1\t5\t.\tT\tG\t.\t.\t.\tGT:MT:MB\t1|0:.|.:.|."),
            "expected variant 1 row, got:\n{s}"
        );
    }
}

#[cfg(test)]
mod reader_error_tests {
    use super::*;
    use crate::vcf::genotype::{Genotype, VariantRecord};

    /// Invalid UTF-8 in the VCF body must surface as `MalformedRecord`, not
    /// be silently coerced to an empty body (which would return an
    /// all-zeros methylation table — indistinguishable from a legitimately
    /// unmethylated contig).
    #[test]
    fn invalid_utf8_body_is_rejected() {
        let buf: &[u8] = &[0xFFu8, 0xFE, b'\n'];
        let err = read_contig_methylation(buf, "chr1", b"ACGT", &[], 2).unwrap_err();
        let msg = format!("{err}");
        assert!(matches!(err, ReadError::MalformedRecord { .. }), "got {err:?}");
        assert!(msg.contains("not valid UTF-8"), "unexpected message: {msg}");
    }

    /// VCF POS is 1-based, so `0` is invalid and must surface as a clean
    /// `MalformedRecord` rather than a `u32` underflow inside the standalone
    /// reader (which does `pos_1based - 1`).
    #[test]
    fn standalone_record_with_pos_zero_is_rejected() {
        let buf = b"chr1\t0\t.\tC\t.\t.\t.\t.\tMT:MB\t0|0:0|0\n";
        let err = read_contig_methylation(buf, "chr1", b"ACGT", &[], 2).unwrap_err();
        let msg = format!("{err}");
        assert!(matches!(err, ReadError::MalformedRecord { .. }), "got {err:?}");
        assert!(msg.contains("POS must be 1-based"), "unexpected message: {msg}");
    }

    #[test]
    fn malformed_record_with_too_few_columns_is_rejected() {
        // Six tab-separated fields, fewer than the required ten.
        let buf = b"chr1\t2\t.\tC\t.\t.\n";
        let err = read_contig_methylation(buf, "chr1", b"ACGT", &[], 2).unwrap_err();
        assert!(
            matches!(err, ReadError::MalformedRecord { .. }),
            "expected MalformedRecord, got {err:?}"
        );
    }

    #[test]
    fn invalid_mtmb_character_in_standalone_record_is_rejected() {
        // MT contains '2' (not a valid 0/1/.).
        let buf = b"chr1\t2\t.\tC\t.\t.\t.\t.\tMT:MB\t2|0:0|0\n";
        let err = read_contig_methylation(buf, "chr1", b"ACGT", &[], 2).unwrap_err();
        assert!(
            matches!(err, ReadError::InvalidMtMbChar { ch: '2', .. }),
            "expected InvalidMtMbChar for '2', got {err:?}"
        );
    }

    #[test]
    fn standalone_record_with_wrong_mt_entry_count_is_rejected() {
        // Ploidy=2 but MT carries 3 |-separated entries — extra entries
        // should not be silently truncated.
        let buf = b"chr1\t2\t.\tC\t.\t.\t.\t.\tMT:MB\t1|0|0:0|0\n";
        let err = read_contig_methylation(buf, "chr1", b"ACGT", &[], 2).unwrap_err();
        assert!(
            matches!(err, ReadError::PloidyEntryCountMismatch { actual: 3, expected: 2, .. }),
            "expected PloidyEntryCountMismatch{{actual:3,expected:2}}, got {err:?}"
        );
    }

    #[test]
    fn standalone_record_with_too_few_mb_entries_is_rejected() {
        // Ploidy=2 but MB carries only 1 |-separated entry — missing
        // trailing entries should not be silently left untouched.
        let buf = b"chr1\t2\t.\tC\t.\t.\t.\t.\tMT:MB\t1|0:0\n";
        let err = read_contig_methylation(buf, "chr1", b"ACGT", &[], 2).unwrap_err();
        assert!(
            matches!(err, ReadError::PloidyEntryCountMismatch { actual: 1, expected: 2, .. }),
            "expected PloidyEntryCountMismatch{{actual:1,expected:2}}, got {err:?}"
        );
    }

    #[test]
    fn variant_record_with_wrong_mt_entry_count_is_rejected() {
        use crate::vcf::genotype::{Genotype, VariantRecord};
        // Ploidy=2 but MT carries 3 |-separated entries on a variant row.
        let variants = vec![VariantRecord {
            position: 1,
            ref_allele: b"T".to_vec(),
            alt_alleles: vec![b"C".to_vec()],
            genotype: Genotype::parse("1|0").unwrap(),
        }];
        let buf = b"chr1\t2\t.\tT\tC\t.\t.\t.\tGT:MT:MB\t1|0:.|.|.:.|.\n";
        let err = read_contig_methylation(buf, "chr1", b"ATTT", &variants, 2).unwrap_err();
        assert!(
            matches!(err, ReadError::PloidyEntryCountMismatch { actual: 3, expected: 2, .. }),
            "expected PloidyEntryCountMismatch{{actual:3,expected:2}}, got {err:?}"
        );
    }

    #[test]
    fn variant_record_with_too_few_mb_entries_is_rejected() {
        use crate::vcf::genotype::{Genotype, VariantRecord};
        let variants = vec![VariantRecord {
            position: 1,
            ref_allele: b"T".to_vec(),
            alt_alleles: vec![b"C".to_vec()],
            genotype: Genotype::parse("1|0").unwrap(),
        }];
        let buf = b"chr1\t2\t.\tT\tC\t.\t.\t.\tGT:MT:MB\t1|0:.|.:.\n";
        let err = read_contig_methylation(buf, "chr1", b"ATTT", &variants, 2).unwrap_err();
        assert!(
            matches!(err, ReadError::PloidyEntryCountMismatch { actual: 1, expected: 2, .. }),
            "expected PloidyEntryCountMismatch{{actual:1,expected:2}}, got {err:?}"
        );
    }

    #[test]
    fn variant_record_mtmb_length_mismatch_is_rejected() {
        // ref ATTT, variant T->C at pos 1 (GT=1|0). Hap0's alt 'C' contains
        // zero CpGs, so the expected MT bitstring length on hap0 is 0; a
        // length-1 string ('1') here is a mismatch.
        let variants = vec![VariantRecord {
            position: 1,
            ref_allele: b"T".to_vec(),
            alt_alleles: vec![b"C".to_vec()],
            genotype: Genotype::parse("1|0").unwrap(),
        }];
        let buf = b"chr1\t2\t.\tT\tC\t.\t.\t.\tGT:MT:MB\t1|0:1|.:.|.\n";
        let err = read_contig_methylation(buf, "chr1", b"ATTT", &variants, 2).unwrap_err();
        assert!(
            matches!(err, ReadError::MtMbLengthMismatch { actual: 1, expected: 0, .. }),
            "expected MtMbLengthMismatch{{actual:1,expected:0}}, got {err:?}"
        );
    }

    /// `"."` is a placeholder for "this haplotype owns zero CpGs at this
    /// variant". When the variants slice says the haplotype actually owns
    /// CpGs (e.g. the `1|0` hap on an `A → ACG` insertion), the parser must
    /// reject a `.` MT/MB entry rather than silently drop the truth bits.
    #[test]
    fn variant_record_dot_when_owned_cpgs_expected_is_rejected() {
        // ref AAAAAAAA + heterozygous `A → ACG` insertion at pos 1 on hap 0.
        // Hap 0's ALT allele inserts a `CG` dinucleotide → 1 owned CpG.
        // The malformed row uses `.` for hap 0 instead of `0` or `1`.
        let variants = vec![VariantRecord {
            position: 1,
            ref_allele: b"A".to_vec(),
            alt_alleles: vec![b"ACG".to_vec()],
            genotype: Genotype::parse("1|0").unwrap(),
        }];
        let buf = b"chr1\t2\t.\tA\tACG\t.\t.\t.\tGT:MT:MB\t1|0:.|.:.|.\n";
        let err = read_contig_methylation(buf, "chr1", b"AAAAAAAA", &variants, 2).unwrap_err();
        assert!(
            matches!(err, ReadError::MtMbLengthMismatch { actual: 0, expected: 1, .. }),
            "expected MtMbLengthMismatch{{actual:0,expected:1}}, got {err:?}"
        );
    }

    /// Two variants share the same POS but live on different haplotypes
    /// (a decomposed multi-allelic site, GT `1|0` and `0|1`). Each one's
    /// alt allele introduces its own CpG on its own haplotype. The
    /// writer/reader index must disambiguate these by full `(POS, REF, ALT,
    /// GT)` identity, not by POS alone.
    ///
    /// Regression: keying `pos_to_vi` by POS alone made the reader resolve
    /// every row at that POS to the last-inserted variant index, then
    /// reject one of the two rows with `MtMbLengthMismatch` because the
    /// expected per-hap CpG count belonged to the other variant. The
    /// round-trip therefore failed altogether, and even when lengths
    /// happened to align, bits would have landed on the wrong haplotype.
    #[test]
    fn round_trip_disambiguates_two_variants_sharing_pos() {
        use crate::meth::{ContigMethylation, MethylationTable};
        use crate::vcf::genotype::Genotype;
        use std::io::Cursor;
        // Reference AAAAAAAA — no CpGs. Two phased insertions at pos 1 that
        // each plant a CpG on a different haplotype.
        let reference = b"AAAAAAAA".to_vec();
        let variants = vec![
            // Variant A: hap 0 gets the insertion `A → ACG` at pos 1.
            VariantRecord {
                position: 1,
                ref_allele: b"A".to_vec(),
                alt_alleles: vec![b"ACG".to_vec()],
                genotype: Genotype::parse("1|0").unwrap(),
            },
            // Variant B: hap 1 gets the insertion `A → ACG` at the SAME pos 1.
            VariantRecord {
                position: 1,
                ref_allele: b"A".to_vec(),
                alt_alleles: vec![b"ACG".to_vec()],
                genotype: Genotype::parse("0|1").unwrap(),
            },
        ];

        // Per-hap methylation tables: hap 0 has its CpG (from variant A)
        // methylated on both strands at hap-coord 2 / 3 (alt span starts at
        // ref pos 1 → hap-coord 1, the inserted `CG` sits at hap-coords 2
        // and 3). Hap 1 has its CpG (from variant B) UNmethylated. The
        // asymmetry is what catches a wrong-variant route on read.
        let mut h0 = MethylationTable::with_len(reference.len() + 2);
        h0.set_top(2, true);
        h0.set_bottom(3, true);
        let h1 = MethylationTable::with_len(reference.len() + 2);
        let cm_in = ContigMethylation::from_tables(vec![h0, h1]);

        let mut buf = Vec::new();
        write_contig(&mut Cursor::new(&mut buf), "chr1", &reference, &variants, &cm_in, 2).unwrap();
        let cm_out = read_contig_methylation(&buf, "chr1", &reference, &variants, 2)
            .expect("round-trip must succeed when records share POS but differ in GT");

        // Hap 0: top + bottom bits of the inserted CpG must survive.
        assert!(cm_out.table_for(0).is_methylated(2, false), "hap0 top should be set");
        assert!(cm_out.table_for(0).is_methylated(3, true), "hap0 bottom should be set");
        // Hap 1: nothing methylated.
        assert!(!cm_out.table_for(1).is_methylated(2, false), "hap1 top must stay unset");
        assert!(!cm_out.table_for(1).is_methylated(3, true), "hap1 bottom must stay unset");
    }

    /// Haploid samples must produce single-entry MT/MB strings on
    /// variant-free contigs too. Before sample-level ploidy resolution, a
    /// per-contig `unwrap_or(2)` fallback would default empty contigs to
    /// diploid and emit `0|0` instead of `0`, breaking round-trip with any
    /// contig that does have a haploid variant.
    #[test]
    fn haploid_variant_free_contig_writes_single_entry_mt_mb() {
        use crate::meth::{ContigMethylation, MethylationTable};
        use std::io::Cursor;
        // No variants → loader can't infer haploidy on its own; we pass
        // `sample_ploidy = 1` explicitly the way `run_simulation` /
        // `methylate::execute` would after the global resolution step.
        let reference = b"ACGT"; // single CpG at pos 1
        let cm = ContigMethylation::from_tables(vec![MethylationTable::with_len(4)]);
        let mut buf = Vec::new();
        write_contig(&mut Cursor::new(&mut buf), "chr1", reference, &[], &cm, 1).unwrap();
        let s = String::from_utf8(buf).unwrap();
        // One standalone row at POS 2 with a single-entry MT and MB field.
        assert!(s.contains("chr1\t2\t.\tC\t.\t.\t.\t.\tMT:MB\t0:0\n"), "got:\n{s}");

        // Round-trip with sample_ploidy=1.
        let cm_out = read_contig_methylation(s.as_bytes(), "chr1", reference, &[], 1).unwrap();
        assert_eq!(cm_out.len(), 1, "haploid round-trip must yield 1 table");
    }
}

#[cfg(test)]
mod record_probe_tests {
    use super::*;
    use std::io::Write as _;

    fn write_temp_vcf(content: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f.flush().unwrap();
        f
    }

    /// Header declares MT/MB AND a real record names them in FORMAT → true.
    #[test]
    fn header_and_record_with_mt_mb_returns_true() {
        let vcf = "##fileformat=VCFv4.4\n\
                   ##FORMAT=<ID=MT,Number=.,Type=String,Description=\"top\">\n\
                   ##FORMAT=<ID=MB,Number=.,Type=String,Description=\"bottom\">\n\
                   #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
                   chr1\t2\t.\tC\t.\t.\t.\t.\tMT:MB\t0:0\n";
        let f = write_temp_vcf(vcf);
        assert!(vcf_has_mt_mb_records(f.path()).unwrap());
    }

    /// Variant-style FORMAT (`GT:MT:MB`) also counts as a real record.
    #[test]
    fn variant_record_with_mt_mb_returns_true() {
        let vcf = "##fileformat=VCFv4.4\n\
                   ##FORMAT=<ID=MT,Number=.,Type=String,Description=\"top\">\n\
                   ##FORMAT=<ID=MB,Number=.,Type=String,Description=\"bottom\">\n\
                   #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
                   chr1\t2\t.\tT\tC\t.\t.\t.\tGT:MT:MB\t1|0:1|.:1|.\n";
        let f = write_temp_vcf(vcf);
        assert!(vcf_has_mt_mb_records(f.path()).unwrap());
    }

    /// Regression: header declares MT/MB but the file has NO records. This
    /// must be `false` so `simulate --methylation-mode` rejects it up front
    /// instead of hitting an internal-invariant failure later.
    #[test]
    fn header_only_no_records_returns_false() {
        let vcf = "##fileformat=VCFv4.4\n\
                   ##FORMAT=<ID=MT,Number=.,Type=String,Description=\"top\">\n\
                   ##FORMAT=<ID=MB,Number=.,Type=String,Description=\"bottom\">\n\
                   #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
        let f = write_temp_vcf(vcf);
        assert!(!vcf_has_mt_mb_records(f.path()).unwrap());
    }

    /// Header declares MT/MB but the only record carries a plain `GT` FORMAT
    /// (no methylation truth) → false.
    #[test]
    fn record_without_mt_mb_format_returns_false() {
        let vcf = "##fileformat=VCFv4.4\n\
                   ##FORMAT=<ID=MT,Number=.,Type=String,Description=\"top\">\n\
                   ##FORMAT=<ID=MB,Number=.,Type=String,Description=\"bottom\">\n\
                   #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
                   chr1\t1\t.\tA\tT\t.\t.\t.\tGT\t0|1\n";
        let f = write_temp_vcf(vcf);
        assert!(!vcf_has_mt_mb_records(f.path()).unwrap());
    }

    #[test]
    fn header_with_only_mt_returns_false() {
        let vcf = "##fileformat=VCFv4.4\n\
                   ##FORMAT=<ID=MT,Number=.,Type=String,Description=\"top\">\n\
                   #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
                   chr1\t2\t.\tC\t.\t.\t.\t.\tMT\t0\n";
        let f = write_temp_vcf(vcf);
        assert!(!vcf_has_mt_mb_records(f.path()).unwrap());
    }

    #[test]
    fn header_with_only_mb_returns_false() {
        let vcf = "##fileformat=VCFv4.4\n\
                   ##FORMAT=<ID=MB,Number=.,Type=String,Description=\"bottom\">\n\
                   #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
                   chr1\t2\t.\tC\t.\t.\t.\t.\tMB\t0\n";
        let f = write_temp_vcf(vcf);
        assert!(!vcf_has_mt_mb_records(f.path()).unwrap());
    }

    #[test]
    fn header_with_neither_returns_false() {
        let vcf = "##fileformat=VCFv4.4\n\
                   #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
        let f = write_temp_vcf(vcf);
        assert!(!vcf_has_mt_mb_records(f.path()).unwrap());
    }

    /// When the header lacks MT/MB, the probe must early-exit at the
    /// `#CHROM` line and never decode the data body. Pin that by appending
    /// invalid UTF-8 after the column header: a body-reading implementation
    /// would fail to decode and surface an `InvalidData` I/O error.
    #[test]
    fn no_mt_mb_header_does_not_read_body() {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(
            b"##fileformat=VCFv4.4\n\
              #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
        )
        .unwrap();
        // Invalid UTF-8 (lone 0xFF bytes) where the data body would be.
        f.write_all(&[0xFFu8; 64]).unwrap();
        f.flush().unwrap();
        // Must return false without erroring: the header lacked MT/MB, so we
        // bail at the `#CHROM` line before touching the 0xFF block.
        assert!(!vcf_has_mt_mb_records(f.path()).unwrap());
    }
}
