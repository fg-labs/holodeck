use coitrees::{COITree, Interval, IntervalTree};

use crate::vcf::genotype::VariantRecord;

/// A single variant assigned to a specific haplotype.
#[derive(Debug, Clone)]
pub struct HaplotypeVariant {
    /// 0-based reference position where the variant starts.
    pub ref_pos: u32,
    /// Length of the reference allele in bases. Uses `u32` to support
    /// structural variants with large reference spans.
    pub ref_len: u32,
    /// Alternate allele bases.
    pub alt_bases: Vec<u8>,
}

/// Sparse representation of one haplotype — a reference overlay of variants.
///
/// Instead of materializing a full haplotype sequence (which would require
/// ~250MB per haplotype for human chr1), this stores only the differences from
/// the reference as a sorted set of variants in a [`COITree`] for efficient
/// range queries.
///
/// Fragment extraction works by walking the reference sequence and
/// substituting alt alleles at variant positions on the fly.
///
/// Because coitrees requires `Copy + Default` for metadata, we store variant
/// indices (as `u32`) in the tree and keep the actual variant data in a
/// separate `Vec`.
pub struct Haplotype {
    /// 0-based haplotype allele index (e.g., 0 or 1 for diploid).
    allele_index: usize,
    /// Variant data, indexed by position in this vec.
    variant_data: Vec<HaplotypeVariant>,
    /// Interval tree mapping genomic ranges to indices in `variant_data`.
    variant_tree: COITree<u32, u32>,
    /// Cumulative-delta index used by [`Self::hap_position_for`]. Each entry
    /// is `(var_end, cumulative_delta)`, sorted by `var_end` ascending.
    /// `cumulative_delta` at index `i` is the sum of `(alt_len - ref_len)`
    /// across all variants whose `var_end` is `<= entries[i].var_end`.
    /// Enables `O(log n)` lookup instead of an `O(n)` scan in the
    /// per-fragment hot path.
    end_prefix_deltas: Vec<(u32, i64)>,
}

impl Haplotype {
    /// Return the allele index of this haplotype.
    #[must_use]
    pub fn allele_index(&self) -> usize {
        self.allele_index
    }

    /// Extract a fragment from this haplotype at the given reference
    /// coordinates.
    ///
    /// Walks the reference from `ref_start` and produces `fragment_len` bases,
    /// substituting alternate alleles where this haplotype has variants.
    /// Returns the fragment bases, a list of reference positions
    /// corresponding to each fragment base (for golden BAM coordinate
    /// mapping), and the haplotype-coordinate position the fragment starts
    /// at (for per-haplotype methylation lookups).
    ///
    /// # Arguments
    /// * `reference` — Full reference sequence for this contig.
    /// * `ref_start` — 0-based start position on the reference.
    /// * `fragment_len` — Desired number of output bases.
    ///
    /// # Returns
    /// A tuple of `(fragment_bases, ref_positions, hap_start)`:
    /// * `fragment_bases[i]` is the base at fragment-internal position `i`.
    /// * `ref_positions[i]` is the reference position corresponding to
    ///   `fragment_bases[i]`. For inserted bases, the reference position is
    ///   that of the base preceding the insertion.
    /// * `hap_start` is the first haplotype-coordinate position covered by
    ///   the fragment. Equals `ref_start` plus net upstream insertion length
    ///   (insertions add bases, deletions remove bases). With no upstream
    ///   indels it equals `ref_start`.
    #[must_use]
    #[expect(clippy::cast_possible_wrap, reason = "genomic coords < i32::MAX")]
    pub fn extract_fragment(
        &self,
        reference: &[u8],
        ref_start: u32,
        fragment_len: usize,
    ) -> (Vec<u8>, Vec<u32>, u32) {
        let mut bases = Vec::with_capacity(fragment_len);
        let mut ref_positions = Vec::with_capacity(fragment_len);

        // Collect overlapping variant indices. We query exactly the range we
        // need: [ref_start, ref_start + fragment_len). Deletions whose ref
        // allele extends past our window are still caught because the tree
        // stores them with their full ref allele span.
        let query_end = (ref_start as usize + fragment_len).min(reference.len());
        #[expect(
            clippy::cast_possible_truncation,
            reason = "query_end bounded by reference.len() which fits i32"
        )]
        let query_end_i32 = query_end.saturating_sub(1) as i32;
        let mut overlapping_indices: Vec<u32> = Vec::new();
        self.variant_tree.query(ref_start as i32, query_end_i32, |node| {
            // clone() rather than *deref because coitrees' query callback
            // yields &T on NEON (ARM) but T directly on nosimd (x86).
            #[allow(clippy::clone_on_copy)]
            overlapping_indices.push(node.metadata.clone());
        });

        // Sort variants by position for sequential processing.
        overlapping_indices.sort_unstable_by_key(|&idx| self.variant_data[idx as usize].ref_pos);

        let mut ref_pos = ref_start as usize;
        let mut var_idx = 0;

        // Advance past variants entirely before our window, and handle
        // variants that start before ref_start but span into it (e.g., a
        // deletion that started upstream).
        while var_idx < overlapping_indices.len() {
            let var = &self.variant_data[overlapping_indices[var_idx] as usize];
            let var_end = var.ref_pos as usize + var.ref_len as usize;
            if var_end <= ref_pos {
                // Variant is entirely before our window — skip it.
                var_idx += 1;
            } else if (var.ref_pos as usize) < ref_pos {
                // Variant starts before our window but spans into it.
                // For a deletion, the ref bases are consumed — skip past them.
                // For an insertion at this position, the alt bases were already
                // partially consumed upstream, so we skip the variant entirely.
                ref_pos = var_end;
                var_idx += 1;
            } else {
                break;
            }
        }

        // Compute the haplotype start position from the post-skip `ref_pos`
        // rather than the original `ref_start`. When `ref_start` lands inside
        // a deletion span, the pre-loop above has already advanced `ref_pos`
        // to the first surviving reference base; deriving `hap_start` from
        // that post-skip coordinate ensures downstream methylation lookups
        // align with the actual first emitted base.
        #[expect(clippy::cast_possible_truncation, reason = "ref positions fit in u32")]
        let hap_start = self.hap_position_for(ref_pos as u32);

        while bases.len() < fragment_len && ref_pos < reference.len() {
            // Check if the current reference position is a variant start.
            if var_idx < overlapping_indices.len() {
                let var = &self.variant_data[overlapping_indices[var_idx] as usize];
                if var.ref_pos as usize == ref_pos {
                    // Emit alt allele bases.
                    for &b in &var.alt_bases {
                        if bases.len() >= fragment_len {
                            break;
                        }
                        bases.push(b);
                        #[expect(
                            clippy::cast_possible_truncation,
                            reason = "ref positions fit in u32"
                        )]
                        ref_positions.push(ref_pos as u32);
                    }
                    // Skip over reference allele bases.
                    ref_pos += var.ref_len as usize;
                    var_idx += 1;
                    continue;
                }
            }

            // Emit reference base.
            bases.push(reference[ref_pos]);
            #[expect(clippy::cast_possible_truncation, reason = "ref positions fit in u32")]
            ref_positions.push(ref_pos as u32);
            ref_pos += 1;
        }

        // Truncate to exact fragment length (alt alleles may have added extra).
        bases.truncate(fragment_len);
        ref_positions.truncate(fragment_len);

        (bases, ref_positions, hap_start)
    }

    /// Map a 0-based reference position to its 0-based haplotype position.
    ///
    /// Sums net length changes (`alt_len - ref_len`) across every variant
    /// strictly upstream of `ref_pos` (i.e. `var.ref_pos + var.ref_len <=
    /// ref_pos`). Variants that straddle `ref_pos` do not yet contribute —
    /// the caller is presumed to be at a position outside any deletion's
    /// interior (otherwise the haplotype position would be ambiguous).
    ///
    /// Used by [`Self::extract_fragment`] to compute `hap_start` and by
    /// methylation-table construction to relate haplotype scans to per-
    /// haplotype indices. Implementation is `O(log n)` via binary search
    /// over the precomputed [`Self::end_prefix_deltas`] index.
    #[must_use]
    pub fn hap_position_for(&self, ref_pos: u32) -> u32 {
        // partition_point returns the first index where the predicate is
        // false; i.e. the first entry with `var_end > ref_pos`. Subtracting
        // one gives the last entry with `var_end <= ref_pos`, whose
        // cumulative delta is the answer.
        let idx = self.end_prefix_deltas.partition_point(|&(end, _)| end <= ref_pos);
        let delta = if idx == 0 { 0 } else { self.end_prefix_deltas[idx - 1].1 };
        let hp = i64::from(ref_pos) + delta;
        debug_assert!(hp >= 0, "hap_position_for produced negative position");
        #[expect(clippy::cast_sign_loss, reason = "haplotype length is non-negative")]
        #[expect(clippy::cast_possible_truncation, reason = "haplotype length fits in u32")]
        let result = hp as u32;
        result
    }
}

/// Build haplotypes from a set of variant records for one contig.
///
/// For each allele index up to `max_ploidy`, constructs a [`Haplotype`]
/// containing only the variants assigned to that allele.
///
/// Phased genotypes assign alleles deterministically. Unphased genotypes
/// assign non-reference alleles to haplotypes using the provided RNG.
///
/// # Arguments
/// * `variants` — Sorted variant records for this contig.
/// * `max_ploidy` — Maximum ploidy across all variants (e.g. 2 for diploid).
/// * `rng` — Random number generator for unphased genotype assignment.
pub fn build_haplotypes(
    variants: &[VariantRecord],
    max_ploidy: usize,
    rng: &mut impl rand::Rng,
) -> Vec<Haplotype> {
    // Collect variant data and tree intervals per haplotype.
    let mut variant_data_per_hap: Vec<Vec<HaplotypeVariant>> =
        (0..max_ploidy).map(|_| Vec::new()).collect();
    let mut intervals_per_hap: Vec<Vec<Interval<u32>>> =
        (0..max_ploidy).map(|_| Vec::new()).collect();

    for vr in variants {
        let gt = &vr.genotype;

        // For unphased genotypes, generate a random permutation mapping
        // allele indices to haplotype indices. This avoids artificial
        // phasing of nearby variants while ensuring each allele goes to
        // exactly one haplotype (unlike independent random draws, which
        // would incorrectly place both alleles of a hom-alt on the same
        // haplotype 25% of the time for diploid).
        let hap_permutation: Vec<usize> = if gt.is_phased() {
            (0..max_ploidy).collect()
        } else {
            let mut perm: Vec<usize> = (0..max_ploidy).collect();
            // Fisher-Yates shuffle.
            for i in (1..perm.len()).rev() {
                let j = rng.random_range(0..=i);
                perm.swap(i, j);
            }
            perm
        };

        for (allele_idx, allele) in gt.alleles().iter().enumerate() {
            if allele_idx >= max_ploidy {
                break;
            }

            let Some(allele_num) = allele else { continue };
            if *allele_num == 0 {
                continue;
            }

            let Some(alt_bases) = vr.allele_bases(*allele_num) else {
                continue;
            };

            let hap_var = HaplotypeVariant {
                ref_pos: vr.position,
                #[expect(clippy::cast_possible_truncation, reason = "ref allele < 4 GB")]
                ref_len: vr.ref_allele.len() as u32,
                alt_bases: alt_bases.to_vec(),
            };

            let target_hap = hap_permutation[allele_idx];

            // Store the variant and create a tree interval pointing to it.
            let data_idx = variant_data_per_hap[target_hap].len();
            variant_data_per_hap[target_hap].push(hap_var);

            let end_pos = (vr.position as usize + vr.ref_allele.len()).saturating_sub(1);
            #[expect(
                clippy::cast_possible_wrap,
                reason = "genomic coords and variant index < i32::MAX / u32::MAX"
            )]
            #[expect(
                clippy::cast_possible_truncation,
                reason = "genomic coords and variant index < i32::MAX / u32::MAX"
            )]
            let iv = Interval::new(vr.position as i32, end_pos as i32, data_idx as u32);
            intervals_per_hap[target_hap].push(iv);
        }
    }

    variant_data_per_hap
        .into_iter()
        .zip(intervals_per_hap)
        .enumerate()
        .map(|(i, (data, ivs))| {
            let end_prefix_deltas = build_end_prefix_deltas(&data);
            Haplotype {
                allele_index: i,
                variant_data: data,
                variant_tree: COITree::new(&ivs),
                end_prefix_deltas,
            }
        })
        .collect()
}

/// Build the cumulative-delta index used by [`Haplotype::hap_position_for`].
/// Sorts by variant end position so binary search can find the largest
/// entry with `var_end <= ref_pos` for any query.
#[expect(clippy::cast_possible_wrap, reason = "alt_bases.len() and ref_len fit in i64")]
fn build_end_prefix_deltas(variants: &[HaplotypeVariant]) -> Vec<(u32, i64)> {
    let mut entries: Vec<(u32, i64)> = variants
        .iter()
        .map(|v| {
            let var_end = v.ref_pos + v.ref_len;
            let delta = v.alt_bases.len() as i64 - i64::from(v.ref_len);
            (var_end, delta)
        })
        .collect();
    entries.sort_by_key(|&(end, _)| end);
    let mut cum: i64 = 0;
    entries
        .into_iter()
        .map(|(end, delta)| {
            cum += delta;
            (end, cum)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::genotype::Genotype;

    /// Build a simple SNP variant record.
    fn snp(pos: u32, ref_base: u8, alt_base: u8, gt: &str) -> VariantRecord {
        VariantRecord {
            position: pos,
            ref_allele: vec![ref_base],
            alt_alleles: vec![vec![alt_base]],
            genotype: Genotype::parse(gt).unwrap(),
        }
    }

    /// Build an indel variant record.
    fn indel(pos: u32, ref_allele: &[u8], alt_allele: &[u8], gt: &str) -> VariantRecord {
        VariantRecord {
            position: pos,
            ref_allele: ref_allele.to_vec(),
            alt_alleles: vec![alt_allele.to_vec()],
            genotype: Genotype::parse(gt).unwrap(),
        }
    }

    #[test]
    fn test_extract_fragment_no_variants() {
        let reference = b"ACGTACGTACGT";
        let haps = build_haplotypes(&[], 2, &mut rand::rng());
        assert_eq!(haps.len(), 2);

        let (bases, positions, hap_start) = haps[0].extract_fragment(reference, 2, 5);
        assert_eq!(&bases, b"GTACG");
        assert_eq!(&positions, &[2, 3, 4, 5, 6]);
        assert_eq!(hap_start, 2);
    }

    #[test]
    fn test_extract_fragment_with_snp() {
        let reference = b"AAAAAAAA";
        let variants = vec![snp(3, b'A', b'T', "0|1")];
        let haps = build_haplotypes(&variants, 2, &mut rand::rng());

        // Haplotype 0 should have reference (allele 0).
        let (bases, _, hs0) = haps[0].extract_fragment(reference, 0, 8);
        assert_eq!(&bases, b"AAAAAAAA");
        assert_eq!(hs0, 0);

        // Haplotype 1 should have the SNP (allele 1).
        let (bases, _, hs1) = haps[1].extract_fragment(reference, 0, 8);
        assert_eq!(&bases, b"AAATAAAA");
        assert_eq!(hs1, 0);
    }

    #[test]
    fn test_extract_fragment_with_insertion() {
        let reference = b"AAAAAAAA";
        // Insertion: A -> ATT at position 3 (ref allele len 1, alt len 3).
        let variants = vec![indel(3, b"A", b"ATT", "0|1")];
        let haps = build_haplotypes(&variants, 2, &mut rand::rng());

        // Haplotype 1 has the insertion.
        let (bases, positions, hap_start) = haps[1].extract_fragment(reference, 0, 10);
        assert_eq!(&bases, b"AAAATTAAAA");
        // Inserted bases all map back to the ref position of the variant (3).
        assert_eq!(&positions, &[0, 1, 2, 3, 3, 3, 4, 5, 6, 7]);
        // Fragment starts at ref 0, no upstream variants → hap_start == 0.
        assert_eq!(hap_start, 0);
    }

    #[test]
    fn test_extract_fragment_with_deletion() {
        let reference = b"ACGTACGTAC";
        // Deletion: ACG -> A at position 4. Replaces 3 ref bases (ACG at
        // positions 4-6) with 1 alt base (A).
        let variants = vec![indel(4, b"ACG", b"A", "0|1")];
        let haps = build_haplotypes(&variants, 2, &mut rand::rng());

        let (bases, _, _) = haps[1].extract_fragment(reference, 0, 8);
        assert_eq!(&bases, b"ACGTATAC");
    }

    #[test]
    fn test_hom_alt_both_haplotypes_affected() {
        let reference = b"AAAA";
        let variants = vec![snp(1, b'A', b'T', "1|1")];
        let haps = build_haplotypes(&variants, 2, &mut rand::rng());

        let (bases0, _, _) = haps[0].extract_fragment(reference, 0, 4);
        let (bases1, _, _) = haps[1].extract_fragment(reference, 0, 4);
        assert_eq!(&bases0, b"ATAA");
        assert_eq!(&bases1, b"ATAA");
    }

    #[test]
    fn test_phased_allele_assignment() {
        let reference = b"AAAA";
        // Phased 1|0: alt on haplotype 0, ref on haplotype 1.
        let variants = vec![snp(1, b'A', b'T', "1|0")];
        let haps = build_haplotypes(&variants, 2, &mut rand::rng());

        let (bases0, _, _) = haps[0].extract_fragment(reference, 0, 4);
        let (bases1, _, _) = haps[1].extract_fragment(reference, 0, 4);
        assert_eq!(&bases0, b"ATAA");
        assert_eq!(&bases1, b"AAAA");
    }

    #[test]
    fn test_fragment_starts_mid_reference() {
        let reference = b"ACGTACGTAC";
        let variants = vec![snp(5, b'C', b'T', "0|1")];
        let haps = build_haplotypes(&variants, 2, &mut rand::rng());

        // Fragment starting at position 3, length 5: covers pos 3-7.
        let (bases, positions, hap_start) = haps[1].extract_fragment(reference, 3, 5);
        assert_eq!(&bases, b"TATGT");
        assert_eq!(&positions, &[3, 4, 5, 6, 7]);
        // SNP at pos 5 doesn't change net length, and the only variant
        // straddles ref_start=3 → hap_start matches ref_start.
        assert_eq!(hap_start, 3);
    }

    #[test]
    fn test_fragment_starts_within_deletion() {
        // Reference: ACGTACGTAC (positions 0-9)
        // Deletion at pos 2: GTA (3 bases) -> G (1 base)
        // After variant: AC + G + CGTAC = ACGCGTAC
        let reference = b"ACGTACGTAC";
        let variants = vec![indel(2, b"GTA", b"G", "0|1")];
        let haps = build_haplotypes(&variants, 2, &mut rand::rng());

        // Fragment starting at position 3 (mid-deletion). The deletion
        // consumes ref positions 2-4, so starting at 3 means we're inside
        // the deletion. The pre-loop handler should skip past the deletion
        // end (position 5) and continue from there.
        let (bases, _, _) = haps[1].extract_fragment(reference, 3, 5);
        assert_eq!(&bases, b"CGTAC");
    }

    #[test]
    fn test_adjacent_variants() {
        // Two adjacent SNPs with no reference gap between them.
        let reference = b"AAAA";
        let variants = vec![snp(1, b'A', b'T', "0|1"), snp(2, b'A', b'C', "0|1")];
        let haps = build_haplotypes(&variants, 2, &mut rand::rng());

        let (bases, _, _) = haps[1].extract_fragment(reference, 0, 4);
        assert_eq!(&bases, b"ATCA");
    }

    #[test]
    fn test_variant_at_position_zero() {
        let reference = b"ACGT";
        let variants = vec![snp(0, b'A', b'T', "0|1")];
        let haps = build_haplotypes(&variants, 2, &mut rand::rng());

        let (bases, _, _) = haps[1].extract_fragment(reference, 0, 4);
        assert_eq!(&bases, b"TCGT");
    }

    #[test]
    fn test_unphased_hom_alt_both_haplotypes() {
        // Unphased hom-alt must place the alt on both haplotypes.
        let reference = b"AAAA";
        let variants = vec![snp(1, b'A', b'T', "1/1")];
        let haps = build_haplotypes(&variants, 2, &mut rand::rng());

        let (bases0, _, _) = haps[0].extract_fragment(reference, 0, 4);
        let (bases1, _, _) = haps[1].extract_fragment(reference, 0, 4);
        assert_eq!(&bases0, b"ATAA");
        assert_eq!(&bases1, b"ATAA");
    }

    #[test]
    fn test_extract_fragment_starts_within_deletion_hap_start() {
        // Reference: ACGTACGTAC (positions 0-9). Deletion at pos 2 spans
        // 5 bases (REF=GTACG → ALT=G), so var_end = 7 and net delta = -4.
        // A fragment starting at ref_start = 4 lands inside the deletion.
        // The pre-loop advances ref_pos to 7, and the *first emitted base*
        // is reference[7]. In haplotype coordinates that base lives at
        // 7 + (-4) = 3, so `hap_start` must equal 3.
        //
        // Regression: prior code derived `hap_start` from the original
        // `ref_start` (= 4 here), which doesn't account for the straddling
        // deletion. That mis-aligned downstream methylation lookups for
        // any fragment whose start fell inside a deletion span.
        let reference = b"ACGTACGTAC";
        let variants = vec![indel(2, b"GTACG", b"G", "0|1")];
        let haps = build_haplotypes(&variants, 2, &mut rand::rng());

        let (bases, _, hap_start) = haps[1].extract_fragment(reference, 4, 3);
        assert_eq!(&bases, b"TAC");
        assert_eq!(hap_start, 3, "hap_start must reflect the post-skip ref_pos");
    }

    #[test]
    fn test_hap_position_for_with_insertion() {
        // Insertion adds 2 bases at ref pos 3 (A -> ATT, alt-ref = +2).
        let variants = vec![indel(3, b"A", b"ATT", "0|1")];
        let haps = build_haplotypes(&variants, 2, &mut rand::rng());

        // Position 3 itself is inside the variant span → no offset.
        assert_eq!(haps[1].hap_position_for(3), 3);
        // After the variant ends (var_end = 3 + 1 = 4), positions shift by +2.
        assert_eq!(haps[1].hap_position_for(4), 6);
        assert_eq!(haps[1].hap_position_for(7), 9);
    }

    #[test]
    fn test_hap_position_for_with_deletion() {
        // Deletion removes 2 bases at ref pos 4 (ACG -> A, alt-ref = -2).
        let variants = vec![indel(4, b"ACG", b"A", "0|1")];
        let haps = build_haplotypes(&variants, 2, &mut rand::rng());

        // Positions before the deletion end: no offset.
        assert_eq!(haps[1].hap_position_for(4), 4);
        assert_eq!(haps[1].hap_position_for(6), 6);
        // After the deletion ends (var_end = 4 + 3 = 7), positions shift by -2.
        assert_eq!(haps[1].hap_position_for(7), 5);
        assert_eq!(haps[1].hap_position_for(9), 7);
    }
}
