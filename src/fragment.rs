//! Fragment extraction and read base generation.
//!
//! Handles extracting genomic fragments from haplotypes, reverse-complement
//! operations, adapter sequence appending for short fragments, and R1/R2
//! base extraction from fragments.
//!
//! # Lowercase-as-marker convention
//!
//! [`crate::fasta::Fasta::load_contig`] stores bases resolved from IUPAC
//! ambiguity codes (including `N`) in lowercase, while real `A`/`C`/`G`/`T`
//! bases remain uppercase. Reads that carry those synthesized bases through
//! from the reference can be detected by counting lowercase bytes via
//! [`lowercase_fraction`], and the buffers are upper-cased in place with
//! [`uppercase_in_place`] before being emitted.

use crate::haplotype::Haplotype;

/// The DNA complement of a base. Preserves the case of the input so the
/// lowercase "synthesized from ambiguity" marker (see module docs) survives
/// reverse-complementing for R2 reads.
#[must_use]
fn complement(base: u8) -> u8 {
    let upper = match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    };
    if base.is_ascii_lowercase() { upper.to_ascii_lowercase() } else { upper }
}

/// Reverse-complement a DNA sequence in place.
pub fn reverse_complement(seq: &mut [u8]) {
    seq.reverse();
    for base in seq.iter_mut() {
        *base = complement(*base);
    }
}

/// A simulated fragment extracted from a haplotype.
///
/// Contains the fragment bases and reference coordinate mapping, plus
/// metadata about the source haplotype and strand.
#[derive(Debug)]
pub struct Fragment {
    /// Fragment bases in the forward orientation (5' to 3').
    pub bases: Vec<u8>,
    /// Reference position for each base in the fragment.
    pub ref_positions: Vec<u32>,
    /// 0-based reference start position of the fragment.
    pub ref_start: u32,
    /// Whether the fragment is on the forward strand.
    pub is_forward: bool,
    /// Index of the haplotype this fragment came from.
    pub haplotype_index: usize,
}

/// Extract a fragment from a haplotype at a given position and strand.
///
/// The fragment is extracted in the forward orientation from the haplotype,
/// then reverse-complemented if on the reverse strand.
///
/// # Arguments
/// * `haplotype` — The haplotype to extract from.
/// * `reference` — The full reference sequence for this contig.
/// * `ref_start` — 0-based start position on the reference.
/// * `fragment_len` — Desired fragment length in bases.
/// * `is_forward` — Whether to read the forward or reverse strand.
#[must_use]
pub fn extract_fragment(
    haplotype: &Haplotype,
    reference: &[u8],
    ref_start: u32,
    fragment_len: usize,
    is_forward: bool,
) -> Fragment {
    let (bases, ref_positions) = haplotype.extract_fragment(reference, ref_start, fragment_len);

    Fragment {
        bases,
        ref_positions,
        ref_start,
        is_forward,
        haplotype_index: haplotype.allele_index(),
    }
}

/// Extract read bases from a fragment, handling adapter sequences for short
/// fragments.
///
/// For R1: takes the first `read_length` bases from the fragment. If the
/// fragment is shorter than `read_length`, appends adapter bases.
///
/// For R2: takes the last `read_length` bases from the fragment and
/// reverse-complements them. If the fragment is shorter than `read_length`,
/// appends adapter bases after the reverse-complemented fragment.
///
/// # Arguments
/// * `fragment` — The source fragment bases.
/// * `read_length` — Desired read length.
/// * `adapter` — Adapter sequence to append for short fragments.
/// * `is_r2` — If true, extract from the end and reverse complement.
#[must_use]
pub fn extract_read_bases(
    fragment_bases: &[u8],
    read_length: usize,
    adapter: &[u8],
    is_r2: bool,
) -> Vec<u8> {
    let frag_len = fragment_bases.len();

    if is_r2 {
        // R2: last `read_length` bases from fragment, reverse-complemented,
        // then adapter.
        let mut bases = Vec::with_capacity(read_length);
        let take = frag_len.min(read_length);
        let start = frag_len.saturating_sub(take);
        bases.extend_from_slice(&fragment_bases[start..]);
        reverse_complement(&mut bases);
        append_adapter_and_pad(&mut bases, read_length, adapter);
        bases
    } else {
        // R1: first `read_length` bases from fragment, then adapter.
        let mut bases = Vec::with_capacity(read_length);
        let take = frag_len.min(read_length);
        bases.extend_from_slice(&fragment_bases[..take]);
        append_adapter_and_pad(&mut bases, read_length, adapter);
        bases
    }
}

/// Fraction of lowercase ASCII bytes in `bases`.
///
/// Lowercase bases mark positions that were synthesized from IUPAC
/// ambiguity codes (including `N`) at reference-load time, so this is a
/// cheap proxy for "how much of this read came from an ambiguous reference
/// region." Returns `0.0` for an empty slice.
#[must_use]
pub fn lowercase_fraction(bases: &[u8]) -> f64 {
    if bases.is_empty() {
        return 0.0;
    }
    // The ASCII 0x20 bit distinguishes lowercase from uppercase letters; any
    // non-alphabetic byte (e.g. adapter pad) has this bit clear for our
    // valid base values (A-Z / a-z) — so counting that bit is equivalent to
    // counting lowercase letters.
    let lower = bases.iter().filter(|&&b| b.is_ascii_lowercase()).count();
    lower as f64 / bases.len() as f64
}

/// Upper-case every ASCII letter in `bases` in place.
pub fn uppercase_in_place(bases: &mut [u8]) {
    for b in bases.iter_mut() {
        b.make_ascii_uppercase();
    }
}

/// Append adapter bases and pad with N to reach `target_len`.
fn append_adapter_and_pad(bases: &mut Vec<u8>, target_len: usize, adapter: &[u8]) {
    if bases.len() < target_len {
        let need = target_len - bases.len();
        let adapter_take = need.min(adapter.len());
        bases.extend_from_slice(&adapter[..adapter_take]);
        // Pad with N if adapter is also too short.
        while bases.len() < target_len {
            bases.push(b'N');
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        let mut seq = b"ACGT".to_vec();
        reverse_complement(&mut seq);
        assert_eq!(&seq, b"ACGT"); // palindrome

        let mut seq2 = b"AACG".to_vec();
        reverse_complement(&mut seq2);
        assert_eq!(&seq2, b"CGTT");
    }

    #[test]
    fn test_reverse_complement_with_n() {
        let mut seq = b"ANGC".to_vec();
        reverse_complement(&mut seq);
        assert_eq!(&seq, b"GCNT");
    }

    #[test]
    fn test_reverse_complement_preserves_lowercase() {
        // Lowercase bases carry an ambiguity-resolved marker that must
        // survive reverse-complementing.
        let mut seq = b"aCgT".to_vec();
        reverse_complement(&mut seq);
        assert_eq!(&seq, b"AcGt");
    }

    #[test]
    fn test_extract_r1_full_fragment() {
        let fragment = b"ACGTACGTAC";
        let bases = extract_read_bases(fragment, 10, b"ADAPTER", false);
        assert_eq!(&bases, b"ACGTACGTAC");
    }

    #[test]
    fn test_extract_r1_fragment_longer_than_read() {
        let fragment = b"ACGTACGTAC";
        let bases = extract_read_bases(fragment, 5, b"ADAPTER", false);
        assert_eq!(&bases, b"ACGTA");
    }

    #[test]
    fn test_extract_r1_short_fragment_with_adapter() {
        let fragment = b"ACG";
        let adapter = b"TTTTTT";
        let bases = extract_read_bases(fragment, 8, adapter, false);
        assert_eq!(&bases, b"ACGTTTTT");
    }

    #[test]
    fn test_extract_r2_full_fragment() {
        let fragment = b"ACGTACGTAC";
        let bases = extract_read_bases(fragment, 10, b"ADAPTER", true);
        // Reverse complement of ACGTACGTAC is GTACGTACGT
        assert_eq!(&bases, b"GTACGTACGT");
    }

    #[test]
    fn test_extract_r2_fragment_longer_than_read() {
        let fragment = b"ACGTACGTAC";
        let bases = extract_read_bases(fragment, 5, b"ADAPTER", true);
        // Last 5 bases: CGTAC, reverse complement: GTACG
        assert_eq!(&bases, b"GTACG");
    }

    #[test]
    fn test_extract_r2_short_fragment_with_adapter() {
        let fragment = b"ACG";
        let adapter = b"TTTTTT";
        let bases = extract_read_bases(fragment, 8, adapter, true);
        // Reverse complement of ACG is CGT (3 bytes), then 5 adapter T's.
        assert_eq!(&bases, b"CGTTTTTT"); // CGT + TTTTT = 8 bytes
        assert_eq!(bases.len(), 8);
    }

    #[test]
    fn test_lowercase_fraction_empty() {
        assert!(lowercase_fraction(b"").abs() < 1e-12);
    }

    #[test]
    fn test_lowercase_fraction_all_upper() {
        assert!(lowercase_fraction(b"ACGTACGT").abs() < 1e-12);
    }

    #[test]
    fn test_lowercase_fraction_all_lower() {
        assert!((lowercase_fraction(b"acgtacgt") - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_lowercase_fraction_mixed() {
        // 3 lowercase out of 10.
        let f = lowercase_fraction(b"ACaGcTAtCA");
        assert!((f - 0.3).abs() < 1e-10, "expected 0.3, got {f}");
    }

    #[test]
    fn test_lowercase_fraction_ignores_non_letters() {
        // N-pad and '-' are non-letter / uppercase — only 'a' should count.
        // "A-Na" is length 4 with 1 lowercase letter.
        let f = lowercase_fraction(b"A-Na");
        assert!((f - 0.25).abs() < 1e-10, "expected 0.25, got {f}");
    }

    #[test]
    fn test_uppercase_in_place() {
        let mut bases = b"aCgTnN".to_vec();
        uppercase_in_place(&mut bases);
        assert_eq!(&bases, b"ACGTNN");
    }

    #[test]
    fn test_uppercase_in_place_empty() {
        let mut bases: Vec<u8> = Vec::new();
        uppercase_in_place(&mut bases);
        assert!(bases.is_empty());
    }

    #[test]
    fn test_extract_empty_fragment_all_adapter() {
        let fragment = b"";
        let adapter = b"AGATCGG";
        let bases = extract_read_bases(fragment, 5, b"AGATCGG", false);
        assert_eq!(&bases, b"AGATC");

        let bases_r2 = extract_read_bases(fragment, 5, adapter, true);
        assert_eq!(&bases_r2, b"AGATC");
    }
}
