//! Bisulfite-aware genomic edit distance for truth-vs-aligner concordance.
//!
//! `NM:i` / `MD:Z` tags are reference- and convention-dependent: a bisulfite
//! aligner may report edits against the original 4-letter reference (every
//! unmethylated `C->T` is a "mismatch"), against a `C->T`-converted reference
//! (conversions match), or in a bisulfite-aware convention. Comparing two
//! aligners' raw tags therefore measures which convention each picked, not
//! whether they aligned correctly.
//!
//! Instead, [`genomic_edits`] recomputes a read's edits directly against the
//! reference and excludes bisulfite conversions using the read's TRUE strand
//! (taken from the golden truth, so it works even for aligners that emit no
//! `XG`). The result — non-conversion mismatches plus indels — is a
//! convention-independent genomic edit distance that is comparable across
//! aligners. [`RefCache`] loads reference contigs on demand to support it.

use std::collections::BTreeSet;
use std::collections::HashMap;

use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record_buf::Cigar;
use rand::SeedableRng;
use rand::rngs::SmallRng;

use super::golden::ConvDir;
use crate::fasta::Fasta;
use crate::seed::compute_seed;

/// A read's genomic edits against the reference: mismatches and indels that are
/// NOT explained by bisulfite conversion. Independent of the aligner's NM/MD
/// tagging convention.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct GenomicEdits {
    /// Edit distance: non-conversion mismatches + inserted + deleted bases
    /// (the bisulfite-aware analogue of `NM:i`).
    pub nm: u32,
    /// Reference positions of non-conversion mismatches and deletions (the
    /// bisulfite-aware analogue of the `MD:Z`-encoded edit set). Insertions
    /// carry no reference position and so contribute to `nm` only.
    pub positions: BTreeSet<u32>,
}

/// Whether `read_base` at a reference `ref_base` is a valid bisulfite
/// conversion under the read's true strand `conv` (and therefore not a genomic
/// edit): `C->T` on the CT strand, `G->A` on the GA strand.
fn is_conversion(ref_base: u8, read_base: u8, conv: Option<ConvDir>) -> bool {
    match conv {
        Some(ConvDir::Ct) => ref_base == b'C' && read_base == b'T',
        Some(ConvDir::Ga) => ref_base == b'G' && read_base == b'A',
        None => false,
    }
}

/// Compute a read's genomic edits against `ref_seq` (the full contig the read
/// aligns to), walking its `cigar` from 0-based `start0`. `conv` is the read's
/// TRUE bisulfite strand; conversions consistent with it are excluded so the
/// result reflects only genuine genomic differences.
#[must_use]
pub fn genomic_edits(
    seq: &[u8],
    cigar: &Cigar,
    start0: u32,
    ref_seq: &[u8],
    conv: Option<ConvDir>,
) -> GenomicEdits {
    let mut edits = GenomicEdits::default();
    let mut ref_pos = start0 as usize;
    let mut read_pos = 0usize;
    for op in cigar.as_ref() {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for k in 0..len {
                    let r = ref_seq.get(ref_pos + k).map(u8::to_ascii_uppercase);
                    let q = seq.get(read_pos + k).map(u8::to_ascii_uppercase);
                    if let (Some(r), Some(q)) = (r, q)
                        && r != q
                        && !is_conversion(r, q, conv)
                    {
                        edits.nm += 1;
                        edits.positions.insert(u32::try_from(ref_pos + k).unwrap_or(0));
                    }
                }
                ref_pos += len;
                read_pos += len;
            }
            Kind::Insertion => {
                edits.nm += u32::try_from(len).unwrap_or(0);
                read_pos += len;
            }
            Kind::Deletion => {
                edits.nm += u32::try_from(len).unwrap_or(0);
                for k in 0..len {
                    edits.positions.insert(u32::try_from(ref_pos + k).unwrap_or(0));
                }
                ref_pos += len;
            }
            Kind::Skip => ref_pos += len,
            Kind::SoftClip => read_pos += len,
            Kind::HardClip | Kind::Pad => {}
        }
    }
    edits
}

/// On-demand cache of uppercased reference contig sequences over an indexed
/// FASTA, so the eval pass loads only the contigs its reads actually touch.
pub struct RefCache {
    fasta: Fasta,
    cache: HashMap<String, Option<Vec<u8>>>,
}

impl RefCache {
    /// Wrap an opened reference FASTA.
    #[must_use]
    pub fn new(fasta: Fasta) -> Self {
        Self { fasta, cache: HashMap::new() }
    }

    /// The uppercased sequence of `contig`, loaded and cached on first use.
    /// Returns `None` (cached) when the contig is absent from the reference.
    pub fn contig(&mut self, contig: &str) -> Option<&[u8]> {
        if !self.cache.contains_key(contig) {
            // IUPAC ambiguity codes (rare, and not at variant sites) are
            // resolved randomly; seed a fresh RNG per contig from a deterministic
            // FNV-1a hash of the contig name so the resolution does not depend on
            // the order contigs are first requested.
            let mut rng = SmallRng::seed_from_u64(compute_seed(contig));
            // Cache the attempt either way (so a failure is not retried). Warn
            // on a genuine load error rather than silently dropping the contig's
            // reads from NM/MD concordance — it is indistinguishable from an
            // absent contig at the call site otherwise.
            let loaded = match self.fasta.load_contig(contig, &mut rng) {
                Ok(seq) => Some(seq),
                Err(e) => {
                    log::warn!(
                        "reference contig {contig:?} unavailable ({e:#}); \
                         its reads are excluded from NM/MD concordance"
                    );
                    None
                }
            };
            self.cache.insert(contig.to_string(), loaded);
        }
        self.cache.get(contig).and_then(Option::as_deref)
    }
}

#[cfg(test)]
mod tests {
    use noodles::sam::alignment::record::cigar::op::Op;

    use super::*;

    fn cigar(ops: &[(Kind, usize)]) -> Cigar {
        Cigar::from(ops.iter().map(|&(k, n)| Op::new(k, n)).collect::<Vec<_>>())
    }

    #[test]
    fn counts_a_plain_mismatch_as_a_genomic_edit() {
        // read T vs ref A at ref pos 12 (offset 2 in a 10M from start 10).
        let edits = genomic_edits(
            b"AATAAAAAAA",
            &cigar(&[(Kind::Match, 10)]),
            10,
            b"AAAAAAAAAAAAAAA",
            None,
        );
        assert_eq!(edits.nm, 1);
        assert_eq!(edits.positions.iter().copied().collect::<Vec<_>>(), vec![12]);
    }

    #[test]
    fn excludes_bisulfite_conversion_on_its_strand() {
        // ref C, read T at every position. On the CT strand these are all
        // conversions (no genomic edits); with no strand they are all edits.
        let seq = b"TTTT";
        let cig = cigar(&[(Kind::Match, 4)]);
        let reference = b"CCCC";
        assert_eq!(genomic_edits(seq, &cig, 0, reference, Some(ConvDir::Ct)).nm, 0);
        assert_eq!(genomic_edits(seq, &cig, 0, reference, None).nm, 4);
        // On the GA strand a C->T is NOT the freed cell, so it stays an edit.
        assert_eq!(genomic_edits(seq, &cig, 0, reference, Some(ConvDir::Ga)).nm, 4);
    }

    #[test]
    fn a_real_variant_survives_conversion_masking() {
        // ref C, read A (a transversion variant) on the CT strand: A is not the
        // conversion product T, so it is a genuine genomic edit.
        let edits =
            genomic_edits(b"A", &cigar(&[(Kind::Match, 1)]), 5, b"CCCCCCC", Some(ConvDir::Ct));
        assert_eq!(edits.nm, 1);
        assert!(edits.positions.contains(&5));
    }

    #[test]
    fn counts_indels_in_nm_and_deletions_in_positions() {
        // 2M1I2M1D2M over ref AAAAAAA: insertion adds to nm only; deletion adds
        // to nm and contributes its reference position.
        let cig = cigar(&[
            (Kind::Match, 2),
            (Kind::Insertion, 1),
            (Kind::Match, 2),
            (Kind::Deletion, 1),
            (Kind::Match, 2),
        ]);
        // read consumes 2+1+2+0+2 = 7 bases, all matching ref where aligned.
        let edits = genomic_edits(b"AAAAAAA", &cig, 0, b"AAAAAAAAAA", None);
        assert_eq!(edits.nm, 2); // 1 inserted + 1 deleted
        // deletion reference position: 2M(0,1) 2M(2,3) D at ref 4.
        assert!(edits.positions.contains(&4));
    }
}
