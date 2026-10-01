//! CIGAR geometry helpers shared by the truth-aware eval metrics.
//!
//! Both operate on a [`Cigar`] paired with the alignment's 0-based reference
//! start. [`reference_len`] gives the reference span consumed;
//! [`ref_pos_to_read_offset`] maps a reference position to the read offset of
//! the base aligned there, or `None` when that position is deleted, skipped,
//! or otherwise not covered by an aligned (`M`/`=`/`X`) operation.

use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record_buf::Cigar;

/// Whether a CIGAR operation consumes reference bases.
fn consumes_reference(kind: Kind) -> bool {
    matches!(
        kind,
        Kind::Match | Kind::Deletion | Kind::Skip | Kind::SequenceMatch | Kind::SequenceMismatch
    )
}

/// Whether a CIGAR operation consumes query (read) bases.
fn consumes_query(kind: Kind) -> bool {
    matches!(
        kind,
        Kind::Match
            | Kind::Insertion
            | Kind::SoftClip
            | Kind::SequenceMatch
            | Kind::SequenceMismatch
    )
}

/// Number of reference bases the alignment spans.
#[must_use]
pub fn reference_len(cigar: &Cigar) -> u32 {
    let mut len: u32 = 0;
    for op in cigar.as_ref() {
        if consumes_reference(op.kind()) {
            len += u32::try_from(op.len()).unwrap_or(0);
        }
    }
    len
}

/// Map a 0-based reference position to the 0-based read offset aligned there.
///
/// `aln_start0` is the alignment's 0-based reference start. Returns `None` when
/// `target_ref0` falls outside the alignment, inside a deletion/skip, or is
/// otherwise not covered by an aligned base.
#[must_use]
pub fn ref_pos_to_read_offset(cigar: &Cigar, aln_start0: u32, target_ref0: u32) -> Option<usize> {
    if target_ref0 < aln_start0 {
        return None;
    }
    let mut ref_pos = aln_start0;
    let mut read_pos: usize = 0;
    for op in cigar.as_ref() {
        let kind = op.kind();
        let span = u32::try_from(op.len()).unwrap_or(0);
        let consumes_ref = consumes_reference(kind);
        let consumes_q = consumes_query(kind);

        if consumes_ref && consumes_q {
            // Aligned run: target may land within it.
            if target_ref0 < ref_pos + span {
                let within = (target_ref0 - ref_pos) as usize;
                return Some(read_pos + within);
            }
            ref_pos += span;
            read_pos += op.len();
        } else if consumes_ref {
            // Deletion/skip: target inside this run is not represented by a base.
            if target_ref0 < ref_pos + span {
                return None;
            }
            ref_pos += span;
        } else if consumes_q {
            // Insertion/soft-clip: advances the read only.
            read_pos += op.len();
        }
    }
    None
}

/// Invoke `f(read_offset, ref_pos0)` for each aligned (`M`/`=`/`X`) base,
/// walking the CIGAR from `aln_start0`. Insertions and soft-clips advance the
/// read only; deletions and skips advance the reference only.
pub fn for_each_aligned(cigar: &Cigar, aln_start0: u32, mut f: impl FnMut(usize, u32)) {
    let mut ref_pos = aln_start0;
    let mut read_pos: usize = 0;
    for op in cigar.as_ref() {
        let kind = op.kind();
        let len = op.len();
        let span = u32::try_from(len).unwrap_or(0);
        match (consumes_reference(kind), consumes_query(kind)) {
            (true, true) => {
                for k in 0..len {
                    f(read_pos + k, ref_pos + u32::try_from(k).unwrap_or(0));
                }
                ref_pos += span;
                read_pos += len;
            }
            (true, false) => ref_pos += span,
            (false, true) => read_pos += len,
            (false, false) => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::record::cigar::op::Op;

    fn cigar(ops: &[(Kind, usize)]) -> Cigar {
        Cigar::from(ops.iter().map(|&(k, n)| Op::new(k, n)).collect::<Vec<_>>())
    }

    #[test]
    fn reference_len_sums_ref_consuming_ops() {
        // 10M2I5M3D4M -> reference span 10 + 5 + 3 + 4 = 22 (insertion excluded).
        let c = cigar(&[
            (Kind::Match, 10),
            (Kind::Insertion, 2),
            (Kind::Match, 5),
            (Kind::Deletion, 3),
            (Kind::Match, 4),
        ]);
        assert_eq!(reference_len(&c), 22);
    }

    #[test]
    fn ref_offset_simple_match() {
        // 100M starting at ref 1000: ref 1005 -> read offset 5.
        let c = cigar(&[(Kind::Match, 100)]);
        assert_eq!(ref_pos_to_read_offset(&c, 1000, 1005), Some(5));
        assert_eq!(ref_pos_to_read_offset(&c, 1000, 1000), Some(0));
    }

    #[test]
    fn ref_offset_after_insertion_shifts_read() {
        // 5M2I5M at ref 1000: ref 1006 is in the second M; read offset = 5 (M) + 2 (I) + 1.
        let c = cigar(&[(Kind::Match, 5), (Kind::Insertion, 2), (Kind::Match, 5)]);
        assert_eq!(ref_pos_to_read_offset(&c, 1000, 1006), Some(8));
    }

    #[test]
    fn ref_offset_after_deletion_shifts_ref() {
        // 5M3D5M at ref 1000: ref 1008 is the first base after the deletion;
        // read offset = 5 (no read bases consumed by D).
        let c = cigar(&[(Kind::Match, 5), (Kind::Deletion, 3), (Kind::Match, 5)]);
        assert_eq!(ref_pos_to_read_offset(&c, 1000, 1008), Some(5));
    }

    #[test]
    fn ref_offset_inside_deletion_is_none() {
        let c = cigar(&[(Kind::Match, 5), (Kind::Deletion, 3), (Kind::Match, 5)]);
        // ref 1006 is inside the 3bp deletion (1005..1008).
        assert_eq!(ref_pos_to_read_offset(&c, 1000, 1006), None);
    }

    #[test]
    fn ref_offset_softclip_offsets_read() {
        // 4S10M at ref 1000: ref 1000 -> read offset 4 (soft-clip consumes read).
        let c = cigar(&[(Kind::SoftClip, 4), (Kind::Match, 10)]);
        assert_eq!(ref_pos_to_read_offset(&c, 1000, 1000), Some(4));
    }

    #[test]
    fn ref_offset_out_of_span_is_none() {
        let c = cigar(&[(Kind::Match, 10)]);
        assert_eq!(ref_pos_to_read_offset(&c, 1000, 1010), None); // 1000..1010 only
        assert_eq!(ref_pos_to_read_offset(&c, 1000, 999), None);
    }
}
