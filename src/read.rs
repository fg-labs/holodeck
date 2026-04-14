//! Simulated read pair generation.
//!
//! Combines fragment extraction, error model application, and read naming
//! into complete [`ReadPair`] objects ready for FASTQ output. Computes
//! proper CIGARs from the haplotype-to-reference coordinate mapping for
//! golden BAM output.

use noodles::sam::alignment::record::cigar::op::{Kind, Op};
use noodles::sam::alignment::record_buf::Cigar;
use rand::Rng;

use crate::error_model::{self, ErrorModel, ReadEnd};
use crate::fragment::{Fragment, extract_read_bases};
use crate::read_naming::{TruthAlignment, encoded_pe_name, encoded_se_name, simple_name};

/// A single simulated read with bases, quality scores, and metadata.
#[derive(Debug, Clone)]
pub struct SimulatedRead {
    /// Read name (shared between R1 and R2 of a pair).
    pub name: String,
    /// Base sequence (possibly with errors applied).
    pub bases: Vec<u8>,
    /// Quality scores (Phred+33 encoded).
    pub qualities: Vec<u8>,
}

/// A simulated read pair (or single read for SE mode).
#[derive(Debug)]
pub struct ReadPair {
    /// First read (always present).
    pub read1: SimulatedRead,
    /// Second read (present only for paired-end mode).
    pub read2: Option<SimulatedRead>,
    /// Truth alignment for R1.
    pub r1_truth: TruthAlignment,
    /// Truth alignment for R2 (present only for paired-end mode).
    pub r2_truth: Option<TruthAlignment>,
    /// Truth CIGAR for R1, reflecting haplotype variants and adapter
    /// soft-clipping.
    pub r1_cigar: Cigar,
    /// Truth CIGAR for R2. `None` for single-end reads.
    pub r2_cigar: Option<Cigar>,
}

/// Generate a simulated read pair from a fragment.
///
/// Extracts R1 and optionally R2 bases from the fragment, applies the error
/// model, constructs read names with truth information, and computes CIGARs
/// from the fragment's reference coordinate mapping.
///
/// # Arguments
/// * `fragment` — The source fragment with bases and reference positions.
/// * `contig_name` — Contig name for read naming.
/// * `read_num` — 1-based read pair number.
/// * `read_length` — Desired read length.
/// * `paired` — Whether to generate both R1 and R2.
/// * `adapter_r1` — Adapter sequence for R1.
/// * `adapter_r2` — Adapter sequence for R2.
/// * `error_model` — Error model to apply.
/// * `simple_names` — Use simple names instead of encoded truth names.
/// * `rng` — Random number generator.
#[allow(clippy::too_many_arguments)] // Simulation parameters are naturally numerous
pub fn generate_read_pair(
    fragment: &Fragment,
    contig_name: &str,
    read_num: u64,
    read_length: usize,
    paired: bool,
    adapter_r1: &[u8],
    adapter_r2: &[u8],
    model: &impl ErrorModel,
    simple_names: bool,
    rng: &mut impl Rng,
) -> ReadPair {
    let frag_len = fragment.bases.len();
    let genomic = frag_len.min(read_length);
    let adapter_bases = read_length.saturating_sub(genomic);
    let right_start = frag_len.saturating_sub(genomic);

    // For a forward (top-strand) fragment, R1 reads from the left end
    // (forward strand) and R2 reads from the right end (reverse strand).
    // For a reverse (bottom-strand) fragment, R1 reads from the right end
    // (reverse strand) and R2 reads from the left end (forward strand).
    // The left read is always forward-strand and the right read is always
    // reverse-strand; what changes is which becomes R1 vs R2.
    let r1_negative_strand = !fragment.is_forward;

    // Extract R1 bases.
    let mut r1_bases =
        extract_read_bases(&fragment.bases, read_length, adapter_r1, r1_negative_strand);
    let (r1_errors, r1_quals) =
        error_model::apply_errors(model, &mut r1_bases, ReadEnd::Read1, rng);

    // R1 ref_positions and CIGAR.  Fragment ref_positions are always in
    // ascending (forward) reference order; negative-strand reads take from
    // the right end of the fragment but the positions are still ascending.
    let r1_positions = if r1_negative_strand {
        &fragment.ref_positions[right_start..frag_len]
    } else {
        &fragment.ref_positions[..genomic]
    };
    let r1_cigar = cigar_from_ref_positions(r1_positions, adapter_bases, r1_negative_strand);

    // R1 truth position: leftmost reference coordinate (1-based).
    let r1_ref_pos = if fragment.ref_positions.is_empty() {
        0
    } else if r1_negative_strand {
        fragment.ref_positions[right_start] + 1
    } else {
        fragment.ref_positions[0] + 1
    };

    let r1_truth = TruthAlignment {
        contig: contig_name.to_string(),
        position: r1_ref_pos,
        is_forward: fragment.is_forward,
        haplotype: fragment.haplotype_index,
        n_errors: r1_errors,
    };

    if !paired {
        let name =
            if simple_names { simple_name(read_num) } else { encoded_se_name(read_num, &r1_truth) };

        return ReadPair {
            read1: SimulatedRead { name, bases: r1_bases, qualities: r1_quals },
            read2: None,
            r1_truth,
            r2_truth: None,
            r1_cigar,
            r2_cigar: None,
        };
    }

    // Paired-end: extract R2 bases from the opposite end of the fragment.
    // R2 is on the negative strand whenever the fragment is forward-strand.
    let r2_negative_strand = fragment.is_forward;
    let mut r2_bases =
        extract_read_bases(&fragment.bases, read_length, adapter_r2, r2_negative_strand);
    let (r2_errors, r2_quals) =
        error_model::apply_errors(model, &mut r2_bases, ReadEnd::Read2, rng);

    // R2 ref_positions and CIGAR.
    let r2_positions = if r2_negative_strand {
        &fragment.ref_positions[right_start..frag_len]
    } else {
        &fragment.ref_positions[..genomic]
    };
    let r2_cigar = cigar_from_ref_positions(r2_positions, adapter_bases, r2_negative_strand);

    // R2 truth position: leftmost reference coordinate (1-based).
    let r2_ref_pos = if fragment.ref_positions.is_empty() {
        0
    } else if r2_negative_strand {
        fragment.ref_positions[right_start] + 1
    } else {
        fragment.ref_positions[0] + 1
    };

    let r2_truth = TruthAlignment {
        contig: contig_name.to_string(),
        position: r2_ref_pos,
        is_forward: !fragment.is_forward,
        haplotype: fragment.haplotype_index,
        n_errors: r2_errors,
    };

    let name = if simple_names {
        simple_name(read_num)
    } else {
        encoded_pe_name(read_num, &r1_truth, &r2_truth)
    };

    ReadPair {
        read1: SimulatedRead { name: name.clone(), bases: r1_bases, qualities: r1_quals },
        read2: Some(SimulatedRead { name, bases: r2_bases, qualities: r2_quals }),
        r1_truth,
        r2_truth: Some(r2_truth),
        r1_cigar,
        r2_cigar: Some(r2_cigar),
    }
}

/// Compute a CIGAR from a slice of ascending reference positions, plus
/// optional adapter soft-clipping.
///
/// Positions must be in ascending order (forward strand). Consecutive
/// positions incrementing by 1 produce M ops, same position produces I ops,
/// and gaps produce D ops.
///
/// This works for both R1 and R2: BAM CIGARs are always expressed in forward
/// reference order from the leftmost aligned position, so even negative-strand
/// reads use ascending positions.
///
/// Adapter bases are appended as a soft-clip (S) operation. For
/// negative-strand reads (`negative_strand = true`) the adapter is sequenced
/// at the 3' end of the read but sits at the left (5') end of the stored BAM
/// record after reverse-complementing, so the S op is placed at the *start*
/// of the CIGAR. For forward-strand reads the S op is placed at the *end*.
#[must_use]
pub fn cigar_from_ref_positions(
    positions: &[u32],
    adapter_bases: usize,
    negative_strand: bool,
) -> Cigar {
    let mut ops: Vec<Op> = Vec::new();

    if positions.is_empty() {
        if adapter_bases > 0 {
            ops.push(Op::new(Kind::SoftClip, adapter_bases));
        }
        return Cigar::from(ops);
    }

    let mut match_run: usize = 1; // First base is always a match.
    let mut ins_run: usize = 0;

    for i in 1..positions.len() {
        let prev = positions[i - 1];
        let curr = positions[i];

        if curr == prev + 1 {
            // Sequential ascending position: flush any insertion, extend match.
            if ins_run > 0 {
                ops.push(Op::new(Kind::Insertion, ins_run));
                ins_run = 0;
            }
            match_run += 1;
        } else if curr == prev {
            // Same position: insertion base. Flush match run, extend insertion.
            if match_run > 0 {
                ops.push(Op::new(Kind::Match, match_run));
                match_run = 0;
            }
            ins_run += 1;
        } else {
            // Gap in positions: deletion. Flush current runs.
            if ins_run > 0 {
                ops.push(Op::new(Kind::Insertion, ins_run));
                ins_run = 0;
            }
            if match_run > 0 {
                ops.push(Op::new(Kind::Match, match_run));
            }
            let gap = curr.saturating_sub(prev).saturating_sub(1);
            if gap > 0 {
                ops.push(Op::new(Kind::Deletion, gap as usize));
            }
            match_run = 1; // Current base starts a new match.
        }
    }

    // Flush remaining runs.
    if ins_run > 0 {
        ops.push(Op::new(Kind::Insertion, ins_run));
    }
    if match_run > 0 {
        ops.push(Op::new(Kind::Match, match_run));
    }

    // Adapter bases as soft-clip. Negative-strand reads are stored
    // reverse-complemented in BAM, so the adapter (at the 3' end in read
    // order) moves to the left (5') end of the record — a leading S op.
    if adapter_bases > 0 {
        if negative_strand {
            ops.insert(0, Op::new(Kind::SoftClip, adapter_bases));
        } else {
            ops.push(Op::new(Kind::SoftClip, adapter_bases));
        }
    }

    Cigar::from(ops)
}

/// Format a CIGAR as a human-readable string (e.g. "50M3I20M2S").
#[must_use]
pub fn cigar_to_string(cigar: &Cigar) -> String {
    use std::fmt::Write;
    cigar.as_ref().iter().fold(String::new(), |mut s, op| {
        let kind_char = match op.kind() {
            Kind::Match => 'M',
            Kind::Insertion => 'I',
            Kind::Deletion => 'D',
            Kind::SoftClip => 'S',
            Kind::HardClip => 'H',
            Kind::Skip => 'N',
            Kind::Pad => 'P',
            Kind::SequenceMatch => '=',
            Kind::SequenceMismatch => 'X',
        };
        let _ = write!(s, "{}{kind_char}", op.len());
        s
    })
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand::rngs::SmallRng;

    use super::*;
    use crate::error_model::illumina::IlluminaErrorModel;

    /// Build a simple fragment for testing.
    fn test_fragment(bases: &[u8], ref_start: u32) -> Fragment {
        #[expect(clippy::cast_possible_truncation, reason = "test data is small")]
        let ref_positions: Vec<u32> = (ref_start..ref_start + bases.len() as u32).collect();
        Fragment {
            bases: bases.to_vec(),
            ref_positions,
            ref_start,
            is_forward: true,
            haplotype_index: 0,
        }
    }

    // --- CIGAR generation tests ---

    #[test]
    fn test_cigar_all_match() {
        let cigar = cigar_from_ref_positions(&[0, 1, 2, 3, 4], 0, false);
        assert_eq!(cigar_to_string(&cigar), "5M");
    }

    #[test]
    fn test_cigar_with_insertion() {
        // Positions 0,1,2,2,2,3,4: two inserted bases at ref pos 2.
        let cigar = cigar_from_ref_positions(&[0, 1, 2, 2, 2, 3, 4], 0, false);
        assert_eq!(cigar_to_string(&cigar), "3M2I2M");
    }

    #[test]
    fn test_cigar_with_deletion() {
        // Gap from 2 to 5: 2 deleted ref bases.
        let cigar = cigar_from_ref_positions(&[0, 1, 2, 5, 6], 0, false);
        assert_eq!(cigar_to_string(&cigar), "3M2D2M");
    }

    #[test]
    fn test_cigar_with_adapter_softclip_forward() {
        // Forward-strand: adapter is a trailing soft-clip.
        let cigar = cigar_from_ref_positions(&[0, 1, 2], 2, false);
        assert_eq!(cigar_to_string(&cigar), "3M2S");
    }

    #[test]
    fn test_cigar_with_adapter_softclip_negative_strand() {
        // Negative-strand: adapter moves to a leading soft-clip after RC.
        let cigar = cigar_from_ref_positions(&[0, 1, 2], 2, true);
        assert_eq!(cigar_to_string(&cigar), "2S3M");
    }

    #[test]
    fn test_cigar_all_adapter() {
        // All-adapter read: placement doesn't matter, but it still round-trips.
        let cigar = cigar_from_ref_positions(&[], 5, false);
        assert_eq!(cigar_to_string(&cigar), "5S");
    }

    #[test]
    fn test_cigar_with_insertion_and_deletion() {
        // Insertion at pos 2 (two extra bases), then deletion of 2 ref bases.
        let cigar = cigar_from_ref_positions(&[0, 1, 2, 2, 5, 6], 0, false);
        assert_eq!(cigar_to_string(&cigar), "3M1I2D2M");
    }

    #[test]
    fn test_cigar_with_adapter_and_deletion() {
        let cigar = cigar_from_ref_positions(&[0, 1, 4, 5], 3, false);
        assert_eq!(cigar_to_string(&cigar), "2M2D2M3S");
    }

    #[test]
    fn test_cigar_single_base() {
        let cigar = cigar_from_ref_positions(&[42], 0, false);
        assert_eq!(cigar_to_string(&cigar), "1M");
    }

    #[test]
    fn test_cigar_high_positions() {
        // Negative-strand R2: positions start from a high offset (ascending),
        // no adapter — CIGAR is identical to forward strand.
        let cigar = cigar_from_ref_positions(&[100, 101, 102, 103, 104], 0, true);
        assert_eq!(cigar_to_string(&cigar), "5M");
    }

    // --- Read pair generation tests ---

    #[test]
    fn test_generate_pe_read_pair() {
        let fragment = test_fragment(b"ACGTACGTACGTACGTACGT", 100);
        let model = IlluminaErrorModel::new(10, 0.0, 0.0);
        let mut rng = SmallRng::seed_from_u64(42);

        let pair = generate_read_pair(
            &fragment, "chr1", 1, 10, true, b"ADAPTER", b"ADAPTER", &model, false, &mut rng,
        );

        assert_eq!(pair.read1.bases, b"ACGTACGTAC");
        assert!(pair.read2.is_some());
        assert_eq!(pair.r1_truth.position, 101);
        assert!(pair.r2_truth.is_some());
        assert_eq!(cigar_to_string(&pair.r1_cigar), "10M");
        assert_eq!(cigar_to_string(pair.r2_cigar.as_ref().unwrap()), "10M");
    }

    #[test]
    fn test_generate_se_read() {
        let fragment = test_fragment(b"ACGTACGTAC", 100);
        let model = IlluminaErrorModel::new(10, 0.0, 0.0);
        let mut rng = SmallRng::seed_from_u64(42);

        let pair = generate_read_pair(
            &fragment, "chr1", 5, 10, false, b"ADAPTER", b"ADAPTER", &model, false, &mut rng,
        );

        assert!(pair.read2.is_none());
        assert!(pair.r2_cigar.is_none());
        assert_eq!(cigar_to_string(&pair.r1_cigar), "10M");
    }

    #[test]
    fn test_adapter_cigar_softclip() {
        // Forward-strand fragment: R1 is positive-strand (trailing S), R2 is
        // negative-strand (leading S after BAM reverse-complement).
        let fragment = test_fragment(b"AC", 0);
        let model = IlluminaErrorModel::new(5, 0.0, 0.0);
        let mut rng = SmallRng::seed_from_u64(42);

        let pair = generate_read_pair(
            &fragment, "chr1", 1, 5, true, b"TTTTT", b"GGGGG", &model, false, &mut rng,
        );

        assert_eq!(cigar_to_string(&pair.r1_cigar), "2M3S");
        assert_eq!(cigar_to_string(pair.r2_cigar.as_ref().unwrap()), "3S2M");
    }

    #[test]
    fn test_simple_name_mode() {
        let fragment = test_fragment(b"ACGT", 0);
        let model = IlluminaErrorModel::new(4, 0.0, 0.0);
        let mut rng = SmallRng::seed_from_u64(42);

        let pair =
            generate_read_pair(&fragment, "chr1", 42, 4, true, b"A", b"A", &model, true, &mut rng);

        assert_eq!(pair.read1.name, "holodeck:42");
    }

    #[test]
    fn test_quality_scores_correct_length() {
        let fragment = test_fragment(b"ACGTACGTAC", 0);
        let model = IlluminaErrorModel::new(10, 0.001, 0.01);
        let mut rng = SmallRng::seed_from_u64(42);

        let pair =
            generate_read_pair(&fragment, "chr1", 1, 10, true, b"A", b"A", &model, false, &mut rng);

        assert_eq!(pair.read1.qualities.len(), 10);
        assert_eq!(pair.read2.as_ref().unwrap().qualities.len(), 10);
    }
}
