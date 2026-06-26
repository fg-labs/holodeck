//! Golden-BAM truth: per-read true alignment loaded for placement and for
//! variant / MD-tag concordance scoring.
//!
//! The golden BAM written by `holodeck simulate --golden-bam` carries the
//! true alignment of every read (MAPQ 60, correct CIGAR, the sequence the read
//! was given) and — for methylation runs — Bismark-style `NM:i` / `MD:Z` call
//! tags. This module indexes those records by read end so the eval pass can
//! look up each mapped read's truth in O(1), including the actual allele the
//! read carries at any reference position (see [`GoldenInfo::base_at`]).

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};
use bstr::ByteSlice;
use noodles::bam;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::Cigar;
use noodles::sam::alignment::record_buf::data::field::Value;

use super::cigar;

/// Key identifying one read end: read name plus whether it is the last segment
/// (R2). R1 and single-end reads use `false`.
pub type ReadKey = (Vec<u8>, bool);

/// True alignment for one read end, taken from the golden BAM.
#[derive(Debug, Clone)]
pub struct GoldenInfo {
    /// True reference contig name.
    pub contig: String,
    /// 0-based true start position.
    pub start0: u32,
    /// Reference bases consumed by the true alignment.
    pub ref_len: u32,
    /// `NM:i` edit distance, if present.
    pub nm: Option<i64>,
    /// `MD:Z` string, if present.
    pub md: Option<String>,
    /// Uppercased read sequence of the true alignment. Paired with `cigar` and
    /// `start0`, this is the per-read oracle for which allele a read actually
    /// carries at a variant site — making variant-representation scoring
    /// independent of whether the truth VCF is phased.
    pub sequence: Vec<u8>,
    /// CIGAR of the true alignment, for mapping a reference position to a read
    /// offset within `sequence`.
    pub cigar: Cigar,
}

impl GoldenInfo {
    /// 0-based exclusive end of the true alignment.
    #[must_use]
    pub fn end0(&self) -> u32 {
        self.start0 + self.ref_len
    }

    /// The uppercased base this read carries at 0-based reference position
    /// `ref_pos0`, or `None` when that position is deleted, clipped, or outside
    /// the alignment. This is read straight from the golden sequence, so it
    /// reflects exactly what the simulator placed on this read (the alt allele
    /// for a read sequenced from the alt copy, the reference base otherwise).
    #[must_use]
    pub fn base_at(&self, ref_pos0: u32) -> Option<u8> {
        let offset = cigar::ref_pos_to_read_offset(&self.cigar, self.start0, ref_pos0)?;
        self.sequence.get(offset).map(u8::to_ascii_uppercase)
    }
}

/// Load every primary, mapped golden record keyed by `(name, is_last_segment)`.
///
/// # Errors
/// Returns an error if the golden BAM cannot be read.
pub fn load(path: &Path) -> Result<HashMap<ReadKey, GoldenInfo>> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(path)
        .with_context(|| format!("Failed to open golden BAM: {}", path.display()))?;
    let header = reader.read_header()?;

    let mut map = HashMap::new();
    for result in reader.record_bufs(&header) {
        let record = result.context("Failed to read golden BAM record")?;
        let flags = record.flags();
        if flags.is_secondary() || flags.is_supplementary() || flags.is_unmapped() {
            continue;
        }

        let Some(name) = record.name() else { continue };
        let Some(ref_id) = record.reference_sequence_id() else { continue };
        let Some((contig_name, _)) = header.reference_sequences().get_index(ref_id) else {
            continue;
        };
        let Some(start) = record.alignment_start() else { continue };

        let start0 = u32::try_from(usize::from(start).saturating_sub(1)).unwrap_or(0);
        let cigar = record.cigar().clone();
        let info = GoldenInfo {
            contig: contig_name.to_str_lossy().into_owned(),
            start0,
            ref_len: cigar::reference_len(&cigar),
            nm: int_tag(&record, b'N', b'M'),
            md: string_tag(&record, b'M', b'D'),
            sequence: record.sequence().as_ref().to_vec(),
            cigar,
        };
        map.insert((name.to_vec(), flags.is_last_segment()), info);
    }

    Ok(map)
}

/// Resolve a record's reference contig name via the header.
pub(super) fn contig_name(header: &noodles::sam::Header, ref_id: usize) -> Option<String> {
    let (name, _) = header.reference_sequences().get_index(ref_id)?;
    Some(name.to_str_lossy().into_owned())
}

/// Read an integer auxiliary tag from a record buffer.
pub(super) fn int_tag(record: &noodles::sam::alignment::RecordBuf, a: u8, b: u8) -> Option<i64> {
    record.data().get(&Tag::new(a, b)).and_then(Value::as_int)
}

/// Read a string (`Z`) auxiliary tag from a record buffer.
pub(super) fn string_tag(
    record: &noodles::sam::alignment::RecordBuf,
    a: u8,
    b: u8,
) -> Option<String> {
    match record.data().get(&Tag::new(a, b)) {
        Some(Value::String(s)) => Some(s.to_str_lossy().into_owned()),
        _ => None,
    }
}
