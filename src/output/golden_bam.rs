//! Golden BAM output for ground-truth alignments.
//!
//! Writes a BAM file containing the true alignment of each simulated read.
//! Each record has the correct reference position, strand, and CIGAR that
//! reflects haplotype variants (insertions, deletions) and adapter
//! soft-clipping. Custom tags encode haplotype and error counts.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::Path;

use anyhow::{Context, Result};
use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::data::field::Tag as DataTag;
use noodles::sam::alignment::record_buf::data::field::Value as DataValue;
use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::program::tag as pg_tag;
use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
use noodles::sam::header::record::value::map::{Program, ReadGroup, ReferenceSequence};

use crate::meth::{MethylationAnnotation, MethylationRecordTags};
use crate::read::ReadPair;
use crate::read_naming::TruthAlignment;
use crate::sequence_dict::SequenceDictionary;

/// Program ID used for the `@PG` header entry and the program name (`PN`).
const PROGRAM_ID: &str = "holodeck";

/// Sequencing platform string used for the `@RG` `PL` field.  Holodeck
/// currently only models an Illumina-style error profile.
const PLATFORM: &str = "ILLUMINA";

/// Read group identifier used for the single `@RG` entry and as the value of
/// the per-record `RG:Z` tag.
pub const READ_GROUP_ID: &str = "A";

/// Metadata used to populate the `@PG` and `@RG` header entries of a golden
/// BAM file.
#[derive(Debug, Clone)]
pub struct GoldenBamMetadata {
    /// Verbatim command-line invocation, written to the `@PG CL` field.
    pub command_line: String,
    /// Holodeck version string, written to the `@PG VN` field.
    pub version: String,
    /// Sample name, used for the `@RG SM` and `@RG LB` fields.
    pub sample: String,
}

/// A writer for golden (ground-truth) BAM files.
///
/// Records are written in the order they are generated (unsorted). Sort
/// downstream if a coordinate-sorted BAM is needed.  Supports both
/// single-threaded (noodles-bgzf) and multi-threaded (pooled-writer)
/// compression depending on how it is constructed.
pub struct GoldenBamWriter {
    /// Underlying BAM writer wrapping a boxed BGZF-producing writer.
    writer: bam::io::Writer<Box<dyn Write>>,
    /// SAM header for the BAM file.
    header: sam::Header,
}

impl GoldenBamWriter {
    /// Create a new single-threaded golden BAM writer at the given path.
    ///
    /// Uses noodles-bgzf with the specified compression level (0-12).
    ///
    /// # Errors
    /// Returns an error if the file cannot be created, the header cannot be
    /// built, or the header cannot be written.
    pub fn new(
        path: &Path,
        dict: &SequenceDictionary,
        compression: u8,
        meta: &GoldenBamMetadata,
    ) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create golden BAM: {}", path.display()))?;
        let level = noodles_bgzf::io::writer::CompressionLevel::new(compression)
            .ok_or_else(|| anyhow::anyhow!("invalid compression level: {compression}"))?;
        let bgzf_writer = noodles_bgzf::io::writer::Builder::default()
            .set_compression_level(level)
            .build_from_writer(BufWriter::new(file));
        Self::from_writer(Box::new(bgzf_writer), dict, meta)
    }

    /// Create a golden BAM writer from an existing writer that handles BGZF
    /// compression (e.g. a [`pooled_writer::PooledWriter`]).
    ///
    /// # Errors
    /// Returns an error if the header cannot be built or written.
    pub fn from_writer(
        writer: Box<dyn Write>,
        dict: &SequenceDictionary,
        meta: &GoldenBamMetadata,
    ) -> Result<Self> {
        let mut writer = bam::io::Writer::from(writer);
        let header = build_header(dict, meta)?;
        writer.write_header(&header)?;
        Ok(Self { writer, header })
    }

    /// Write a read pair to the golden BAM as one or two alignment records.
    ///
    /// For paired-end data, writes two records (R1 and R2) with proper paired
    /// flags, mate information, and CIGARs reflecting haplotype variants.
    /// For single-end data, writes one record.
    ///
    /// # Errors
    /// Returns an error if writing fails.
    pub fn write_pair(&mut self, pair: &ReadPair) -> Result<()> {
        let r1_truth = &pair.r1_truth;
        let r1_contig_idx = self.resolve_contig_index(&r1_truth.contig);

        let r1_tags = pair.methylation.as_ref().and_then(MethylationAnnotation::r1_tags);
        let r2_tags = pair.methylation.as_ref().and_then(MethylationAnnotation::r2_tags);

        if let (Some(r2_truth), Some(r2), Some(r2_cigar)) =
            (&pair.r2_truth, &pair.read2, &pair.r2_cigar)
        {
            let r2_contig_idx = self.resolve_contig_index(&r2_truth.contig);

            let r1_record = build_record(
                &pair.read1.name,
                &pair.read1.bases,
                &pair.read1.qualities,
                r1_truth,
                r1_contig_idx,
                pe_r1_flags(r1_truth, r2_truth),
                r2_contig_idx,
                Some(r2_truth.position),
                &pair.r1_cigar,
                r1_tags.as_ref(),
            );
            self.writer.write_alignment_record(&self.header, &r1_record)?;

            let r2_record = build_record(
                &pair.read1.name,
                &r2.bases,
                &r2.qualities,
                r2_truth,
                r2_contig_idx,
                pe_r2_flags(r1_truth, r2_truth),
                r1_contig_idx,
                Some(r1_truth.position),
                r2_cigar,
                r2_tags.as_ref(),
            );
            self.writer.write_alignment_record(&self.header, &r2_record)?;
        } else {
            let r1_record = build_record(
                &pair.read1.name,
                &pair.read1.bases,
                &pair.read1.qualities,
                r1_truth,
                r1_contig_idx,
                se_flags(r1_truth),
                None,
                None,
                &pair.r1_cigar,
                r1_tags.as_ref(),
            );
            self.writer.write_alignment_record(&self.header, &r1_record)?;
        }

        Ok(())
    }

    /// Finalize the BAM file.
    pub fn close(self) {
        drop(self.writer);
    }

    /// Resolve a contig name to its header index.
    fn resolve_contig_index(&self, contig: &str) -> Option<usize> {
        self.header.reference_sequences().get_index_of(contig.as_bytes())
    }
}

/// Build a SAM header from a sequence dictionary plus program/read-group
/// metadata.  Emits one `@SQ` line per contig, a single `@PG` entry describing
/// the holodeck invocation, and a single `@RG` entry tying every record in
/// the file to the simulated sample.
fn build_header(dict: &SequenceDictionary, meta: &GoldenBamMetadata) -> Result<sam::Header> {
    let mut builder = sam::Header::builder();
    for contig_meta in dict.iter() {
        let ref_seq = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(contig_meta.length()).expect("contig len > 0"),
        );
        builder = builder.add_reference_sequence(contig_meta.name().as_bytes(), ref_seq);
    }

    let program = Map::<Program>::builder()
        .insert(pg_tag::NAME, PROGRAM_ID)
        .insert(pg_tag::VERSION, meta.version.as_str())
        .insert(pg_tag::COMMAND_LINE, meta.command_line.as_str())
        .build()
        .map_err(|e| anyhow::anyhow!("failed to build @PG entry: {e}"))?;
    builder = builder.add_program(PROGRAM_ID, program);

    let read_group = Map::<ReadGroup>::builder()
        .insert(rg_tag::SAMPLE, meta.sample.as_str())
        .insert(rg_tag::LIBRARY, meta.sample.as_str())
        .insert(rg_tag::PLATFORM, PLATFORM)
        .build()
        .map_err(|e| anyhow::anyhow!("failed to build @RG entry: {e}"))?;
    builder = builder.add_read_group(READ_GROUP_ID, read_group);

    Ok(builder.build())
}

/// Build a single BAM `RecordBuf` from truth alignment data and a
/// pre-computed CIGAR.
///
/// BAM requires that all bases, qualities, and CIGARs are stored in the
/// orientation of the reference (left to right).  When the
/// `REVERSE_COMPLEMENTED` flag is set the supplied `bases` are in the 5'→3'
/// read orientation (i.e. the FASTQ sequence); this function reverse-
/// complements the sequence and reverses the quality scores before storing
/// them so that the record conforms to the BAM convention.
#[allow(clippy::too_many_arguments)]
fn build_record(
    name: &str,
    bases: &[u8],
    qualities: &[u8],
    truth: &TruthAlignment,
    contig_idx: Option<usize>,
    flags: Flags,
    mate_contig_idx: Option<usize>,
    mate_position: Option<u32>,
    cigar: &Cigar,
    methylation: Option<&MethylationRecordTags<'_>>,
) -> RecordBuf {
    // Normalise to reference orientation when the read is reverse-complemented.
    let is_rc = flags.is_reverse_complemented();
    let (seq, quals): (Vec<u8>, Vec<u8>) = if is_rc {
        let mut rc_bases = bases.to_vec();
        crate::fragment::reverse_complement(&mut rc_bases);
        let rev_quals: Vec<u8> = qualities.iter().rev().map(|&q| q.saturating_sub(33)).collect();
        (rc_bases, rev_quals)
    } else {
        let fwd_quals: Vec<u8> = qualities.iter().map(|&q| q.saturating_sub(33)).collect();
        (bases.to_vec(), fwd_quals)
    };

    let mut builder = RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_sequence(Sequence::from(seq.as_slice()))
        .set_quality_scores(QualityScores::from(quals))
        .set_cigar(cigar.clone());

    // Set reference and position.
    if let Some(idx) = contig_idx {
        builder = builder.set_reference_sequence_id(idx);
        if let Some(pos) = noodles::core::Position::new(truth.position as usize) {
            builder = builder
                .set_alignment_start(pos)
                .set_mapping_quality(sam::alignment::record::MappingQuality::new(60).unwrap());
        }
    }

    // Set mate info for paired-end reads.
    if let Some(idx) = mate_contig_idx {
        builder = builder.set_mate_reference_sequence_id(idx);
    }
    if let Some(mate_pos) = mate_position
        && let Some(pos) = noodles::core::Position::new(mate_pos as usize)
    {
        builder = builder.set_mate_alignment_start(pos);
    }

    // Add standard RG tag plus custom tags: hp (haplotype index), ne
    // (number of errors).
    let mut data = sam::alignment::record_buf::Data::default();
    data.insert(DataTag::READ_GROUP, DataValue::from(READ_GROUP_ID));
    #[expect(clippy::cast_possible_truncation, reason = "haplotype index is small")]
    #[expect(clippy::cast_possible_wrap, reason = "haplotype index is small")]
    let hp_val = truth.haplotype as i32;
    data.insert(DataTag::new(b'h', b'p'), DataValue::from(hp_val));
    #[expect(clippy::cast_possible_truncation, reason = "error count is small")]
    #[expect(clippy::cast_possible_wrap, reason = "error count is small")]
    let ne_val = truth.n_errors as i32;
    data.insert(DataTag::new(b'n', b'e'), DataValue::from(ne_val));

    if let Some(mt) = methylation {
        // XG:Z (Bismark) -- genome-strand indicator. `CT` for reads
        // derived from the top genome strand (OT/CTOT in Bismark's
        // four-strand model), `GA` for the bottom (OB/CTOB). Fragment-
        // level: R1 and R2 of a pair share the value. Identical for
        // em-seq and TAPS -- the chemistry's *biological* meaning of an
        // observed C→T differs between modes, but the genome-strand
        // indicator does not.
        let strand_tag = mt.conversion_type.as_tag_str();
        data.insert(DataTag::new(b'X', b'G'), DataValue::from(strand_tag));
        // cf:i (holodeck) -- ground-truth conversion-failure flag. `1` when
        // the source molecule was drawn as a whole-molecule conversion
        // failure, `0` otherwise. A molecule property, so R1 and R2 share
        // the value. Not recoverable from SEQ/YS alone (a molecule that
        // converted normally can coincidentally retain C's), so it is
        // stamped here as the truth an evaluator scores against.
        data.insert(DataTag::new(b'c', b'f'), DataValue::from(i32::from(mt.conversion_failed)));
        // XR:Z (Bismark) -- read-conversion direction in the read's
        // own orientation. Under a directional library the read shows
        // a C→T conversion pattern when read 5'→3' if it is R1 (or SE),
        // and a G→A pattern if it is R2 (which is the PCR-synthesized
        // complement of the source strand). So `XR:Z:CT` for R1/SE and
        // `XR:Z:GA` for R2, regardless of XG.
        let xr_tag = if flags.is_segmented() && flags.is_last_segment() { "GA" } else { "CT" };
        data.insert(DataTag::new(b'X', b'R'), DataValue::from(xr_tag));

        // XM:Z, YM:Z, NM:i, MD:Z (Bismark-style call tags) -- emitted
        // only when the simulator computed them (which it does when
        // both methylation chemistry and a golden BAM are requested).
        if let Some(call) = mt.call_tags {
            // XM and YM are stored in genomic orientation alongside SEQ.
            let xm_str = std::str::from_utf8(&call.xm).expect("XM is ASCII");
            data.insert(DataTag::new(b'X', b'M'), DataValue::from(xm_str));
            let ym_str = std::str::from_utf8(&call.ym).expect("YM is ASCII");
            data.insert(DataTag::new(b'Y', b'M'), DataValue::from(ym_str));
            #[expect(clippy::cast_possible_wrap, reason = "edit distance fits in i32")]
            let edit_distance = call.nm as i32;
            data.insert(DataTag::new(b'N', b'M'), DataValue::from(edit_distance));
            data.insert(DataTag::new(b'M', b'D'), DataValue::from(call.md.as_str()));
        }
        // YS:Z must be reference-oriented to match the SEQ field. The
        // captured pre-conversion bases are in read 5'->3' (FASTQ)
        // orientation; reverse-complement when the record is reverse.
        let ys_bytes: Vec<u8> = if is_rc {
            let mut tmp = mt.pre_conversion_bases.to_vec();
            crate::fragment::reverse_complement(&mut tmp);
            tmp
        } else {
            mt.pre_conversion_bases.to_vec()
        };
        // Pre-conversion bases are valid ASCII (uppercased ACGTN) so the
        // UTF-8 conversion is infallible, and DataValue::from(&str) yields
        // a String-typed BAM tag.
        let ys_str = std::str::from_utf8(&ys_bytes).expect("pre-conversion bases are valid ASCII");
        data.insert(DataTag::new(b'Y', b'S'), DataValue::from(ys_str));
    }

    builder = builder.set_data(data);

    builder.build()
}

/// Compute SAM flags for R1 of a paired-end read.
fn pe_r1_flags(r1: &TruthAlignment, r2: &TruthAlignment) -> Flags {
    let mut flags = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
    if !r1.is_forward {
        flags |= Flags::REVERSE_COMPLEMENTED;
    }
    if !r2.is_forward {
        flags |= Flags::MATE_REVERSE_COMPLEMENTED;
    }
    if r1.contig == r2.contig {
        flags |= Flags::PROPERLY_SEGMENTED;
    }
    flags
}

/// Compute SAM flags for R2 of a paired-end read.
fn pe_r2_flags(r1: &TruthAlignment, r2: &TruthAlignment) -> Flags {
    let mut flags = Flags::SEGMENTED | Flags::LAST_SEGMENT;
    if !r2.is_forward {
        flags |= Flags::REVERSE_COMPLEMENTED;
    }
    if !r1.is_forward {
        flags |= Flags::MATE_REVERSE_COMPLEMENTED;
    }
    if r1.contig == r2.contig {
        flags |= Flags::PROPERLY_SEGMENTED;
    }
    flags
}

/// Compute SAM flags for a single-end read.
fn se_flags(truth: &TruthAlignment) -> Flags {
    let mut flags = Flags::empty();
    if !truth.is_forward {
        flags |= Flags::REVERSE_COMPLEMENTED;
    }
    flags
}
