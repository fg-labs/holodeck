//! Shared test helpers for integration tests.
//!
//! Provides builders for creating temporary FASTA, VCF, and BED files
//! programmatically, avoiding any committed test data files.

// Each integration test binary compiles this module independently, so not all
// functions are used in every binary.
#![allow(dead_code)]

use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

use std::num::NonZeroUsize;
use std::process::Command;

use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::header::record::value::map::header::{sort_order, tag};
use noodles::sam::header::record::value::{Map, map};
use noodles::vcf;
use tempfile::TempDir;

/// A VCF variant to be written by [`TestEnv::write_vcf`].
#[allow(dead_code)]
pub struct VcfVariant<'a> {
    /// Contig/chromosome name.
    pub chrom: &'a str,
    /// 1-based VCF position.
    pub pos_1based: u32,
    /// Reference allele bases.
    pub ref_allele: &'a str,
    /// One or more alternate allele strings.
    pub alt_alleles: &'a [&'a str],
    /// Genotype string e.g. "1|1", "0|1", "0/1".
    pub gt: &'a str,
}

/// A test environment with a temporary directory and a reference FASTA.
pub struct TestEnv {
    /// The temporary directory (cleaned up on drop).
    pub dir: TempDir,
    /// Path to the reference FASTA file.
    pub fasta_path: PathBuf,
}

impl TestEnv {
    /// Create a new test environment with a FASTA file containing the given
    /// contigs.
    ///
    /// Each entry is `(contig_name, sequence_bytes)`.
    ///
    /// # Panics
    /// Panics on I/O errors.
    #[must_use]
    pub fn new(contigs: &[(&str, &[u8])]) -> Self {
        let dir = TempDir::new().unwrap();
        let fasta_path = write_fasta(dir.path(), contigs);
        Self { dir, fasta_path }
    }

    /// Return the output prefix path for writing test output files.
    #[must_use]
    pub fn output_prefix(&self) -> PathBuf {
        self.dir.path().join("output")
    }

    /// Write a VCF file with the given contigs, sample, and variants.
    ///
    /// Returns the path to the written `.vcf` file. Uses the noodles builder
    /// API throughout -- no raw VCF string assembly.
    ///
    /// # Panics
    /// Panics on I/O errors or if a genotype string cannot be parsed.
    #[allow(dead_code)]
    #[must_use]
    pub fn write_vcf(
        &self,
        sample_name: &str,
        contig_lengths: &[(&str, usize)],
        variants: &[VcfVariant<'_>],
    ) -> PathBuf {
        use noodles::core::Position;
        use vcf::header::record::value::{
            Map as VcfMap,
            map::{Contig, Format},
        };
        use vcf::variant::io::Write as VcfWrite;
        use vcf::variant::record::samples::keys::key;
        use vcf::variant::record_buf::{
            AlternateBases, Filters as VcfFilters, Samples,
            samples::{Keys, sample::Value as SampleValue},
        };

        let path = self.dir.path().join("variants.vcf");

        // Build contig map entries with lengths.
        let mut header_builder = vcf::Header::builder();
        for &(name, len) in contig_lengths {
            let contig = VcfMap::<Contig>::builder().set_length(len).build().unwrap();
            header_builder = header_builder.add_contig(name, contig);
        }

        // Add FORMAT=GT definition and the sample column.
        let gt_format = VcfMap::<Format>::from(key::GENOTYPE);
        header_builder =
            header_builder.add_format(key::GENOTYPE, gt_format).add_sample_name(sample_name);

        let header = header_builder.build();

        let file = std::fs::File::create(&path).unwrap();
        let mut writer = vcf::io::Writer::new(std::io::BufWriter::new(file));
        writer.write_header(&header).unwrap();

        // Write each variant record.
        let gt_keys: Keys = [String::from(key::GENOTYPE)].into_iter().collect();

        for v in variants {
            let alt = AlternateBases::from(
                v.alt_alleles.iter().map(|a| String::from(*a)).collect::<Vec<_>>(),
            );

            let genotype: vcf::variant::record_buf::samples::sample::value::Genotype =
                v.gt.parse().unwrap();

            let samples =
                Samples::new(gt_keys.clone(), vec![vec![Some(SampleValue::Genotype(genotype))]]);

            let record = vcf::variant::RecordBuf::builder()
                .set_reference_sequence_name(v.chrom)
                .set_variant_start(Position::try_from(v.pos_1based as usize).unwrap())
                .set_reference_bases(v.ref_allele)
                .set_alternate_bases(alt)
                .set_filters(VcfFilters::pass())
                .set_samples(samples)
                .build();

            writer.write_variant_record(&header, &record).unwrap();
        }

        path
    }

    /// Write a VCF file with the given sample names and contig definitions but
    /// no variant records.
    ///
    /// Useful for testing sample-level validation logic without needing actual
    /// variants.  Uses raw VCF text to support arbitrary sample configurations
    /// (single-sample, multi-sample, zero-sample).
    ///
    /// # Panics
    /// Panics on I/O errors.
    #[must_use]
    pub fn write_vcf_header_only(
        &self,
        sample_names: &[&str],
        contig_lengths: &[(&str, usize)],
    ) -> PathBuf {
        let path = self.dir.path().join("samples.vcf");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "##fileformat=VCFv4.3").unwrap();
        for &(name, len) in contig_lengths {
            writeln!(f, "##contig=<ID={name},length={len}>").unwrap();
        }
        writeln!(f, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">").unwrap();
        write!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT").unwrap();
        for name in sample_names {
            write!(f, "\t{name}").unwrap();
        }
        writeln!(f).unwrap();
        f.flush().unwrap();
        path
    }

    /// Write a BED file with the given intervals and return the path.
    ///
    /// Each entry is `(contig, start, end)` in 0-based half-open coordinates.
    ///
    /// # Panics
    /// Panics on I/O errors.
    #[must_use]
    pub fn write_bed(&self, intervals: &[(&str, u32, u32)]) -> PathBuf {
        let path = self.dir.path().join("targets.bed");
        let mut f = std::fs::File::create(&path).unwrap();
        for (contig, start, end) in intervals {
            writeln!(f, "{contig}\t{start}\t{end}").unwrap();
        }
        f.flush().unwrap();
        path
    }
}

/// Write a FASTA file and its `.fai` index from programmatic data.
///
/// # Panics
/// Panics on I/O errors.
fn write_fasta(dir: &Path, contigs: &[(&str, &[u8])]) -> PathBuf {
    let fasta_path = dir.join("ref.fa");
    let fai_path = dir.join("ref.fa.fai");

    let mut fa = std::fs::File::create(&fasta_path).unwrap();
    let mut fai = std::fs::File::create(&fai_path).unwrap();

    let line_len: usize = 80;
    let mut offset: u64 = 0;

    for &(name, seq) in contigs {
        let header = format!(">{name}\n");
        fa.write_all(header.as_bytes()).unwrap();
        offset += header.len() as u64;

        let seq_offset = offset;
        for chunk in seq.chunks(line_len) {
            fa.write_all(chunk).unwrap();
            fa.write_all(b"\n").unwrap();
        }

        let seq_len = seq.len();
        let line_width = line_len + 1;
        writeln!(fai, "{name}\t{seq_len}\t{seq_offset}\t{line_len}\t{line_width}").unwrap();

        let full_lines = seq_len / line_len;
        let remainder = seq_len % line_len;
        offset += (full_lines * line_width) as u64;
        if remainder > 0 {
            offset += (remainder + 1) as u64;
        }
    }

    fasta_path
}

/// Generate a deterministic non-repetitive DNA sequence of the given length.
///
/// Uses a simple LCG to avoid repetitive patterns that would confuse pileup
/// reasoning. The same seed always produces the same sequence.
pub fn non_repetitive_seq(len: usize) -> Vec<u8> {
    let bases = b"ACGT";
    let mut x: u64 = 0xdead_beef_cafe_1234;
    (0..len)
        .map(|_| {
            x = x.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
            bases[((x >> 33) & 3) as usize]
        })
        .collect()
}

/// Per-position summary produced by `pileup_bases`.
#[derive(Debug, Default, Clone)]
#[allow(dead_code)]
pub struct PileupColumn {
    /// 0-based reference position.
    pub pos: u32,
    pub a: u32,
    pub c: u32,
    pub g: u32,
    pub t: u32,
    pub n: u32,
    /// Number of reads with an insertion immediately after this position.
    pub ins: u32,
    /// Counts of each distinct inserted sequence at this anchor position.
    pub ins_seqs: HashMap<Vec<u8>, u32>,
    /// Total reference-consuming depth = a + c + g + t + n.
    pub total: u32,
}

/// Read the given BAM, sort it by coordinate, write a sorted copy with a BAI
/// index, then return a per-position pileup for every contig.
///
/// Returns `Vec<Vec<PileupColumn>>` where the outer `Vec` is indexed by contig
/// order in the BAM header and the inner `Vec` is indexed by 0-based reference
/// position within that contig.  Positions with no coverage have all-zero
/// counts.
///
/// Uses noodles for all BAM reading, sorting, writing, indexing, and querying.
/// The per-base accumulation is a straightforward CIGAR walk over the decoded
/// `RecordBuf` values returned by the indexed query.
///
/// # Panics
/// Panics on any I/O or indexing error — appropriate for test helpers.
pub fn pileup_bases(bam_path: &Path) -> Vec<Vec<PileupColumn>> {
    // ── 1. Read header and all records ──────────────────────────────────────
    let mut reader = bam::io::reader::Builder.build_from_path(bam_path).unwrap();
    let header = reader.read_header().unwrap();

    let mut records: Vec<RecordBuf> = reader.record_bufs(&header).map(|r| r.unwrap()).collect();

    // ── 2. Sort by (reference_sequence_id, alignment_start) ─────────────────
    records.sort_by_key(|r| {
        let ref_id = r.reference_sequence_id().unwrap_or(usize::MAX);
        let start = r.alignment_start().map_or(0usize, usize::from);
        (ref_id, start)
    });

    // ── 3. Write sorted BAM — bam::fs::index requires SO:coordinate ─────────
    let sorted_bam_path = {
        let mut p = bam_path.as_os_str().to_owned();
        p.push(".sorted.bam");
        PathBuf::from(p)
    };
    {
        let mut sorted_header = header.clone();
        let hd = Map::<map::Header>::builder()
            .insert(tag::SORT_ORDER, sort_order::COORDINATE)
            .build()
            .expect("build HD header map");
        *sorted_header.header_mut() = Some(hd);

        let file = std::fs::File::create(&sorted_bam_path).unwrap();
        let mut writer = bam::io::Writer::new(std::io::BufWriter::new(file));
        writer.write_header(&sorted_header).unwrap();
        for record in &records {
            writer.write_alignment_record(&sorted_header, record).unwrap();
        }
        writer.try_finish().unwrap();
    }

    // ── 4. Build and write BAI ───────────────────────────────────────────────
    let index = bam::fs::index(&sorted_bam_path).unwrap();
    let bai_path = {
        let mut p = sorted_bam_path.as_os_str().to_owned();
        p.push(".bai");
        PathBuf::from(p)
    };
    bam::bai::fs::write(&bai_path, &index).unwrap();

    // ── 5. Allocate pileup result sized by contig lengths ────────────────────
    let contig_count = header.reference_sequences().len();
    let contig_lengths: Vec<usize> =
        header.reference_sequences().iter().map(|(_, rs)| usize::from(rs.length())).collect();

    let mut pileup: Vec<Vec<PileupColumn>> = contig_lengths
        .iter()
        .map(|&len| {
            #[expect(clippy::cast_possible_truncation, reason = "contig positions fit in u32")]
            (0..len).map(|pos| PileupColumn { pos: pos as u32, ..Default::default() }).collect()
        })
        .collect();

    // ── 6. Open indexed reader and accumulate pileup per contig ─────────────
    let mut indexed_reader =
        bam::io::indexed_reader::Builder::default().build_from_path(&sorted_bam_path).unwrap();
    let sorted_header = indexed_reader.read_header().unwrap();

    for contig_idx in 0..contig_count {
        let (name_bytes, _) = sorted_header.reference_sequences().get_index(contig_idx).unwrap();
        let contig_name = std::str::from_utf8(name_bytes).unwrap();
        let region: noodles::core::Region = contig_name.parse().unwrap();

        let query = indexed_reader.query(&sorted_header, &region).unwrap();
        for result in query.records() {
            let lazy_record = result.unwrap();
            let record =
                RecordBuf::try_from_alignment_record(&sorted_header, &lazy_record).unwrap();
            accumulate_record(&record, contig_idx, &mut pileup);
        }
    }

    pileup
}

// ── Command runners ─────────────────────────────────────────────────────────

/// Run `holodeck simulate` with the given arguments and return
/// `(success, stdout, stderr)`.
pub fn run_simulate(args: &[&str]) -> (bool, String, String) {
    run_holodeck(args)
}

/// Run `holodeck mutate` with the given arguments and return
/// `(success, stdout, stderr)`.
pub fn run_mutate(args: &[&str]) -> (bool, String, String) {
    run_holodeck(args)
}

/// Run `holodeck eval` with the given arguments and return
/// `(success, stdout, stderr)`.
pub fn run_eval(args: &[&str]) -> (bool, String, String) {
    run_holodeck(args)
}

/// Run the holodeck binary with arbitrary arguments.
fn run_holodeck(args: &[&str]) -> (bool, String, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_holodeck"))
        .args(args)
        .output()
        .expect("Failed to run holodeck");
    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();
    (output.status.success(), stdout, stderr)
}

// ── VCF reading ─────────────────────────────────────────────────────────────

/// A parsed VCF variant record for test assertions.
#[derive(Debug)]
#[allow(dead_code)]
pub struct ParsedVcfRecord {
    /// Contig name.
    pub chrom: String,
    /// 1-based position.
    pub pos: usize,
    /// Reference allele.
    pub ref_allele: String,
    /// Alternate allele strings.
    pub alt_alleles: Vec<String>,
    /// Genotype string (e.g. "0/1", "1|1").
    pub gt: String,
}

/// Format a noodles Genotype value back to a VCF GT string (e.g. "0/1").
fn format_gt(gt: &vcf::variant::record_buf::samples::sample::value::Genotype) -> String {
    use vcf::variant::record::samples::series::value::genotype::Phasing;

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

/// Read all variant records from a VCF file, returning parsed records.
///
/// Extracts chrom, pos, ref, alts, and the GT for the first sample.
///
/// # Panics
/// Panics on I/O or parse errors — appropriate for tests.
#[allow(dead_code)]
pub fn read_vcf_records(path: &Path) -> Vec<ParsedVcfRecord> {
    use vcf::variant::record::samples::keys::key;
    use vcf::variant::record_buf::samples::sample::Value as SampleValue;

    let mut reader = vcf::io::reader::Builder::default().build_from_path(path).unwrap();
    let header = reader.read_header().unwrap();

    let mut records = Vec::new();
    for result in reader.record_bufs(&header) {
        let record: vcf::variant::RecordBuf = result.unwrap();

        let chrom = record.reference_sequence_name().to_string();
        let pos = record.variant_start().unwrap().get();
        let ref_allele = record.reference_bases().to_string();
        let alt_alleles: Vec<String> = record.alternate_bases().as_ref().to_vec();

        // Extract GT from first sample.
        let samples = record.samples();
        let sample = samples.get_index(0).unwrap();
        let gt_value = sample.get(key::GENOTYPE).unwrap();
        let gt = match gt_value {
            Some(SampleValue::Genotype(g)) => format_gt(g),
            _ => panic!("Expected Genotype value"),
        };

        records.push(ParsedVcfRecord { chrom, pos, ref_allele, alt_alleles, gt });
    }
    records
}

// ── BAM writing ─────────────────────────────────────────────────────────────

/// Write a minimal BAM file with the given contig info and records.
///
/// Each record is specified as `(name, contig_idx, pos_1based, flags, mapq)`.
/// Sequence and qualities are set to a single 'A' base with Q30.
///
/// # Panics
/// Panics on I/O errors — appropriate for tests.
/// A single BAM record specification for [`write_bam`].
pub type BamRecordSpec<'a> = (&'a str, Option<usize>, Option<u32>, Flags, u8);

pub fn write_bam(path: &Path, contigs: &[(&str, usize)], records: &[BamRecordSpec<'_>]) {
    let mut header_builder = sam::Header::builder();
    for &(name, len) in contigs {
        let ref_seq = Map::<map::ReferenceSequence>::new(NonZeroUsize::new(len).unwrap());
        header_builder = header_builder.add_reference_sequence(name.as_bytes(), ref_seq);
    }
    let header = header_builder.build();

    let file = std::fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(std::io::BufWriter::new(file));
    writer.write_header(&header).unwrap();

    for &(name, contig_idx, pos, flags, mapq) in records {
        let mut builder = RecordBuf::builder()
            .set_name(name)
            .set_flags(flags)
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"A".as_ref()))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![30]))
            .set_cigar(
                [sam::alignment::record::cigar::Op::new(Kind::Match, 1)].into_iter().collect(),
            );

        if let Some(idx) = contig_idx {
            builder = builder.set_reference_sequence_id(idx);
        }
        if let Some(p) = pos {
            if let Some(noodles_pos) = noodles::core::Position::new(p as usize) {
                builder = builder.set_alignment_start(noodles_pos);
            }
            if mapq > 0 {
                builder = builder.set_mapping_quality(
                    sam::alignment::record::MappingQuality::new(mapq).unwrap(),
                );
            }
        }

        let record = builder.build();
        writer.write_alignment_record(&header, &record).unwrap();
    }

    writer.try_finish().unwrap();
}

// ── Misc helpers ────────────────────────────────────────────────────────────

/// Read a gzipped file and return its contents as a string.
pub fn read_gzipped(path: &Path) -> String {
    use std::io::Read;
    let file = std::fs::File::open(path).unwrap();
    let mut decoder = flate2::read::MultiGzDecoder::new(file);
    let mut contents = String::new();
    decoder.read_to_string(&mut contents).unwrap();
    contents
}

/// Count the number of FASTQ records in a string (4 lines per record).
pub fn count_fastq_records(contents: &str) -> usize {
    let lines: Vec<&str> = contents.lines().collect();
    assert_eq!(lines.len() % 4, 0, "FASTQ file has non-multiple-of-4 lines");
    lines.len() / 4
}

/// Build a simple 1000bp reference for testing.
pub fn simple_env() -> TestEnv {
    let seq = b"ACGTACGTAC".repeat(100); // 1000 bp
    TestEnv::new(&[("chr1", &seq)])
}

/// Walk a single record's CIGAR and accumulate per-base counts into `pileup`.
fn accumulate_record(record: &RecordBuf, contig_idx: usize, pileup: &mut [Vec<PileupColumn>]) {
    let Some(alignment_start) = record.alignment_start() else { return };
    let mut ref_pos = usize::from(alignment_start) - 1; // convert 1-based → 0-based

    let seq = record.sequence().as_ref(); // &[u8] of ASCII bases
    let mut read_pos = 0usize;
    let contig_len = pileup[contig_idx].len();

    for op in record.cigar().as_ref() {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in 0..len {
                    let rp = ref_pos + i;
                    let rd = read_pos + i;
                    if rp >= contig_len {
                        break;
                    }
                    if rd >= seq.len() {
                        break;
                    }
                    let base = seq[rd];
                    let col = &mut pileup[contig_idx][rp];
                    match base {
                        b'A' => col.a += 1,
                        b'C' => col.c += 1,
                        b'G' => col.g += 1,
                        b'T' => col.t += 1,
                        _ => col.n += 1,
                    }
                    col.total += 1;
                }
                ref_pos += len;
                read_pos += len;
            }
            Kind::Insertion => {
                // Insertion is anchored at the last aligned ref base (ref_pos - 1).
                if ref_pos > 0 {
                    let anchor = ref_pos - 1;
                    if anchor < contig_len && read_pos + len <= seq.len() {
                        let ins_seq = seq[read_pos..read_pos + len].to_vec();
                        let col = &mut pileup[contig_idx][anchor];
                        col.ins += 1;
                        *col.ins_seqs.entry(ins_seq).or_insert(0) += 1;
                    }
                }
                read_pos += len;
            }
            Kind::Deletion | Kind::Skip => {
                ref_pos += len;
            }
            Kind::SoftClip => {
                read_pos += len;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
}
