use std::fs::File;
use std::io::{Seek, SeekFrom};
use std::path::Path;

use anyhow::{Context, Result, anyhow};
use noodles::core::Region;
use noodles::fasta;

use crate::sequence_dict::SequenceDictionary;

/// A loaded FASTA reference with random-access sequence retrieval.
///
/// Wraps a noodles `IndexedReader` for efficient seek-based access.  Contig
/// metadata (names, lengths) is pre-computed from the `.fai` index at
/// construction time.
pub struct Fasta {
    /// Noodles indexed FASTA reader.
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    /// Sequence dictionary built from the `.fai` index.
    dict: SequenceDictionary,
}

impl Fasta {
    /// Open a FASTA file and its `.fai` index.
    ///
    /// The `.fai` file must exist alongside the FASTA (e.g. `ref.fa.fai`).
    ///
    /// # Errors
    /// Returns an error if the FASTA or its index cannot be read.
    pub fn from_path(path: &Path) -> Result<Self> {
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .with_context(|| format!("Failed to open indexed FASTA: {}", path.display()))?;

        let dict = SequenceDictionary::from(reader.index());

        Ok(Self { reader, dict })
    }

    /// Return a reference to the underlying sequence dictionary.
    #[must_use]
    pub fn dict(&self) -> &SequenceDictionary {
        &self.dict
    }

    /// Load and return the full sequence of `contig_name` as uppercase bytes.
    ///
    /// Uses seek-based access via the FAI index for efficiency, avoiding the
    /// intermediate `Record` allocation that `query()` performs.
    ///
    /// # Errors
    /// Returns an error if the contig is unknown or the FASTA cannot be read.
    pub fn load_contig(&mut self, contig_name: &str) -> Result<Vec<u8>> {
        let meta = self
            .dict
            .get_by_name(contig_name)
            .ok_or_else(|| anyhow!("Contig '{contig_name}' not found in reference"))?;
        let expected_len = meta.length();

        let region = Region::new(contig_name, ..);
        let pos =
            self.reader.index().query(&region).with_context(|| {
                format!("Failed to look up contig '{contig_name}' in FASTA index")
            })?;
        self.reader
            .get_mut()
            .seek(SeekFrom::Start(pos))
            .with_context(|| format!("Failed to seek to contig '{contig_name}' in FASTA"))?;

        let mut seq = Vec::with_capacity(expected_len);
        self.reader
            .read_sequence(&mut seq)
            .with_context(|| format!("Failed to read contig '{contig_name}' from FASTA"))?;

        // Uppercase all bases for consistent matching.
        for b in &mut seq {
            *b = b.to_ascii_uppercase();
        }
        Ok(seq)
    }

    /// Return contig names in index order.
    #[must_use]
    pub fn contig_names(&self) -> Vec<&str> {
        self.dict.names()
    }
}

/// Write a FASTA file and its `.fai` index from programmatic data.
///
/// Each entry in `contigs` is `(name, sequence)`. Returns the path to the
/// FASTA file.
///
/// # Panics
/// Panics on any I/O error. This is a test-only helper.
#[cfg(test)]
#[must_use]
pub fn write_test_fasta(dir: &std::path::Path, contigs: &[(&str, &[u8])]) -> std::path::PathBuf {
    use std::io::Write;

    let fasta_path = dir.join("ref.fa");
    let fai_path = dir.join("ref.fa.fai");

    let mut fa = std::fs::File::create(&fasta_path).unwrap();
    let mut fai = std::fs::File::create(&fai_path).unwrap();

    let mut offset: u64 = 0;
    for &(name, seq) in contigs {
        // Write FASTA: >name\nsequence\n
        let header = format!(">{name}\n");
        fa.write_all(header.as_bytes()).unwrap();
        offset += header.len() as u64;

        let seq_offset = offset;
        // Write sequence in 80-char lines.
        let line_len = 80;
        for chunk in seq.chunks(line_len) {
            fa.write_all(chunk).unwrap();
            fa.write_all(b"\n").unwrap();
        }

        let seq_len = seq.len();
        let line_bases = line_len;
        let line_width = line_len + 1; // bases + newline
        // FAI format: name\tlength\toffset\tline_bases\tline_width
        writeln!(fai, "{name}\t{seq_len}\t{seq_offset}\t{line_bases}\t{line_width}").unwrap();

        // Update offset: seq bytes + newlines
        let full_lines = seq_len / line_len;
        let remainder = seq_len % line_len;
        offset += (full_lines * line_width) as u64;
        if remainder > 0 {
            offset += (remainder + 1) as u64;
        }
    }

    fasta_path
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_contig_single() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_test_fasta(dir.path(), &[("chr1", b"ACGTACGT")]);
        let mut fasta = Fasta::from_path(&path).unwrap();

        assert_eq!(fasta.contig_names(), vec!["chr1"]);
        assert_eq!(fasta.dict().len(), 1);

        let seq = fasta.load_contig("chr1").unwrap();
        assert_eq!(&seq, b"ACGTACGT");
    }

    #[test]
    fn test_load_contig_multiple() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_test_fasta(
            dir.path(),
            &[("chr1", b"AAAA"), ("chr2", b"CCCCGGGG"), ("chr3", b"TT")],
        );
        let mut fasta = Fasta::from_path(&path).unwrap();

        assert_eq!(fasta.contig_names(), vec!["chr1", "chr2", "chr3"]);

        let seq1 = fasta.load_contig("chr1").unwrap();
        assert_eq!(&seq1, b"AAAA");

        let seq2 = fasta.load_contig("chr2").unwrap();
        assert_eq!(&seq2, b"CCCCGGGG");

        let seq3 = fasta.load_contig("chr3").unwrap();
        assert_eq!(&seq3, b"TT");
    }

    #[test]
    fn test_load_contig_out_of_order() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_test_fasta(dir.path(), &[("chr1", b"AAAA"), ("chr2", b"CCCC")]);
        let mut fasta = Fasta::from_path(&path).unwrap();

        // Load chr2 first, then chr1 — verify seeking works correctly.
        let seq2 = fasta.load_contig("chr2").unwrap();
        assert_eq!(&seq2, b"CCCC");

        let seq1 = fasta.load_contig("chr1").unwrap();
        assert_eq!(&seq1, b"AAAA");
    }

    #[test]
    fn test_load_contig_uppercase() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_test_fasta(dir.path(), &[("chr1", b"acgtNn")]);
        let mut fasta = Fasta::from_path(&path).unwrap();

        let seq = fasta.load_contig("chr1").unwrap();
        assert_eq!(&seq, b"ACGTNN");
    }

    #[test]
    fn test_load_contig_unknown() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_test_fasta(dir.path(), &[("chr1", b"ACGT")]);
        let mut fasta = Fasta::from_path(&path).unwrap();

        let result = fasta.load_contig("chrZ");
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("not found"));
    }
}
