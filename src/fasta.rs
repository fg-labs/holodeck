use std::fs::File;
use std::io::{Seek, SeekFrom};
use std::path::Path;

use anyhow::{Context, Result, anyhow};
use noodles::core::Region;
use noodles::fasta;
use rand::Rng;

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

    /// Load and return the full sequence of `contig_name`.
    ///
    /// Uses seek-based access via the FAI index for efficiency, avoiding the
    /// intermediate `Record` allocation that `query()` performs.
    ///
    /// Bases are normalized during loading:
    /// - `A`, `C`, `G`, `T` (case-insensitive) are stored as uppercase.
    /// - `U` (RNA) is converted to `T`.
    /// - Any other IUPAC ambiguity code (`R`, `Y`, `S`, `W`, `K`, `M`, `B`,
    ///   `D`, `H`, `V`, `N`) is resolved to a uniformly-random base drawn
    ///   from that code's ambiguity set and stored in **lowercase**. The
    ///   lowercase marker allows downstream sampling to detect and optionally
    ///   reject reads that contain too many synthesized bases (see
    ///   [`crate::fragment::lowercase_fraction`]).
    /// - Any other byte is rejected with an error.
    ///
    /// Passing the same `rng` state produces a deterministic resolution, so
    /// two runs seeded identically produce identical reference buffers.
    ///
    /// # Errors
    /// Returns an error if the contig is unknown, the FASTA cannot be read,
    /// or the contig contains an unrecognized base.
    pub fn load_contig(&mut self, contig_name: &str, rng: &mut impl Rng) -> Result<Vec<u8>> {
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

        normalize_bases(&mut seq, contig_name, rng)?;
        Ok(seq)
    }

    /// Return contig names in index order.
    #[must_use]
    pub fn contig_names(&self) -> Vec<&str> {
        self.dict.names()
    }
}

/// Resolve a non-ACGT IUPAC code to one uniformly-chosen valid base.
///
/// The input byte must already be uppercase. Returns `None` if the input
/// is not a recognized IUPAC code (including plain `A`/`C`/`G`/`T`, which
/// the caller should handle before calling this function).
fn resolve_iupac(code: u8, rng: &mut impl Rng) -> Option<u8> {
    // Ambiguity sets per the IUPAC nucleotide code table.
    let alts: &[u8] = match code {
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"CG",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' => b"ACGT",
        _ => return None,
    };
    Some(alts[rng.random_range(0..alts.len())])
}

/// Normalize a contig's bases in place per the rules documented on
/// [`Fasta::load_contig`].
///
/// Invalid bytes produce an error identifying the offending byte and its
/// 0-based position within the contig.
fn normalize_bases(seq: &mut [u8], contig_name: &str, rng: &mut impl Rng) -> Result<()> {
    for (i, b) in seq.iter_mut().enumerate() {
        let upper = b.to_ascii_uppercase();
        match upper {
            b'A' | b'C' | b'G' | b'T' => *b = upper,
            b'U' => *b = b'T',
            b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' | b'N' => {
                // Resolve to lowercase (marker for "synthesized from ambiguity").
                let picked = resolve_iupac(upper, rng).expect("code is a valid IUPAC code");
                *b = picked.to_ascii_lowercase();
            }
            _ => {
                return Err(anyhow!(
                    "Contig '{contig_name}' contains unrecognized base {:?} (0x{:02X}) at position {i}",
                    char::from(*b),
                    *b
                ));
            }
        }
    }
    Ok(())
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
    use rand::SeedableRng;
    use rand::rngs::SmallRng;

    /// Build a fresh deterministic RNG for tests.
    fn test_rng() -> SmallRng {
        SmallRng::seed_from_u64(42)
    }

    #[test]
    fn test_load_contig_single() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_test_fasta(dir.path(), &[("chr1", b"ACGTACGT")]);
        let mut fasta = Fasta::from_path(&path).unwrap();

        assert_eq!(fasta.contig_names(), vec!["chr1"]);
        assert_eq!(fasta.dict().len(), 1);

        let seq = fasta.load_contig("chr1", &mut test_rng()).unwrap();
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

        let seq1 = fasta.load_contig("chr1", &mut test_rng()).unwrap();
        assert_eq!(&seq1, b"AAAA");

        let seq2 = fasta.load_contig("chr2", &mut test_rng()).unwrap();
        assert_eq!(&seq2, b"CCCCGGGG");

        let seq3 = fasta.load_contig("chr3", &mut test_rng()).unwrap();
        assert_eq!(&seq3, b"TT");
    }

    #[test]
    fn test_load_contig_out_of_order() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_test_fasta(dir.path(), &[("chr1", b"AAAA"), ("chr2", b"CCCC")]);
        let mut fasta = Fasta::from_path(&path).unwrap();

        // Load chr2 first, then chr1 — verify seeking works correctly.
        let seq2 = fasta.load_contig("chr2", &mut test_rng()).unwrap();
        assert_eq!(&seq2, b"CCCC");

        let seq1 = fasta.load_contig("chr1", &mut test_rng()).unwrap();
        assert_eq!(&seq1, b"AAAA");
    }

    #[test]
    fn test_load_contig_uppercases_acgt() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_test_fasta(dir.path(), &[("chr1", b"acgtACGT")]);
        let mut fasta = Fasta::from_path(&path).unwrap();

        let seq = fasta.load_contig("chr1", &mut test_rng()).unwrap();
        assert_eq!(&seq, b"ACGTACGT");
    }

    #[test]
    fn test_load_contig_u_to_t() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_test_fasta(dir.path(), &[("chr1", b"AUGCUu")]);
        let mut fasta = Fasta::from_path(&path).unwrap();

        let seq = fasta.load_contig("chr1", &mut test_rng()).unwrap();
        assert_eq!(&seq, b"ATGCTT");
    }

    #[test]
    fn test_load_contig_unknown() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_test_fasta(dir.path(), &[("chr1", b"ACGT")]);
        let mut fasta = Fasta::from_path(&path).unwrap();

        let result = fasta.load_contig("chrZ", &mut test_rng());
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("not found"));
    }

    #[test]
    fn test_load_contig_rejects_unknown_byte() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_test_fasta(dir.path(), &[("chr1", b"ACZT")]);
        let mut fasta = Fasta::from_path(&path).unwrap();

        let err = fasta.load_contig("chr1", &mut test_rng()).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("unrecognized base"), "message was: {msg}");
        assert!(msg.contains("position 2"), "message was: {msg}");
    }

    #[test]
    fn test_load_contig_iupac_resolved_to_lowercase() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_test_fasta(dir.path(), &[("chr1", b"ARYSWKMBDHVN")]);
        let mut fasta = Fasta::from_path(&path).unwrap();

        let seq = fasta.load_contig("chr1", &mut test_rng()).unwrap();
        assert_eq!(seq.len(), 12);
        // First base is a real A.
        assert_eq!(seq[0], b'A');
        // All ambiguity-resolved bases should be lowercase ACGT.
        for (i, &b) in seq.iter().enumerate().skip(1) {
            assert!(matches!(b, b'a' | b'c' | b'g' | b't'), "position {i} got {b:?}");
        }
    }

    #[test]
    fn test_load_contig_resolution_respects_ambiguity_set() {
        let dir = tempfile::tempdir().unwrap();
        // Long runs so we observe multiple draws per code.
        let path = write_test_fasta(
            dir.path(),
            &[("chr1", &[b'R'; 200]), ("chr2", &[b'Y'; 200]), ("chr3", &[b'K'; 200])],
        );
        let mut fasta = Fasta::from_path(&path).unwrap();

        let r = fasta.load_contig("chr1", &mut test_rng()).unwrap();
        let y = fasta.load_contig("chr2", &mut test_rng()).unwrap();
        let k = fasta.load_contig("chr3", &mut test_rng()).unwrap();

        // R ∈ {A,G} (lowercase)
        for &b in &r {
            assert!(matches!(b, b'a' | b'g'), "R resolved to {b:?}");
        }
        // Y ∈ {C,T}
        for &b in &y {
            assert!(matches!(b, b'c' | b't'), "Y resolved to {b:?}");
        }
        // K ∈ {G,T}
        for &b in &k {
            assert!(matches!(b, b'g' | b't'), "K resolved to {b:?}");
        }
    }

    #[test]
    fn test_load_contig_resolution_reproducible() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_test_fasta(dir.path(), &[("chr1", b"NNNNNNNNNNRRYYSSWWKKMMBBDDHHVV")]);
        let mut fasta_a = Fasta::from_path(&path).unwrap();
        let mut fasta_b = Fasta::from_path(&path).unwrap();

        let a = fasta_a.load_contig("chr1", &mut SmallRng::seed_from_u64(1234)).unwrap();
        let b = fasta_b.load_contig("chr1", &mut SmallRng::seed_from_u64(1234)).unwrap();
        assert_eq!(a, b);
    }
}
