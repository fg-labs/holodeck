//! FASTQ output for simulated reads.
//!
//! Writes bgzf-compressed FASTQ files. Supports both single-threaded
//! compression (via noodles-bgzf) and multi-threaded compression (via
//! pooled-writer) depending on how the writer is constructed.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use noodles_bgzf as bgzf;

use crate::read::SimulatedRead;

/// A FASTQ writer that writes bgzf-compressed output.
///
/// The inner writer handles BGZF block framing and compression.  Use
/// [`new`](Self::new) for single-threaded compression or
/// [`from_writer`](Self::from_writer) to supply a pre-configured writer
/// (e.g. a [`pooled_writer::PooledWriter`] for multi-threaded compression).
pub struct FastqWriter {
    writer: Box<dyn Write>,
}

impl FastqWriter {
    /// Create a new single-threaded FASTQ writer at the given path.
    ///
    /// Uses noodles-bgzf with the specified compression level (0-12).
    ///
    /// # Errors
    /// Returns an error if the file cannot be created or the compression
    /// level is invalid.
    pub fn new(path: &Path, compression: u8) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create FASTQ file: {}", path.display()))?;
        let level = bgzf::io::writer::CompressionLevel::new(compression)
            .ok_or_else(|| anyhow::anyhow!("invalid compression level: {compression}"))?;
        let writer = bgzf::io::writer::Builder::default()
            .set_compression_level(level)
            .build_from_writer(BufWriter::new(file));
        Ok(Self { writer: Box::new(writer) })
    }

    /// Create a FASTQ writer from an existing writer that handles BGZF
    /// compression (e.g. a [`pooled_writer::PooledWriter`]).
    pub fn from_writer(writer: impl Write + 'static) -> Self {
        Self { writer: Box::new(writer) }
    }

    /// Write a single read as a FASTQ record.
    ///
    /// Writes the four-line FASTQ format:
    /// ```text
    /// @read_name
    /// BASES
    /// +
    /// QUALITIES
    /// ```
    ///
    /// # Errors
    /// Returns an error if writing fails.
    pub fn write_read(&mut self, read: &SimulatedRead) -> Result<()> {
        self.writer.write_all(b"@")?;
        self.writer.write_all(read.name.as_bytes())?;
        self.writer.write_all(b"\n")?;
        self.writer.write_all(&read.bases)?;
        self.writer.write_all(b"\n+\n")?;
        self.writer.write_all(&read.qualities)?;
        self.writer.write_all(b"\n")?;
        Ok(())
    }

    /// Finalize the FASTQ file.
    ///
    /// Drops the underlying writer, which flushes any buffered data and
    /// (for BGZF writers) writes the EOF marker.  For pooled writers, this
    /// sends remaining data to the compression pool.
    pub fn close(self) {
        drop(self.writer);
    }
}

#[cfg(test)]
mod tests {
    use std::io::Read;

    use super::*;
    use crate::read::SimulatedRead;

    #[test]
    fn test_fastq_write_and_read_back() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.fastq.gz");

        // Write two reads.
        {
            let mut w = FastqWriter::new(&path, 1).unwrap();
            w.write_read(&SimulatedRead {
                name: "read1".to_string(),
                bases: b"ACGT".to_vec(),
                qualities: b"IIII".to_vec(),
            })
            .unwrap();
            w.write_read(&SimulatedRead {
                name: "read2".to_string(),
                bases: b"TTAA".to_vec(),
                qualities: b"????".to_vec(),
            })
            .unwrap();
            w.close();
        }

        // Read back and verify.
        let file = File::open(&path).unwrap();
        let mut decoder = flate2::read::MultiGzDecoder::new(file);
        let mut contents = String::new();
        decoder.read_to_string(&mut contents).unwrap();

        assert_eq!(contents, "@read1\nACGT\n+\nIIII\n@read2\nTTAA\n+\n????\n");
    }

    #[test]
    fn test_fastq_write_empty_read() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.fastq.gz");

        {
            let mut w = FastqWriter::new(&path, 1).unwrap();
            w.write_read(&SimulatedRead {
                name: "r1".to_string(),
                bases: Vec::new(),
                qualities: Vec::new(),
            })
            .unwrap();
            w.close();
        }

        let file = File::open(&path).unwrap();
        let mut decoder = flate2::read::MultiGzDecoder::new(file);
        let mut contents = String::new();
        decoder.read_to_string(&mut contents).unwrap();

        assert_eq!(contents, "@r1\n\n+\n\n");
    }
}
