//! Extension-aware VCF output.
//!
//! A VCF's extension is a contract about its codec: `.gz`/`.bgz` means
//! BGZF-compressed, anything else means plain text. Downstream tools rely on
//! it — `tabix`/`bcftools` expect a `*.vcf.gz` to be a real BGZF stream, and an
//! extension-based reader rejects a plain-text file named `*.vcf.gz` with an
//! opaque BGZF error. [`VcfWriter`] honours that contract: it BGZF-compresses
//! when (and only when) the output path ends in `.gz`/`.bgz`, so every VCF
//! holodeck writes is named truthfully and round-trips through
//! `simulate`/`methylate` regardless of name. (holodeck's own VCF reader sniffs
//! the leading magic bytes, so it would tolerate a mislabelled file, but other
//! tools — and users — should not have to.)

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use noodles_bgzf as bgzf;

/// A VCF writer whose output is plain text or BGZF-compressed, chosen by file
/// extension so the written file is named truthfully (`.gz`/`.bgz` => BGZF,
/// otherwise plain text).
///
/// The two codecs are distinct variants (rather than a single `Box<dyn Write>`)
/// so that [`close`](Self::close) can call BGZF's consuming `finish` and
/// surface a failed EOF-block write — a `Drop`-based finalize would silently
/// discard that error and leave a truncated file. Callers write records with
/// `write!`/`writeln!` via the [`Write`] impl, then [`close`](Self::close) the
/// writer to flush and (for BGZF) emit the EOF block.
pub enum VcfWriter {
    /// Uncompressed text output (path does not end in `.gz`/`.bgz`).
    Plain(BufWriter<File>),
    /// BGZF-compressed output (path ends in `.gz`/`.bgz`).
    Bgzf(bgzf::io::Writer<File>),
}

impl VcfWriter {
    /// Open a VCF writer at `path`, BGZF-compressing when the extension is
    /// `.gz` or `.bgz` (case-insensitive) and writing plain text otherwise.
    ///
    /// # Errors
    /// Returns an error if the output file cannot be created.
    pub fn new(path: &Path) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create VCF file: {}", path.display()))?;
        Ok(if is_bgzf_path(path) {
            VcfWriter::Bgzf(bgzf::io::Writer::new(file))
        } else {
            VcfWriter::Plain(BufWriter::new(file))
        })
    }

    /// Finalize the VCF file: flush buffered data and, for BGZF, write and
    /// verify the EOF block.
    ///
    /// # Errors
    /// Returns an error if the final flush or BGZF EOF-block write fails —
    /// surfacing close-time truncation that a `Drop`-based finalize would
    /// silently swallow.
    pub fn close(self) -> Result<()> {
        match self {
            VcfWriter::Plain(mut w) => w.flush()?,
            VcfWriter::Bgzf(w) => {
                w.finish()?;
            }
        }
        Ok(())
    }
}

impl Write for VcfWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            VcfWriter::Plain(w) => w.write(buf),
            VcfWriter::Bgzf(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            VcfWriter::Plain(w) => w.flush(),
            VcfWriter::Bgzf(w) => w.flush(),
        }
    }
}

/// Whether `path` should be BGZF-compressed, i.e. its final extension is `.gz`
/// or `.bgz` (case-insensitive) — the conventional contract that `tabix`/
/// `bcftools` and other ecosystem tools assume for a `*.vcf.gz` name.
fn is_bgzf_path(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz") || ext.eq_ignore_ascii_case("bgz"))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// First two bytes of any gzip/BGZF member.
    const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

    /// The canonical 28-byte BGZF EOF marker that closes a well-formed BGZF
    /// stream (mirrors `noodles_bgzf`'s internal `BGZF_EOF`).
    const BGZF_EOF: [u8; 28] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];

    #[test]
    fn is_bgzf_path_detects_gz_and_bgz_case_insensitively() {
        assert!(is_bgzf_path(Path::new("muts.vcf.gz")));
        assert!(is_bgzf_path(Path::new("muts.vcf.bgz")));
        assert!(is_bgzf_path(Path::new("MUTS.VCF.GZ")));
        assert!(!is_bgzf_path(Path::new("muts.vcf")));
        assert!(!is_bgzf_path(Path::new("muts")));
        // Only the final extension matters, matching the reader.
        assert!(!is_bgzf_path(Path::new("muts.gz.vcf")));
    }

    /// A `.gz` path must produce a real BGZF stream: gzip magic up front and a
    /// trailing BGZF EOF block once closed.
    #[test]
    fn gz_path_writes_bgzf_with_eof_block() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.vcf.gz");
        let mut writer = VcfWriter::new(&path).unwrap();
        writeln!(writer, "##fileformat=VCFv4.3").unwrap();
        writer.close().unwrap();

        let bytes = std::fs::read(&path).unwrap();
        assert_eq!(&bytes[..2], &GZIP_MAGIC, "expected gzip/BGZF magic bytes");
        assert!(bytes.ends_with(&BGZF_EOF), "expected trailing BGZF EOF block");
    }

    /// A non-`.gz` path must produce uncompressed text -- no gzip magic.
    #[test]
    fn plain_path_writes_uncompressed_text() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.vcf");
        let mut writer = VcfWriter::new(&path).unwrap();
        writeln!(writer, "##fileformat=VCFv4.3").unwrap();
        writer.close().unwrap();

        let bytes = std::fs::read(&path).unwrap();
        assert_ne!(&bytes[..2], &GZIP_MAGIC, "plain output must not be gzip");
        assert_eq!(bytes, b"##fileformat=VCFv4.3\n");
    }
}
