# Changelog

All notable changes to holodeck are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- `methylate` no longer assigns methylation to CpGs in ambiguous reference
  sequence. `N`s and other ambiguity codes are resolved to random bases, which
  contain a `CG` about every 16 bases, and each of these was written as a CpG
  record; on hs38DH they were 26% of all records. A `CG` now counts as a CpG
  only when both bases are real reference or variant bases. `simulate` applies
  the same rule to its CpG truth bedGraph and to the `XM`/`YM` tags in the
  golden BAM. For a given `--seed`, methylation values differ from earlier
  versions on any reference that contains ambiguity codes. Methylation VCFs
  written by earlier versions should be regenerated: most still load, but one
  with a variant next to an ambiguity code is rejected with an MT/MB length
  mismatch.
- `methylate` is faster and, without an input VCF, uses less memory. On hs38DH,
  single-threaded: 87 s to 24 s and 2.1 GB to 0.6 GB without a VCF; 126 s to
  49 s with a 3.5-million-variant VCF. Part of this comes from the change
  above; the rest leaves the output unchanged. `simulate` shares the faster
  haplotype fragment extraction.
- Lowercase bases in a VCF alternate allele are treated as real bases, not as
  ambiguous reference positions, when `simulate` filters reads with
  `--max-n-frac`.

### Fixed

- `methylate --vcf` no longer slows down in proportion to the number of
  variants on a contig. Each CpG's variant lookup scanned every variant on the
  contig; it now takes constant time, so run time is close to that of
  `methylate` without a VCF (on human chr22 with 42,000 variants, 87 s before
  and under 3 s after). Output is unchanged.
- `simulate` and `methylate` now tolerate VCFs that redeclare a header ID
  (e.g. `duplicate INFO ID: BREAKSIMLENGTH`). Such duplicates are common in
  files from upstream tools and are accepted by bcftools; holodeck drops the
  repeated definitions (keeping the first) instead of erroring. Compression is
  now detected from the file's magic bytes rather than its extension, so the
  variant reader also accepts plain-gzip and extension-less inputs.
- `mutate` and `methylate` now choose VCF output compression from the file
  extension: `.gz`/`.bgz` paths are BGZF-compressed, everything else is plain
  text. Previously `mutate` always wrote uncompressed text (so a `*.vcf.gz`
  output was rejected by `simulate` with an opaque BGZF error) and `methylate`
  always wrote BGZF (so a `*.vcf` output was a BGZF stream the reader treated as
  plain text). Both now round-trip through `simulate`/`methylate` regardless of
  name.
- A reference FASTA without a `.fai` index now fails with an error that names
  the missing index and the `samtools faidx` command that creates it, instead
  of `Failed to open indexed FASTA: ref.fa` followed by `No such file or
  directory`, which read as though the FASTA itself were missing. A missing
  FASTA is likewise reported by name.

## [0.3.0] - 2026-06-13

### Added

- End-to-end methylation simulation as two composable steps. The new
  `methylate` subcommand scans every CpG in the reference (after applying any
  VCF variants) and writes per-haplotype, per-strand methylation truth as
  `MT`/`MB` FORMAT fields, using a context-aware, spatially correlated model:
  CpGs are classified de novo into island / shore / open-sea (Gardiner-Garden)
  with per-context target rates and correlation lengths, walked by a two-state
  Markov chain so methylation forms realistic autocorrelated runs. Methylation
  is symmetric by default with a low sporadic hemimethylation rate, and
  per-haplotype draws yield allele-specific methylation.
- `simulate --methylation-mode` applies EM-seq/bisulfite (unmethylated C->T) or
  TAPS (methylated C->T) conversion when the input VCF carries `MT`/`MB` truth,
  including a bimodal per-molecule conversion-failure model.
- Golden BAM now emits Bismark-compatible methylation tags (`XG`/`XR`/`XM`/`NM`/`MD`)
  plus holodeck truth tags (`YM`/`YS`/`cf`), so the perfect-truth BAM drops
  straight into Bismark / MethylDackel / IGV.
- Methylation truth outputs: `--cpg-truth-bedgraph` (coverage-weighted, from the
  reads that covered each CpG) and a closed-form population-fraction bedGraph,
  both in MethylDackel `extract` format.

## [0.2.1] - 2026-05-01

### Changed

- Non-ACGT bases in the reference FASTA are now resolved to a concrete base
  (e.g. `N`->`ACGT`, `R`->`AG`) when simulating reads, in a way that keeps the
  base identifiable, rather than passed through unchecked.
- Lowered the default `simulate --max-n-frac` to `0.02`; reads or read pairs are
  rejected when any read exceeds this fraction of ambiguous bases.

### Fixed

- Use absolute URLs for logo images in the README.

## [0.2.0] - 2026-04-22

### Added

- Source fragment length is now encoded in read names.

### Changed

- The read-name field separator changed from `:` to `::` to cleanly handle
  contig names containing single colons (e.g. HLA alleles).

## [0.1.0] - 2026-04-22

### Added

- First release of holodeck, an NGS read simulator.

[unreleased]: https://github.com/fg-labs/holodeck/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/fg-labs/holodeck/compare/v0.2.1...v0.3.0
[0.2.1]: https://github.com/fg-labs/holodeck/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/fg-labs/holodeck/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/fg-labs/holodeck/releases/tag/v0.1.0
