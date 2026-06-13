# Changelog

All notable changes to holodeck are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

- First release of holodeck, an NGS read simulator.

[unreleased]: https://github.com/fg-labs/holodeck/compare/v0.2.1...HEAD
[0.2.1]: https://github.com/fg-labs/holodeck/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/fg-labs/holodeck/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/fg-labs/holodeck/releases/tag/v0.1.0
