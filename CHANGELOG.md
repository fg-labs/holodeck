# Changelog

All notable changes to holodeck are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `eval` now scores accuracy against holodeck's own truth beyond placement.
  With `--truth` (the golden BAM) it takes per-read true positions, spans, and
  sequences from the golden alignment rather than only the encoded read name,
  and `--variants` reports how faithfully aligned reads represent the simulated
  substitutions. The allele a read truly carries at each truth site is read
  from the golden read's own sequence — the per-read oracle — so scoring is
  correct whether or not the truth VCF is phased (a read sequenced from the
  reference copy shows the reference base and is not expected to carry the
  alt). For every such expected substitution it walks the mapped read's CIGAR
  to the variant position and checks the observed base, accumulating the
  represented fraction with the read's MAPQ and alignment score. With
  `--reference` it also reports per-read NM/MD concordance as a **bisulfite-aware
  genomic edit distance**: rather than comparing raw `NM:i`/`MD:Z` tags (which
  are convention-dependent — a bisulfite aligner may score against the original
  or the converted reference, so the tags differ even when both placed the read
  correctly), it recomputes each read's edits against the reference and excludes
  conversions using the read's TRUE strand (taken from the golden truth, so it
  works even for aligners such as bwameth that emit no `XG`). The result is
  comparable across aligners; without `--reference` NM/MD concordance is `NA`.
  `--meth` breaks the variant results down by bisulfite substitution class,
  labelling the conversion-confounded `C->T`/`G->A` cell as such. `--cpg-truth`
  correlates the aligner's Bismark `XM` calls against the simulated cpg-truth
  bedGraph (Pearson r and RMSE; `NA` for aligners that emit no `XM`). Results
  are written to `<prefix>.variants.tsv` and `<prefix>.meth.tsv` alongside the
  existing `<prefix>.eval.txt`.

### Changed

- `methylate` now methylates contigs in parallel. Each contig is independent —
  its RNGs are seeded purely from `(seed, contig)` with no state carried between
  contigs (the methylation Markov chain runs within a single contig) — so the
  per-contig loop runs as a work-stealing parallel map and the output is
  byte-identical to the previous single-threaded version regardless of thread
  count. Single-item jobs spread the wildly-uneven per-contig cost (chr1 ≫ a
  50 kb alt) evenly across the pool, and one reused FASTA handle per worker
  avoids re-parsing the sequence dictionary per contig. On a whole human genome
  this is roughly 5× faster on a 12-core host (~110 s → ~25 s including BGZF
  output); the thread count honors `RAYON_NUM_THREADS`.

### Fixed

- `methylate --vcf` (allele-specific methylation) no longer slows down
  quadratically with variant density. The per-haplotype CpG classifier looked
  up each CpG's variant/reference source with a linear scan over every variant
  on the contig, making it O(CpGs × variants) per haplotype — on a whole human
  genome with a few-million-variant VCF this was ~100× the work of the
  reference-only path (≈43 min vs ≈25 s). Because the alt spans are sorted and
  disjoint and the CpG scan is ascending, a single monotonic cursor now resolves
  each lookup in O(1) amortized, making classification O(haplotype length +
  variants); output is byte-identical. The reference-only path was already fast
  and is unchanged.
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
