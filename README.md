[![Build](https://github.com/fg-labs/holodeck/actions/workflows/check.yml/badge.svg)](https://github.com/fg-labs/holodeck/actions/workflows/check.yml)
[![Version at crates.io](https://img.shields.io/crates/v/holodeck)](https://crates.io/crates/holodeck)
[![Documentation at docs.rs](https://img.shields.io/docsrs/holodeck)](https://docs.rs/holodeck)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/holodeck.svg?label=bioconda)](https://bioconda.github.io/recipes/holodeck/README.html)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fg-labs/holodeck/blob/main/LICENSE)

# Holodeck

Modern NGS read simulator written in Rust.

<p>
<a href="https://fulcrumgenomics.com">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset=".github/logos/fulcrumgenomics-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset=".github/logos/fulcrumgenomics-light.svg">
  <img alt="Fulcrum Genomics" src=".github/logos/fulcrumgenomics-light.svg" height="100">
</picture>
</a>
</p>

[Visit us at Fulcrum Genomics](https://www.fulcrumgenomics.com) to learn more about how we can power your Bioinformatics with holodeck and beyond.

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-%2338b44a.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-%2326a8e0.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

## Overview

Holodeck generates simulated sequencing reads from a reference genome with optional variants from a VCF file.  It supports Illumina-style paired-end and single-end reads with a position-dependent error model, targeted sequencing with BED files, and arbitrary ploidy with per-contig overrides.

Simulated reads include ground-truth information encoded in read names, and optionally in a golden BAM file, enabling downstream evaluation of alignment and variant-calling accuracy.

## Installation

### From source

Requires Rust 1.94.0 or later.

```bash
git clone https://github.com/fulcrumgenomics/holodeck.git
cd holodeck
cargo build --release
# Binary is at target/release/holodeck
```

## Commands

| Command | Description |
|---------|-------------|
| `simulate` | Generate simulated reads from a reference genome with optional VCF variants |
| `mutate` | Generate a random VCF of mutations from a reference genome |
| `eval` | Evaluate alignment accuracy of simulated reads against truth |

## Quick Start

```bash
# Index your reference (if not already done)
samtools faidx ref.fa

# Generate random mutations
holodeck mutate -r ref.fa -o mutations.vcf --snp-rate 0.001

# Simulate 30x paired-end reads with ground-truth BAM
holodeck simulate -r ref.fa -v mutations.vcf -o sim --coverage 30 --golden-bam

# Outputs: sim.r1.fastq.gz, sim.r2.fastq.gz, sim.golden.bam

# Align reads and evaluate accuracy
minimap2 -a ref.fa sim.r1.fastq.gz sim.r2.fastq.gz | samtools sort -o mapped.bam
samtools index mapped.bam
holodeck eval --mapped mapped.bam -o eval_results
# Output: eval_results.eval.txt
```

## Simulate

Generate paired-end or single-end FASTQ files from a reference genome with optional variants from a VCF.  Reads are sampled with a position-dependent Illumina error model where the error rate ramps from a minimum at the start of reads to a maximum at the end.

```bash
# Whole-genome 30x with variants
holodeck simulate -r ref.fa -v variants.vcf -o output --coverage 30

# Targeted sequencing (exome/panel)
holodeck simulate -r ref.fa -v variants.vcf -b targets.bed -o output --coverage 100

# Single-end, no errors, with golden BAM
holodeck simulate -r ref.fa -o output --single-end --min-error-rate 0 --max-error-rate 0 --golden-bam

# Custom fragment size distribution and compression threads
holodeck simulate -r ref.fa -o output --fragment-mean 400 --fragment-stddev 80 -t 8
```

**Key options:**

| Option | Default | Description |
|--------|---------|-------------|
| `-r, --reference` | required | Indexed FASTA reference |
| `-v, --vcf` | none | VCF with variants to apply |
| `-b, --targets` | none | BED file of target regions |
| `-o, --output` | required | Output prefix |
| `-c, --coverage` | 30 | Mean coverage depth |
| `-l, --read-length` | 150 | Read length in bases |
| `-d, --fragment-mean` | 300 | Mean fragment size |
| `-s, --fragment-stddev` | 50 | Fragment size standard deviation |
| `--min-error-rate` | 0.001 | Error rate at start of reads |
| `--max-error-rate` | 0.01 | Error rate at end of reads |
| `--max-n-frac` | 0.5 | Reject reads with >this fraction of bases from ambiguous reference positions (see [Ambiguous reference bases](#ambiguous-reference-bases)) |
| `--golden-bam` | off | Write ground-truth BAM |
| `--single-end` | off | Generate SE instead of PE reads |
| `--simple-names` | off | Use `holodeck::N` names instead of encoded truth |
| `--compression` | 1 | BGZF compression level (0-12) |
| `-t, --threads` | 4 | Threads for BGZF compression |
| `--seed` | auto | Random seed (deterministic by default) |

### Ambiguous reference bases

Real references contain a mix of `A`/`C`/`G`/`T`, large stretches of `N` (assembly gaps), and — rarely — IUPAC ambiguity codes (`R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`) for known-ambiguous positions. Holodeck reads these and turns them into ACGT reads without leaking non-ACGT characters into emitted FASTQ, using a two-step strategy:

1. **Reference normalization at load time.** When a contig is loaded, every byte is classified:
   - `A`/`C`/`G`/`T` (either case) → stored as uppercase.
   - `U` → converted to `T` (so RNA-style references work).
   - Any IUPAC ambiguity code (including `N`) → resolved to a uniformly-random base drawn from that code's ambiguity set (e.g. `R` → `A` or `G`, `N` → `A`/`C`/`G`/`T`). The drawn base is stored in **lowercase** as a marker that this position was synthesized from ambiguity.
   - Anything else → hard error with the offending byte and position.

   The RNG used for this normalization is seeded deterministically from `--seed` plus the contig name, so two runs with the same seed produce byte-identical outputs even when the reference has ambiguous positions.

2. **Read rejection at sampling time.** Each generated read counts how many of its bases are lowercase (i.e. came from ambiguous positions). If either R1 or R2 has a lowercase fraction above `--max-n-frac` (default `0.5`), the pair is rejected and resampled. Accepted reads are upper-cased in place before the error model runs, so emitted FASTQ and BAM contain only `A`/`C`/`G`/`T` (plus a rare `N` when the configured adapter is shorter than the bases needed past the insert — a separate, pre-existing behavior of the adapter-padding code).

   Set `--max-n-frac 1.0` to disable the filter (accept reads from any region). Set `--max-n-frac 0.0` to require every base in every read to come from an unambiguous reference position.

**Known limitation:** requested `--coverage` is computed from raw contig/BED lengths, not from the non-ambiguous territory. For a reference like hs38DH (~5% N), rejection is noise and coverage lands where you'd expect. For simulations targeted at heavily-N contigs (or with `--max-n-frac 0.0` in N-dense regions), effective coverage will be slightly below the requested value; a warning is logged if the resampling budget is exhausted.

## Mutate

Generate a VCF of random mutations from a reference genome.  The output can be fed directly to `holodeck simulate`.  Supports independent control of SNP, indel, and MNP rates with configurable ploidy.

```bash
# Default rates
holodeck mutate -r ref.fa -o mutations.vcf

# Realistic human-like rates (~1 variant per 1000bp, 10:1 SNP:indel)
holodeck mutate -r ref.fa -o mutations.vcf --snp-rate 0.000727 --indel-rate 0.0000727 --mnp-rate 0

# Male sample with haploid chrX/chrY
holodeck mutate -r ref.fa -o mutations.vcf \
  --ploidy 2 --ploidy-override chrX=1 --ploidy-override chrY=1

# Restrict mutations to target regions
holodeck mutate -r ref.fa -o mutations.vcf -b targets.bed
```

**Key options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--snp-rate` | 0.001 | SNP rate per base |
| `--indel-rate` | 0.0001 | Indel rate per base |
| `--mnp-rate` | 0.00005 | MNP rate per base |
| `--het-hom-ratio` | 2.0 | Ratio of heterozygous to homozygous variants |
| `--ploidy` | 2 | Default ploidy |
| `--ploidy-override` | none | Per-contig/region ploidy (e.g. `chrX=1`) |
| `--indel-length-param` | 0.7 | Geometric distribution parameter for indel lengths |

## Eval

Evaluate alignment accuracy by comparing mapped positions against truth positions encoded in read names.  Reports accuracy stratified by MAPQ bin.

```bash
holodeck eval --mapped aligned.bam -o eval_results
holodeck eval --mapped aligned.bam -o eval_results --wiggle 10
```

**Key options:**

| Option | Default | Description |
|--------|---------|-------------|
| `-m, --mapped` | required | BAM file of mapped reads |
| `-o, --output` | required | Output prefix (writes `.eval.txt`) |
| `--wiggle` | 5 | Max distance (bp) for a correct mapping |

## Features

- **Position-dependent error model** -- error rate ramps across the read, with R2 having higher rates than R1 (configurable multiplier)
- **Multi-sample VCF support** -- select a sample with `--sample` from multi-sample VCFs
- **Arbitrary ploidy** -- phased and unphased genotypes with any ploidy; per-contig and per-region overrides for sex chromosomes and PAR regions
- **BED target regions** -- efficient padded-interval sampling for exome/panel simulations
- **Adapter simulation** -- configurable adapter sequences appended when fragment < read length
- **Encoded read names** -- truth coordinates embedded in read names for downstream evaluation without needing the golden BAM
- **Golden BAM output** -- ground-truth alignments with correct positions, CIGARs (including variant-induced indels), and adapter soft-clips
- **Multi-threaded compression** -- BGZF output compression parallelized across configurable threads via `pooled-writer`
- **Deterministic by default** -- same parameters always produce the same output via hash-derived seeding
- **Sparse haplotype representation** -- variants stored in interval trees overlaid on the reference; no full-genome copies

## Read Name Format

Read names encode ground-truth alignment information for use by `holodeck eval` and other tools.

**Encoded format** (default):
```text
PE: @holodeck::READ_NUM::FRAG_LEN::CONTIG::POS1+STRAND::POS2+STRAND::HAP::ERRS1::ERRS2
SE: @holodeck::READ_NUM::FRAG_LEN::CONTIG::POS+STRAND::HAP::ERRS
```

Example:
```text
@holodeck::42::600::chr1::10000F::10450R::0::2::1
```

Fields: read number, source fragment length (bp), contig, R1 position (1-based) + strand (F/R), R2 position + strand, haplotype index, R1 errors, R2 errors.

`FRAG_LEN` is the length of the originating template. When it is smaller than the read length, the remainder of each read (`read_length - FRAG_LEN` bases) is adapter sequence (possibly padded with `N` if the configured adapter is shorter than that remainder). This lets consumers locate the adapter boundary in both PE and SE reads without the golden BAM.

Fields are separated by `::` (double colon) so contig names that legally contain single `:` characters (e.g. HLA alleles like `HLA-A*01:01:01:01`) parse unambiguously. Contig names must not contain `@` (the FASTQ header prefix) or `::`; both are rejected by a debug assertion in the read-name formatter.

**Simple format** (`--simple-names`):
```text
@holodeck::42
```

## Performance

Holodeck is designed for high throughput.  Typical performance on Apple M-series hardware:

| Scenario | Reads/sec |
|----------|-----------|
| FASTQ only | ~330K read pairs/sec |
| FASTQ + golden BAM | ~170K read pairs/sec |

Performance scales linearly with read count.  Adding realistic variant density (~1 per 1,000bp) has no measurable impact on throughput.

## License

[MIT](LICENSE)
