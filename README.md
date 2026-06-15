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
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/fg-labs/holodeck/main/.github/logos/fulcrumgenomics-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/fg-labs/holodeck/main/.github/logos/fulcrumgenomics-light.svg">
  <img alt="Fulcrum Genomics" src="https://raw.githubusercontent.com/fg-labs/holodeck/main/.github/logos/fulcrumgenomics-light.svg" height="100">
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
git clone https://github.com/fg-labs/holodeck.git
cd holodeck
cargo build --release
# Binary is at target/release/holodeck
```

## Commands

| Command | Description |
|---------|-------------|
| `mutate` | Generate a random VCF of mutations from a reference genome |
| `methylate` | Annotate a VCF with per-haplotype CpG methylation state |
| `simulate` | Simulate sequencing reads from a reference (and optional methylation-annotated VCF) |
| `eval` | Evaluate alignment accuracy of simulated reads against truth |

## Quick Start

```bash
# Index your reference (if not already done)
samtools faidx ref.fa

# Generate random mutations
holodeck mutate -r ref.fa -o mutations.vcf --snp-rate 0.001

# Annotate the VCF with per-haplotype CpG methylation truth
holodeck methylate -r ref.fa -v mutations.vcf -o methylated.vcf.gz --seed 42

# Simulate 30x paired-end reads with ground-truth BAM
# (methylation chemistry applies because the input VCF has MT/MB fields from holodeck methylate;
#  for a non-methylation run, drop --methylation-mode and pass -v mutations.vcf instead)
holodeck simulate -r ref.fa -v methylated.vcf.gz -o sim --coverage 30 \
    --methylation-mode em-seq --methylation-conversion-rate 0.999 \
    --methylation-failure-rate 0.01 --golden-bam

# Outputs: sim.r1.fastq.gz, sim.r2.fastq.gz, sim.golden.bam (records tagged cf:i for conversion failures)

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

# Bisulfite/em-seq (two-step: methylate + simulate). Methylate uses
# context-aware defaults (islands hypo-, open-sea hyper-methylated).
holodeck methylate -r ref.fa -o methylated.vcf.gz --seed 42
holodeck simulate -r ref.fa -v methylated.vcf.gz -o output \
    --methylation-mode em-seq --methylation-conversion-rate 0.99 --golden-bam

# TAPS chemistry on the same methylated genome
holodeck simulate -r ref.fa -v methylated.vcf.gz -o output_taps \
    --methylation-mode taps --methylation-conversion-rate 0.99 --golden-bam
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
| `--min-fragment-length` | 20 | Minimum fragment length; sampled lengths below this are clamped up. Must be ≥1 |
| `--adapter-r1` | TruSeq R1 | Adapter appended to read 1 when the fragment is shorter than the read length |
| `--adapter-r2` | TruSeq R2 | Adapter appended to read 2 when the fragment is shorter than the read length |
| `--min-error-rate` | 0.001 | Error rate at start of reads |
| `--max-error-rate` | 0.01 | Error rate at end of reads |
| `--max-n-frac` | 0.02 | Reject reads with >this fraction of bases from ambiguous reference positions (see [Ambiguous reference bases](#ambiguous-reference-bases)) |
| `--clip-5p-rate` | 0.0 | Probability a read's 5' end carries a terminal soft-clip artifact (see [Terminal soft-clip artifacts](#terminal-soft-clip-artifacts)) |
| `--clip-3p-rate` | 0.0 | Probability a read's 3' end carries a terminal soft-clip artifact, independent of 3' adapter read-through |
| `--clip-length-mean` | 8 | Mean injected clip length (truncated geometric); only used when a clip rate is non-zero |
| `--clip-length-max` | 20 | Maximum injected clip length; only used when a clip rate is non-zero |
| `--methylation-mode` | none | Methylation chemistry: `em-seq` (or `bisulfite`) or `taps`. Presence enables methylation simulation |
| `--methylation-conversion-rate` | 0.999 | Chemistry efficiency for normally-converting molecules (probability that the converting class of C converts to T) |
| `--methylation-failure-rate` | 0.01 | Fraction of molecules that are whole-molecule conversion failures (convert at `1 − conversion-rate`). Requires `--methylation-mode` |
| `--cpg-truth-bedgraph` | none | Write per-CpG ground-truth methylation tally in MethylDackel `extract` bedGraph format. Requires `--methylation-mode` |
| `--golden-bam` | off | Write ground-truth BAM |
| `--golden-vcf` | off | Reserved for a coverage-annotated ground-truth VCF; **not yet implemented** (logs a warning and is otherwise a no-op) |
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

2. **Read rejection at sampling time.** Each generated read counts how many of its bases are lowercase (i.e. came from ambiguous positions). If either R1 or R2 has a lowercase fraction above `--max-n-frac` (default `0.02`), the pair is rejected and resampled. Accepted reads are upper-cased in place before the error model runs, so emitted FASTQ and BAM contain only `A`/`C`/`G`/`T` (plus a rare `N` when the configured adapter is shorter than the bases needed past the insert — a separate, pre-existing behavior of the adapter-padding code).

   Set `--max-n-frac 1.0` to disable the filter (accept reads from any region). Set `--max-n-frac 0.0` to require every base in every read to come from an unambiguous reference position.

**Known limitation:** requested `--coverage` is computed from raw contig/BED lengths, not from the non-ambiguous territory. For a reference like hs38DH (~5% N), rejection is noise and coverage lands where you'd expect. For simulations targeted at heavily-N contigs (or with `--max-n-frac 0.0` in N-dense regions), effective coverage will be slightly below the requested value; a warning is logged if the resampling budget is exhausted.

### Terminal soft-clip artifacts

By default holodeck's reads align essentially end-to-end: aside from VCF-derived indels and 3' adapter read-through, the simulated CIGAR is one long match block. Real libraries are messier — end-repair fill-in, damaged or non-templated read ends, and similar effects make aligners soft-clip a few bases off one or both ends of a meaningful fraction of reads. Twist EM-seq, for example, shows ~8% of reads carrying a soft clip, strongly 5'-biased and modal at 6–10 bp. The mechanism is protocol-agnostic, so the model is independent of `--methylation-mode` and useful for both methylation and non-methylation simulation.

`--clip-5p-rate` and `--clip-3p-rate` turn this on (both default to `0.0`, i.e. disabled). For each read, the 5' and 3' ends independently receive a clip with their configured probabilities; a clipped end's length is drawn from a truncated geometric distribution with mean `--clip-length-mean`, capped at `--clip-length-max`. The 5' and 3' rates are separate because real libraries are markedly asymmetric (the 5' end clips several times more often), while the length distribution is shared across both ends.

When an end is clipped, holodeck replaces those terminal bases with random sequence — so a downstream aligner re-derives the soft-clip rather than forcing a match — and records the soft-clip in the golden-BAM CIGAR (and advances the alignment start when a forward read's 5' end is clipped). Ground truth therefore stays exact: the aligned block still maps to the correct reference span. A disabled (or omitted) model draws no randomness, so existing seeded runs remain byte-identical.

### Golden BAM

Pass `--golden-bam` to write a perfect-truth BAM alongside the FASTQ output. Every record carries the alignment the simulator generated it from — contig, position, strand, CIGAR — at MAPQ 60, in unsorted (generation) order. Sort downstream if a coordinate-sorted BAM is needed.

The same truth coordinates are also encoded in the FASTQ read names (see [Read Name Format](#read-name-format)), so the golden BAM is optional for evaluation. It exists to plug holodeck output directly into BAM-consuming tooling (samtools, IGV, methylation extractors) without re-aligning.

#### Header

| Line | Fields | Notes |
| ---  | ---    | --- |
| `@SQ` | `SN` (contig name), `LN` (length) | One per contig in the reference, in dictionary order. |
| `@PG` | `ID:holodeck`, `PN:holodeck`, `VN:<crate version>`, `CL:<full command line>` | Single entry. `CL` records the verbatim invocation for reproducibility. |
| `@RG` | `ID:A`, `SM:<sample>`, `LB:<sample>`, `PL:ILLUMINA` | Single entry; every record's `RG:Z` tag points to `ID:A`. `PL` is hard-coded to `ILLUMINA` (the only error model holodeck ships). `SM` and `LB` are both set to the sample name — taken from the VCF when `--vcf` is set (use `--sample` to disambiguate multi-sample VCFs), otherwise defaulting to `holodeck-simulation`. |

#### Per-record tags

| Tag | Type | Description |
| --- | ---  | --- |
| `RG:Z` | string | Read-group identifier; always `A`. Ties the record to the single `@RG` header entry. |
| `hp:i` | integer | 0-based haplotype index the read was sampled from. `0` for haploid contigs and the first haplotype of polyploid contigs; `1` for the second haplotype, etc. Useful for restricting evaluation to a single haplotype, or for measuring allele-specific behavior. |
| `ne:i` | integer | Number of substituted bases the simulator injected into the record. Holodeck's error model is substitution-only (no indels), so this is a per-base substitution count. Lets you stratify alignment-accuracy or methylation-call evaluation by per-read error load without re-running with `--max-error-rate 0`. |

`SEQ` and `QUAL` are stored in reference (forward-strand) orientation, with the `REVERSE_COMPLEMENTED` flag set for reverse-strand records — i.e. the BAM convention, not the FASTQ orientation.

When `--methylation-mode` is also set, an additional eight methylation tags are emitted per record. See [Methylation simulation → Interpreting the methylation-simulated golden BAM](#interpreting-the-methylation-simulated-golden-bam) for the tag set.

### Methylation simulation

Holodeck models methylation biology and sequencing chemistry as two independent, composable steps.

#### Two-command flow

1. **`holodeck methylate`** annotates a VCF (or the bare reference, if no VCF is given) with per-haplotype, per-strand CpG methylation truth. It writes a BGZF-compressed VCF containing every variant from the input VCF plus new `MT`/`MB` FORMAT fields that encode which CpG sites on which haplotype × strand are methylated.
2. **`holodeck simulate`** reads that methylation-annotated VCF and, when `--methylation-mode` is set, applies the appropriate chemistry conversion to each read based on the `MT`/`MB` truth embedded in the VCF.

#### Why two commands?

Separating biology from chemistry lets you methylate a genome once and then re-sequence it under different chemistries — for example, compare em-seq vs TAPS coverage bias — without re-running the methylation assignment step. It also makes the methylation truth an explicit, inspectable artifact rather than a transient simulation parameter.

#### `methylate` flags

| Option | Default | Description |
|--------|---------|-------------|
| `-r, --reference` | required | Indexed FASTA reference |
| `-v, --vcf` | none | Input VCF with variants to methylate; methylates the unmodified reference if omitted |
| `--sample` | none | Sample name to select from a multi-sample input VCF (requires `--vcf`); the selected name is also used for the output sample column. When no VCF sample is in play, the output column is named `METHYLATE` |
| `--methylation-rate-island` | 0.1 | Target methylation fraction for CpG-island-interior CpGs (hypomethylated) |
| `--methylation-rate-shore` | 0.5 | Target methylation fraction for island-shore CpGs (within 2 kb of an island) |
| `--methylation-rate-open-sea` | 0.85 | Target methylation fraction for open-sea CpGs (the hypermethylated bulk) |
| `--methylation-correlation-length-island` | 1000 | Spatial correlation length (bp) for island CpGs; larger → longer like-methylated runs |
| `--methylation-correlation-length-shore` | 1000 | Spatial correlation length (bp) for shore CpGs |
| `--methylation-correlation-length-open-sea` | 1000 | Spatial correlation length (bp) for open-sea CpGs |
| `--hemimethylation-rate` | 0.01 | Probability a methylated CpG is made hemimethylated (one strand left unmethylated) |
| `--seed` | auto | Random seed for deterministic methylation draws |
| `-o, --output` | required | Output BGZF-compressed VCF path |
| `--bedgraph` | none | Write a MethylDackel-format population-fraction bedGraph from the methylation bitmap |

Methylation is assigned with a context-aware Markov model rather than an independent per-CpG coin flip: each CpG is classified into island / shore / open-sea (detected de novo from the sequence via Gardiner-Garden criteria), and a two-state chain walks the CpG list so that each context's mean methylation matches its target rate while neighbouring CpGs are spatially correlated (runs whose length scales with the correlation length). Methylation is symmetric across strands by default, with a low sporadic hemimethylation rate. For **uniform** methylation with no island structure (e.g. for controlled tests), set the three `--methylation-rate-*` flags equal. Note the no-flags default is now biologically realistic (islands hypo-, open-sea hyper-methylated), **not** fully methylated.

#### MT/MB FORMAT schema

Methylation state is stored in two per-sample FORMAT fields:

- **`MT`** (top-strand methylation) and **`MB`** (bottom-strand methylation), both `Type=String`, `Number=.`. Per-haplotype values are joined with `|`; each per-haplotype entry is either a bit string of `0`/`1` characters (one bit per CpG owned at this record, in 5'→3' order) or `.` to indicate "this haplotype carries REF" or "no CpGs owned here for this haplotype."

There are two record kinds:

- **Standalone methylation-only records** at reference CpGs outside any variant span: `REF=C`, `ALT=.`, no GT. FORMAT is `MT:MB`. Each per-haplotype entry is a single `0`/`1` character (one CpG per record). Example for a heterozygous methylation state on a diploid sample:

  ```text
  chr1  100200  .  C  .  .  .  .  MT:MB  1|0:0|1
  ```

  Top-C at `chr1:100200` is methylated on the top strand of haplotype 0 and on the bottom strand of haplotype 1.

- **Variant records** with methylation annotated on the variant: standard variant fields plus a per-haplotype bit string per FORMAT field. FORMAT is `GT:MT:MB`. Example for an insertion `T→ACGTACG` that introduces two CpGs on the variant haplotype, both methylated on both strands:

  ```text
  chr1  100500  .  T  ACGTACG  .  .  .  GT:MT:MB  1|0:11|.:11|.
  ```

  Hap0 carries the insertion (`GT=1`), its alt contains two CpGs with top-strand bits `11` and bottom-strand bits `11`. Hap1 carries REF (`GT=0`), so MT and MB are `.`.

  CpGs straddling a variant boundary (one base in the alt, one in the flanking reference) are owned by the variant. When two adjacent variants jointly form a CpG, the upstream variant owns it (deterministic tiebreaker).

#### `simulate` flags for methylation

| Option | Default | Description |
|--------|---------|-------------|
| `--methylation-mode` | none | Methylation chemistry: `em-seq` (or `bisulfite`) or `taps`. Presence enables methylation simulation |
| `--methylation-conversion-rate` | 0.999 | Chemistry efficiency for normally-converting molecules: probability that the converting class of C converts to T |
| `--methylation-failure-rate` | 0.01 | Fraction of molecules that are whole-molecule conversion failures. Requires `--methylation-mode` |
| `--cpg-truth-bedgraph` | none | Write per-CpG ground-truth methylation tally in MethylDackel `extract` bedGraph format. Requires `--methylation-mode` |

Note: the methylation-*rate* flags are **not** `simulate` flags — they live on `holodeck methylate`. At `simulate` time the methylation truth comes exclusively from the `MT`/`MB` fields in the input VCF produced by `holodeck methylate`.

#### Methylation chemistry modes

- **em-seq** (alias `bisulfite`): unmethylated cytosines convert to thymine; methylated cytosines are preserved. Matches both classical bisulfite and enzymatic methyl-seq (em-seq, NEBNext) — the conversion patterns are identical.
- **taps**: methylated cytosines convert to thymine; unmethylated cytosines are preserved. The inverse of em-seq — a `C→T` event at a CpG signals methylation rather than the absence of it.

The "converting class" of cytosines (unmethylated for em-seq, methylated for taps) converts to thymine with probability `--methylation-conversion-rate` (default `0.999`) in a molecule that converted normally. Non-CpG cytosines are always treated as unmethylated.

#### Conversion failure (per-molecule)

Real bisulfite/EM-seq conversion is effectively **bimodal** at the molecule level: most molecules convert near-completely, while a small fraction (fragments that fail to denature, or re-anneal too fast) escape conversion as a coherent unit. The dataset-wide "conversion rate" you see quoted (~0.98–0.99) is the *mean of this mixture*, not a rate any individual molecule sits at.

Holodeck models this with `--methylation-failure-rate` (default `0.01`): each fragment is independently drawn as a *conversion failure* with that probability. A failed molecule converts its should-convert cytosines at `1 − --methylation-conversion-rate` (near-zero), so it coherently retains almost all of them as C — and because the draw happens once per fragment, both mates of a pair agree. Pass `--methylation-failure-rate 0.0` (with `--methylation-conversion-rate 1.0`) to recover perfectly deterministic, lossless conversion.

When `--golden-bam` is set, every record is stamped with `cf:i:{0|1}` recording whether its source molecule was drawn as a failure (see the golden-BAM tag section) — the one piece of conversion-failure truth that is *not* recoverable from the read sequence alone.

#### Behavior under different conversion rates

| Chemistry | Conversion rate | Expected outcome |
|-----------|-----------------|------------------|
| em-seq | 1.0 | Unmethylated CpGs read as T; methylated CpGs read as C |
| em-seq | < 1.0 | Some unmethylated CpGs escape conversion and still read as C |
| taps | 1.0 | Methylated CpGs read as T; unmethylated CpGs read as C |
| taps | < 1.0 | Some methylated CpGs escape conversion and still read as C |

#### `simulate` validation matrix

How `simulate` reacts to the four combinations of input-VCF methylation content and the `--methylation-mode` flag:

| Input VCF has `MT`/`MB`? | `--methylation-mode` set? | Behavior |
|--------------------------|---------------------------|----------|
| yes | yes | Run methylation chemistry simulation against the embedded truth. |
| yes | no  | Log a single warning and continue with variants-only output (no chemistry applied). |
| no  | yes | Hard error: "--methylation-mode requires a methylation-annotated VCF (MT/MB FORMAT fields); run `holodeck methylate` first". |
| no  | no  | Variants-only output (unchanged from non-methylation runs). |

The "warn and continue" case lets you experimentally drop chemistry from a methylated VCF without rebuilding the input file; the "hard error" case catches the common mistake of forgetting to run `holodeck methylate` before `holodeck simulate`.

#### `methylate --bedgraph` vs `simulate --cpg-truth-bedgraph`

Both flags write MethylDackel-compatible CpG bedGraphs, but they represent different quantities:

- **`methylate --bedgraph`**: computed directly from the methylation bitmap stored in the VCF. It reports the closed-form population fraction — the fraction of haplotype × strand draws at each CpG that were set to methylated — without simulating any reads. Use this to inspect the methylation truth before running simulate.
- **`simulate --cpg-truth-bedgraph`**: computed from the reads that actually covered each CpG during simulation. It is coverage-weighted: CpGs covered by more reads contribute more observations. Use this as the ground truth for comparing against a MethylDackel extract run on aligned simulated reads.

Both outputs have value: `methylate --bedgraph` validates the methylation assignment; `simulate --cpg-truth-bedgraph` validates the full simulate-align-extract pipeline.

#### Complete two-command worked example

```bash
holodeck methylate \
  --reference ref.fa \
  --vcf variants.vcf.gz --sample HG00100 \
  --methylation-rate-open-sea 0.7 \
  --seed 42 \
  --output methylated.vcf.gz

holodeck simulate \
  --reference ref.fa \
  --vcf methylated.vcf.gz --sample HG00100 \
  --coverage 30 \
  --methylation-mode em-seq \
  --methylation-conversion-rate 0.99 \
  --output sim
```

When `--cpg-truth-bedgraph <PATH>` is set on `simulate`, holodeck writes a per-CpG ground-truth methylation tally in MethylDackel's `extract` CpG bedGraph format (`chrom  start  end  rate(0–100)  n_methylated  n_unmethylated`, with a `track` header line). For every reference CpG covered by at least one simulated read, it counts how many read mates called the site as methylated vs unmethylated according to the simulator's per-haplotype, per-strand methylation bitmap. Sites are emitted in genomic order. The output is format-identical to MethylDackel's so a downstream concordance script can compare aligner-derived methylation calls (from real or simulated BAMs through MethylDackel) against ground truth using the same code paths.

#### Interpreting the methylation-simulated golden BAM

When `--golden-bam` is enabled together with `--methylation-mode`, holodeck emits a per-record tag set so that downstream methylation tools (`bismark_methylation_extractor`, MethylDackel `extract`, IGV's bisulfite track) read the perfect-truth BAM directly. Tag values are biology-faithful under both em-seq and TAPS — `Z` always means methylated, regardless of chemistry — but the *parity invariants* that fall out (Bismark/MethylDackel agreement with the truth bedGraph) are em-seq only, because those tools assume bisulfite chemistry when re-deriving calls.

##### Tag reference

A non-methylation `--golden-bam` run emits **none** of the tags below; the BAM only carries `RG`, `hp` (haplotype index), and `ne` (number of error events). All eight methylation tags appear together when `--methylation-mode` is set.

| Tag    | Source   | One-line summary |
| ---    | ---      | --- |
| `XG:Z` | Bismark  | Genome-strand indicator (`CT` = top, `GA` = bottom). |
| `XR:Z` | Bismark  | Read-conversion direction in read orientation. |
| `XM:Z` | Bismark  | Observation-derived per-base methylation call string. |
| `YM:Z` | holodeck | Truth-derived per-base methylation call string. |
| `NM:i` | Bismark  | Edit distance against the unconverted reference, chemistry events suppressed. |
| `MD:Z` | Bismark  | Match/mismatch description against the unconverted reference, chemistry events suppressed. |
| `YS:Z` | holodeck | Pre-conversion read sequence in reference-forward orientation. |
| `cf:i` | holodeck | Conversion-failure flag: `1` if the source molecule was drawn as a whole-molecule conversion failure, else `0`. |

###### `XG:Z` — genome-strand indicator

`CT` = read derived from the top genome strand (OT or CTOT in Bismark's four-strand model). `GA` = bottom genome strand (OB or CTOB). Fragment-level: R1 and R2 of a pair carry the same value. Identical between em-seq and TAPS — strand info is chemistry-agnostic.

###### `XR:Z` — read-conversion direction

Conversion direction in *read* (5'→3') orientation. `CT` for R1 / SE reads (the read 5'→3' shows a `C→T` conversion pattern). `GA` for R2 (the PCR-synthesized complement of the source strand, which displays `G→A`). Fixed by mate index under a directional library, independent of `XG`.

###### `XM:Z` — observation-derived methylation call

Per-base methylation call string in `SEQ` orientation, derived from the observed read base versus the unconverted reference. Same length as `SEQ`.

Alphabet: `Z`/`z` for methylated/unmethylated CpG, `X`/`x` for CHG, `H`/`h` for CHH, `.` for any position where no call applies (non-cytosine on the source strand, indels, soft-clips, or an observed nucleotide that's neither the cytosine nor its chemistry-converted form).

The alphabet is biology-faithful: `Z` always means methylated, regardless of chemistry. The mapping from observed base to call therefore inverts between modes:

| Mode   | Preserved C/G reads as | Converted (C→T / G→A) reads as |
| ---    | ---                    | --- |
| em-seq | methylated (`Z`)       | unmethylated (`z`) |
| TAPS   | unmethylated (`z`)     | methylated (`Z`)   |

This matches Bismark's spec for the alphabet but deviates from Bismark's *implementation*, which hard-codes bisulfite when re-deriving calls from `SEQ` vs reference. Consequence: a bisulfite-only extractor run against a TAPS BAM will return the inverse of holodeck's `XM` — see "Invariant checks" below.

###### `YM:Z` — truth-derived methylation call

Holodeck-specific. Same shape and alphabet as `XM:Z`, but the methylated/unmethylated decision comes from the simulator's per-haplotype methylation bitmap rather than from the observed read base. Identical to `XM:Z` for zero-error simulations under either chemistry; diverges where errors corrupt cytosines (e.g. `C → A` at a methylated CpG: `XM = '.'`, `YM = 'Z'`). Holodeck's truth model is CpG-only, so `YM` always emits lowercase `x`/`h` for CHG/CHH calls.

###### `NM:i` — edit distance, chemistry-aware

Edit distance against the **unconverted** reference, with chemistry events suppressed. For `XG:Z:CT` reads, ref-`C` → seq-`T` is treated as a chemistry event and not counted; for `XG:Z:GA` reads, ref-`G` → seq-`A` is treated likewise. Insertions and deletions count one per base. The suppression rule is the same under em-seq and TAPS because both chemistries convert in the same direction (C→T / G→A); only the biological *meaning* of the conversion flips.

###### `MD:Z` — match/mismatch description, chemistry-aware

Match/mismatch description against the unconverted reference, same suppression rule as `NM:i`. Chemistry-allowed mismatches fold into the surrounding match runs. Insertions are not represented (per the MD spec); deletions appear as `^<bases>`.

###### `YS:Z` — pre-conversion read sequence

Holodeck-specific. The pre-chemistry read sequence in reference-forward orientation. Diff `SEQ` against `YS` base-for-base to recover the ground-truth chemistry events the simulator applied to this record. Identical mechanics under em-seq and TAPS; the recovered events carry the mode's biological meaning.

###### `cf:i` — conversion-failure flag

Holodeck-specific. `1` if the source molecule was drawn as a whole-molecule conversion failure (see [Conversion failure](#conversion-failure-per-molecule)), else `0`. A molecule property, so R1 and R2 of a pair always carry the same value. This is the one piece of conversion-failure ground truth that *cannot* be recovered from the read alone: a normally-converted molecule can coincidentally retain cytosines and look failed, and vice versa, so a detector scored on `SEQ`/`YS` needs `cf` as the gold label. The *consequences* of failure are still observable without it — a failed read has `SEQ ≈ YS` and its `XM:Z` diverges from `YM:Z` at retained CpH cytosines.

##### Truth vs observation: when do `XM` and `YM` diverge?

`XM:Z` asks "what does the observed read base say about methylation at this reference cytosine?" That answer is corrupted by sequencing errors — an error that converts a methylated `C` to `A` reads as `XM:Z:.`, indistinguishable from a genuine mismatch.

`YM:Z` instead asks "what *was* the methylation state, according to the simulator's truth bitmap?" It's invariant to errors: a methylated CpG always reports `Z` regardless of what the read base ended up being.

For a zero-error simulation (`--max-error-rate 0`), `XM == YM` byte-for-byte under either chemistry. With a non-zero error rate, the diff `XM != YM` is exactly the set of cytosines whose methylation call would be miscalled by an observation-based extractor — useful for stress-testing aligners and methylation callers under error pressure.

##### Invariant checks the tag set unlocks

Two correctness checks that fall out of running existing tools on the perfect-truth golden BAM. **Both are em-seq only** — Bismark and MethylDackel re-derive methylation calls from `SEQ` vs reference using bisulfite formulas, so they will return inverted calls on a TAPS BAM regardless of what `XM:Z` says. For TAPS, validate against `--cpg-truth-bedgraph` directly using `XM`/`YM` (which are biology-faithful in both modes) or via a TAPS-aware extractor.

1. **Bismark methylation extractor on the em-seq golden BAM ≈ truth bedGraph.** Run `bismark_methylation_extractor --comprehensive --bedGraph` against the golden BAM. The resulting per-CpG bedGraph should match `--cpg-truth-bedgraph` exactly under zero errors. With errors, the divergence is bounded by the rate at which the error model corrupts CpG cytosines.
2. **MethylDackel extract on the em-seq golden BAM ≈ truth bedGraph.** `MethylDackel extract --mergeContext` reads `XG:Z` and emits a CpG bedGraph in the same format as `--cpg-truth-bedgraph`. Identical-comparison invariant under zero errors.

Note: `samtools calmd` is **not** a valid validator for the methylation-tagged golden BAM. It recomputes standard `NM`/`MD` from the unconverted reference and does not preserve the suppression of `C→T` (top) and `G→A` (bottom) events, so running it against a correctly-tagged record will produce spurious diffs.

##### Limitations

- **Directional libraries only** — TruSeq Methyl, NEBNext EM-seq, Twist Methyl, KAPA HyperMeth. R2 is simulated as the PCR-synthesized complement of the converted source strand (`revcomp(c2t(top))` for top-strand-derived fragments), matching what bwameth/Bismark expect from `--directional` libraries. Non-directional / PBAT protocols (scBS-seq, sci-MET, single-cell methods) have different strand semantics and are not supported; use Bismark's `--pbat`/`--non_directional` modes for those.
- **TAPS downstream tooling is non-standard** — Bismark and MethylDackel assume bisulfite chemistry and re-derive methylation calls by treating preserved cytosines as methylated; on a TAPS BAM they will report every CpG with the methylated/unmethylated call inverted. Holodeck's `XM:Z` and `YM:Z` are biology-faithful under both chemistries, but downstream extraction tools must be TAPS-aware (e.g. `asTair`, custom pipelines from the TAPS authors) for the bedGraph output to mean what it says. The Bismark/MethylDackel parity invariants in the previous section apply to em-seq only.
- **CpG context only** — Non-CpG cytosines (CHG, CHH) are always treated as unmethylated in holodeck's truth model, so `YM:Z` emits lowercase `x` and `h` for CHG/CHH calls. `XM:Z` is observation-derived and may still show uppercase `X`/`H` when the sequenced base remains unconverted. If you need realistic non-CpG methylation (plant genomes, embryonic stem cells), holodeck's truth model does not represent it.
- **Context-class methylation, not per-position** — `holodeck methylate` assigns methylation from three CpG context classes (island / shore / open-sea) with a spatial-correlation model, not from per-position rates (e.g. an imported methylation BED or beta values). Methylation is symmetric across strands by default with only sporadic hemimethylation; the strands are not modelled as independent. Allele-specific methylation arises because each haplotype is drawn independently.
- **Variant-driven CpG changes are correctly handled** — SNPs and indels in `--vcf` haplotypes pass through the per-haplotype methylation table by construction. SNPs that create or destroy a CpG, and indels that shift CpG positions, all get the correct chemistry on the appropriate haplotype. (This is a feature, listed here for completeness.)
- **`MD:Z` and `NM:i` are emitted only with `--methylation-mode`** — Without methylation chemistry, holodeck's golden BAM does not emit `NM`/`MD` at all. `samtools calmd` will fill them in with standard semantics if you need them on a non-methylation run.

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
