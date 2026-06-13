# Holodeck - Claude Code Guide

## Build & Test Commands

```bash
cargo build                    # Debug build
cargo build --release          # Release build
cargo ci-test                  # Run tests via nextest
cargo ci-fmt                   # Check formatting
cargo ci-lint                  # Run clippy with pedantic warnings
cargo fmt                      # Auto-format code
```

## Architecture

Holodeck is a single-crate Rust project with a binary (`holodeck`) and library (`holodeck_lib`).

### Subcommands

- **simulate** -- Generate reads from reference + optional VCF. Core simulation engine. Optional methylation chemistry via `--methylation-mode`.
- **mutate** -- Generate a random VCF from a reference. Independent of simulator.
- **methylate** -- Annotate a VCF with per-CpG methylation truth (MT/MB FORMAT fields) for downstream `simulate --vcf` consumption.
- **eval** -- Evaluate alignment accuracy by comparing truth vs mapped positions.

### Module Overview

| Module | Purpose |
|--------|---------|
| `commands/simulate.rs` | Full simulation pipeline: load ref/VCF/BED, build haplotypes, sample fragments, generate reads |
| `commands/mutate.rs` | Random VCF generation with SNP/indel/MNP rates and ploidy overrides |
| `commands/methylate.rs` | Annotate a VCF with per-CpG MT/MB methylation truth |
| `commands/eval.rs` | Alignment accuracy evaluation from encoded read names |
| `commands/common.rs` | Shared CLI option groups (reference, output, VCF, BED, seed) |
| `bed.rs` | BED file loading with `coitrees` for overlap queries |
| `vcf/mod.rs` | VCF reading with noodles, sample selection |
| `vcf/genotype.rs` | GT field parsing supporting arbitrary ploidy and phasing |
| `haplotype.rs` | Sparse haplotype variant overlay (reference + COITree of variants) |
| `meth.rs` | Methylation primitives: per-haplotype/per-strand bitmaps, chemistry conversion, conversion-failure model, reference CpG scanning |
| `methylation_tags.rs` | Bismark-style methylation call tags (XM/YM/NM/MD) for the golden BAM |
| `vcf/methylation.rs` | MT/MB FORMAT field read/write; CpG classification across haplotypes |
| `output/cpg_truth.rs` | Coverage-weighted per-CpG truth bedGraph (from simulated reads) |
| `output/methylation_bedgraph.rs` | Closed-form population-fraction methylation bedGraph |
| `fragment.rs` | Fragment extraction, reverse complement, adapter padding |
| `read.rs` | Read pair generation combining fragments + error model + naming |
| `error_model/mod.rs` | ErrorModel trait + apply_errors free function |
| `error_model/illumina.rs` | Position-dependent Illumina error model with precomputed lookup tables |
| `read_naming.rs` | Encoded and simple read name formatting + parsing |
| `ploidy.rs` | PloidyMap with per-contig/per-region overrides |
| `seed.rs` | Deterministic FNV-1a seed computation |
| `sequence_dict.rs` | Sequence dictionary (name/index/length lookups) |
| `fasta.rs` | Indexed FASTA reader |
| `output/fastq.rs` | BGZF-compressed FASTQ writer (single-threaded or pooled) |
| `output/golden_bam.rs` | Ground-truth BAM writer (single-threaded or pooled) |

### Key Patterns

- **CLI**: `clap` derive with styled help. Shared option groups via `#[command(flatten)]`.
- **Commands**: `Command` trait with `execute()` via `enum_dispatch`.
- **Intervals**: `coitrees` crate for all overlap queries (BED targets, variant lookup).
- **Haplotypes**: Sparse variant overlay on reference -- NOT full sequence copies. Variants stored in COITree with index metadata (because COITree requires `Copy + Default`).
- **Output**: Multi-threaded BGZF compression via `pooled-writer` when `--threads > 1`; single-threaded via `noodles-bgzf` otherwise. Both FASTQ and BAM writers accept `Box<dyn Write>`.
- **Errors**: `anyhow` for application errors, `thiserror` for library error types.
- **Allocator**: `mimalloc` as global allocator.
- **RNG**: `SmallRng` (Xoshiro256++) with deterministic FNV-1a seed for reproducibility.
- **Error model**: Per-cycle error probabilities and base quality scores are precomputed into lookup tables at model construction time; the per-base hot loop does table lookups, not floating-point math.

### Conventions

- `#![deny(unsafe_code)]` in both lib.rs and main.rs.
- Doc comments on all public and non-trivial private items.
- Module-level `//!` documentation on all modules.
- Generate all test data programmatically -- never commit test data files.
- Many small individual tests over parameterized/table-driven tests.
- Encoded read names use `::` (double colon) as the field separator so contig names may legally contain single `:` characters (e.g. HLA alleles).
- Contig names must not contain `@` (FASTQ header prefix) or `::` (field separator); both are rejected by a debug assertion in the read-name formatter.
