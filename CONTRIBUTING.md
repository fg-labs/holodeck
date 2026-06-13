# Contributing to Holodeck

## Getting Started

```bash
git clone https://github.com/fg-labs/holodeck.git
cd holodeck
cargo build
cargo ci-test
```

Requires Rust 1.94.0 or later (edition 2024).

## Development Workflow

### Build and Test

```bash
cargo build                    # Debug build
cargo build --release          # Release build
cargo ci-test                  # Run tests via nextest
cargo ci-fmt                   # Check formatting
cargo ci-lint                  # Run clippy with pedantic warnings
cargo fmt                      # Auto-format code
```

All three checks (`ci-test`, `ci-fmt`, `ci-lint`) must pass before submitting a PR.  These same checks run in CI on every push and pull request.

### Running a Single Test

```bash
cargo test --test test_simulate -- test_name_here
cargo test --lib -- module::tests::test_name
```

## Code Conventions

### General

- `#![deny(unsafe_code)]` in both `lib.rs` and `main.rs` -- no unsafe code.
- Doc comments (`///`) on all public items and non-trivial private items.
- Module-level documentation (`//!`) on all modules.
- Comments explain **why**, not what.  No `// increment x` above `x += 1`.
- Prefer small/medium functions with clear inputs and outputs over large monolithic functions.

### Error Handling

- `anyhow::Result` for application-level errors (commands, I/O).
- `thiserror` for library error types where callers need to match on variants.
- No `unwrap()` in production code paths; `expect()` only where invariants are documented.

### CLI

- `clap` derive API with `#[command(flatten)]` for shared option groups.
- Styled help text via `clap::builder::styling`.
- Shared options in `commands/common.rs`: `ReferenceOptions`, `OutputPrefixOptions`, `VcfOptions`, `BedOptions`, `SeedOptions`.

### Testing

- **Generate all test data programmatically** -- never commit test data files.  Build references, VCFs, BEDs, and BAMs inline so they're visible to reviewers.
- **Many small individual tests** over parameterized/table-driven tests.
- **Test function, not implementation** -- tests should survive a significant refactor.
- Test expected results, error conditions, and boundary cases.
- Integration tests in `tests/` exercise the compiled binary end-to-end via `run_simulate()`, `run_mutate()`, and `run_eval()` helpers in `tests/helpers/mod.rs`.
- Unit tests in `#[cfg(test)]` modules within each source file.

### Naming

- Encoded read names use `::` (double colon) as the field separator so contig names may legally contain single `:` characters (e.g. HLA alleles).
- Contig names must not contain `@` (FASTQ header delimiter) or `::` (field separator); both are rejected by a debug assertion in the read-name formatter.

## Architecture

Holodeck is a single-crate project with a binary (`holodeck`) and library (`holodeck_lib`).

### Module Overview

| Module | Purpose |
|--------|---------|
| `commands/simulate.rs` | Simulation pipeline: load ref/VCF/BED, build haplotypes, sample fragments, generate reads |
| `commands/mutate.rs` | Random VCF generation with configurable rates and ploidy |
| `commands/eval.rs` | Alignment accuracy evaluation from encoded read names |
| `commands/common.rs` | Shared CLI option groups |
| `haplotype.rs` | Sparse haplotype variant overlay (reference + COITree) |
| `fragment.rs` | Fragment extraction, reverse complement, adapter padding |
| `read.rs` | Read pair generation: fragments + error model + naming + CIGARs |
| `error_model/` | ErrorModel trait + Illumina position-dependent implementation |
| `vcf/` | VCF reading, sample selection, genotype parsing |
| `bed.rs` | BED file loading, overlap queries, padded interval sampling |
| `read_naming.rs` | Encoded/simple read name formatting and parsing |
| `ploidy.rs` | PloidyMap with per-contig/per-region overrides |
| `output/fastq.rs` | BGZF-compressed FASTQ writer |
| `output/golden_bam.rs` | Ground-truth BAM writer |

### Key Design Decisions

- **Sparse haplotypes**: Variants are stored in `COITree` interval trees overlaid on the reference, not as full-sequence copies.  Fragment extraction walks the reference and substitutes alt alleles on-the-fly.
- **Padded interval sampling**: For targeted sequencing, fragment start positions are drawn from padded target regions rather than rejection-sampled across the whole contig.
- **Multi-threaded compression**: BGZF compression is offloaded to a `pooled-writer` thread pool shared across all output files (FASTQ R1, R2, and golden BAM).
- **Precomputed error model**: Per-cycle error probabilities and quality scores are computed once at model construction and stored in lookup tables.
- **Deterministic seeding**: When no explicit seed is given, a deterministic seed is derived from all simulation parameters via FNV-1a hashing, so identical parameters always produce identical output.

## CI

GitHub Actions runs these jobs on every push to `main` and on pull requests:

1. **test** -- `cargo ci-test` (nextest)
2. **lint** -- `cargo ci-lint` (clippy with pedantic warnings)
3. **format** -- `cargo ci-fmt` (rustfmt check)
4. **security** -- `cargo audit` (advisory-database scan)

## Releasing

Releases are cut **locally** from `main` with [`cargo release`](https://github.com/crate-ci/cargo-release)
(`cargo install cargo-release`).  Configuration lives in `release.toml`.

1. Make sure `main` is up to date, green in CI, and `[Unreleased]` in
   `CHANGELOG.md` describes everything in the release.
2. Preview the release (this is the default -- no changes are made):

   ```bash
   cargo release patch    # or: minor / major
   ```

3. When the preview looks right, execute it:

   ```bash
   cargo release patch --execute
   ```

   This runs the full gate (`cargo ci-fmt && cargo ci-lint && cargo ci-test`),
   bumps the version in `Cargo.toml`/`Cargo.lock`, rolls the `[Unreleased]`
   CHANGELOG section into a dated version section, commits, creates an annotated
   `vX.Y.Z` tag, pushes the commit and tag to `origin`, and publishes to
   crates.io.

4. Create the GitHub release from the new tag (not yet automated):

   ```bash
   gh release create vX.Y.Z --generate-notes
   ```

Publishing to crates.io requires a one-time `cargo login` with a token from
<https://crates.io/me>.
