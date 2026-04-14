use std::path::PathBuf;

/// Help heading for core options shared by most commands.
pub const CORE_OPTIONS: &str = "Core Options";

/// Options for specifying the reference FASTA file.
#[derive(clap::Args, Debug, Clone)]
pub struct ReferenceOptions {
    /// Reference FASTA file (must be indexed with .fai).
    #[arg(short = 'r', long, value_name = "FASTA")]
    pub reference: PathBuf,
}

/// Options for specifying the output path prefix (for commands that produce
/// multiple output files with shared prefix, e.g. `simulate`).
#[derive(clap::Args, Debug, Clone)]
pub struct OutputPrefixOptions {
    /// Output path prefix for generated files.
    #[arg(short = 'o', long, value_name = "PREFIX")]
    pub output: PathBuf,
}

/// Options for specifying an input VCF file and optional sample selection.
#[derive(clap::Args, Debug, Clone)]
pub struct VcfOptions {
    /// VCF file with variants to apply. Must have GT field in genotypes.
    #[arg(short = 'v', long, value_name = "VCF")]
    pub vcf: Option<PathBuf>,

    /// Sample name to select from a multi-sample VCF. Required if the VCF
    /// contains more than one sample; ignored for single-sample VCFs.
    #[arg(long, value_name = "SAMPLE")]
    pub sample: Option<String>,
}

/// Options for specifying a BED file of target regions.
#[derive(clap::Args, Debug, Clone)]
pub struct BedOptions {
    /// BED file of target regions. Only fragments overlapping these regions
    /// will be emitted.
    #[arg(short = 'b', long = "targets", value_name = "BED")]
    pub targets: Option<PathBuf>,
}

/// Options for controlling the random seed.
#[derive(clap::Args, Debug, Clone)]
pub struct SeedOptions {
    /// Random seed for reproducibility. If omitted, a deterministic seed is
    /// derived from all other options so that identical parameters produce
    /// identical output.
    #[arg(long, value_name = "SEED")]
    pub seed: Option<u64>,
}
