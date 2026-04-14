use std::path::{Path, PathBuf};

use anyhow::Result;

/// Trait implemented by every holodeck subcommand.
#[enum_dispatch::enum_dispatch]
pub trait Command {
    /// Execute the command.
    ///
    /// # Errors
    /// Returns an error if the command fails.
    fn execute(&self) -> Result<()>;
}

/// Build an output path by appending a suffix to a prefix path.
///
/// E.g. `output_path(Path::new("out/sample"), ".r1.fastq.gz")` → `out/sample.r1.fastq.gz`
#[must_use]
pub fn output_path(prefix: &Path, suffix: &str) -> PathBuf {
    PathBuf::from(format!("{}{suffix}", prefix.to_string_lossy()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_output_path_simple() {
        let result = output_path(Path::new("out/sample"), ".r1.fastq.gz");
        assert_eq!(result, PathBuf::from("out/sample.r1.fastq.gz"));
    }

    #[test]
    fn test_output_path_no_directory() {
        let result = output_path(Path::new("prefix"), ".txt");
        assert_eq!(result, PathBuf::from("prefix.txt"));
    }
}
