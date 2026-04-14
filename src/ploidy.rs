use anyhow::{Result, bail};

/// A ploidy entry specifying the ploidy for a genomic region.
#[derive(Debug, Clone)]
struct PloidyEntry {
    /// Contig name (e.g. "chrX"). None means all contigs.
    contig: Option<String>,
    /// Genomic range within the contig. None means the entire contig.
    /// Coordinates are 0-based half-open `[start, end)`.
    range: Option<(u32, u32)>,
    /// Ploidy value.
    ploidy: u8,
}

/// Maps genomic positions to ploidy values, supporting a global default with
/// per-contig and per-region overrides.
///
/// Overrides are processed in order and use last-writer-wins semantics: later
/// entries override earlier ones for overlapping positions. This allows
/// patterns like:
///
/// ```text
/// --ploidy 2 --ploidy-override chrX=1 --ploidy-override chrX:10001-2781479=2
/// ```
///
/// which sets diploid genome-wide, then haploid chrX, then diploid PAR1 on
/// chrX.
#[derive(Debug, Clone)]
pub struct PloidyMap {
    /// Default ploidy for positions not covered by any override.
    default_ploidy: u8,
    /// Ordered list of overrides; later entries take precedence.
    overrides: Vec<PloidyEntry>,
}

impl PloidyMap {
    /// Create a new ploidy map with the given default and overrides.
    ///
    /// Override strings have the format `CONTIG=PLOIDY` (whole contig) or
    /// `CONTIG:START-END=PLOIDY` (1-based inclusive coordinates, converted
    /// internally to 0-based half-open).
    ///
    /// # Errors
    /// Returns an error if any override string is malformed.
    pub fn new(default_ploidy: u8, override_specs: &[String]) -> Result<Self> {
        let mut overrides = Vec::with_capacity(override_specs.len());

        for spec in override_specs {
            let entry = parse_ploidy_override(spec)?;
            overrides.push(entry);
        }

        Ok(Self { default_ploidy, overrides })
    }

    /// Look up the ploidy for a specific contig and 0-based position.
    ///
    /// Returns the last matching override, or the default if no override
    /// applies.
    #[must_use]
    pub fn ploidy_at(&self, contig: &str, position: u32) -> u8 {
        let mut result = self.default_ploidy;

        for entry in &self.overrides {
            let contig_matches = entry.contig.as_ref().is_none_or(|c| c == contig);

            if !contig_matches {
                continue;
            }

            let range_matches = entry
                .range
                .as_ref()
                .is_none_or(|&(start, end)| position >= start && position < end);

            if range_matches {
                result = entry.ploidy;
            }
        }

        result
    }

    /// Return the default ploidy.
    #[must_use]
    pub fn default_ploidy(&self) -> u8 {
        self.default_ploidy
    }
}

/// Parse a single ploidy override specification.
///
/// Formats:
///   - `CONTIG=PLOIDY` — applies to entire contig
///   - `CONTIG:START-END=PLOIDY` — applies to 1-based inclusive range
fn parse_ploidy_override(spec: &str) -> Result<PloidyEntry> {
    let (location, ploidy_str) = spec
        .rsplit_once('=')
        .ok_or_else(|| anyhow::anyhow!("Invalid ploidy override format: '{spec}'. Expected CONTIG=PLOIDY or CONTIG:START-END=PLOIDY"))?;

    let ploidy: u8 = ploidy_str
        .parse()
        .map_err(|_| anyhow::anyhow!("Invalid ploidy value in override: '{ploidy_str}'"))?;

    if let Some((contig, range_str)) = location.split_once(':') {
        // CONTIG:START-END format (1-based inclusive input).
        let (start_str, end_str) = range_str.split_once('-').ok_or_else(|| {
            anyhow::anyhow!(
                "Invalid range format in ploidy override: '{range_str}'. Expected START-END"
            )
        })?;

        let start_1based: u32 = start_str.parse().map_err(|_| {
            anyhow::anyhow!("Invalid start position in ploidy override: '{start_str}'")
        })?;
        let end_1based: u32 = end_str
            .parse()
            .map_err(|_| anyhow::anyhow!("Invalid end position in ploidy override: '{end_str}'"))?;

        if start_1based == 0 {
            bail!("Ploidy override start position must be >= 1: '{spec}'");
        }
        if end_1based < start_1based {
            bail!("Ploidy override end < start: '{spec}'");
        }

        // Convert to 0-based half-open.
        let start_0based = start_1based - 1;
        let end_0based = end_1based; // 1-based inclusive end == 0-based exclusive end

        Ok(PloidyEntry {
            contig: Some(contig.to_string()),
            range: Some((start_0based, end_0based)),
            ploidy,
        })
    } else {
        // CONTIG=PLOIDY format (whole contig).
        Ok(PloidyEntry { contig: Some(location.to_string()), range: None, ploidy })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_ploidy() {
        let map = PloidyMap::new(2, &[]).unwrap();
        assert_eq!(map.ploidy_at("chr1", 0), 2);
        assert_eq!(map.ploidy_at("chrX", 1000), 2);
    }

    #[test]
    fn test_whole_contig_override() {
        let map = PloidyMap::new(2, &["chrX=1".to_string()]).unwrap();
        assert_eq!(map.ploidy_at("chr1", 0), 2);
        assert_eq!(map.ploidy_at("chrX", 0), 1);
        assert_eq!(map.ploidy_at("chrX", 999_999), 1);
    }

    #[test]
    fn test_range_override() {
        let overrides = vec![
            "chrX=1".to_string(),
            "chrX:10001-2781479=2".to_string(), // PAR1
        ];
        let map = PloidyMap::new(2, &overrides).unwrap();

        // Autosome: default 2
        assert_eq!(map.ploidy_at("chr1", 100), 2);
        // chrX before PAR1: haploid
        assert_eq!(map.ploidy_at("chrX", 5000), 1);
        // chrX in PAR1 (1-based 10001-2781479 → 0-based 10000-2781479): diploid
        assert_eq!(map.ploidy_at("chrX", 10000), 2);
        assert_eq!(map.ploidy_at("chrX", 2_781_478), 2);
        // chrX just past PAR1: haploid again
        assert_eq!(map.ploidy_at("chrX", 2_781_479), 1);
    }

    #[test]
    fn test_last_writer_wins() {
        // Two overlapping overrides for chrX.
        let overrides = vec!["chrX=1".to_string(), "chrX=3".to_string()];
        let map = PloidyMap::new(2, &overrides).unwrap();
        // Last writer (3) wins.
        assert_eq!(map.ploidy_at("chrX", 0), 3);
    }

    #[test]
    fn test_parse_errors() {
        // Missing equals sign.
        assert!(PloidyMap::new(2, &["chrX".to_string()]).is_err());
        // Non-numeric ploidy.
        assert!(PloidyMap::new(2, &["chrX=abc".to_string()]).is_err());
        // Bad range format.
        assert!(PloidyMap::new(2, &["chrX:100=1".to_string()]).is_err());
        // Zero start position.
        assert!(PloidyMap::new(2, &["chrX:0-100=1".to_string()]).is_err());
        // End < start.
        assert!(PloidyMap::new(2, &["chrX:200-100=1".to_string()]).is_err());
    }

    #[test]
    fn test_multiple_contigs_and_ranges() {
        let overrides = vec![
            "chrX=1".to_string(),
            "chrY=1".to_string(),
            "chrX:10001-2781479=2".to_string(),
            "chrX:155701383-156030895=2".to_string(),
        ];
        let map = PloidyMap::new(2, &overrides).unwrap();

        assert_eq!(map.ploidy_at("chr1", 100), 2);
        assert_eq!(map.ploidy_at("chrX", 5000), 1);
        assert_eq!(map.ploidy_at("chrX", 50000), 2); // PAR1
        assert_eq!(map.ploidy_at("chrX", 100_000_000), 1); // middle of chrX
        assert_eq!(map.ploidy_at("chrX", 155_800_000), 2); // PAR2
        assert_eq!(map.ploidy_at("chrY", 100), 1);
    }
}
