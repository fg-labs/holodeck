//! Read naming schemes for simulated reads.
//!
//! Supports two modes: **encoded** names that embed truth coordinates
//! (contig, position, strand, haplotype, error count) for downstream
//! evaluation, and **simple** names that are just sequential identifiers.
//!
//! Contig names may contain `:` characters (e.g. HLA contigs). The parser
//! handles this by splitting from the right where the field count is fixed.

/// Truth alignment data for a single read, used to encode position information
/// into the read name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TruthAlignment {
    /// Reference contig name.
    pub contig: String,
    /// 1-based reference start position.
    pub position: u32,
    /// Whether the read is on the forward strand.
    pub is_forward: bool,
    /// 0-based haplotype index.
    pub haplotype: usize,
    /// Number of sequencing errors introduced.
    pub n_errors: usize,
}

/// Format a read name in encoded mode for a paired-end read.
///
/// Format: `holodeck:READ_NUM:CONTIG:POS1+STRAND:POS2+STRAND:HAP:ERRS1:ERRS2`
///
/// Example: `holodeck:42:chr1:10000F:10450R:0:2:1`
///
/// Both R1 and R2 must be on the same contig; `r2.contig` is not encoded
/// separately.
#[must_use]
pub fn encoded_pe_name(read_num: u64, r1: &TruthAlignment, r2: &TruthAlignment) -> String {
    debug_assert_eq!(
        r1.contig, r2.contig,
        "R1 and R2 must be on the same contig for encoded PE names"
    );
    format!(
        "holodeck:{}:{}:{}{}:{}{}:{}:{}:{}",
        read_num,
        r1.contig,
        r1.position,
        strand_char(r1.is_forward),
        r2.position,
        strand_char(r2.is_forward),
        r1.haplotype,
        r1.n_errors,
        r2.n_errors,
    )
}

/// Format a read name in encoded mode for a single-end read.
///
/// Format: `holodeck:READ_NUM:CONTIG:POS+STRAND:HAP:ERRS`
///
/// Example: `holodeck:42:chr1:10000F:0:2`
#[must_use]
pub fn encoded_se_name(read_num: u64, r1: &TruthAlignment) -> String {
    format!(
        "holodeck:{}:{}:{}{}:{}:{}",
        read_num,
        r1.contig,
        r1.position,
        strand_char(r1.is_forward),
        r1.haplotype,
        r1.n_errors,
    )
}

/// Format a simple read name with no truth information.
///
/// Format: `holodeck:READ_NUM`
#[must_use]
pub fn simple_name(read_num: u64) -> String {
    format!("holodeck:{read_num}")
}

/// Return the strand character for a boolean forward flag.
fn strand_char(is_forward: bool) -> char {
    if is_forward { 'F' } else { 'R' }
}

/// Parse an encoded paired-end read name back into truth alignments.
///
/// Parses from the right to handle contig names that contain `:` characters
/// (e.g. HLA contigs like `HLA-A*01:01:01:01`).
///
/// Returns `(read_num, r1_truth, r2_truth)` or `None` if the name doesn't
/// match the expected format.
#[must_use]
pub fn parse_encoded_pe_name(name: &str) -> Option<(u64, TruthAlignment, TruthAlignment)> {
    // Split from the right: the last 5 fields are always
    // pos2+strand, hap, errs_r1, errs_r2, and pos1+strand.
    // The prefix is "holodeck:read_num:contig" (contig may contain colons).
    let mut rev_parts: Vec<&str> = name.rsplitn(6, ':').collect();
    if rev_parts.len() != 6 {
        return None;
    }
    // rsplitn gives parts in reverse order:
    // [n_errors_r2, n_errors_r1, hap, pos2+strand, pos1+strand, "holodeck:read_num:contig"]
    rev_parts.reverse();
    // Now: ["holodeck:read_num:contig", pos1+strand, pos2+strand, hap, n_errors_r1, n_errors_r2]

    let prefix = rev_parts[0];
    let (pos1, fwd1) = parse_pos_strand(rev_parts[1])?;
    let (pos2, fwd2) = parse_pos_strand(rev_parts[2])?;
    let haplotype: usize = rev_parts[3].parse().ok()?;
    let n_errors_r1: usize = rev_parts[4].parse().ok()?;
    let n_errors_r2: usize = rev_parts[5].parse().ok()?;

    // Parse prefix: "holodeck:read_num:contig_with_possible_colons"
    let prefix = prefix.strip_prefix("holodeck:")?;
    let (read_num_str, contig) = prefix.split_once(':')?;
    let read_num: u64 = read_num_str.parse().ok()?;

    let r1 = TruthAlignment {
        contig: contig.to_string(),
        position: pos1,
        is_forward: fwd1,
        haplotype,
        n_errors: n_errors_r1,
    };
    let r2 = TruthAlignment {
        contig: contig.to_string(),
        position: pos2,
        is_forward: fwd2,
        haplotype,
        n_errors: n_errors_r2,
    };

    Some((read_num, r1, r2))
}

/// Parse an encoded single-end read name back into a truth alignment.
///
/// Parses from the right to handle contig names that contain `:` characters.
///
/// Returns `(read_num, truth)` or `None` if the name doesn't match.
#[must_use]
pub fn parse_encoded_se_name(name: &str) -> Option<(u64, TruthAlignment)> {
    // Last 3 fields from the right: pos+strand, hap, n_errors
    // Prefix: "holodeck:read_num:contig"
    let mut rev_parts: Vec<&str> = name.rsplitn(4, ':').collect();
    if rev_parts.len() != 4 {
        return None;
    }
    rev_parts.reverse();
    // ["holodeck:read_num:contig", pos+strand, hap, n_errors]

    let prefix = rev_parts[0];
    let (pos, fwd) = parse_pos_strand(rev_parts[1])?;
    let haplotype: usize = rev_parts[2].parse().ok()?;
    let n_errors: usize = rev_parts[3].parse().ok()?;

    let prefix = prefix.strip_prefix("holodeck:")?;
    let (read_num_str, contig) = prefix.split_once(':')?;
    let read_num: u64 = read_num_str.parse().ok()?;

    let truth = TruthAlignment {
        contig: contig.to_string(),
        position: pos,
        is_forward: fwd,
        haplotype,
        n_errors,
    };

    Some((read_num, truth))
}

/// Parse a position+strand field like "10000F" or "450R".
fn parse_pos_strand(field: &str) -> Option<(u32, bool)> {
    if field.is_empty() {
        return None;
    }
    let strand_char = field.as_bytes()[field.len() - 1];
    let is_forward = match strand_char {
        b'F' => true,
        b'R' => false,
        _ => return None,
    };
    let pos: u32 = field[..field.len() - 1].parse().ok()?;
    Some((pos, is_forward))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_truth(contig: &str, pos: u32, fwd: bool, hap: usize, errs: usize) -> TruthAlignment {
        TruthAlignment {
            contig: contig.to_string(),
            position: pos,
            is_forward: fwd,
            haplotype: hap,
            n_errors: errs,
        }
    }

    #[test]
    fn test_encoded_pe_name() {
        let r1 = make_truth("chr1", 10000, true, 0, 2);
        let r2 = make_truth("chr1", 10450, false, 0, 1);
        let name = encoded_pe_name(42, &r1, &r2);
        assert_eq!(name, "holodeck:42:chr1:10000F:10450R:0:2:1");
    }

    #[test]
    fn test_encoded_se_name() {
        let r1 = make_truth("chr1", 10000, true, 0, 2);
        let name = encoded_se_name(42, &r1);
        assert_eq!(name, "holodeck:42:chr1:10000F:0:2");
    }

    #[test]
    fn test_simple_name() {
        assert_eq!(simple_name(42), "holodeck:42");
        assert_eq!(simple_name(1), "holodeck:1");
    }

    #[test]
    fn test_parse_pe_roundtrip() {
        let r1 = make_truth("chr1", 10000, true, 0, 2);
        let r2 = make_truth("chr1", 10450, false, 0, 1);
        let name = encoded_pe_name(42, &r1, &r2);

        let (num, parsed_r1, parsed_r2) = parse_encoded_pe_name(&name).unwrap();
        assert_eq!(num, 42);
        assert_eq!(parsed_r1, r1);
        assert_eq!(parsed_r2, r2);
    }

    #[test]
    fn test_parse_se_roundtrip() {
        let r1 = make_truth("chrX", 500, false, 1, 0);
        let name = encoded_se_name(99, &r1);

        let (num, parsed) = parse_encoded_se_name(&name).unwrap();
        assert_eq!(num, 99);
        assert_eq!(parsed, r1);
    }

    #[test]
    fn test_parse_pe_with_colon_in_contig() {
        // HLA contig names contain colons.
        let r1 = make_truth("HLA-A*01:01:01:01", 100, true, 0, 0);
        let r2 = make_truth("HLA-A*01:01:01:01", 400, false, 0, 1);
        let name = encoded_pe_name(1, &r1, &r2);

        let (num, parsed_r1, parsed_r2) = parse_encoded_pe_name(&name).unwrap();
        assert_eq!(num, 1);
        assert_eq!(parsed_r1.contig, "HLA-A*01:01:01:01");
        assert_eq!(parsed_r2.contig, "HLA-A*01:01:01:01");
        assert_eq!(parsed_r1.position, 100);
        assert_eq!(parsed_r2.position, 400);
    }

    #[test]
    fn test_parse_se_with_colon_in_contig() {
        let r1 = make_truth("HLA-B*07:02", 50, false, 1, 3);
        let name = encoded_se_name(5, &r1);

        let (num, parsed) = parse_encoded_se_name(&name).unwrap();
        assert_eq!(num, 5);
        assert_eq!(parsed, r1);
    }

    #[test]
    fn test_cross_format_rejection() {
        // SE name should not parse as PE.
        let se_name = encoded_se_name(1, &make_truth("chr1", 10, true, 0, 0));
        assert!(parse_encoded_pe_name(&se_name).is_none());

        // PE name should not parse as SE (wrong field structure from the right).
        let pe_name = encoded_pe_name(
            1,
            &make_truth("chr1", 10, true, 0, 0),
            &make_truth("chr1", 20, false, 0, 0),
        );
        // SE parser expects pos+strand as the 3rd-from-right field, but for PE
        // that field is the haplotype index (a plain number, no F/R suffix).
        assert!(parse_encoded_se_name(&pe_name).is_none());
    }

    #[test]
    fn test_parse_invalid_names() {
        assert!(parse_encoded_pe_name("not_holodeck:1:chr1:10F:20R:0:0:0").is_none());
        assert!(parse_encoded_pe_name("").is_none());
        assert!(parse_encoded_se_name("holodeck:1").is_none());
    }

    #[test]
    fn test_parse_pos_strand() {
        assert_eq!(parse_pos_strand("10000F"), Some((10000, true)));
        assert_eq!(parse_pos_strand("450R"), Some((450, false)));
        assert_eq!(parse_pos_strand("0F"), Some((0, true)));
        assert_eq!(parse_pos_strand("X"), None);
        assert_eq!(parse_pos_strand(""), None);
    }
}
