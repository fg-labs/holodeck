//! Read naming schemes for simulated reads.
//!
//! Supports two modes: **encoded** names that embed truth coordinates
//! (contig, position, strand, haplotype, error count) plus the source
//! fragment length for downstream evaluation, and **simple** names that are
//! just sequential identifiers.
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
/// Format: `holodeck:READ_NUM:FRAG_LEN:CONTIG:POS1+STRAND:POS2+STRAND:HAP:ERRS1:ERRS2`
///
/// Example: `holodeck:42:450:chr1:10000F:10450R:0:2:1`
///
/// `FRAG_LEN` is the length in bases of the source fragment. When it is
/// shorter than the read length, the remainder of each read is adapter (plus
/// `N`-padding if the configured adapter is too short).
///
/// Both R1 and R2 must be on the same contig; `r2.contig` is not encoded
/// separately.
#[must_use]
pub fn encoded_pe_name(
    read_num: u64,
    fragment_length: u32,
    r1: &TruthAlignment,
    r2: &TruthAlignment,
) -> String {
    debug_assert_eq!(
        r1.contig, r2.contig,
        "R1 and R2 must be on the same contig for encoded PE names"
    );
    format!(
        "holodeck:{}:{}:{}:{}{}:{}{}:{}:{}:{}",
        read_num,
        fragment_length,
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
/// Format: `holodeck:READ_NUM:FRAG_LEN:CONTIG:POS+STRAND:HAP:ERRS`
///
/// Example: `holodeck:42:275:chr1:10000F:0:2`
///
/// `FRAG_LEN` is the length in bases of the source fragment. When it is
/// shorter than the read length, the remainder of the read is adapter (plus
/// `N`-padding if the configured adapter is too short).
#[must_use]
pub fn encoded_se_name(read_num: u64, fragment_length: u32, r1: &TruthAlignment) -> String {
    format!(
        "holodeck:{}:{}:{}:{}{}:{}:{}",
        read_num,
        fragment_length,
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
/// Returns `(read_num, fragment_length, r1_truth, r2_truth)` or `None` if the
/// name doesn't match the expected format.
#[must_use]
pub fn parse_encoded_pe_name(name: &str) -> Option<(u64, u32, TruthAlignment, TruthAlignment)> {
    // Split from the right: the last 5 fields are always
    // pos1+strand, pos2+strand, hap, errs_r1, errs_r2.
    // The prefix is "holodeck:read_num:frag_len:contig" (contig may contain colons).
    let mut rev_parts: Vec<&str> = name.rsplitn(6, ':').collect();
    if rev_parts.len() != 6 {
        return None;
    }
    rev_parts.reverse();
    // ["holodeck:read_num:frag_len:contig", pos1+strand, pos2+strand, hap, n_errors_r1, n_errors_r2]

    let prefix = rev_parts[0];
    let (pos1, fwd1) = parse_pos_strand(rev_parts[1])?;
    let (pos2, fwd2) = parse_pos_strand(rev_parts[2])?;
    let haplotype: usize = rev_parts[3].parse().ok()?;
    let n_errors_r1: usize = rev_parts[4].parse().ok()?;
    let n_errors_r2: usize = rev_parts[5].parse().ok()?;

    let (read_num, fragment_length, contig) = parse_prefix(prefix)?;

    let r1 = TruthAlignment {
        contig: contig.clone(),
        position: pos1,
        is_forward: fwd1,
        haplotype,
        n_errors: n_errors_r1,
    };
    let r2 = TruthAlignment {
        contig,
        position: pos2,
        is_forward: fwd2,
        haplotype,
        n_errors: n_errors_r2,
    };

    Some((read_num, fragment_length, r1, r2))
}

/// Parse an encoded single-end read name back into a truth alignment.
///
/// Parses from the right to handle contig names that contain `:` characters.
///
/// Returns `(read_num, fragment_length, truth)` or `None` if the name doesn't
/// match.
#[must_use]
pub fn parse_encoded_se_name(name: &str) -> Option<(u64, u32, TruthAlignment)> {
    // Last 3 fields from the right: pos+strand, hap, n_errors.
    // Prefix: "holodeck:read_num:frag_len:contig".
    let mut rev_parts: Vec<&str> = name.rsplitn(4, ':').collect();
    if rev_parts.len() != 4 {
        return None;
    }
    rev_parts.reverse();
    // ["holodeck:read_num:frag_len:contig", pos+strand, hap, n_errors]

    let prefix = rev_parts[0];
    let (pos, fwd) = parse_pos_strand(rev_parts[1])?;
    let haplotype: usize = rev_parts[2].parse().ok()?;
    let n_errors: usize = rev_parts[3].parse().ok()?;

    let (read_num, fragment_length, contig) = parse_prefix(prefix)?;

    let truth = TruthAlignment { contig, position: pos, is_forward: fwd, haplotype, n_errors };

    Some((read_num, fragment_length, truth))
}

/// Parse the fixed prefix `"holodeck:READ_NUM:FRAG_LEN:CONTIG"` into its
/// components. The contig may contain `:` characters, so everything after the
/// third `:` is treated as the contig name.
fn parse_prefix(prefix: &str) -> Option<(u64, u32, String)> {
    let rest = prefix.strip_prefix("holodeck:")?;
    let (read_num_str, rest) = rest.split_once(':')?;
    let (frag_len_str, contig) = rest.split_once(':')?;
    let read_num: u64 = read_num_str.parse().ok()?;
    let fragment_length: u32 = frag_len_str.parse().ok()?;
    if contig.is_empty() {
        return None;
    }
    Some((read_num, fragment_length, contig.to_string()))
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
        let name = encoded_pe_name(42, 600, &r1, &r2);
        assert_eq!(name, "holodeck:42:600:chr1:10000F:10450R:0:2:1");
    }

    #[test]
    fn test_encoded_se_name() {
        let r1 = make_truth("chr1", 10000, true, 0, 2);
        let name = encoded_se_name(42, 275, &r1);
        assert_eq!(name, "holodeck:42:275:chr1:10000F:0:2");
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
        let name = encoded_pe_name(42, 600, &r1, &r2);

        let (num, frag_len, parsed_r1, parsed_r2) = parse_encoded_pe_name(&name).unwrap();
        assert_eq!(num, 42);
        assert_eq!(frag_len, 600);
        assert_eq!(parsed_r1, r1);
        assert_eq!(parsed_r2, r2);
    }

    #[test]
    fn test_parse_se_roundtrip() {
        let r1 = make_truth("chrX", 500, false, 1, 0);
        let name = encoded_se_name(99, 275, &r1);

        let (num, frag_len, parsed) = parse_encoded_se_name(&name).unwrap();
        assert_eq!(num, 99);
        assert_eq!(frag_len, 275);
        assert_eq!(parsed, r1);
    }

    #[test]
    fn test_parse_pe_with_colon_in_contig() {
        // HLA contig names contain colons.
        let r1 = make_truth("HLA-A*01:01:01:01", 100, true, 0, 0);
        let r2 = make_truth("HLA-A*01:01:01:01", 400, false, 0, 1);
        let name = encoded_pe_name(1, 450, &r1, &r2);

        let (num, frag_len, parsed_r1, parsed_r2) = parse_encoded_pe_name(&name).unwrap();
        assert_eq!(num, 1);
        assert_eq!(frag_len, 450);
        assert_eq!(parsed_r1.contig, "HLA-A*01:01:01:01");
        assert_eq!(parsed_r2.contig, "HLA-A*01:01:01:01");
        assert_eq!(parsed_r1.position, 100);
        assert_eq!(parsed_r2.position, 400);
    }

    #[test]
    fn test_parse_se_with_colon_in_contig() {
        let r1 = make_truth("HLA-B*07:02", 50, false, 1, 3);
        let name = encoded_se_name(5, 120, &r1);

        let (num, frag_len, parsed) = parse_encoded_se_name(&name).unwrap();
        assert_eq!(num, 5);
        assert_eq!(frag_len, 120);
        assert_eq!(parsed, r1);
    }

    #[test]
    fn test_short_fragment_encodes_adapter_boundary() {
        // A 40-base fragment read at 150bp read length -- the trailing 110 bases
        // are adapter. Consumers can recover this from fragment_length alone.
        let r1 = make_truth("chr1", 10000, true, 0, 0);
        let r2 = make_truth("chr1", 10000, false, 0, 0);
        let name = encoded_pe_name(7, 40, &r1, &r2);

        let (_num, frag_len, _r1, _r2) = parse_encoded_pe_name(&name).unwrap();
        assert_eq!(frag_len, 40);
    }

    #[test]
    fn test_cross_format_rejection() {
        // SE name should not parse as PE.
        let se_name = encoded_se_name(1, 150, &make_truth("chr1", 10, true, 0, 0));
        assert!(parse_encoded_pe_name(&se_name).is_none());

        // PE name should not parse as SE (wrong field structure from the right).
        let pe_name = encoded_pe_name(
            1,
            200,
            &make_truth("chr1", 10, true, 0, 0),
            &make_truth("chr1", 20, false, 0, 0),
        );
        // SE parser expects pos+strand as the 3rd-from-right field, but for PE
        // that field is the haplotype index (a plain number, no F/R suffix).
        assert!(parse_encoded_se_name(&pe_name).is_none());
    }

    #[test]
    fn test_parse_invalid_names() {
        assert!(parse_encoded_pe_name("not_holodeck:1:200:chr1:10F:20R:0:0:0").is_none());
        assert!(parse_encoded_pe_name("").is_none());
        assert!(parse_encoded_se_name("holodeck:1").is_none());
        // Missing fragment-length field (pre-existing format) should not parse.
        assert!(parse_encoded_pe_name("holodeck:1:chr1:10F:20R:0:0:0").is_none());
        assert!(parse_encoded_se_name("holodeck:1:chr1:10F:0:0").is_none());
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
