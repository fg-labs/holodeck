//! Read naming schemes for simulated reads.
//!
//! Supports two modes: **encoded** names that embed truth coordinates
//! (contig, position, strand, haplotype, fragment length, error count) for
//! downstream evaluation, and **simple** names that are just sequential
//! identifiers.
//!
//! Fields are separated by `::` (double colon). Contig names may legally
//! contain single `:` characters (e.g. HLA contigs like `HLA-A*01:01:01:01`);
//! using `::` as the separator keeps parsing unambiguous without requiring
//! right-to-left tricks.

/// Truth alignment data for a single read, used to encode position
/// information into the read name.
///
/// In paired-end mode both reads of a pair carry the same `fragment_length`
/// and `haplotype`.
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
    /// Length in bases of the source fragment. When shorter than the read
    /// length, the remainder of the read is adapter (optionally N-padded).
    pub fragment_length: u32,
    /// Number of sequencing errors introduced.
    pub n_errors: usize,
}

/// Field separator for encoded read names. Chosen so that contig names
/// containing single `:` characters remain unambiguous.
const SEP: &str = "::";

/// Prefix that identifies an encoded holodeck read name.
const PREFIX: &str = "holodeck";

/// Format a read name in encoded mode for a paired-end read.
///
/// Format: `holodeck::READ_NUM::FRAG_LEN::CONTIG::POS1+STRAND::POS2+STRAND::HAP::ERRS1::ERRS2`
///
/// Example: `holodeck::42::600::chr1::10000F::10450R::0::2::1`
///
/// Both R1 and R2 must share the same `contig`, `haplotype`, and
/// `fragment_length` — these fields are encoded once from `r1` only.
#[must_use]
pub fn encoded_pe_name(read_num: u64, r1: &TruthAlignment, r2: &TruthAlignment) -> String {
    debug_assert_eq!(
        r1.contig, r2.contig,
        "R1 and R2 must be on the same contig for encoded PE names"
    );
    debug_assert_eq!(
        r1.fragment_length, r2.fragment_length,
        "R1 and R2 must share the same fragment length"
    );
    debug_assert_eq!(r1.haplotype, r2.haplotype, "R1 and R2 must share the same haplotype");
    assert_contig_is_safe(&r1.contig);
    format!(
        "{PREFIX}{SEP}{read_num}{SEP}{frag_len}{SEP}{contig}{SEP}{p1}{s1}{SEP}{p2}{s2}{SEP}{hap}{SEP}{e1}{SEP}{e2}",
        frag_len = r1.fragment_length,
        contig = r1.contig,
        p1 = r1.position,
        s1 = strand_char(r1.is_forward),
        p2 = r2.position,
        s2 = strand_char(r2.is_forward),
        hap = r1.haplotype,
        e1 = r1.n_errors,
        e2 = r2.n_errors,
    )
}

/// Format a read name in encoded mode for a single-end read.
///
/// Format: `holodeck::READ_NUM::FRAG_LEN::CONTIG::POS+STRAND::HAP::ERRS`
///
/// Example: `holodeck::42::275::chr1::10000F::0::2`
#[must_use]
pub fn encoded_se_name(read_num: u64, r1: &TruthAlignment) -> String {
    assert_contig_is_safe(&r1.contig);
    format!(
        "{PREFIX}{SEP}{read_num}{SEP}{frag_len}{SEP}{contig}{SEP}{pos}{strand}{SEP}{hap}{SEP}{errs}",
        frag_len = r1.fragment_length,
        contig = r1.contig,
        pos = r1.position,
        strand = strand_char(r1.is_forward),
        hap = r1.haplotype,
        errs = r1.n_errors,
    )
}

/// Format a simple read name with no truth information.
///
/// Format: `holodeck::READ_NUM`
#[must_use]
pub fn simple_name(read_num: u64) -> String {
    format!("{PREFIX}{SEP}{read_num}")
}

/// Panic in debug builds if `contig` contains characters that would corrupt
/// the encoded read name: `@` (FASTQ header prefix) or `::` (field
/// separator).
fn assert_contig_is_safe(contig: &str) {
    debug_assert!(
        !contig.contains('@'),
        "contig name must not contain '@' (FASTQ header prefix): {contig:?}",
    );
    debug_assert!(
        !contig.contains(SEP),
        "contig name must not contain '::' (field separator): {contig:?}",
    );
}

/// Return the strand character for a boolean forward flag.
fn strand_char(is_forward: bool) -> char {
    if is_forward { 'F' } else { 'R' }
}

/// Parse an encoded paired-end read name back into truth alignments.
///
/// Returns `(read_num, r1_truth, r2_truth)` or `None` if the name doesn't
/// match the expected format. The returned alignments share `contig`,
/// `haplotype`, and `fragment_length`.
#[must_use]
pub fn parse_encoded_pe_name(name: &str) -> Option<(u64, TruthAlignment, TruthAlignment)> {
    let parts: Vec<&str> = name.split(SEP).collect();
    // holodeck, read_num, frag_len, contig, pos1+strand, pos2+strand, hap, errs1, errs2.
    if parts.len() != 9 || parts[0] != PREFIX {
        return None;
    }

    let read_num: u64 = parts[1].parse().ok()?;
    let fragment_length: u32 = parts[2].parse().ok()?;
    let contig = parts[3];
    if contig.is_empty() {
        return None;
    }
    let (pos1, fwd1) = parse_pos_strand(parts[4])?;
    let (pos2, fwd2) = parse_pos_strand(parts[5])?;
    let haplotype: usize = parts[6].parse().ok()?;
    let n_errors_r1: usize = parts[7].parse().ok()?;
    let n_errors_r2: usize = parts[8].parse().ok()?;

    let r1 = TruthAlignment {
        contig: contig.to_string(),
        position: pos1,
        is_forward: fwd1,
        haplotype,
        fragment_length,
        n_errors: n_errors_r1,
    };
    let r2 = TruthAlignment {
        contig: contig.to_string(),
        position: pos2,
        is_forward: fwd2,
        haplotype,
        fragment_length,
        n_errors: n_errors_r2,
    };

    Some((read_num, r1, r2))
}

/// Parse an encoded single-end read name back into a truth alignment.
///
/// Returns `(read_num, truth)` or `None` if the name doesn't match.
#[must_use]
pub fn parse_encoded_se_name(name: &str) -> Option<(u64, TruthAlignment)> {
    let parts: Vec<&str> = name.split(SEP).collect();
    // holodeck, read_num, frag_len, contig, pos+strand, hap, errs.
    if parts.len() != 7 || parts[0] != PREFIX {
        return None;
    }

    let read_num: u64 = parts[1].parse().ok()?;
    let fragment_length: u32 = parts[2].parse().ok()?;
    let contig = parts[3];
    if contig.is_empty() {
        return None;
    }
    let (pos, fwd) = parse_pos_strand(parts[4])?;
    let haplotype: usize = parts[5].parse().ok()?;
    let n_errors: usize = parts[6].parse().ok()?;

    let truth = TruthAlignment {
        contig: contig.to_string(),
        position: pos,
        is_forward: fwd,
        haplotype,
        fragment_length,
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

    fn make_truth(
        contig: &str,
        pos: u32,
        fwd: bool,
        hap: usize,
        frag_len: u32,
        errs: usize,
    ) -> TruthAlignment {
        TruthAlignment {
            contig: contig.to_string(),
            position: pos,
            is_forward: fwd,
            haplotype: hap,
            fragment_length: frag_len,
            n_errors: errs,
        }
    }

    #[test]
    fn test_encoded_pe_name() {
        let r1 = make_truth("chr1", 10000, true, 0, 600, 2);
        let r2 = make_truth("chr1", 10450, false, 0, 600, 1);
        let name = encoded_pe_name(42, &r1, &r2);
        assert_eq!(name, "holodeck::42::600::chr1::10000F::10450R::0::2::1");
    }

    #[test]
    fn test_encoded_se_name() {
        let r1 = make_truth("chr1", 10000, true, 0, 275, 2);
        let name = encoded_se_name(42, &r1);
        assert_eq!(name, "holodeck::42::275::chr1::10000F::0::2");
    }

    #[test]
    fn test_simple_name() {
        assert_eq!(simple_name(42), "holodeck::42");
        assert_eq!(simple_name(1), "holodeck::1");
    }

    #[test]
    fn test_parse_pe_roundtrip() {
        let r1 = make_truth("chr1", 10000, true, 0, 600, 2);
        let r2 = make_truth("chr1", 10450, false, 0, 600, 1);
        let name = encoded_pe_name(42, &r1, &r2);

        let (num, parsed_r1, parsed_r2) = parse_encoded_pe_name(&name).unwrap();
        assert_eq!(num, 42);
        assert_eq!(parsed_r1, r1);
        assert_eq!(parsed_r2, r2);
    }

    #[test]
    fn test_parse_se_roundtrip() {
        let r1 = make_truth("chrX", 500, false, 1, 275, 0);
        let name = encoded_se_name(99, &r1);

        let (num, parsed) = parse_encoded_se_name(&name).unwrap();
        assert_eq!(num, 99);
        assert_eq!(parsed, r1);
    }

    #[test]
    fn test_parse_pe_with_colon_in_contig() {
        // HLA contig names contain single colons, which is the whole point of
        // using `::` as the separator.
        let r1 = make_truth("HLA-A*01:01:01:01", 100, true, 0, 450, 0);
        let r2 = make_truth("HLA-A*01:01:01:01", 400, false, 0, 450, 1);
        let name = encoded_pe_name(1, &r1, &r2);

        let (num, parsed_r1, parsed_r2) = parse_encoded_pe_name(&name).unwrap();
        assert_eq!(num, 1);
        assert_eq!(parsed_r1, r1);
        assert_eq!(parsed_r2, r2);
    }

    #[test]
    fn test_parse_se_with_colon_in_contig() {
        let r1 = make_truth("HLA-B*07:02", 50, false, 1, 120, 3);
        let name = encoded_se_name(5, &r1);

        let (num, parsed) = parse_encoded_se_name(&name).unwrap();
        assert_eq!(num, 5);
        assert_eq!(parsed, r1);
    }

    #[test]
    fn test_short_fragment_encodes_adapter_boundary() {
        // A 40-base fragment read at 150bp read length: bases[40..] are adapter.
        let r1 = make_truth("chr1", 10000, true, 0, 40, 0);
        let r2 = make_truth("chr1", 10000, false, 0, 40, 0);
        let name = encoded_pe_name(7, &r1, &r2);

        let (_num, parsed_r1, _parsed_r2) = parse_encoded_pe_name(&name).unwrap();
        assert_eq!(parsed_r1.fragment_length, 40);
    }

    #[test]
    fn test_cross_format_rejection() {
        // SE name should not parse as PE and vice versa — field counts differ.
        let se_name = encoded_se_name(1, &make_truth("chr1", 10, true, 0, 150, 0));
        assert!(parse_encoded_pe_name(&se_name).is_none());

        let pe_name = encoded_pe_name(
            1,
            &make_truth("chr1", 10, true, 0, 200, 0),
            &make_truth("chr1", 20, false, 0, 200, 0),
        );
        assert!(parse_encoded_se_name(&pe_name).is_none());
    }

    #[test]
    fn test_parse_invalid_names() {
        // Wrong prefix.
        assert!(parse_encoded_pe_name("not_holodeck::1::200::chr1::10F::20R::0::0::0").is_none());
        // Empty.
        assert!(parse_encoded_pe_name("").is_none());
        assert!(parse_encoded_se_name("holodeck::1").is_none());
        // Pre-existing single-colon format must not parse under the new scheme.
        assert!(parse_encoded_pe_name("holodeck:1:200:chr1:10F:20R:0:0:0").is_none());
        assert!(parse_encoded_se_name("holodeck:1:150:chr1:10F:0:0").is_none());
        // Empty contig.
        assert!(parse_encoded_se_name("holodeck::1::150::::10F::0::0").is_none());
    }

    #[test]
    fn test_parse_pos_strand() {
        assert_eq!(parse_pos_strand("10000F"), Some((10000, true)));
        assert_eq!(parse_pos_strand("450R"), Some((450, false)));
        assert_eq!(parse_pos_strand("0F"), Some((0, true)));
        assert_eq!(parse_pos_strand("X"), None);
        assert_eq!(parse_pos_strand(""), None);
    }

    #[test]
    #[should_panic(expected = "contig name must not contain '@'")]
    fn test_formatter_rejects_at_in_contig_se() {
        let r1 = make_truth("bad@contig", 1, true, 0, 100, 0);
        let _ = encoded_se_name(1, &r1);
    }

    #[test]
    #[should_panic(expected = "contig name must not contain '::'")]
    fn test_formatter_rejects_double_colon_in_contig_pe() {
        let r1 = make_truth("bad::contig", 1, true, 0, 100, 0);
        let r2 = make_truth("bad::contig", 20, false, 0, 100, 0);
        let _ = encoded_pe_name(1, &r1, &r2);
    }
}
