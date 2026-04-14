use anyhow::{Result, bail};

/// A parsed genotype for a single sample at a variant site.
///
/// Alleles are 0-indexed where 0 means reference and 1+ means the Nth
/// alternate allele.  `None` represents a missing allele (`.` in VCF).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Genotype {
    /// Allele indices for each haplotype.  Length equals ploidy at this site.
    alleles: Vec<Option<usize>>,
    /// Whether the genotype is phased (`|` separator) or unphased (`/`).
    phased: bool,
}

impl Genotype {
    /// Parse a VCF GT field string (e.g. `"0/1"`, `"1|0"`, `"0|1|2"`,
    /// `"./."`, `"0"`).
    ///
    /// Returns the parsed genotype with allele indices and phasing information.
    ///
    /// # Errors
    /// Returns an error if the GT string is empty or contains invalid allele
    /// values.
    pub fn parse(gt: &str) -> Result<Self> {
        if gt.is_empty() {
            bail!("Empty GT field");
        }

        // Determine separator: phased (|) or unphased (/).
        // A haploid genotype (e.g. "0") has no separator; treat as phased.
        let (sep, phased) = if gt.contains('|') {
            ('|', true)
        } else if gt.contains('/') {
            ('/', false)
        } else {
            // Haploid: single allele, no separator.
            let allele = parse_allele(gt)?;
            return Ok(Self { alleles: vec![allele], phased: true });
        };

        let alleles: Result<Vec<Option<usize>>> = gt.split(sep).map(parse_allele).collect();
        let alleles = alleles?;

        if alleles.is_empty() {
            bail!("GT field has no alleles: {gt}");
        }

        Ok(Self { alleles, phased })
    }

    /// Return the allele indices (one per haplotype).
    #[must_use]
    pub fn alleles(&self) -> &[Option<usize>] {
        &self.alleles
    }

    /// Return whether this genotype is phased.
    #[must_use]
    pub fn is_phased(&self) -> bool {
        self.phased
    }

    /// Return the ploidy (number of allele slots) at this site.
    #[must_use]
    pub fn ploidy(&self) -> usize {
        self.alleles.len()
    }

    /// Return `true` if any allele is non-reference (index > 0).
    #[must_use]
    pub fn has_alt(&self) -> bool {
        self.alleles.iter().any(|a| matches!(a, Some(idx) if *idx > 0))
    }

    /// Return `true` if all alleles are the same non-reference allele.
    #[must_use]
    pub fn is_hom_alt(&self) -> bool {
        if !self.has_alt() {
            return false;
        }
        let first = self.alleles[0];
        first.is_some() && self.alleles.iter().all(|a| *a == first)
    }

    /// Return `true` if all alleles are reference (index 0).
    #[must_use]
    pub fn is_hom_ref(&self) -> bool {
        self.alleles.iter().all(|a| *a == Some(0))
    }
}

/// Parse a single allele field: `"."` → `None`, digits → `Some(n)`.
fn parse_allele(s: &str) -> Result<Option<usize>> {
    if s == "." {
        Ok(None)
    } else {
        let idx: usize = s.parse().map_err(|_| anyhow::anyhow!("Invalid allele value: '{s}'"))?;
        Ok(Some(idx))
    }
}

/// A variant record with its position, alleles, and genotype for the selected
/// sample.
#[derive(Debug, Clone)]
pub struct VariantRecord {
    /// 0-based reference position.
    pub position: u32,
    /// Reference allele bases.
    pub ref_allele: Vec<u8>,
    /// Alternate allele sequences, indexed by alt allele number (0-based among
    /// alternates, so alt allele 1 in VCF is index 0 here).
    pub alt_alleles: Vec<Vec<u8>>,
    /// Parsed genotype for the selected sample.
    pub genotype: Genotype,
}

impl VariantRecord {
    /// Return the allele sequence for a given allele index (0 = ref, 1+ = alt).
    ///
    /// Returns `None` if the index is out of range.
    #[must_use]
    pub fn allele_bases(&self, allele_index: usize) -> Option<&[u8]> {
        if allele_index == 0 {
            Some(&self.ref_allele)
        } else {
            self.alt_alleles.get(allele_index - 1).map(Vec::as_slice)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_diploid_unphased_het() {
        let gt = Genotype::parse("0/1").unwrap();
        assert_eq!(gt.alleles(), &[Some(0), Some(1)]);
        assert!(!gt.is_phased());
        assert_eq!(gt.ploidy(), 2);
        assert!(gt.has_alt());
        assert!(!gt.is_hom_alt());
        assert!(!gt.is_hom_ref());
    }

    #[test]
    fn test_parse_diploid_phased_het() {
        let gt = Genotype::parse("1|0").unwrap();
        assert_eq!(gt.alleles(), &[Some(1), Some(0)]);
        assert!(gt.is_phased());
        assert!(gt.has_alt());
    }

    #[test]
    fn test_parse_diploid_hom_ref() {
        let gt = Genotype::parse("0/0").unwrap();
        assert_eq!(gt.alleles(), &[Some(0), Some(0)]);
        assert!(!gt.has_alt());
        assert!(gt.is_hom_ref());
        assert!(!gt.is_hom_alt());
    }

    #[test]
    fn test_parse_diploid_hom_alt() {
        let gt = Genotype::parse("1/1").unwrap();
        assert_eq!(gt.alleles(), &[Some(1), Some(1)]);
        assert!(gt.has_alt());
        assert!(gt.is_hom_alt());
        assert!(!gt.is_hom_ref());
    }

    #[test]
    fn test_parse_triploid() {
        let gt = Genotype::parse("0|1|2").unwrap();
        assert_eq!(gt.alleles(), &[Some(0), Some(1), Some(2)]);
        assert!(gt.is_phased());
        assert_eq!(gt.ploidy(), 3);
    }

    #[test]
    fn test_parse_missing_alleles() {
        let gt = Genotype::parse("./.").unwrap();
        assert_eq!(gt.alleles(), &[None, None]);
        assert!(!gt.has_alt());
        assert!(!gt.is_hom_ref());
    }

    #[test]
    fn test_parse_haploid() {
        let gt = Genotype::parse("0").unwrap();
        assert_eq!(gt.alleles(), &[Some(0)]);
        assert!(gt.is_phased()); // haploid treated as phased
        assert_eq!(gt.ploidy(), 1);

        let gt = Genotype::parse("1").unwrap();
        assert!(gt.has_alt());
    }

    #[test]
    fn test_parse_empty_errors() {
        assert!(Genotype::parse("").is_err());
    }

    #[test]
    fn test_parse_invalid_allele_errors() {
        assert!(Genotype::parse("0/abc").is_err());
    }

    #[test]
    fn test_variant_record_allele_bases() {
        let vr = VariantRecord {
            position: 100,
            ref_allele: b"A".to_vec(),
            alt_alleles: vec![b"T".to_vec(), b"GCC".to_vec()],
            genotype: Genotype::parse("0/1").unwrap(),
        };
        assert_eq!(vr.allele_bases(0), Some(b"A".as_slice()));
        assert_eq!(vr.allele_bases(1), Some(b"T".as_slice()));
        assert_eq!(vr.allele_bases(2), Some(b"GCC".as_slice()));
        assert_eq!(vr.allele_bases(3), None);
    }
}
