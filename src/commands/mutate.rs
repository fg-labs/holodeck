use std::io::Write;
use std::path::PathBuf;

use anyhow::{Result, bail};
use clap::Parser;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, Geometric};

use super::command::Command;
use super::common::{BedOptions, ReferenceOptions, SeedOptions};
use crate::bed::TargetRegions;
use crate::fasta::Fasta;
use crate::ploidy::PloidyMap;
use crate::seed::{derive_seed, resolve_seed};
use crate::vcf::writer::VcfWriter;

/// DNA bases for random mutation generation.
const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Generate a VCF of random mutations from a reference genome.
///
/// Produces a standard VCF with proper GT fields that can be fed directly to
/// `holodeck simulate`.  Supports independent control of SNP, indel, and MNP
/// rates, with configurable ploidy including per-contig and per-region
/// overrides for handling sex chromosomes and pseudo-autosomal regions.
#[derive(Parser, Debug)]
#[command(after_long_help = "EXAMPLES:\n  \
    holodeck mutate -r ref.fa -o muts.vcf --snp-rate 0.001\n  \
    holodeck mutate -r ref.fa -o muts.vcf --snp-rate 0.001 --ploidy 2 \
    --ploidy-override chrX=1 --ploidy-override chrX:10001-2781479=2")]
pub struct Mutate {
    #[command(flatten)]
    pub reference: ReferenceOptions,

    #[command(flatten)]
    pub bed: BedOptions,

    #[command(flatten)]
    pub seed: SeedOptions,

    /// Output VCF file path.
    #[arg(short = 'o', long, value_name = "VCF")]
    pub output: PathBuf,

    /// Rate of SNP mutations per base.
    #[arg(long, default_value_t = 0.001, value_name = "FLOAT")]
    pub snp_rate: f64,

    /// Rate of indel mutations per base.
    #[arg(long, default_value_t = 0.0001, value_name = "FLOAT")]
    pub indel_rate: f64,

    /// Rate of MNP (multi-nucleotide polymorphism) mutations per base.
    #[arg(long, default_value_t = 0.00005, value_name = "FLOAT")]
    pub mnp_rate: f64,

    /// Parameter for the geometric distribution of indel lengths. Larger values
    /// produce shorter indels on average.
    #[arg(long, default_value_t = 0.7, value_name = "FLOAT")]
    pub indel_length_param: f64,

    /// Ratio of heterozygous to homozygous variants. A value of 2.0 means
    /// twice as many het as hom variants.
    #[arg(long, default_value_t = 2.0, value_name = "FLOAT")]
    pub het_hom_ratio: f64,

    /// Default ploidy for all contigs.
    #[arg(long, default_value_t = 2, value_name = "INT")]
    pub ploidy: u8,

    /// Per-contig or per-region ploidy override. Format: `CONTIG=PLOIDY` or
    /// `CONTIG:START-END=PLOIDY` (1-based coordinates). Later values override
    /// earlier ones for the same positions (last-writer-wins).
    ///
    /// Example: `--ploidy-override chrX=1 --ploidy-override chrX:10001-2781479=2`
    #[arg(long, value_name = "SPEC")]
    pub ploidy_override: Vec<String>,
}

impl Command for Mutate {
    fn execute(&self) -> Result<()> {
        self.validate()?;
        self.run_mutation()
    }
}

impl Mutate {
    /// Validate command-line arguments.
    fn validate(&self) -> Result<()> {
        if self.snp_rate < 0.0 || self.indel_rate < 0.0 || self.mnp_rate < 0.0 {
            bail!("Individual mutation rates must be >= 0");
        }
        let total_rate = self.snp_rate + self.indel_rate + self.mnp_rate;
        if !(0.0..=1.0).contains(&total_rate) {
            bail!("Combined mutation rate must be in [0, 1], got {total_rate}");
        }
        if self.indel_length_param <= 0.0 || self.indel_length_param >= 1.0 {
            bail!("--indel-length-param must be in (0, 1)");
        }
        if self.het_hom_ratio < 0.0 {
            bail!("--het-hom-ratio must be >= 0");
        }
        if self.ploidy == 0 {
            bail!("--ploidy must be >= 1");
        }
        Ok(())
    }

    /// Run the mutation generation pipeline.
    fn run_mutation(&self) -> Result<()> {
        let seed_desc = format!(
            "mutate:{}:{}:{}:{}:{}",
            self.reference.reference.display(),
            self.snp_rate,
            self.indel_rate,
            self.mnp_rate,
            self.ploidy,
        );
        let seed = resolve_seed(self.seed.seed, &seed_desc);
        let mut rng = SmallRng::seed_from_u64(seed);
        log::info!("Using random seed: {seed}");

        let mut fasta = Fasta::from_path(&self.reference.reference)?;
        let dict = fasta.dict().clone();
        log::info!(
            "Loaded reference with {} contigs, total {} bp",
            dict.len(),
            dict.total_length()
        );

        let targets = match &self.bed.targets {
            Some(bed_path) => {
                let t = TargetRegions::from_path(bed_path, &dict)?;
                log::info!(
                    "Restricting mutations to {} bp of target territory",
                    t.total_territory()
                );
                Some(t)
            }
            None => None,
        };

        let ploidy_map = PloidyMap::new(self.ploidy, &self.ploidy_override)?;
        let indel_dist = Geometric::new(self.indel_length_param)
            .map_err(|e| anyhow::anyhow!("Invalid indel length distribution: {e}"))?;

        // Open output VCF. Compression follows the file extension (`.gz`/`.bgz`
        // → BGZF, else plain text) so the file is named truthfully for
        // `holodeck simulate` and any external tool that keys codec off the
        // extension (e.g. `tabix`/`bcftools`).
        let mut vcf_out = VcfWriter::new(&self.output)?;
        Self::write_vcf_header(&mut vcf_out, &dict)?;

        let total_rate = self.snp_rate + self.indel_rate + self.mnp_rate;
        let mut total_variants = 0u64;
        let contig_names: Vec<String> = dict.names().into_iter().map(String::from).collect();

        for contig_name in &contig_names {
            // Use a per-contig RNG for reference normalization so ambiguity
            // resolution is reproducible and independent of the mutation RNG.
            let contig_seed = derive_seed(seed, contig_name);
            let mut ref_rng = SmallRng::seed_from_u64(contig_seed);
            let reference = fasta.load_contig(contig_name, &mut ref_rng)?;
            let contig_idx = dict.get_by_name(contig_name).unwrap().index();

            let mut pos = 0u32;
            #[expect(clippy::cast_possible_truncation, reason = "contig length fits u32")]
            let contig_len = reference.len() as u32;

            while pos < contig_len {
                // Skip non-ACGT bases (Ns).
                if !BASES.contains(&reference[pos as usize]) {
                    pos += 1;
                    continue;
                }

                // Check if this position is within targets (if specified).
                if let Some(tgt) = &targets
                    && !tgt.overlaps(contig_idx, pos, pos + 1)
                {
                    pos += 1;
                    continue;
                }

                // Decide if a mutation occurs at this position.
                if rng.random::<f64>() >= total_rate {
                    pos += 1;
                    continue;
                }

                // Determine mutation type by drawing a second random value
                // and partitioning [0, total_rate) into SNP/indel/MNP ranges.
                let ploidy = ploidy_map.ploidy_at(contig_name, pos);
                let type_roll: f64 = rng.random::<f64>() * total_rate;
                let (ref_allele, alt_allele, advance) = if type_roll < self.snp_rate {
                    generate_snp(reference[pos as usize], &mut rng)
                } else if type_roll < self.snp_rate + self.indel_rate {
                    generate_indel(&reference, pos, contig_len, &indel_dist, &mut rng)
                } else {
                    // MNP requires at least 2 bases remaining.
                    if pos + 2 > contig_len {
                        pos += 1;
                        continue;
                    }
                    generate_mnp(&reference, pos, contig_len, &mut rng)
                };

                // The anchor-position check above only covers position `pos`.
                // Multi-base REF alleles (MNPs and deletions) may still span
                // an ambiguity-resolved lowercase byte at `pos+1` or later,
                // which would emit a noncanonical REF into the VCF — skip
                // those variants.
                if !ref_allele.iter().all(|b| BASES.contains(b)) {
                    pos += 1;
                    continue;
                }

                // Assign genotype (het vs hom).
                let gt = generate_genotype(ploidy, self.het_hom_ratio, &mut rng);

                // Write VCF record (1-based position).
                writeln!(
                    vcf_out,
                    "{contig_name}\t{}\t.\t{}\t{}\t100\tPASS\t.\tGT\t{gt}",
                    pos + 1,
                    String::from_utf8_lossy(&ref_allele),
                    String::from_utf8_lossy(&alt_allele),
                )?;

                total_variants += 1;
                pos += advance;
            }
        }

        // Finalize the VCF: flush buffered data and (when BGZF) write the EOF
        // block.
        vcf_out.close()?;

        log::info!("Generated {total_variants} variants");
        Ok(())
    }

    /// Write the VCF header to the output file.
    fn write_vcf_header<W: Write>(
        out: &mut W,
        dict: &crate::sequence_dict::SequenceDictionary,
    ) -> Result<()> {
        writeln!(out, "##fileformat=VCFv4.3")?;
        writeln!(out, "##source=holodeck-mutate")?;
        for meta in dict.iter() {
            writeln!(out, "##contig=<ID={},length={}>", meta.name(), meta.length())?;
        }
        writeln!(out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
        writeln!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE")?;
        Ok(())
    }
}

/// Generate a SNP: single base substitution.
/// Returns `(ref_allele, alt_allele, advance)`.
fn generate_snp(ref_base: u8, rng: &mut impl Rng) -> (Vec<u8>, Vec<u8>, u32) {
    let alt = loop {
        let candidate = BASES[rng.random_range(0..4)];
        if candidate != ref_base {
            break candidate;
        }
    };
    (vec![ref_base], vec![alt], 1)
}

/// Generate an indel (insertion or deletion).
/// Returns `(ref_allele, alt_allele, advance)`.
fn generate_indel(
    reference: &[u8],
    pos: u32,
    contig_len: u32,
    indel_dist: &Geometric,
    rng: &mut impl Rng,
) -> (Vec<u8>, Vec<u8>, u32) {
    // Draw indel length from geometric distribution (minimum 1).
    #[expect(
        clippy::cast_possible_truncation,
        reason = "geometric distribution produces small values"
    )]
    let length = (indel_dist.sample(rng) as u32).max(1);
    let is_insertion: bool = rng.random();
    let anchor = reference[pos as usize];

    if is_insertion {
        // Insertion: REF = anchor base, ALT = anchor + random bases.
        let mut alt = vec![anchor];
        for _ in 0..length {
            alt.push(BASES[rng.random_range(0..4)]);
        }
        (vec![anchor], alt, 1)
    } else {
        // Deletion: REF = anchor + deleted bases, ALT = anchor.
        // Need at least 1 base after the anchor to delete.
        let available = contig_len.saturating_sub(pos + 1);
        if available == 0 {
            // Fall back to insertion at this position.
            let mut alt = vec![anchor];
            for _ in 0..length {
                alt.push(BASES[rng.random_range(0..4)]);
            }
            return (vec![anchor], alt, 1);
        }
        let del_len = length.min(available);
        let del_end = pos + 1 + del_len;
        let ref_allele: Vec<u8> = reference[pos as usize..del_end as usize].to_vec();
        #[expect(clippy::cast_possible_truncation, reason = "allele length fits u32")]
        let advance = ref_allele.len() as u32;
        (ref_allele, vec![anchor], advance)
    }
}

/// Generate an MNP (multi-nucleotide polymorphism, 2-3 bases).
/// Returns `(ref_allele, alt_allele, advance)`.
fn generate_mnp(
    reference: &[u8],
    pos: u32,
    contig_len: u32,
    rng: &mut impl Rng,
) -> (Vec<u8>, Vec<u8>, u32) {
    let length = rng.random_range(2..=3u32).min(contig_len - pos);
    let ref_allele: Vec<u8> = reference[pos as usize..(pos + length) as usize].to_vec();

    // Generate alt allele that differs in at least one position.
    let mut alt_allele = ref_allele.clone();
    for base in &mut alt_allele {
        if rng.random::<f64>() < 0.8 {
            // Mutate this position.
            let original = *base;
            *base = loop {
                let candidate = BASES[rng.random_range(0..4)];
                if candidate != original {
                    break candidate;
                }
            };
        }
    }
    // Ensure at least one base differs.
    if alt_allele == ref_allele {
        let idx = rng.random_range(0..alt_allele.len());
        let original = alt_allele[idx];
        alt_allele[idx] = loop {
            let candidate = BASES[rng.random_range(0..4)];
            if candidate != original {
                break candidate;
            }
        };
    }

    (ref_allele, alt_allele, length)
}

/// Generate a genotype string for the given ploidy and het/hom ratio.
///
/// The het/hom ratio controls the probability of heterozygous vs homozygous
/// genotypes. For a ratio of 2.0, ~67% of variants will be het and ~33% hom.
fn generate_genotype(ploidy: u8, het_hom_ratio: f64, rng: &mut impl Rng) -> String {
    let p_het = het_hom_ratio / (1.0 + het_hom_ratio);
    let is_het = ploidy > 1 && rng.random::<f64>() < p_het;

    if ploidy == 1 {
        // Haploid: always the alt allele.
        "1".to_string()
    } else if is_het {
        // Heterozygous: randomly choose which haplotype(s) get the alt.
        let mut alleles: Vec<&str> = vec!["0"; ploidy as usize];
        let alt_hap = rng.random_range(0..ploidy as usize);
        alleles[alt_hap] = "1";
        alleles.join("/")
    } else {
        // Homozygous alt: all haplotypes get the alt.
        vec!["1"; ploidy as usize].join("/")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::common::{BedOptions, ReferenceOptions, SeedOptions};

    /// Build a `Mutate` for a given reference and output path with a high SNP
    /// rate (so even a tiny contig yields several variants) and a fixed seed.
    fn mutate_for(reference: std::path::PathBuf, output: std::path::PathBuf) -> Mutate {
        Mutate {
            reference: ReferenceOptions { reference },
            bed: BedOptions { targets: None },
            seed: SeedOptions { seed: Some(42) },
            output,
            snp_rate: 0.05,
            indel_rate: 0.0,
            mnp_rate: 0.0,
            indel_length_param: 0.7,
            het_hom_ratio: 2.0,
            ploidy: 2,
            ploidy_override: Vec::new(),
        }
    }

    /// A `.vcf.gz` output path must produce a real BGZF stream (gzip magic up
    /// front) — not plain text mislabelled with a `.gz` name — that round-trips
    /// through `parse_variants_by_contig`, the reader `holodeck simulate` uses.
    #[test]
    fn mutate_gz_output_is_bgzf_and_reads_back() {
        let dir = tempfile::tempdir().unwrap();
        let seq: Vec<u8> = b"ACGT".iter().copied().cycle().take(2000).collect();
        let ref_path = crate::fasta::write_test_fasta(dir.path(), &[("chr1", &seq)]);
        let out_path = dir.path().join("muts.vcf.gz");

        mutate_for(ref_path, out_path.clone()).execute().unwrap();

        // Output must be gzip/BGZF, not plain text.
        let bytes = std::fs::read(&out_path).unwrap();
        assert_eq!(&bytes[..2], &[0x1f, 0x8b], "expected BGZF/gzip magic bytes");

        // And it must read back through the simulate consumption path.
        let dict = crate::sequence_dict::SequenceDictionary::from_entries(vec![
            crate::sequence_dict::SequenceMetadata::new(0, "chr1".to_string(), 2000),
        ]);
        let parsed = crate::vcf::parse_variants_by_contig(&out_path, None, &dict).unwrap();
        assert!(
            parsed.by_contig.get("chr1").is_some_and(|v| !v.is_empty()),
            "expected at least one variant to round-trip from the .gz output"
        );
    }

    /// A plain `.vcf` output path must produce uncompressed text that also
    /// reads back through the simulate consumption path.
    #[test]
    fn mutate_plain_output_is_text_and_reads_back() {
        let dir = tempfile::tempdir().unwrap();
        let seq: Vec<u8> = b"ACGT".iter().copied().cycle().take(2000).collect();
        let ref_path = crate::fasta::write_test_fasta(dir.path(), &[("chr1", &seq)]);
        let out_path = dir.path().join("muts.vcf");

        mutate_for(ref_path, out_path.clone()).execute().unwrap();

        let bytes = std::fs::read(&out_path).unwrap();
        assert_ne!(&bytes[..2], &[0x1f, 0x8b], "plain .vcf output must not be gzip");
        assert!(bytes.starts_with(b"##fileformat=VCFv4.3"), "expected plain VCF header");

        let dict = crate::sequence_dict::SequenceDictionary::from_entries(vec![
            crate::sequence_dict::SequenceMetadata::new(0, "chr1".to_string(), 2000),
        ]);
        let parsed = crate::vcf::parse_variants_by_contig(&out_path, None, &dict).unwrap();
        assert!(parsed.by_contig.get("chr1").is_some_and(|v| !v.is_empty()));
    }

    #[test]
    fn test_generate_snp() {
        let mut rng = rand::rng();
        for _ in 0..100 {
            let (ref_a, alt_a, advance) = generate_snp(b'A', &mut rng);
            assert_eq!(ref_a, vec![b'A']);
            assert_eq!(alt_a.len(), 1);
            assert_ne!(alt_a[0], b'A');
            assert_eq!(advance, 1);
        }
    }

    #[test]
    fn test_generate_indel_insertion() {
        let reference = b"ACGTACGT";
        let indel_dist = Geometric::new(0.7).unwrap();
        let mut rng = SmallRng::seed_from_u64(42);

        let mut found_insertion = false;
        for _ in 0..100 {
            let (ref_a, alt_a, advance) = generate_indel(reference, 2, 8, &indel_dist, &mut rng);
            if alt_a.len() > ref_a.len() {
                found_insertion = true;
                assert_eq!(ref_a.len(), 1); // Anchor base only
                assert!(alt_a.len() >= 2); // Anchor + at least 1 inserted base
                assert_eq!(ref_a[0], alt_a[0]); // Same anchor
                assert_eq!(advance, 1);
            }
        }
        assert!(found_insertion, "Should have generated at least one insertion");
    }

    #[test]
    fn test_generate_indel_deletion() {
        let reference = b"ACGTACGT";
        let indel_dist = Geometric::new(0.7).unwrap();
        let mut rng = SmallRng::seed_from_u64(99);

        let mut found_deletion = false;
        for _ in 0..100 {
            let (ref_a, alt_a, _advance) = generate_indel(reference, 2, 8, &indel_dist, &mut rng);
            if ref_a.len() > alt_a.len() {
                found_deletion = true;
                assert_eq!(alt_a.len(), 1); // Anchor base only
                assert!(ref_a.len() >= 2); // Anchor + at least 1 deleted base
                assert_eq!(ref_a[0], alt_a[0]); // Same anchor
            }
        }
        assert!(found_deletion, "Should have generated at least one deletion");
    }

    #[test]
    fn test_generate_mnp() {
        let reference = b"ACGTACGT";
        let mut rng = rand::rng();
        for _ in 0..100 {
            let (ref_a, alt_a, advance) = generate_mnp(reference, 1, 8, &mut rng);
            assert!(ref_a.len() >= 2 && ref_a.len() <= 3);
            assert_eq!(ref_a.len(), alt_a.len());
            assert_ne!(ref_a, alt_a);
            assert_eq!(advance as usize, ref_a.len());
        }
    }

    #[test]
    fn test_generate_genotype_haploid() {
        let mut rng = rand::rng();
        let gt = generate_genotype(1, 2.0, &mut rng);
        assert_eq!(gt, "1");
    }

    #[test]
    fn test_generate_genotype_diploid_het() {
        let mut rng = rand::rng();
        let mut het_count = 0;
        let mut hom_count = 0;
        for _ in 0..1000 {
            let gt = generate_genotype(2, 2.0, &mut rng);
            if gt == "0/1" || gt == "1/0" {
                het_count += 1;
            } else if gt == "1/1" {
                hom_count += 1;
            }
        }
        // With het_hom_ratio=2.0, expect ~67% het, ~33% hom.
        assert!(het_count > 500, "expected majority het, got {het_count}");
        assert!(hom_count > 200, "expected some hom, got {hom_count}");
    }

    #[test]
    fn test_generate_genotype_triploid() {
        let mut rng = rand::rng();
        let gt = generate_genotype(3, 0.0, &mut rng);
        // hom-alt for triploid.
        assert_eq!(gt, "1/1/1");
    }
}
