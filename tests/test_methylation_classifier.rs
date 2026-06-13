//! Tests for the per-CpG ownership classifier in `src/vcf/methylation.rs`.
//! Each test constructs a (reference, variants, per-haplotype bitmap) input
//! and asserts the classifier's record-level output matches expectation.

use holodeck_lib::vcf::methylation::{ClassifyError, CpgPlacement, classify_cpgs};

#[test]
fn empty_reference_produces_no_records() {
    let reference: &[u8] = b"";
    let variants = vec![];
    let placements = classify_cpgs(reference, &variants).unwrap();
    assert!(placements.is_empty());
}

#[test]
fn reference_cpgs_no_variants_all_standalone() {
    // ref: ACGTACG     positions 0..7
    //       ^   ^      CpGs at top-C positions 1 and 5
    let reference = b"ACGTACG";
    let variants = vec![];
    let placements = classify_cpgs(reference, &variants).unwrap();
    assert_eq!(
        placements,
        vec![CpgPlacement::Standalone { ref_pos: 1 }, CpgPlacement::Standalone { ref_pos: 5 },],
    );
}

#[test]
fn case_04_snp_destroys_reference_cpg() {
    use holodeck_lib::vcf::genotype::{Genotype, VariantRecord};
    // ref: ACGTACG; SNP at position 1: C→T on hap0, hap1 carries REF.
    // Hap0 has no CpG at pos 1; hap1 keeps the CpG.
    // Hap0 still has the CpG at pos 5; hap1 also.
    let reference = b"ACGTACG";
    let variants = vec![VariantRecord {
        position: 1,
        ref_allele: b"C".to_vec(),
        alt_alleles: vec![b"T".to_vec()],
        genotype: Genotype::parse("1|0").unwrap(),
    }];
    let placements = classify_cpgs(reference, &variants).unwrap();
    // Expectation: the ref CpG at pos 1 is still a standalone record
    // (it exists on hap1); the per-haplotype MT/MB will be set to
    // "0|.,." style by the writer, but classification just lists the
    // record. The variant record itself carries no CpG annotation
    // because the alt (T) has no CpG context on hap0.
    // The ref CpG at pos 5 is unaffected.
    assert_eq!(
        placements,
        vec![CpgPlacement::Standalone { ref_pos: 1 }, CpgPlacement::Standalone { ref_pos: 5 },],
    );
}

#[test]
fn standalone_and_on_variant_at_same_position_both_emitted_in_order() {
    use holodeck_lib::vcf::genotype::{Genotype, VariantRecord};
    // ref: ACGT (CpG at top-C pos 1). Insertion at pos 1: REF=C, ALT=CG.
    // Variant hap0 materializes as A,C,G,G,T → CpG at hap-coords 1-2 (both from variant alt).
    //   OnVariant { variant_index: 0, hap_offset: 0 }.
    // Reference hap1 keeps the ref CpG at pos 1 → Standalone { ref_pos: 1 }.
    // Both placements share effective ref position 1; verify they're both emitted
    // with Standalone first (canonical order).
    let reference = b"ACGT";
    let variants = vec![VariantRecord {
        position: 1,
        ref_allele: b"C".to_vec(),
        alt_alleles: vec![b"CG".to_vec()],
        genotype: Genotype::parse("1|0").unwrap(),
    }];
    let placements = classify_cpgs(reference, &variants).unwrap();
    assert_eq!(
        placements,
        vec![
            CpgPlacement::Standalone { ref_pos: 1 },
            CpgPlacement::OnVariant { variant_index: 0, hap_offset: 0 },
        ],
    );
}

#[test]
fn case_07_adjacent_snps_jointly_form_cpg_upstream_owns() {
    use holodeck_lib::vcf::genotype::{Genotype, VariantRecord};
    // ref: AAATTAA; two adjacent SNPs on hap0: pos 3 T→C, pos 4 T→G.
    // Hap0 sequence becomes AAACGAA → CpG at hap-coord 3. Hap1 keeps REF: no CpG.
    let reference = b"AAATTAA";
    let variants = vec![
        VariantRecord {
            position: 3,
            ref_allele: b"T".to_vec(),
            alt_alleles: vec![b"C".to_vec()],
            genotype: Genotype::parse("1|0").unwrap(),
        },
        VariantRecord {
            position: 4,
            ref_allele: b"T".to_vec(),
            alt_alleles: vec![b"G".to_vec()],
            genotype: Genotype::parse("1|0").unwrap(),
        },
    ];
    let placements = classify_cpgs(reference, &variants).unwrap();
    assert_eq!(placements, vec![CpgPlacement::OnVariant { variant_index: 0, hap_offset: 0 }],);
}

// ── Edge-case batch (Task 6) ────────────────────────────────────────────────

#[test]
fn case_01_insertion_creates_cpg_at_downstream_boundary() {
    use holodeck_lib::vcf::genotype::{Genotype, VariantRecord};
    // ref: TTGG (no CpG). Insertion at pos 1: REF=T, ALT=TC. GT=1|0.
    //
    // Hap0 carries ALT: materialized as T(pos0), T(pos1), C(pos1), G(pos2), G(pos3) = TTCGG.
    //   var_hap_ranges for vi=0: hap_start = hap_position_for(1) = 1 (var_end=2 > 1 → delta=0).
    //   hap_end = 1 + len("TC") = 3.
    //   Scan TTCGG for CpG:
    //     h=2: C at hap-coord 2, G at hap-coord 3.
    //       src_h(2):  2 ∈ [1,3) → Variant{vi=0, hap_start=1}.
    //       src_h1(3): 3 ∈ [1,3)? 3 < 3 is false → Reference.
    //       hap_offset = hpos - hap_start = 2 - 1 = 1.
    //       → OnVariant { variant_index: 0, hap_offset: 1 }.
    // Hap1 keeps REF (TTGG): no CpG.
    let reference = b"TTGG";
    let variants = vec![VariantRecord {
        position: 1,
        ref_allele: b"T".to_vec(),
        alt_alleles: vec![b"TC".to_vec()],
        genotype: Genotype::parse("1|0").unwrap(),
    }];
    let placements = classify_cpgs(reference, &variants).unwrap();
    assert_eq!(placements, vec![CpgPlacement::OnVariant { variant_index: 0, hap_offset: 1 }],);
}

#[test]
fn case_02_snp_creates_g_at_upstream_boundary_makes_variant_own_cpg() {
    use holodeck_lib::vcf::genotype::{Genotype, VariantRecord};
    // ref: ATCAA. SNP at pos 3: A→G. GT=1|0.
    //
    // Hap0 carries ALT: A,T,C,G,A = ATCGA.
    //   var_hap_ranges for vi=0: hap_start = hap_position_for(3) = 3. hap_end = 4.
    //   Scan ATCGA:
    //     h=2: C at hap-coord 2, G at hap-coord 3.
    //       src_h(2):  2 ∈ [3,4)? No → Reference.
    //       src_h1(3): 3 ∈ [3,4)? Yes → Variant{vi=0, hap_start=3}.
    //       hap_offset = (hpos+1) - hap_start = 3 - 3 = 0.
    //       → OnVariant { variant_index: 0, hap_offset: 0 }.
    // Hap1 keeps REF (ATCAA): no CpG.
    let reference = b"ATCAA";
    let variants = vec![VariantRecord {
        position: 3,
        ref_allele: b"A".to_vec(),
        alt_alleles: vec![b"G".to_vec()],
        genotype: Genotype::parse("1|0").unwrap(),
    }];
    let placements = classify_cpgs(reference, &variants).unwrap();
    assert_eq!(placements, vec![CpgPlacement::OnVariant { variant_index: 0, hap_offset: 0 }],);
}

#[test]
fn case_03_deletion_juxtaposes_ref_c_and_ref_g() {
    use holodeck_lib::vcf::genotype::{Genotype, VariantRecord};
    // ref: CATTG. Deletion at pos 0: REF=CATT, ALT=C (anchor C, delete ATT). GT=1|0.
    //
    // Hap0 carries ALT: cap = hap_position_for(5). var_end=4, delta=1-4=-3. hp=5-3=2. cap=2.
    //   At ref_pos=0: variant matches → emit C; ref_pos += 4 = 4.
    //   At ref_pos=4: emit G (ref[4]).
    //   Hap0 bases: C,G. ref_positions: 0,4.
    //   var_hap_ranges for vi=0: hap_start = hap_position_for(0) = 0. alt_bases=C, len=1. hap_end=1.
    //   Scan CG:
    //     h=0: C,G → CpG!
    //       src_h(0):  0 ∈ [0,1) → Variant{vi=0, hap_start=0}.
    //       src_h1(1): 1 ∈ [0,1)? 1 < 1 is false → Reference.
    //       hap_offset = 0 - 0 = 0.
    //       → OnVariant { variant_index: 0, hap_offset: 0 }.
    // Hap1 keeps REF (CATTG): no CpG (C,A,T,T,G).
    let reference = b"CATTG";
    let variants = vec![VariantRecord {
        position: 0,
        ref_allele: b"CATT".to_vec(),
        alt_alleles: vec![b"C".to_vec()],
        genotype: Genotype::parse("1|0").unwrap(),
    }];
    let placements = classify_cpgs(reference, &variants).unwrap();
    assert_eq!(placements, vec![CpgPlacement::OnVariant { variant_index: 0, hap_offset: 0 }],);
}

#[test]
fn case_05_insertion_splits_reference_cpg() {
    use holodeck_lib::vcf::genotype::{Genotype, VariantRecord};
    // ref: ACGT. CpG at top-C pos 1. Insertion at pos 1: REF=C, ALT=CA. GT=1|0.
    //
    // Hap0 carries ALT: A,C,A,G,T = ACAGT. No CpG (C at hap-coord 1 is followed by A, not G).
    // Hap1 keeps REF (ACGT): CpG at pos 1 → Standalone{1}.
    let reference = b"ACGT";
    let variants = vec![VariantRecord {
        position: 1,
        ref_allele: b"C".to_vec(),
        alt_alleles: vec![b"CA".to_vec()],
        genotype: Genotype::parse("1|0").unwrap(),
    }];
    let placements = classify_cpgs(reference, &variants).unwrap();
    assert_eq!(placements, vec![CpgPlacement::Standalone { ref_pos: 1 }]);
}

#[test]
fn case_06_deletion_removes_reference_cpg_het() {
    use holodeck_lib::vcf::genotype::{Genotype, VariantRecord};
    // ref: ACGT. CpG at pos 1. Deletion at pos 1: REF=CG, ALT=C (anchor C, delete G). GT=1|0.
    //
    // Hap0 carries ALT: A,C,T = ACT. No CpG.
    // Hap1 keeps REF (ACGT): CpG at pos 1 → Standalone{1}.
    let reference = b"ACGT";
    let variants = vec![VariantRecord {
        position: 1,
        ref_allele: b"CG".to_vec(),
        alt_alleles: vec![b"C".to_vec()],
        genotype: Genotype::parse("1|0").unwrap(),
    }];
    let placements = classify_cpgs(reference, &variants).unwrap();
    assert_eq!(placements, vec![CpgPlacement::Standalone { ref_pos: 1 }]);
}

#[test]
fn case_06_deletion_removes_reference_cpg_hom_alt() {
    use holodeck_lib::vcf::genotype::{Genotype, VariantRecord};
    // Same as case_06_het but GT=1|1 — both haps lose the CpG.
    let reference = b"ACGT";
    let variants = vec![VariantRecord {
        position: 1,
        ref_allele: b"CG".to_vec(),
        alt_alleles: vec![b"C".to_vec()],
        genotype: Genotype::parse("1|1").unwrap(),
    }];
    let placements = classify_cpgs(reference, &variants).unwrap();
    assert!(placements.is_empty());
}

#[test]
fn case_08_variant_destroys_then_recreates_cpg() {
    use holodeck_lib::vcf::genotype::{Genotype, VariantRecord};
    // ref: AAGT. SNP1 pos 0 A→C, SNP2 pos 1 A→G. Both GT=1|0.
    //
    // Hap0 carries both ALTs: C(pos0), G(pos1), G(pos2), T(pos3) = CGGT.
    //   var_hap_ranges:
    //     vi=0 (pos 0, alt=C, ref_len=1): hap_start=0, hap_end=1.
    //     vi=1 (pos 1, alt=G, ref_len=1): hap_start=hap_position_for(1)=1 (var_end=1 of SNP1
    //       satisfies end<=1 → delta=0). hap_start=1. hap_end=2.
    //   Scan CGGT:
    //     h=0: C,G → CpG!
    //       src_h(0):  0 ∈ [0,1) → Variant{vi=0, hap_start=0}.
    //       src_h1(1): 1 ∈ [1,2) → Variant{vi=1, hap_start=1}.
    //       Different variants: vi0=0 < vi1=1 → upstream wins.
    //       hap_offset = 0 - 0 = 0. → OnVariant { variant_index: 0, hap_offset: 0 }.
    // Hap1 keeps REF (AAGT): no CpG.
    let reference = b"AAGT";
    let variants = vec![
        VariantRecord {
            position: 0,
            ref_allele: b"A".to_vec(),
            alt_alleles: vec![b"C".to_vec()],
            genotype: Genotype::parse("1|0").unwrap(),
        },
        VariantRecord {
            position: 1,
            ref_allele: b"A".to_vec(),
            alt_alleles: vec![b"G".to_vec()],
            genotype: Genotype::parse("1|0").unwrap(),
        },
    ];
    let placements = classify_cpgs(reference, &variants).unwrap();
    assert_eq!(placements, vec![CpgPlacement::OnVariant { variant_index: 0, hap_offset: 0 }],);
}

#[test]
fn case_09_multi_allelic_each_haplotype_carries_different_alt() {
    use holodeck_lib::vcf::genotype::{Genotype, VariantRecord};
    // ref: ATG. Variant at pos 1 with ALTs=[C, A]. GT=1|2.
    //
    // Hap0 (allele_index=0) carries allele 1 → alt_alleles[0] = C.
    //   Materialized: A,C,G = ACG. CpG at h=1.
    //     src_h(1):  1 ∈ [1,2) → Variant{vi=0, hap_start=1}.
    //     src_h1(2): 2 ∈ [1,2)? No → Reference.
    //     hap_offset = 1 - 1 = 0. → OnVariant { variant_index: 0, hap_offset: 0 }.
    // Hap1 (allele_index=1) carries allele 2 → alt_alleles[1] = A.
    //   Materialized: A,A,G = AAG. No CpG.
    let reference = b"ATG";
    let variants = vec![VariantRecord {
        position: 1,
        ref_allele: b"T".to_vec(),
        alt_alleles: vec![b"C".to_vec(), b"A".to_vec()],
        genotype: Genotype::parse("1|2").unwrap(),
    }];
    let placements = classify_cpgs(reference, &variants).unwrap();
    assert_eq!(placements, vec![CpgPlacement::OnVariant { variant_index: 0, hap_offset: 0 }],);
}

#[test]
fn case_10_homozygous_alt_symmetric_state() {
    use holodeck_lib::vcf::genotype::{Genotype, VariantRecord};
    // ref: ATG. Variant at pos 1: T→C, GT=1|1 (hom-alt).
    //
    // Both haps carry C at pos 1: A,C,G = ACG. CpG at h=1.
    //   OnVariant { variant_index: 0, hap_offset: 0 } from both haps.
    // After dedup: single OnVariant{0, 0}.
    let reference = b"ATG";
    let variants = vec![VariantRecord {
        position: 1,
        ref_allele: b"T".to_vec(),
        alt_alleles: vec![b"C".to_vec()],
        genotype: Genotype::parse("1|1").unwrap(),
    }];
    let placements = classify_cpgs(reference, &variants).unwrap();
    assert_eq!(placements, vec![CpgPlacement::OnVariant { variant_index: 0, hap_offset: 0 }],);
}

#[test]
fn case_13_hemizygous_single_haplotype() {
    use holodeck_lib::vcf::genotype::{Genotype, VariantRecord};
    // ref: ACG (CpG at pos 1). Variant at pos 1: C→G, GT=1 (haploid).
    //
    // max_ploidy=1 → one haplotype built.
    // Hap0 carries allele 1 (G at pos 1): A,G,G = AGG. No CpG.
    // No placements.
    let reference = b"ACG";
    let variants = vec![VariantRecord {
        position: 1,
        ref_allele: b"C".to_vec(),
        alt_alleles: vec![b"G".to_vec()],
        genotype: Genotype::parse("1").unwrap(),
    }];
    let placements = classify_cpgs(reference, &variants).unwrap();
    assert!(placements.is_empty());
}

// ── Task 7: rejection paths ──────────────────────────────────────────────────

#[test]
fn case_11_unphased_gt_rejected() {
    use holodeck_lib::vcf::genotype::{Genotype, VariantRecord};
    // ref: AAATTAA; SNP at pos 3: T→C with unphased GT=0/1.
    // Classifier must reject because methylation truth is haplotype-specific
    // and unphased genotypes have no defined haplotype assignment.
    let reference = b"AAATTAA";
    let variants = vec![VariantRecord {
        position: 3,
        ref_allele: b"T".to_vec(),
        alt_alleles: vec![b"C".to_vec()],
        genotype: Genotype::parse("0/1").unwrap(), // unphased
    }];
    let err = classify_cpgs(reference, &variants).unwrap_err();
    assert!(matches!(err, ClassifyError::UnphasedGenotype { .. }));
}

#[test]
fn case_12_overlapping_variants_rejected() {
    use holodeck_lib::vcf::genotype::{Genotype, VariantRecord};
    // ref: ACGTACGT; two variants with overlapping REF spans.
    // Variant 0: pos 1, REF=CGT (spans [1,4)); Variant 1: pos 2, REF=GT (starts at 2, inside [1,4)).
    // Classifier must reject to avoid ambiguous methylation ownership.
    let reference = b"ACGTACGT";
    let variants = vec![
        VariantRecord {
            position: 1,
            ref_allele: b"CGT".to_vec(),
            alt_alleles: vec![b"C".to_vec()],
            genotype: Genotype::parse("1|0").unwrap(),
        },
        VariantRecord {
            position: 2, // overlaps the previous variant's REF span [1, 4)
            ref_allele: b"GT".to_vec(),
            alt_alleles: vec![b"G".to_vec()],
            genotype: Genotype::parse("1|0").unwrap(),
        },
    ];
    let err = classify_cpgs(reference, &variants).unwrap_err();
    assert!(matches!(err, ClassifyError::OverlappingVariants { .. }));
}

#[test]
fn phased_ref_span_overlap_on_disjoint_haplotypes_is_allowed() {
    use holodeck_lib::vcf::genotype::{Genotype, VariantRecord};
    // Two phased variants whose REF spans overlap in reference coordinates,
    // but the alt alleles are on disjoint haplotypes: GT 1|0 + GT 0|1.
    // Hap0 carries only the first variant's alt; hap1 carries only the
    // second's. Neither haplotype materializes both alt spans, so the
    // overlap is benign and the classifier should accept it.
    let reference = b"ACGTACGT";
    let variants = vec![
        VariantRecord {
            position: 1,
            ref_allele: b"CGT".to_vec(),
            alt_alleles: vec![b"C".to_vec()],
            genotype: Genotype::parse("1|0").unwrap(),
        },
        VariantRecord {
            position: 2, // REF span [2, 4) overlaps the previous [1, 4)
            ref_allele: b"GT".to_vec(),
            alt_alleles: vec![b"G".to_vec()],
            genotype: Genotype::parse("0|1").unwrap(),
        },
    ];
    let _placements = classify_cpgs(reference, &variants)
        .expect("phased overlap on disjoint haplotypes should be accepted");
}

#[test]
fn non_adjacent_overlap_on_shared_haplotype_is_rejected() {
    use holodeck_lib::vcf::genotype::{Genotype, VariantRecord};
    // Adjacent-pair scanning misses this: var0 has a long REF span on hap0
    // that overlaps var2 (also on hap0), but the intervening var1 sits on
    // hap1 only. Window(var0, var1) sees disjoint haplotypes; window(var1,
    // var2) is fully disjoint in reference coordinates and gets skipped.
    // The genuine conflict — hap0 carries both var0 and var2 alts whose
    // REF spans overlap — must still be reported.
    let reference = b"AAAAAAAAAAAA";
    let variants = vec![
        VariantRecord {
            position: 1,
            ref_allele: b"AAAAAAAA".to_vec(), // REF span [1, 9) on hap0
            alt_alleles: vec![b"A".to_vec()],
            genotype: Genotype::parse("1|0").unwrap(),
        },
        VariantRecord {
            position: 3,
            ref_allele: b"A".to_vec(), // REF span [3, 4) on hap1 only
            alt_alleles: vec![b"C".to_vec()],
            genotype: Genotype::parse("0|1").unwrap(),
        },
        VariantRecord {
            position: 5,
            ref_allele: b"A".to_vec(), // REF span [5, 6) on hap0, inside var0's hap0 span
            alt_alleles: vec![b"C".to_vec()],
            genotype: Genotype::parse("1|0").unwrap(),
        },
    ];
    let err = classify_cpgs(reference, &variants).unwrap_err();
    assert!(matches!(err, ClassifyError::OverlappingVariants { .. }));
}

#[test]
fn phased_ref_span_overlap_sharing_a_haplotype_is_rejected() {
    use holodeck_lib::vcf::genotype::{Genotype, VariantRecord};
    // Same overlap geometry as the disjoint-haplotypes case, but both
    // variants are phased so that hap0 carries both alts. That re-creates
    // the genuine conflict the classifier is meant to reject.
    let reference = b"ACGTACGT";
    let variants = vec![
        VariantRecord {
            position: 1,
            ref_allele: b"CGT".to_vec(),
            alt_alleles: vec![b"C".to_vec()],
            genotype: Genotype::parse("1|0").unwrap(),
        },
        VariantRecord {
            position: 2,
            ref_allele: b"GT".to_vec(),
            alt_alleles: vec![b"G".to_vec()],
            genotype: Genotype::parse("1|0").unwrap(), // shares hap0 with var 0
        },
    ];
    let err = classify_cpgs(reference, &variants).unwrap_err();
    assert!(matches!(err, ClassifyError::OverlappingVariants { .. }));
}
