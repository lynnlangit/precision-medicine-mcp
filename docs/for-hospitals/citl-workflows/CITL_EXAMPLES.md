# Clinician-in-the-Loop (CitL) Review Examples

**Purpose:** Complete example reviews demonstrating APPROVE, REVISE, and REJECT scenarios
**Target Audience:** Oncologists, pathologists, radiologists learning the CitL workflow
**Version:** 1.0
**Last Updated:** 2026-01-13

---

## Table of Contents

1. [Example 1: APPROVE - PatientOne Clean Approval](#example-1-approve---patientone-clean-approval)
2. [Example 2: REVISE - Quality Flag Requires Action](#example-2-revise---quality-flag-requires-action)
3. [Example 3: REJECT - Critical Data Error](#example-3-reject---critical-data-error)
4. [Comparison Table](#comparison-table)
5. [Key Takeaways](#key-takeaways)

---

## Example 1: APPROVE - PatientOne Clean Approval

**Scenario:** High-grade serous ovarian carcinoma (HGSOC), Stage IV, platinum-resistant recurrence. All quality checks passed, findings consistent with clinical presentation, NCCN-aligned treatment recommendations.

**Review Time:** 25 minutes
**Expected Outcome:** Final approved report ready for tumor board

---

### Complete Review Form

```json
{
  "patient_id": "PAT001-OVC-2025",
  "report_date": "2026-01-13T14:00:00Z",
  "reviewer": {
    "name": "Dr. Sarah Johnson",
    "email": "sarah.johnson@hospital.org",
    "credentials": "MD, Gynecologic Oncology",
    "role": "oncologist"
  },
  "review_date": "2026-01-13T14:30:00Z",

  "decision": {
    "status": "APPROVE",
    "rationale": "All findings are consistent with clinical presentation and imaging. Molecular results align with observed platinum resistance (ABCB1/MDR1 overexpression, BCL2L1 anti-apoptotic signaling). Treatment recommendations follow NCCN guidelines for recurrent HGSOC. Quality checks passed with one minor warning (tumor_necrotic region 35 spots) that does not affect clinical validity given strong effect sizes and corroborating imaging."
  },

  "per_finding_validation": [
    {
      "finding_id": "DEG_1",
      "gene": "TP53",
      "validation_status": "CONFIRMED",
      "comments": "TP53 R175H is a well-characterized hotspot mutation present in >95% of HGSOC cases. Consistent with p53 IHC staining pattern (diffuse nuclear positivity). Very high confidence (FDR = 5.04e-20). CONFIRMED."
    },
    {
      "finding_id": "DEG_2",
      "gene": "PIK3CA",
      "validation_status": "CONFIRMED",
      "comments": "PIK3CA E545K activation mutation (helical domain). Found in ~10-15% of HGSOC cases. Consistent with PI3K/AKT pathway upregulation observed in pathway analysis. High confidence (FDR = 2.31e-18). Actionable target for PI3K inhibitor therapy. CONFIRMED."
    },
    {
      "finding_id": "DEG_3",
      "gene": "BRCA1",
      "validation_status": "CONFIRMED",
      "comments": "BRCA1 germline mutation confirmed by prior genetic testing (GeneDx report 2024-03-15, variant c.5266dupC, p.Gln1756Profs). Consistent with homologous recombination deficiency signature. Rationale for PARP inhibitor therapy. CONFIRMED."
    },
    {
      "finding_id": "DEG_4",
      "gene": "ABCB1",
      "validation_status": "CONFIRMED",
      "comments": "ABCB1 (MDR1) overexpression consistent with platinum resistance. Patient has platinum-free interval of 4 months (platinum-resistant per GOG criteria). Log2FC = 3.8, FDR = 1.2e-15. Strong biological rationale. CONFIRMED."
    },
    {
      "finding_id": "DEG_5",
      "gene": "BCL2L1",
      "validation_status": "CONFIRMED",
      "comments": "BCL2L1 overexpression indicates anti-apoptotic signaling. Mechanism of therapy resistance. Log2FC = 3.2, FDR = 8.7e-14. Consistent with apoptosis resistance phenotype. Potential target for BCL2 inhibitors. CONFIRMED."
    },
    {
      "finding_id": "DEG_6",
      "gene": "CD8A",
      "validation_status": "CONFIRMED",
      "comments": "CD8A downregulation in tumor core regions consistent with immune exclusion phenotype observed on imaging (lack of T-cell infiltration on IHC). Log2FC = -2.9, FDR = 4.3e-12. Explains poor response to prior checkpoint inhibitor. CONFIRMED."
    },
    {
      "finding_id": "DEG_7",
      "gene": "VEGFA",
      "validation_status": "CONFIRMED",
      "comments": "VEGFA overexpression consistent with angiogenic phenotype. Patient responded to bevacizumab previously. Log2FC = 2.7, FDR = 1.8e-11. Rationale for continued anti-angiogenic therapy. CONFIRMED."
    },
    {
      "finding_id": "DEG_8",
      "gene": "MKI67",
      "validation_status": "CONFIRMED",
      "comments": "MKI67 overexpression indicates high proliferation rate. Consistent with Ki67 IHC (75% positivity). Log2FC = 3.5, FDR = 9.2e-13. Aggressive phenotype. CONFIRMED."
    },
    {
      "finding_id": "DEG_9",
      "gene": "AKT1",
      "validation_status": "CONFIRMED",
      "comments": "AKT1 upregulation consistent with PI3K/AKT pathway activation (PIK3CA E545K mutation). Log2FC = 2.3, FDR = 7.4e-10. Mechanistically linked to PIK3CA mutation. CONFIRMED."
    },
    {
      "finding_id": "DEG_10",
      "gene": "MTOR",
      "validation_status": "CONFIRMED",
      "comments": "MTOR upregulation downstream of PI3K/AKT activation. Log2FC = 1.9, FDR = 2.1e-9. Consistent with pathway analysis showing mTOR signaling activation. Potential combination target (PI3K + mTOR inhibitor). CONFIRMED."
    }
  ],

  "guideline_compliance": {
    "nccn_aligned": "ALIGNED",
    "nccn_deviations": [],
    "institutional_aligned": "ALIGNED",
    "institutional_deviations": []
  },

  "quality_flags_assessment": [
    {
      "flag_id": "sample_size_warning",
      "severity": "warning",
      "reviewer_assessment": "ACCEPTABLE",
      "comments": "Tumor_necrotic region has 35 spots (below ideal 50 threshold). However, strong effect sizes (log2FC > 3) and very low FDR values (< 1e-10) provide sufficient confidence. Additionally, necrotic regions corroborated by CT imaging (hypodense areas consistent with necrosis). Findings from this region are not critical to treatment decisions. ACCEPTABLE."
    }
  ],

  "treatment_recommendations_review": [
    {
      "recommendation_id": "REC_1",
      "therapy_name": "PI3K inhibitor (Alpelisib) + PARP inhibitor (Olaparib)",
      "agreement": "AGREE",
      "comments": "Combination is rational given PIK3CA E545K mutation + BRCA1 germline deficiency. Patient meets eligibility for clinical trial NCT04729387 (alpelisib + olaparib in PIK3CA-mutant, BRCA-deficient ovarian cancer). Platinum-free interval is 4 months (platinum-resistant), making PARP inhibitor monotherapy less effective. Combination addresses both PI3K activation and HR deficiency. AGREE - recommend discussing trial enrollment."
    },
    {
      "recommendation_id": "REC_2",
      "therapy_name": "Bevacizumab continuation",
      "agreement": "AGREE",
      "comments": "Patient responded to bevacizumab in prior line (PFS 9 months, above median). VEGFA overexpression provides biological rationale. No contraindications (no GI perforation history, no recent surgery). Duration is 8 cycles (within institutional protocol limit of 15 cycles). AGREE."
    },
    {
      "recommendation_id": "REC_3",
      "therapy_name": "Liposomal doxorubicin (alternative option)",
      "agreement": "AGREE",
      "comments": "NCCN Category 1 recommendation for platinum-resistant ovarian cancer. Appropriate as alternative if patient declines trial or is ineligible. Cardiac function adequate (LVEF 60%). AGREE as backup option."
    }
  ],

  "attestation": {
    "reviewed_all_findings": true,
    "assessed_compliance": true,
    "clinical_judgment": true,
    "medical_record_acknowledgment": true,
    "signature_hash": "a1b2c3d4e5f6789012345678901234567890123456789012345678901234abcd",
    "timestamp": "2026-01-13T14:55:00Z"
  },

  "revision_count": 0
}
```

---

### Key Decision Points

**Why APPROVE?**
1. ✅ All 10 findings CONFIRMED with strong clinical rationale
2. ✅ Quality checks passed (one minor warning acceptable)
3. ✅ NCCN-aligned treatment recommendations
4. ✅ Findings consistent with clinical presentation (platinum resistance, imaging, IHC)
5. ✅ Confident presenting at tumor board

**Review Time Breakdown:**
- Clinical summary review: 5 minutes
- Per-finding validation (10 findings): 12 minutes
- Guideline compliance: 3 minutes
- Quality flags: 2 minutes
- Treatment recommendations: 3 minutes
- **Total: 25 minutes**

**Next Steps:**
1. Run finalization script:
   ```bash
   python tools/reports/finalize_patient_report.py --patient-id PAT001-OVC-2025
   ```
2. Present at tumor board (Tuesday)
3. Discuss clinical trial enrollment with patient
4. Document decision in EHR

---

## Example 2: REVISE - Quality Flag Requires Action

**Scenario:** Breast cancer, ER+/HER2-, Stage III. Quality flag raised for inadequate sample size in one region. Findings from that region have marginal FDR values. Clinician determines re-analysis is needed with region excluded.

**Review Time:** 30 minutes
**Expected Outcome:** Re-analysis with adjusted parameters, resubmission for review

---

### Complete Review Form

```json
{
  "patient_id": "PAT002-BRC-2025",
  "report_date": "2026-01-12T10:00:00Z",
  "reviewer": {
    "name": "Dr. Michael Chen",
    "email": "michael.chen@hospital.org",
    "credentials": "MD, PhD, Medical Oncology",
    "role": "oncologist"
  },
  "review_date": "2026-01-12T15:45:00Z",

  "decision": {
    "status": "REVISE",
    "rationale": "Quality flag for tumor_edge region (28 spots) requires action. Sample size is below minimum threshold of 30 spots, and FDR values from this region are marginal (0.02-0.04), indicating insufficient statistical power. Findings from adequate regions (tumor_core: 420 spots, stroma: 380 spots) are valid. Request re-analysis excluding tumor_edge region and increasing minimum spots threshold to 50 to prevent similar issues."
  },

  "per_finding_validation": [
    {
      "finding_id": "DEG_1",
      "gene": "ESR1",
      "validation_status": "CONFIRMED",
      "comments": "ESR1 (estrogen receptor) high expression consistent with ER+ status (IHC: 95% positive). From tumor_core region (420 spots). FDR = 1.2e-25. CONFIRMED."
    },
    {
      "finding_id": "DEG_2",
      "gene": "PGR",
      "validation_status": "CONFIRMED",
      "comments": "PGR (progesterone receptor) moderate expression consistent with PR+ status (IHC: 60% positive). From tumor_core region. FDR = 3.4e-18. CONFIRMED."
    },
    {
      "finding_id": "DEG_3",
      "gene": "ERBB2",
      "validation_status": "CONFIRMED",
      "comments": "ERBB2 (HER2) low expression consistent with HER2- status (IHC 0). From tumor_core region. FDR = 5.1e-20. CONFIRMED."
    },
    {
      "finding_id": "DEG_4",
      "gene": "MKI67",
      "validation_status": "UNCERTAIN",
      "comments": "MKI67 upregulation from tumor_edge region (28 spots). FDR = 0.024 (marginal). Ki67 IHC shows heterogeneity (15-30% across tumor). Insufficient statistical power in small region. Mark UNCERTAIN pending re-analysis with larger sample."
    },
    {
      "finding_id": "DEG_5",
      "gene": "CCND1",
      "validation_status": "UNCERTAIN",
      "comments": "CCND1 overexpression from tumor_edge region (28 spots). FDR = 0.037 (marginal). Cyclin D1 IHC inconclusive. Insufficient power. Mark UNCERTAIN pending re-analysis."
    },
    {
      "finding_id": "DEG_6",
      "gene": "CDK4",
      "validation_status": "UNCERTAIN",
      "comments": "CDK4 upregulation from tumor_edge region (28 spots). FDR = 0.041 (marginal). Related to CCND1 finding. Both require validation with larger sample. Mark UNCERTAIN."
    },
    {
      "finding_id": "DEG_7",
      "gene": "FOXC1",
      "validation_status": "CONFIRMED",
      "comments": "FOXC1 downregulation in tumor_core (luminal phenotype). FDR = 2.8e-12. Consistent with ER+ luminal subtype. CONFIRMED."
    },
    {
      "finding_id": "DEG_8",
      "gene": "GATA3",
      "validation_status": "CONFIRMED",
      "comments": "GATA3 high expression in tumor_core. Luminal lineage marker. FDR = 1.5e-16. CONFIRMED."
    },
    {
      "finding_id": "DEG_9",
      "gene": "BCL2",
      "validation_status": "CONFIRMED",
      "comments": "BCL2 overexpression in tumor_core. Associated with ER+ tumors. FDR = 4.2e-14. CONFIRMED."
    },
    {
      "finding_id": "DEG_10",
      "gene": "PTEN",
      "validation_status": "CONFIRMED",
      "comments": "PTEN expression retained in tumor_core. No loss of function. FDR = 8.9e-11. PI3K pathway not constitutively activated. CONFIRMED."
    }
  ],

  "guideline_compliance": {
    "nccn_aligned": "ALIGNED",
    "nccn_deviations": [],
    "institutional_aligned": "ALIGNED",
    "institutional_deviations": []
  },

  "quality_flags_assessment": [
    {
      "flag_id": "sample_size_warning",
      "severity": "warning",
      "reviewer_assessment": "REQUIRES_ACTION",
      "comments": "Tumor_edge region has only 28 spots (below minimum 30). Three findings (MKI67, CCND1, CDK4) from this region have marginal FDR values (0.024-0.041), indicating insufficient statistical power. Cannot confidently use these findings for treatment decisions. REQUIRES_ACTION: Exclude tumor_edge region and re-run differential expression analysis."
    }
  ],

  "treatment_recommendations_review": [
    {
      "recommendation_id": "REC_1",
      "therapy_name": "Endocrine therapy (Aromatase inhibitor + CDK4/6 inhibitor)",
      "agreement": "AGREE",
      "comments": "Standard NCCN Category 1 recommendation for ER+/HER2- advanced breast cancer. ESR1, PGR expression confirmed. However, CDK4 overexpression (rationale for CDK4/6 inhibitor) comes from unreliable tumor_edge region. Still AGREE with recommendation as CDK4/6 inhibitors are standard regardless of CDK4 expression. Re-analysis will strengthen molecular rationale."
    },
    {
      "recommendation_id": "REC_2",
      "therapy_name": "Chemotherapy (alternative if endocrine-resistant)",
      "agreement": "AGREE",
      "comments": "Appropriate alternative if endocrine therapy fails. MKI67 status from tumor_edge region is uncertain, but clinical features (node-positive, large tumor) suggest intermediate-high risk. AGREE."
    }
  ],

  "attestation": {
    "reviewed_all_findings": true,
    "assessed_compliance": true,
    "clinical_judgment": true,
    "medical_record_acknowledgment": true,
    "signature_hash": "b2c3d4e5f678901234567890123456789012345678901234567890123456bcde",
    "timestamp": "2026-01-12T16:15:00Z"
  },

  "revision_instructions": {
    "issues_to_address": [
      "Exclude tumor_edge region (28 spots) due to insufficient sample size. FDR values from this region (0.024-0.041) are unreliable.",
      "Re-run differential expression analysis with minimum spot threshold increased to 50 to prevent similar issues in future analyses.",
      "Three findings (MKI67, CCND1, CDK4) from tumor_edge region should be re-evaluated after excluding this region. These may still be detected in adjacent regions if truly present."
    ],
    "reanalysis_parameters": {
      "regions_to_exclude": ["tumor_edge"],
      "min_spots_per_region": 50,
      "fdr_threshold": 0.05,
      "additional_tests_required": []
    },
    "resubmission_date": "2026-01-15"
  },

  "revision_count": 0
}
```

---

### Key Decision Points

**Why REVISE (not APPROVE)?**
1. ⚠️ Quality flag requires action (sample size < 30)
2. ⚠️ Three findings UNCERTAIN due to marginal FDR values
3. ⚠️ Insufficient statistical power in one region
4. ✅ Other findings valid, no NCCN deviations
5. ⚠️ Cannot confidently use uncertain findings for treatment decisions

**Why REVISE (not REJECT)?**
- Majority of findings (7/10) are CONFIRMED from adequate regions
- Issue is localized to one small region
- Clear path to resolution (exclude problematic region)
- Treatment recommendations still appropriate
- No critical errors or data integrity issues

**Review Time Breakdown:**
- Clinical summary review: 6 minutes
- Per-finding validation: 14 minutes (extra time assessing tumor_edge findings)
- Guideline compliance: 3 minutes
- Quality flags: 3 minutes (assessing impact on validity)
- Treatment recommendations: 2 minutes
- Revision instructions (Section 7): 4 minutes
- **Total: 32 minutes**

**Next Steps:**
1. Bioinformatics team receives revision instructions
2. Re-run analysis excluding tumor_edge region, min spots = 50
3. Generate new draft report (version 2)
4. Clinician reviews version 2 (expected time: 15-20 minutes, faster than initial)
5. Likely outcome: APPROVE on version 2

---

## Example 3: REJECT - Critical Data Error

**Scenario:** Lung adenocarcinoma, Stage IV. Report shows EGFR exon 19 deletion, but patient's prior molecular testing (6 months ago, FoundationOne) showed EGFR wild-type. Critical inconsistency suggests sample mix-up or data error.

**Review Time:** 35 minutes
**Expected Outcome:** Escalation to PI/bioinformatics team, investigation, complete re-analysis

---

### Complete Review Form

```json
{
  "patient_id": "PAT003-NSCLC-2025",
  "report_date": "2026-01-11T09:00:00Z",
  "reviewer": {
    "name": "Dr. Jennifer Martinez",
    "email": "jennifer.martinez@hospital.org",
    "credentials": "MD, Thoracic Oncology",
    "role": "oncologist"
  },
  "review_date": "2026-01-11T14:20:00Z",

  "decision": {
    "status": "REJECT",
    "rationale": "Critical inconsistency: Report identifies EGFR exon 19 deletion (p.E746_A750del), but patient's prior molecular testing (FoundationOne, 2025-07-15, report #FMI-12345) showed EGFR wild-type (no mutations, amplifications, or deletions). Patient has not received EGFR TKI therapy that could select for resistant clones. Suspect sample mix-up (wrong patient specimen) or critical data processing error. Cannot proceed with clinical decision-making until discrepancy is resolved. Recommend verifying sample identity (SNP fingerprinting) and investigating EGFR call. If sample identity confirmed, re-run from raw sequencing data."
  },

  "per_finding_validation": [
    {
      "finding_id": "MUT_1",
      "gene": "EGFR",
      "validation_status": "INCORRECT",
      "comments": "EGFR exon 19 deletion (p.E746_A750del) contradicts prior FoundationOne testing (2025-07-15) which showed EGFR wild-type. Patient has not received EGFR TKI therapy. No biological explanation for acquiring EGFR mutation in 6 months without TKI selection pressure. Suspect sample mix-up or data error. Mark INCORRECT pending investigation."
    },
    {
      "finding_id": "MUT_2",
      "gene": "TP53",
      "validation_status": "UNCERTAIN",
      "comments": "TP53 R273H reported. Prior FoundationOne testing showed TP53 R248Q (different hotspot). Both are gain-of-function mutations, but different amino acid positions. If sample identity is correct, this suggests tumor evolution or sampling from different clone. If sample mix-up, this is from different patient. Mark UNCERTAIN pending sample verification."
    },
    {
      "finding_id": "MUT_3",
      "gene": "KRAS",
      "validation_status": "UNCERTAIN",
      "comments": "KRAS wild-type reported. Consistent with prior testing (KRAS WT). However, EGFR and KRAS mutations are typically mutually exclusive. If EGFR deletion is real, KRAS WT is expected. But EGFR deletion contradicts prior data. Overall validity UNCERTAIN pending investigation."
    },
    {
      "finding_id": "DEG_4",
      "gene": "ALK",
      "validation_status": "UNCERTAIN",
      "comments": "ALK not rearranged (consistent with prior testing). Cannot validate other findings until sample identity confirmed. Mark UNCERTAIN."
    },
    {
      "finding_id": "DEG_5",
      "gene": "ROS1",
      "validation_status": "UNCERTAIN",
      "comments": "ROS1 not rearranged (consistent with prior testing). Cannot validate until sample verified. Mark UNCERTAIN."
    },
    {
      "finding_id": "DEG_6",
      "gene": "BRAF",
      "validation_status": "UNCERTAIN",
      "comments": "BRAF wild-type. Consistent with prior testing. Mark UNCERTAIN pending resolution."
    },
    {
      "finding_id": "DEG_7",
      "gene": "MET",
      "validation_status": "UNCERTAIN",
      "comments": "MET amplification not detected. Consistent with prior testing. Mark UNCERTAIN."
    },
    {
      "finding_id": "DEG_8",
      "gene": "RET",
      "validation_status": "UNCERTAIN",
      "comments": "RET not rearranged. Consistent with prior. Mark UNCERTAIN."
    },
    {
      "finding_id": "DEG_9",
      "gene": "ERBB2",
      "validation_status": "UNCERTAIN",
      "comments": "ERBB2 (HER2) no mutations detected. Consistent with prior. Mark UNCERTAIN."
    },
    {
      "finding_id": "DEG_10",
      "gene": "NTRK1",
      "validation_status": "UNCERTAIN",
      "comments": "NTRK not rearranged. Consistent with prior. Mark UNCERTAIN."
    }
  ],

  "guideline_compliance": {
    "nccn_aligned": "NOT_ALIGNED",
    "nccn_deviations": [
      "Treatment recommendation for EGFR TKI (osimertinib) is based on EGFR exon 19 deletion that contradicts prior molecular testing. Cannot assess guideline compliance until EGFR discrepancy is resolved."
    ],
    "institutional_aligned": "NOT_ALIGNED",
    "institutional_deviations": [
      "EGFR TKI therapy contradicts institutional protocol requiring confirmed EGFR mutation. Prior testing showed wild-type."
    ]
  },

  "quality_flags_assessment": [
    {
      "flag_id": "mutation_concordance_error",
      "severity": "critical",
      "reviewer_assessment": "REQUIRES_ACTION",
      "comments": "CRITICAL: EGFR mutation call discordant with prior molecular testing. This is not a quality flag from automated checks, but a critical finding discovered during clinical review. REQUIRES_ACTION: Immediate investigation required."
    }
  ],

  "treatment_recommendations_review": [
    {
      "recommendation_id": "REC_1",
      "therapy_name": "EGFR TKI (Osimertinib)",
      "agreement": "DISAGREE",
      "comments": "Recommendation is based on EGFR exon 19 deletion that contradicts prior testing. Cannot administer EGFR TKI based on unverified molecular finding. DISAGREE - patient should NOT receive osimertinib until EGFR status is definitively confirmed."
    },
    {
      "recommendation_id": "REC_2",
      "therapy_name": "Platinum-based chemotherapy (alternative)",
      "agreement": "AGREE",
      "comments": "Standard NCCN recommendation for EGFR wild-type NSCLC (per prior FoundationOne testing). This remains appropriate treatment if EGFR mutation is ruled out. AGREE as appropriate alternative."
    }
  ],

  "attestation": {
    "reviewed_all_findings": true,
    "assessed_compliance": true,
    "clinical_judgment": true,
    "medical_record_acknowledgment": true,
    "signature_hash": "c3d4e5f6789012345678901234567890123456789012345678901234567cdef",
    "timestamp": "2026-01-11T14:55:00Z"
  },

  "revision_instructions": {
    "issues_to_address": [
      "CRITICAL: Verify sample identity. Confirm correct patient specimen was analyzed. Perform SNP fingerprinting to match sample to patient germline.",
      "Investigate EGFR exon 19 deletion call. Reconcile discrepancy with prior FoundationOne testing (2025-07-15, EGFR wild-type).",
      "If sample identity confirmed correct: Investigate technical reasons for EGFR discrepancy (sequencing error, variant caller settings, reference genome mismatch). Consider orthogonal validation (Sanger sequencing, ddPCR) for EGFR locus.",
      "If sample mix-up confirmed: Identify correct patient specimen, re-run complete analysis pipeline with correct sample.",
      "Review TP53 discrepancy (R273H vs R248Q from prior testing). If sample correct, this may indicate tumor evolution or multi-clonal disease.",
      "Do NOT proceed with EGFR TKI therapy until EGFR status definitively confirmed. Patient harm could result from inappropriate targeted therapy."
    ],
    "reanalysis_parameters": {
      "regions_to_exclude": [],
      "additional_tests_required": [
        "SNP fingerprinting for sample identity verification",
        "Orthogonal EGFR sequencing (Sanger or ddPCR) for exon 19 region",
        "Re-run from raw FASTQ files if data processing error suspected",
        "Review variant caller logs for EGFR locus",
        "Compare with prior FoundationOne VCF file if available"
      ]
    },
    "resubmission_date": "2026-01-25"
  },

  "revision_count": 0
}
```

---

### Key Decision Points

**Why REJECT (not REVISE)?**
1. ❌ **Critical error:** EGFR mutation contradicts prior definitive testing
2. ❌ **Patient safety:** Recommending inappropriate targeted therapy (EGFR TKI)
3. ❌ **Data integrity concern:** Suspect sample mix-up (most serious error)
4. ❌ **Cannot use ANY findings:** Validity of entire report is questionable
5. ❌ **Requires investigation:** Not just parameter adjustment, need root cause analysis

**Why REJECT (not just flag and continue)?**
- Acting on incorrect molecular data could harm patient (inappropriate therapy, missed correct therapy)
- Sample mix-up affects ALL findings, not just one gene
- Regulatory/ethical requirement to ensure correct patient-sample matching
- EGFR status is actionable and treatment-determining in NSCLC

**Review Time Breakdown:**
- Clinical summary review: 7 minutes
- Identifying EGFR discrepancy: 5 minutes (reviewing prior reports)
- Per-finding validation: 10 minutes (marking as UNCERTAIN due to sample concern)
- Guideline compliance: 3 minutes
- Treatment recommendations: 3 minutes
- Revision instructions (Section 7): 10 minutes (detailed investigation plan)
- **Total: 38 minutes**

**Next Steps:**
1. **Immediate:** Alert bioinformatics team and PI (critical error)
2. **Day 1:** Schedule emergency meeting (clinician, bioinformatics, lab director)
3. **Day 2-3:** Sample verification (SNP fingerprinting)
4. **Day 4-7:** Investigation and root cause analysis
   - If sample mix-up → Identify correct sample, re-run analysis
   - If data error → Debug pipeline, validate EGFR call, re-run from raw data
5. **Day 10-14:** Complete re-analysis and resubmit for review
6. **Patient care:** Continue with standard platinum chemotherapy (per prior FoundationOne results) while investigation ongoing. Do NOT start EGFR TKI.

---

## Comparison Table

| Aspect | APPROVE | REVISE | REJECT |
|--------|---------|---------|---------|
| **Decision Status** | All findings valid, ready for clinical use | Findings mostly valid, needs corrections | Critical error, cannot use report |
| **Confidence Level** | High - comfortable presenting at tumor board | Medium - need adjustments before clinical use | Low - data validity in question |
| **Findings Validation** | 8-10 CONFIRMED, 0-2 UNCERTAIN, 0 INCORRECT | 5-8 CONFIRMED, 2-5 UNCERTAIN, 0-2 INCORRECT | 0-3 CONFIRMED, 5-10 UNCERTAIN, 1-3 INCORRECT |
| **Quality Flags** | All ACCEPTABLE | 1-2 REQUIRES_ACTION | Critical flags or data integrity concern |
| **Guideline Compliance** | ALIGNED or PARTIAL (justified) | ALIGNED (but molecular rationale needs strengthening) | NOT_ALIGNED (based on invalid data) |
| **Treatment Recs** | AGREE with all/most | AGREE with most (pending molecular confirmation) | DISAGREE with key recommendations |
| **Next Step** | Finalize report → Tumor board | Re-analysis → Resubmit for review | Investigation → Root cause analysis → Complete re-analysis |
| **Timeline** | Immediate (same day) | 3-5 days (re-analysis) | 1-2 weeks (investigation + re-analysis) |
| **Severity** | None - routine approval | Minor - needs refinement | Critical - patient safety concern |
| **Section 7 Required?** | No | Yes - specific re-analysis instructions | Yes - detailed investigation plan |
| **Review Time** | 20-25 minutes | 25-35 minutes | 30-40 minutes |
| **Typical Frequency** | 70-80% of reviews | 15-25% of reviews | 5-10% of reviews |

---

## Key Takeaways

### When to APPROVE
✅ All findings consistent with clinical presentation
✅ Quality checks passed (or minor warnings that don't affect validity)
✅ Treatment recommendations follow guidelines
✅ You're confident using this report for clinical decision-making
✅ Ready to present findings at tumor board

### When to REVISE
⚠️ Findings are mostly valid but have issues requiring correction
⚠️ Quality flags that need addressing (e.g., small sample size, borderline FDR values)
⚠️ Treatment recommendations need refinement
⚠️ Clear, actionable path to resolution (exclude regions, adjust parameters, add validation)
⚠️ Not ready for tumor board YET, but will be after corrections

### When to REJECT
❌ Critical errors in findings (wrong patient data, major inconsistencies)
❌ Data integrity concerns (sample mix-up, severe quality issues)
❌ Patient safety at risk (inappropriate therapy based on invalid data)
❌ Requires investigation and root cause analysis (not just parameter tweaking)
❌ Cannot use ANY findings from report until resolved

### Decision-Making Framework

**Ask yourself:**
1. **Would I present these findings at tumor board?**
   - Yes → APPROVE
   - Not yet, but close → REVISE
   - No, something is wrong → REJECT

2. **Would I start treatment based on these findings?**
   - Yes → APPROVE
   - Yes, but with additional validation → REVISE
   - No, data may be incorrect → REJECT

3. **What's the primary issue?**
   - No significant issues → APPROVE
   - Statistical power, minor quality flags → REVISE
   - Data validity, patient safety → REJECT

4. **How confident am I in the findings?**
   - High (>90%) → APPROVE
   - Medium (60-90%) → REVISE
   - Low (<60%) → REJECT

### Best Practices

1. **Always review prior molecular testing reports** - Catch discrepancies early
2. **Don't ignore quality flags** - Even "minor" warnings can indicate real issues
3. **Be specific in revision instructions** - "Improve quality" is not actionable
4. **Document your reasoning** - Future reviewers and auditors will read your comments
5. **Err on the side of caution** - When in doubt, REVISE rather than APPROVE
6. **Patient safety first** - REJECT if any concern about data validity or patient harm

---

**Document Version:** 1.0
**Last Updated:** 2026-01-13
**Part of:** Precision Medicine MCP - Clinician-in-the-Loop Validation Workflow
**See Also:**
- CITL_WORKFLOW_GUIDE.md (step-by-step instructions)
- CITL_REVIEW_TEMPLATE.md (blank review form)
- citl_review_schema.json (JSON schema for validation)
