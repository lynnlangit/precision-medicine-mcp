# PatientOne Full Test Prompt
## Complete End-to-End Multi-Modality Workflow

**Version:** 2.0 (Enhanced Multi-Omics)
**Patient:** PAT001-OVC-2025
**Last Updated:** December 26, 2025

---

## Overview

This document contains the complete test prompts for running the PatientOne precision medicine workflow using all 9 MCP servers and 40 tools.

**Changes in Version 2.0:**
- ⭐ TEST_2_MULTIOMICS updated to TEST_2_MULTIOMICS_ENHANCED.txt
- ⭐ Added preprocessing pipeline (3 new tools)
- ⭐ Added upstream regulator prediction (1 new tool)
- Multi-omics workflow now 8 steps instead of 4

---

## TEST 1: Clinical-Genomic Analysis

**File:** `TEST_1_CLINICAL_GENOMIC.txt`
**Servers Used:** mcp-fgbio (4 tools)
**Expected Runtime:** 15-20 minutes (full analysis mode) / 2-3 minutes (DRY_RUN mode)

### Test Prompt

```
[Paste complete content of TEST_1_CLINICAL_GENOMIC.txt]
```

### Expected Outputs

1. **Variant Call File (VCF)**
   - PIK3CA amplification (copy number = 4)
   - TP53 mutation (p.R273H)
   - PTEN deletion (homozygous loss)

2. **Quality Control Report**
   - Alignment rate: >90%
   - Mean coverage: >100x
   - Variant quality scores: PASS

3. **Clinical Interpretation**
   - High-risk genomic profile
   - PI3K pathway activation predicted
   - TP53 loss → compromised tumor suppression

---

## TEST 2: Multi-Omics Resistance Analysis ⭐ UPDATED

**File:** `TEST_2_MULTIOMICS_ENHANCED.txt` (⚠️ UPDATED - was TEST_2_MULTIOMICS.txt)
**Servers Used:** mcp-multiomics (9 tools - was 5 tools)
**Expected Runtime:** 60-90 minutes (full analysis) / 5-10 minutes (DRY_RUN mode)

### Test Prompt

```
[Paste complete content of TEST_2_MULTIOMICS_ENHANCED.txt]
```

### Workflow Steps (Enhanced)

**STEP 0: PREPROCESSING** ⭐ NEW (CRITICAL for real data)

1. **validate_multiomics_data**
   - Detect batch effects
   - Identify outliers
   - Check sample consistency
   - **Expected:** PC1-batch correlation = 0.82 (CRITICAL)

2. **preprocess_multiomics_data**
   - Apply ComBat batch correction
   - KNN imputation (k=5)
   - Outlier removal
   - Quantile normalization
   - **Expected:** PC1-batch correlation reduced to 0.12

3. **visualize_data_quality**
   - Generate before/after PCA plots
   - Verify batch correction effectiveness
   - **Expected:** PC1-batch r < 0.3 (PASS)

**STEP 1: INTEGRATION**

4. **integrate_omics_data**
   - Load PREPROCESSED data (not raw!)
   - Align 13 samples across 3 modalities
   - **Expected:** 19.5K RNA, 6.8K protein, 4.9K phospho features

**STEP 2: ASSOCIATION & META-ANALYSIS**

5. **run_halla_analysis** (optional)
   - Test RNA-protein associations
   - Chunking strategy (1000 features/chunk)
   - **Expected:** NOMINAL p-values returned

6. **calculate_stouffer_meta**
   - Combine evidence from 3 modalities
   - FDR correction applied AFTER combination
   - **Expected:** 7 resistance genes with q < 0.05

**STEP 3: UPSTREAM REGULATORS** ⭐ NEW

7. **predict_upstream_regulators**
   - Identify activated kinases
   - Predict therapeutic targets
   - Map to FDA-approved drugs
   - **Expected:** AKT1, MTOR, PI3K activated; TP53 inhibited

**STEP 4: VISUALIZATION**

8. **create_multiomics_heatmap** (optional)
9. **run_multiomics_pca** (optional)

### Expected Outputs (Enhanced)

**Preprocessing Results:**
- ✅ Batch effects detected: PC1-batch r=0.82
- ✅ Batch correction applied: PC1-batch r → 0.12
- ✅ Missing values imputed: 2000 protein + 1500 phospho
- ✅ Outliers removed: 2 samples (Sample_07, Sample_12)
- ✅ Final sample count: 13 (7 resistant, 6 sensitive)
- ✅ QC plots generated: 4 PNG files

**Gene-Level Results:**

| Gene   | RNA FC | Prot FC | Phos FC | Z-score | q-value | Direction |
|--------|--------|---------|---------|---------|---------|-----------|
| AKT1   | +2.1   | +1.9    | +2.3    | 4.5     | <0.0001 | UP ↑      |
| PIK3CA | +2.3   | +2.0    | +1.8    | 4.2     | 0.0001  | UP ↑      |
| ABCB1  | +2.5   | +2.2    | +1.9    | 4.1     | 0.0001  | UP ↑      |
| PTEN   | -2.1   | -1.9    | -1.7    | -3.9    | 0.0002  | DOWN ↓    |
| MTOR   | +1.9   | +1.7    | +1.5    | 3.8     | 0.0003  | UP ↑      |
| BCL2L1 | +1.8   | +1.6    | +1.4    | 3.2     | 0.002   | UP ↑      |
| TP53   | -1.5   | -1.3    | -1.1    | -2.8    | 0.005   | DOWN ↓    |

**Upstream Regulator Results:** ⭐ NEW

- **Activated Kinases:**
  - AKT1: Z=3.2, q=0.001
  - MTOR: Z=2.8, q=0.003
  - PI3K: Z=3.0, q=0.002

- **Inhibited Transcription Factors:**
  - TP53: Z=-3.5, q=0.0001 (loss of tumor suppression)

- **Drug Recommendations:**
  - Alpelisib (PI3K inhibitor) - FDA approved
  - Capivasertib (AKT inhibitor) - Phase III
  - Everolimus (mTOR inhibitor) - FDA approved

- **Clinical Trial Match:**
  - NCT03602859: Alpelisib + Capivasertib in PTEN-deficient tumors

**Pathway Summary:**
- ✅ PI3K/AKT/mTOR pathway ACTIVATED
- ✅ Evidence: PIK3CA, AKT1, MTOR upregulated; PTEN downregulated
- ✅ Mechanism: PTEN loss → PI3K hyperactivation → AKT/mTOR signaling
- ✅ Therapeutic strategy: Dual PI3K/AKT inhibition

---

## TEST 3: Spatial Transcriptomics Analysis

**File:** `TEST_3_SPATIAL_TRANSCRIPTOMICS.txt`
**Servers Used:** mcp-spatialtools (8 tools)
**Expected Runtime:** 30-45 minutes (full analysis) / 5-8 minutes (DRY_RUN mode)

### Test Prompt

```
[Paste complete content of TEST_3_SPATIAL_TRANSCRIPTOMICS.txt]
```

### Expected Outputs

1. **Spatial Clusters**
   - Tumor core (epithelial cells)
   - Immune infiltrate (T cells, macrophages)
   - Cancer-associated fibroblasts (CAFs)
   - Stromal regions

2. **Spatial Features**
   - Immune exclusion pattern (low T cell infiltration)
   - CAF abundance correlates with resistance
   - Spatial heterogeneity in PI3K pathway activation

3. **Deconvolution**
   - Cell type proportions by spatial region
   - Immune landscape characterization

---

## TEST 4: Histology Imaging Analysis

**File:** `TEST_4_HISTOLOGY_IMAGING.txt`
**Servers Used:** mcp-openimagedata (3 tools), mcp-deepcell (2 tools)
**Expected Runtime:** 20-30 minutes (full analysis) / 3-5 minutes (DRY_RUN mode)

### Test Prompt

```
[Paste complete content of TEST_4_HISTOLOGY_IMAGING.txt]
```

### Expected Outputs

1. **Nuclear Segmentation**
   - ~50,000 nuclei segmented
   - High-grade morphology confirmed
   - Pleomorphic nuclei identified

2. **Tissue Architecture**
   - Solid tumor architecture
   - Necrotic regions (10-15% of area)
   - Minimal desmoplastic reaction

3. **Quantitative Features**
   - Nuclear area, perimeter, circularity
   - Chromatin texture features
   - Spatial distribution metrics

---

## TEST 5: TCGA Comparison & Survival Analysis

**File:** `TEST_5_TCGA_SURVIVAL.txt`
**Servers Used:** mcp-tcga (5 tools)
**Expected Runtime:** 15-20 minutes (full analysis) / 3-5 minutes (DRY_RUN mode)

### Test Prompt

```
[Paste complete content of TEST_5_TCGA_SURVIVAL.txt]
```

### Expected Outputs

1. **TCGA Subtype Match**
   - C2 (Immunoreactive) subtype most similar
   - Median survival: 18 months
   - PIK3CA amplification common in this subtype

2. **Survival Curves**
   - Kaplan-Meier analysis
   - Patients with PI3K pathway activation: median PFS 6 months
   - Patients without: median PFS 12 months

3. **Expression Comparison**
   - Patient gene expression vs TCGA-OV cohort
   - Genes higher in patient: PIK3CA, AKT1, ABCB1
   - Genes lower in patient: PTEN, TP53

---

## Complete Workflow Integration

### All 5 Tests Combined

**Order of Execution:**
1. TEST 1 (Clinical-Genomic) → Genomic alterations
2. TEST 2 (Multi-Omics) → Resistance mechanisms ⭐ ENHANCED
3. TEST 3 (Spatial) → Tumor microenvironment
4. TEST 4 (Imaging) → Tissue architecture
5. TEST 5 (TCGA) → Prognostic context

**Total Runtime:**
- Full analysis mode: ~2-3 hours
- DRY_RUN mode: ~20-30 minutes

### Expected Final Report

**Genomic Findings (TEST 1):**
- PIK3CA amplification (CN=4)
- TP53 mutation (p.R273H)
- PTEN deletion (homozygous)

**Multi-Omics Findings (TEST 2):** ⭐ ENHANCED
- **Preprocessing:** Batch correction successful (0.82 → 0.12)
- **Pathway:** PI3K/AKT/mTOR activated (all Z > 2.5)
- **Therapeutic Targets:** AKT1, MTOR, PI3K
- **Drug Recommendations:** Alpelisib + Capivasertib
- **Clinical Trial:** NCT03602859

**Spatial Findings (TEST 3):**
- Immune exclusion phenotype
- CAF-rich microenvironment
- Spatial heterogeneity in resistance markers

**Imaging Findings (TEST 4):**
- High-grade serous architecture
- Pleomorphic nuclei
- Necrotic regions present

**TCGA Findings (TEST 5):**
- C2 (Immunoreactive) subtype
- Median survival: 18 months
- PI3K pathway activation common

### Integrated Clinical Interpretation

**Diagnosis:**
- Stage IV HGSOC, platinum-resistant
- High-risk genomic profile (PTEN loss, TP53 mutation, PIK3CA amplification)
- PI3K/AKT/mTOR pathway hyperactivation confirmed across genomic, transcriptomic, proteomic, and spatial data

**Treatment Recommendation:**
- **Primary:** Dual PI3K/AKT inhibition (Alpelisib + Capivasertib)
- **Evidence:** Multi-omics data shows pathway activation, drug targets identified
- **Clinical Trial:** NCT03602859 (patient eligible)
- **Monitoring:** Phospho-AKT/S6 levels, CA-125, imaging every 8 weeks

**Prognosis:**
- Platinum-resistant disease (inherently challenging)
- TCGA comparison: 18-month median survival
- With targeted therapy: Potential for extended PFS (6-9 months)

---

## Validation Checkpoints

### TEST 2 Validation (Multi-Omics) ⭐ UPDATED

**Preprocessing Validation:**
- [ ] Batch effects detected (PC1-batch r > 0.7)
- [ ] Batch correction effective (PC1-batch r < 0.3)
- [ ] Imputation completed (~2000 protein values)
- [ ] QC plots generated (4 PNG files)
- [ ] Final sample count: 13 (7R + 6S)

**Analysis Validation:**
- [ ] All 7 resistance genes analyzed
- [ ] Stouffer's Z-scores > 3 for top genes
- [ ] FDR correction applied AFTER combination
- [ ] All q-values < 0.05

**Upstream Regulator Validation:**
- [ ] Kinases identified: AKT1, MTOR, PI3K
- [ ] TF identified: TP53 (inhibited)
- [ ] Drug targets: Alpelisib, Capivasertib, Everolimus
- [ ] Clinical trial matched: NCT03602859

### All Tests Validation

- [ ] TEST 1: Variants called successfully
- [ ] TEST 2: Multi-omics analysis complete (with preprocessing)
- [ ] TEST 3: Spatial clusters identified
- [ ] TEST 4: Nuclei segmented successfully
- [ ] TEST 5: TCGA comparison complete

---

## Troubleshooting

### Common Issues

**Issue 1: TEST_2 batch correction doesn't work**
- **Symptom:** PC1-batch correlation still high after preprocessing
- **Cause:** Batches confounded with phenotype
- **Solution:** Check batch-phenotype contingency table
- **Expected:** Both batches should have mix of resistant and sensitive samples

**Issue 2: No significant genes after Stouffer's**
- **Symptom:** All q-values > 0.05
- **Cause:** Insufficient statistical power or incorrect FDR workflow
- **Solution:** Verify FDR applied AFTER combination (not before)
- **Check:** Raw p-values should be combined first

**Issue 3: Upstream regulators not identified**
- **Symptom:** No kinases or TFs with q < 0.05
- **Cause:** Insufficient differentially expressed target genes
- **Solution:** Lower initial gene selection threshold (p < 0.05 instead of q < 0.05)

**Issue 4: DRY_RUN mode not working**
- **Symptom:** Tools trying to access real files
- **Cause:** Environment variable not set
- **Solution:** Set `MULTIOMICS_DRY_RUN=true` before running

---

## File Locations

**Test Prompts:**
```
/manual_testing/PatientOne-OvarianCancer/implementation/
├── TEST_1_CLINICAL_GENOMIC.txt
├── TEST_2_MULTIOMICS_ENHANCED.txt  ⭐ UPDATED (was TEST_2_MULTIOMICS.txt)
├── TEST_3_SPATIAL_TRANSCRIPTOMICS.txt
├── TEST_4_HISTOLOGY_IMAGING.txt
└── TEST_5_TCGA_SURVIVAL.txt
```

**Expected Outputs:**
```
/architecture/patient-one/patient-one-outputs/
├── for-developer/
│   ├── MCP_Servers_Reference_Guide.pdf (updated)
│   ├── MCP_Report_PAT001.pdf (updated)
│   └── Full_Test_Prompt.pdf (this document)
├── for-care-team/
│   ├── multiomics_resistance_analysis.png (updated)
│   ├── spatial_transcriptomics_analysis.png
│   ├── histology_imaging_analysis.png
│   └── MCP_Report_PAT001.pdf (updated)
└── for-patient/
    ├── patient_summary.html (updated)
    ├── medication_guide.html (updated)
    └── patient_infographic.png
```

---

## Appendix: Tool Counts by Server

| Server | Tools (v1.0) | Tools (v2.0) | Change |
|--------|--------------|--------------|--------|
| mcp-fgbio | 4 | 4 | - |
| mcp-spatialtools | 8 | 8 | - |
| mcp-openimagedata | 3 | 3 | - |
| mcp-seqera | 3 | 3 | - |
| mcp-huggingface | 3 | 3 | - |
| mcp-deepcell | 2 | 2 | - |
| mcp-mockepic | 3 | 3 | - |
| mcp-tcga | 5 | 5 | - |
| **mcp-multiomics** | **5** | **9** | **+4** ⭐ |
| **TOTAL** | **36** | **40** | **+4** |

**New Tools in v2.0:**
1. validate_multiomics_data (preprocessing)
2. preprocess_multiomics_data (preprocessing)
3. visualize_data_quality (preprocessing)
4. predict_upstream_regulators (therapeutic targets)

---

**Document Status:** ✅ Updated for Version 2.0
**Last Validated:** December 26, 2025
**Ready for Use:** Yes - all test prompts and expected outputs updated
