# Patient One - Testing Guide

**Quick Start:** Use this synthetic dataset to test all 9 MCP servers

⚠️ **IMPORTANT:** The comprehensive end-to-end test may hit Claude Desktop's context limit. Use the **5 smaller focused tests** instead (see below).

---

## Testing Approaches

### Option 1: Focused Tests (Recommended) ✅
**Use these 5 smaller tests to avoid context limits:**
- `TEST_1_CLINICAL_GENOMIC.txt` - Clinical data + genomics (3 servers)
- `TEST_2_MULTIOMICS.txt` - Multi-omics analysis (1 server)
- `TEST_3_SPATIAL.txt` - Spatial transcriptomics (1 server)
- `TEST_4_IMAGING.txt` - Histology and imaging (2 servers)
- `TEST_5_INTEGRATION.txt` - Synthesis and recommendations

**See:** `TESTING_STRATEGY.md` for detailed instructions
**Quick Ref:** `QUICK_TEST_REFERENCE.md` for fast lookup

### Option 2: Comprehensive Test ⚠️
**May hit context limit - use with caution:**
- `COPY_PASTE_PROMPT.txt` - All-in-one test (all 9 servers)
- `END_TO_END_TEST_PROMPT.md` - Full documentation

---

## Test Status

### All Files Created ✅

**Clinical Data:**
- `README.md` - Complete patient case description
- `clinical/patient_demographics.json` - Full demographics, family history
- `clinical/lab_results.json` - CA-125 trends showing progression

**Multi-Omics Data:**
- `multiomics/sample_metadata.csv` - 15 PDX samples (7 resistant, 8 sensitive)
- `multiomics/pdx_rna_seq.csv` - RNA-seq expression data (1000 genes × 15 samples)
- `multiomics/pdx_proteomics.csv` - Proteomics data (500 proteins × 15 samples)
- `multiomics/pdx_phosphoproteomics.csv` - Phosphoproteomics data (300 sites × 15 samples)

**Genomics Data:**
- `genomics/somatic_variants.vcf` - VCF with TP53, PIK3CA, PTEN mutations

**Spatial Transcriptomics Data:**
- `spatial/visium_spatial_coordinates.csv` - 900 spot coordinates
- `spatial/visium_gene_expression.csv` - 31 genes × 900 spots
- `spatial/visium_region_annotations.csv` - Region annotations (tumor core/interface/stroma/immune/necrotic)

**Imaging Data:**
- `imaging/PAT001_tumor_HE_20x.tiff` - H&E histology (512×512)
- `imaging/PAT001_tumor_IF_DAPI.tiff` - DAPI nuclear stain
- `imaging/PAT001_tumor_IF_CD3.tiff` - CD3 T cell marker
- `imaging/PAT001_tumor_IF_CD8.tiff` - CD8 cytotoxic T cell marker
- `imaging/PAT001_tumor_IF_KI67.tiff` - Ki67 proliferation marker
- `imaging/PAT001_tumor_IF_PanCK.tiff` - Pan-cytokeratin epithelial marker
- `imaging/PAT001_tumor_multiplex_IF_TP53_KI67_DAPI.tiff` - 3-channel multiplex IF

---

## End-to-End Testing Prompts

### Test 1: EHR Data Retrieval (mcp-mockepic)
```
For patient PAT001-OVC-2025, please:
1. Retrieve patient demographics
2. Show CA-125 lab results trend
3. Summarize treatment history

Data location: /manual_testing/PatientOne-OvarianCancer/implementation/clinical/
```

**Expected:** Demographics, CA-125 trend from 1456→22→389→289, platinum-resistant progression

### Test 2: Multi-Omics Resistance Analysis (mcp-multiomics)
```
Analyze platinum resistance in PDX samples:

Data: /manual_testing/PatientOne-OvarianCancer/implementation/multiomics/
Metadata: sample_metadata.csv (7 resistant, 8 sensitive)

Key resistance genes to analyze:
- PI3K/AKT pathway: AKT1, PIK3CA, mTOR, PTEN
- Drug resistance: ABCB1 (MDR1), ABCC1
- Apoptosis: BCL2L1, BAX, TP53
- DNA repair: BRCA1, BRCA2

Tasks:
1. Integrate RNA/Protein/Phospho data
2. Run Stouffer's meta-analysis with directionality
3. Apply FDR correction (α=0.05)
4. Identify genes significant across all modalities
5. Which pathways are activated in resistant tumors?
```

**Expected:** AKT1, PIK3CA, ABCB1, BCL2L1 highly significant; PI3K/AKT pathway activated

### Test 3: TCGA Comparison (mcp-tcga)
```
Compare patient PAT001-OVC-2025 to TCGA-OV cohort:

Molecular profile:
- TP53: R175H mutation (hotspot)
- BRCA1: Germline pathogenic variant
- HRD score: 42 (positive)
- Stage: IV
- Histology: High-grade serous

Questions:
1. What TCGA molecular subtype does this patient match?
2. What is typical survival for similar patients?
3. Are there other BRCA1-mutant cases with platinum resistance?
4. What treatments showed benefit in similar cases?
```

**Expected:** C1 (immunoreactive) subtype, BRCA-mutant cohort, poor prognosis with platinum resistance

### Test 4: Complete Precision Medicine Workflow
```
Comprehensive analysis for patient PAT001-OVC-2025:

1. Retrieve EHR data (demographics, labs, treatments)
2. Analyze multi-omics PDX resistance signatures
3. Compare molecular profile to TCGA-OV
4. Identify activated pathways
5. Recommend targeted therapies

Data directory: /manual_testing/PatientOne-OvarianCancer/implementation/

Expected outcome:
- Identify PI3K/AKT/mTOR pathway as key resistance mechanism
- Recommend PI3K inhibitors (alpelisib), AKT inhibitors (capivasertib)
- Note platinum resistance and immunosuppressive tumor microenvironment
```

---

## Key Testing Expectations

### Clinical Data (mcp-mockepic)
- ✅ Patient demographics parsed correctly
- ✅ CA-125 trend visualized: diagnosis→normalized→progression
- ✅ Family history shows BRCA1 mutation context
- ✅ Treatment timeline clear

### Multi-Omics Integration (mcp-multiomics)
- ✅ 15 samples loaded (7R + 8S)
- ✅ Sample alignment across modalities
- ✅ Stouffer's method with directionality
- ✅ FDR correction applied
- ✅ Top genes: AKT1, PIK3CA, ABCB1, BCL2L1
- ✅ Pathway: PI3K/AKT/mTOR activation

### TCGA Comparison (mcp-tcga)
- ✅ Subtype: C1 (immunoreactive) or C2 (differentiated)
- ✅ Cohort: BRCA-mutant HGSOC
- ✅ Survival: Poor with Stage IV + platinum resistance
- ✅ Similar cases identified

### Treatment Recommendations
- ✅ PI3K inhibitors (alpelisib + fulvestrant)
- ✅ AKT inhibitors (capivasertib + paclitaxel)
- ✅ mTOR inhibitors (everolimus)
- ✅ Clinical trials for BRCA-mutant resistant HGSOC

---

## Molecular Profile Summary

**Key Mutations:**
- TP53: R175H (chr17:7,578,406 C>A) - hotspot, loss of function
- BRCA1: Germline c.5266dupC - frameshift, pathogenic
- PIK3CA: E545K (chr3:178,936,091 G>A) - activating, resistance

**Copy Number:**
- Amplified: MYC (8q24), CCNE1 (19q12), AKT2 (19q13), KRAS (12p12)
- Deleted: PTEN (10q23), RB1 (13q14), CDKN2A (9p21), TP53 (17p13)

**HRD Score:** 42 (positive, >33 threshold)

**Resistance Mechanisms:**
1. PI3K/AKT/mTOR pathway activation (PIK3CA mutation, AKT2 amplification, PTEN loss)
2. Drug efflux pumps (ABCB1/MDR1, ABCC1/MRP1 upregulation)
3. Anti-apoptotic signaling (BCL2L1/Bcl-xL upregulation)
4. DNA repair alterations (acquired BRCA reversion, RAD51 upregulation)

---

## File Locations

**Note:** Patient data files have been moved to MCP-accessible location:
```
/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/
├── clinical/
│   ├── patient_demographics.json                ✅ Demographics, family history
│   └── lab_results.json                         ✅ CA-125 trends, CBC, metabolic panel
├── genomics/
│   └── somatic_variants.vcf                     ✅ TP53, PIK3CA, PTEN mutations
├── multiomics/
│   ├── sample_metadata.csv                      ✅ 15 PDX samples
│   ├── pdx_rna_seq.csv                          ✅ 1000 genes × 15 samples
│   ├── pdx_proteomics.csv                       ✅ 500 proteins × 15 samples
│   └── pdx_phosphoproteomics.csv                ✅ 300 sites × 15 samples
├── spatial/
│   ├── visium_spatial_coordinates.csv           ✅ 900 spot coordinates
│   ├── visium_gene_expression.csv               ✅ 31 genes × 900 spots
│   └── visium_region_annotations.csv            ✅ Region annotations
└── imaging/
    ├── PAT001_tumor_HE_20x.tiff                 ✅ H&E histology
    ├── PAT001_tumor_IF_DAPI.tiff                ✅ DAPI nuclear stain
    ├── PAT001_tumor_IF_CD3.tiff                 ✅ T cell marker
    ├── PAT001_tumor_IF_CD8.tiff                 ✅ Cytotoxic T cell marker
    ├── PAT001_tumor_IF_KI67.tiff                ✅ Proliferation marker
    ├── PAT001_tumor_IF_PanCK.tiff               ✅ Epithelial marker
    └── PAT001_tumor_multiplex_IF_TP53_KI67_DAPI.tiff  ✅ 3-channel multiplex
```

**This directory (implementation/) contains:** Test prompts and documentation only.

---

## Ready for Testing! 🚀

All synthetic data files have been generated and are ready for end-to-end testing.

### Quick Start Testing

1. **Test EHR Data (mcp-mockepic):**
   Use Test 1 above to retrieve patient demographics and lab results

2. **Test Multi-Omics Analysis (mcp-multiomics):**
   Use Test 2 above to analyze platinum resistance signatures

3. **Test TCGA Comparison (mcp-tcga):**
   Use Test 3 above to compare patient to TCGA-OV cohort

4. **Test Complete Workflow:**
   Use Test 4 above for comprehensive precision medicine analysis

### Data Summary

- **Clinical:** 2 JSON files (demographics, labs)
- **Multi-Omics:** 4 CSV files (metadata + RNA/Protein/Phospho)
- **Genomics:** 1 VCF file (key mutations)
- **Spatial:** 3 CSV files (coordinates, expression, annotations)
- **Imaging:** 7 TIFF files (H&E, IF, multiplex)

**Total:** 17 data files covering all modalities

---

**Status:** ✅ **ALL FILES COMPLETE AND READY FOR TESTING**

**Test Coverage:** All 9 MCP servers (mockepic, multiomics, tcga, spatialtools, openimagedata, deepcell, fgbio, seqera, huggingface)

**Estimated Testing Time:** 30-45 minutes for complete end-to-end workflows

**Created:** November 12, 2025
**Last Updated:** November 12, 2025
