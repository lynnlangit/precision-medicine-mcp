# Patient One - Stage IV Ovarian Cancer Synthetic Dataset

**Purpose:** Comprehensive synthetic dataset for end-to-end testing of all 9 MCP servers

**Patient ID:** PAT001-OVC-2025
**Diagnosis:** Stage IV High-Grade Serous Ovarian Carcinoma (HGSOC)
**Status:** Treatment-resistant, PDX model generated

---

## 📁 Data Location

**IMPORTANT:** The actual patient data files are located in:
```
/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/
```

This location is accessible to Claude Desktop's MCP servers. The data structure is:
```
data/patient-data/PAT001-OVC-2025/
├── clinical/          (2 files: patient_demographics.json, lab_results.json)
├── genomics/          (1 file: somatic_variants.vcf)
├── multiomics/        (4 files: metadata, RNA, protein, phospho data)
├── spatial/           (3 files: coordinates, expression, annotations)
└── imaging/           (7 files: H&E, IF, multiplex images)
```

**This directory contains:** Test prompts and documentation only (no data files).

**See:** `CLAUDE_DESKTOP_FILE_ACCESS_GUIDE.md` for detailed explanation of why data is in `/data/` location.

---

## Clinical Summary

### Demographics
- Age: 58 years, Female, Caucasian
- Family History: Mother with breast cancer (age 52)

### Diagnosis (January 2024)
- **Stage:** IV (FIGO) - Metastases to omentum, peritoneum, liver, pleura
- **Histology:** High-grade serous carcinoma
- **Molecular:** TP53 R175H mutation, BRCA1 germline mutation, HRD+ (score 42)

### Treatment Course
1. Surgery: Suboptimal debulking (February 2024)
2. First-line: Carboplatin/Paclitaxel + Olaparib maintenance
3. Progression: October 2024 (platinum-resistant, <6 months)
4. Second-line: Pegylated liposomal doxorubicin (ongoing)
5. Research: PDX model from pleural effusion

---

## Dataset Contents

**Note:** All data files are located in `../../../data/patient-data/PAT001-OVC-2025/`

### Clinical Data (2 files)
- `patient_demographics.json` - Demographics, insurance, BRCA1 status
- `lab_results.json` - CA-125 trends showing platinum resistance (1456→22→389→289 U/mL)

### Genomic Data (1 file)
- `somatic_variants.vcf` - TP53 R175H, PIK3CA E545K, PTEN LOH mutations

### Multi-Omics PDX Data (4 files)
**15 samples:** 7 resistant + 8 sensitive to carboplatin
- `sample_metadata.csv` - Treatment response annotations
- `pdx_rna_seq.csv` - 1,000 genes × 15 samples
- `pdx_proteomics.csv` - 500 proteins × 15 samples
- `pdx_phosphoproteomics.csv` - 300 phosphosites × 15 samples

### Spatial Transcriptomics (3 files)
- `visium_spatial_coordinates.csv` - 900 spots across 6 regions
- `visium_gene_expression.csv` - 31 genes × 900 spots
- `visium_region_annotations.csv` - Region labels (tumor_core, tumor_proliferative, etc.)

### Histology Imaging (7 files)
- `PAT001_tumor_HE_20x.tiff` - H&E stained tissue (20x magnification)
- `PAT001_tumor_IF_DAPI.tiff` - Nuclear stain
- `PAT001_tumor_IF_CD3.tiff` - T cell marker
- `PAT001_tumor_IF_CD8.tiff` - Cytotoxic T cell marker
- `PAT001_tumor_IF_KI67.tiff` - Proliferation marker
- `PAT001_tumor_IF_PanCK.tiff` - Epithelial marker
- `PAT001_tumor_multiplex_IF_TP53_KI67_DAPI.tiff` - 3-channel multiplex image

**Total:** 17 files (~3.2 MB)

---

## Key Molecular Features

**Mutations:**
- TP53: R175H (hotspot, chr17:7,578,406)
- BRCA1: Germline c.5266dupC
- PIK3CA: E545K (resistance mechanism)

**Copy Number:**
- Amplified: MYC (8q24), CCNE1, AKT2
- Deleted: PTEN, RB1, CDKN2A

**Resistance Signatures (Multi-Omics):**
- Upregulated: AKT1, PIK3CA, ABCB1 (MDR1), BCL2L1
- Downregulated: BAX, PTEN, FOXO3
- Pathway: PI3K/AKT/mTOR activation

---

## Testing Workflows

### Workflow 1: EHR → TCGA Comparison
```
1. Get patient data (mcp-mockepic)
2. Compare to TCGA-OV cohort (mcp-tcga)
3. Identify similar cases and outcomes
```

### Workflow 2: Genomic Analysis
```
1. Validate FASTQ (mcp-fgbio)
2. Query TP53, BRCA1 mutations
3. Analyze CNVs
```

### Workflow 3: Spatial Analysis
```
1. Load Visium data (mcp-spatialtools)
2. Segment cells (mcp-deepcell)
3. Identify spatial domains
4. Calculate tumor/stroma/immune regions
```

### Workflow 4: Multi-Omics Integration
```
1. Integrate RNA/Protein/Phospho (mcp-multiomics)
2. Run Stouffer's meta-analysis
3. Identify resistance pathways
4. Generate heatmaps
```

### Workflow 5: Complete Precision Medicine
```
1-5. All above workflows
6. Generate treatment recommendations
7. Launch validation (mcp-seqera)
```

---

## Expected Results

**TCGA Comparison:** Clusters with C1 (immunoreactive), BRCA-mutant cohort  
**Multi-Omics:** AKT1, PIK3CA, ABCB1 highly significant (FDR<0.001)  
**Spatial:** Tumor core vs invasive front, immune exclusion phenotype  
**Treatment Predictions:** High sensitivity to PI3K/AKT inhibitors  

---

## Quick Start

**To test with Claude Desktop:**
1. Verify data is in MCP-accessible location:
   ```bash
   ls -la ../../../data/patient-data/PAT001-OVC-2025/
   ```
2. Copy a test prompt:
   ```bash
   cat TEST_1_CLINICAL_GENOMIC.txt
   ```
3. Paste into Claude Desktop and run

**See:** `QUICK_TEST_REFERENCE.md` for detailed testing instructions

---

**Data Location:** `../../../data/patient-data/PAT001-OVC-2025/`
**Total Size:** ~3.2 MB (17 files)
**Format:** JSON (clinical), VCF (variants), CSV (omics), TIFF (images)
**Status:** ✅ Ready for Claude Desktop testing

**Created:** November 11, 2025
**Updated:** November 13, 2025 (removed redundant data copies, clarified location)
