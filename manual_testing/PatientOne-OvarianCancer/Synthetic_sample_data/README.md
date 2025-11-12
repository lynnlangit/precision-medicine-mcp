# Patient One - Stage IV Ovarian Cancer Synthetic Dataset

**Purpose:** Comprehensive synthetic dataset for end-to-end testing of all 9 MCP servers

**Patient ID:** PAT001-OVC-2025
**Diagnosis:** Stage IV High-Grade Serous Ovarian Carcinoma (HGSOC)
**Status:** Treatment-resistant, PDX model generated

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

### Clinical Data (`clinical/`)
- `patient_demographics.json` - Demographics, insurance (mcp-mockepic)
- `diagnosis_history.json` - Staging, pathology reports
- `treatment_history.json` - Surgeries, chemotherapy
- `lab_results.json` - CA-125 trends (35→389 U/mL progression)
- `medications.json` - Current/past medications
- `imaging_reports.json` - CT/MRI findings

### Genomic Data (`genomics/`)
- `tumor_WES_R1/R2.fastq.gz` - Whole exome sequencing (headers only)
- `somatic_variants.vcf` - TP53, PIK3CA mutations
- `germline_variants.vcf` - BRCA1 germline
- `copy_number_alterations.seg` - CNVs (MYC amp, PTEN del)

### Spatial Transcriptomics (`spatial/`)
- `visium_filtered_feature_bc_matrix.h5` - 10x Visium gene counts
- `visium_tissue_positions.csv` - Spatial coordinates
- `tissue_image_hires.png` - H&E image
- `cell_phenotypes.csv` - Cell type annotations

### Multi-Omics PDX Data (`multiomics/`)
**15 samples:** 7 resistant + 8 sensitive to carboplatin
- `pdx_rna_seq.csv` - 15,000 genes × 15 samples
- `pdx_proteomics.csv` - 8,500 proteins × 15 samples
- `pdx_phosphoproteomics.csv` - 4,200 phosphosites × 15 samples
- `sample_metadata.csv` - Treatment response annotations
- `differential_analysis_results.csv` - Pre-computed DE

### Histology Imaging (`imaging/`)
- `primary_tumor_HE.tiff` - Ovarian tumor H&E
- `omentum_metastasis_HE.tiff` - Metastatic site
- `PAX8_IHC.tiff`, `WT1_IHC.tiff` - Positive markers
- `p53_IHC.tiff` - Aberrant expression (mutant)

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

**Total Size:** ~155 MB
**Format:** JSON (clinical), VCF (variants), H5/CSV (omics), TIFF (images)
**Status:** ✅ Ready for end-to-end testing

**Created:** November 11, 2025
