# End-to-End Test Prompt: Patient One - Ovarian Cancer

**Purpose:** Comprehensive test of all 9 MCP servers using realistic synthetic patient data

**Patient ID:** PAT001-OVC-2025

**Data Location:** `/Users/lynnlangit/Documents/GitHub/spatial-mcp/manual_testing/PatientOne-OvarianCancer/Synthetic_sample_data/`

---

## ⚠️ CRITICAL: Use Existing Data Files

**All 17 data files ALREADY EXIST** in the directory above. When testing:
- ✅ **DO:** Read and analyze the existing files
- ❌ **DON'T:** Generate new synthetic data
- ❌ **DON'T:** Create new patient data
- ✅ **DO:** Use the pre-generated files in clinical/, genomics/, multiomics/, spatial/, and imaging/ directories

The prompt has been updated to explicitly instruct Claude Desktop to use existing files.

---

## Complete Test Prompt

Copy and paste this entire prompt into Claude Desktop:

```
I need your help with a comprehensive precision oncology analysis for a 58-year-old female patient with Stage IV high-grade serous ovarian carcinoma who has developed platinum resistance.

⚠️ IMPORTANT: All patient data files ALREADY EXIST in this directory. Please USE THE EXISTING FILES - do NOT generate new synthetic data. Read the actual files from:

/Users/lynnlangit/Documents/GitHub/spatial-mcp/manual_testing/PatientOne-OvarianCancer/Synthetic_sample_data/

## Patient Background

Patient ID: PAT001-OVC-2025 (this patient data is already created and saved in the above directory)

Key Clinical Features:
- Diagnosis: Stage IV HGSOC (High-Grade Serous Ovarian Cancer)
- BRCA1 germline mutation (pathogenic)
- TP53 R175H somatic mutation (hotspot)
- HRD positive (score: 42)
- Treatment status: Platinum-resistant progression

## Analysis Request - Please Complete All Steps:

### PART 1: Clinical Data Review (mcp-mockepic)

1. **Retrieve patient demographics and clinical history:**
   - Load patient_demographics.json from the clinical/ directory
   - Summarize key demographics, family history, and genetic risk factors

2. **Analyze lab results and disease progression:**
   - Load lab_results.json from the clinical/ directory
   - Create a CA-125 trend visualization showing:
     * Initial diagnosis (should show ~1456 U/mL)
     * Response to first-line therapy (normalized to ~22 U/mL)
     * Progression/resistance (rising to ~389 U/mL)
     * Current status on second-line therapy (~289 U/mL)
   - Interpret the trend in context of platinum resistance

### PART 2: Genomic Analysis (mcp-fgbio, mcp-tcga)

3. **Parse somatic variants:**
   - Load somatic_variants.vcf from the genomics/ directory
   - Identify key mutations: TP53, PIK3CA, PTEN
   - Note copy number alterations (MYC, CCNE1, AKT2 amplifications)

4. **Compare to TCGA-OV cohort (mcp-tcga):**
   - What TCGA molecular subtype does this patient match? (C1 immunoreactive vs C2 differentiated vs C4 mesenchymal)
   - What is the expected survival for BRCA1-mutant, TP53-mutant, Stage IV HGSOC?
   - Are there other platinum-resistant cases with similar molecular profiles in TCGA?
   - What pathways are commonly activated in similar cases?

### PART 3: Multi-Omics PDX Resistance Analysis (mcp-multiomics)

5. **Integrate multi-omics data from PDX models:**
   - Load data from multiomics/ directory:
     * sample_metadata.csv (15 samples: 7 resistant, 8 sensitive)
     * pdx_rna_seq.csv (1000 genes × 15 samples)
     * pdx_proteomics.csv (500 proteins × 15 samples)
     * pdx_phosphoproteomics.csv (300 sites × 15 samples)

6. **Run Stouffer's meta-analysis:**
   - Focus on resistance pathway genes:
     * PI3K/AKT pathway: PIK3CA, AKT1, AKT2, MTOR, PTEN
     * Drug resistance: ABCB1 (MDR1), ABCC1 (MRP1)
     * Anti-apoptotic: BCL2L1, MCL1
     * DNA repair: BRCA1, BRCA2, RAD51
   - Use directionality from effect sizes
   - Apply FDR correction (α = 0.05)
   - Which genes are significantly dysregulated across ALL three modalities?
   - What pathways are activated in resistant vs. sensitive PDX models?

### PART 4: Spatial Tumor Analysis (mcp-spatialtools, mcp-openimagedata, mcp-deepcell)

7. **Analyze spatial transcriptomics data:**
   - Load spatial data from spatial/ directory:
     * visium_spatial_coordinates.csv (900 spots)
     * visium_gene_expression.csv (31 genes × 900 spots)
     * visium_region_annotations.csv (6 spatial regions)
   - Identify spatial regions: tumor core, tumor-stroma interface, immune-infiltrated areas
   - Which resistance markers (PIK3CA, AKT1, ABCB1) show spatial heterogeneity?
   - Where are immune cells (CD3, CD8, CD68) located relative to tumor?

8. **Process histology images:**
   - Load H&E image: imaging/PAT001_tumor_HE_20x.tiff
   - Analyze tissue architecture and identify necrotic regions
   - Load immunofluorescence images: IF_DAPI.tiff, IF_CD8.tiff, IF_KI67.tiff
   - Quantify:
     * Tumor cellularity
     * CD8+ T cell infiltration
     * Ki67 proliferation index
     * Spatial distribution patterns

9. **Cell segmentation and phenotyping (if possible with mcp-deepcell):**
   - Segment cells from multiplex IF image: PAT001_tumor_multiplex_IF_TP53_KI67_DAPI.tiff
   - Identify cell types and phenotypes
   - Calculate cell density in different tumor regions

### PART 5: Integrated Analysis & Recommendations

10. **Synthesize findings across all modalities:**
    - What are the PRIMARY resistance mechanisms based on:
      * Genomics (TP53, PIK3CA mutations)
      * Multi-omics (activated pathways across RNA/Protein/Phospho)
      * Spatial patterns (heterogeneity, immune exclusion)
      * Imaging (proliferation, immune infiltration)

11. **Therapeutic recommendations:**
    Based on the integrated analysis, recommend:
    - Targeted therapies (PI3K/AKT/mTOR inhibitors?)
    - Combination strategies
    - Clinical trials for BRCA-mutant platinum-resistant HGSOC
    - Immunotherapy considerations based on immune landscape

12. **Biomarker identification:**
    - Which molecular features could be monitored for treatment response?
    - Are there spatial biomarkers predictive of therapy resistance?
    - What should be tracked in serial biopsies?

## Output Format

Please provide:

1. **Executive Summary** (1-2 paragraphs)
   - Patient status, key findings, primary resistance mechanisms

2. **Detailed Analysis by Modality**
   - Clinical trajectory and CA-125 trends
   - Genomic alterations and TCGA comparison
   - Multi-omics resistance signatures
   - Spatial tumor heterogeneity
   - Imaging and cellular phenotypes

3. **Integrated Pathway Analysis**
   - Which pathways are consistently activated across all data types?
   - Confidence level for each finding (genomics + proteomics + spatial = high confidence)

4. **Clinical Recommendations**
   - Ranked list of treatment options with rationale
   - Clinical trial opportunities
   - Monitoring strategy

5. **Technical Summary**
   - Which MCP servers were used for each analysis step
   - Data quality metrics
   - Any limitations or caveats

## Expected MCP Server Usage

This analysis should invoke:
- ✅ mcp-mockepic (clinical data retrieval)
- ✅ mcp-fgbio (VCF parsing, gene annotations)
- ✅ mcp-tcga (TCGA-OV cohort comparison)
- ✅ mcp-multiomics (Stouffer's meta-analysis, multi-omics integration)
- ✅ mcp-spatialtools (spatial transcriptomics analysis)
- ✅ mcp-openimagedata (histology image processing)
- ✅ mcp-deepcell (cell segmentation, if available)
- ✅ mcp-huggingface (optional: ML model predictions)
- ✅ mcp-seqera (optional: workflow orchestration)

Please proceed with the analysis and let me know if you need any clarification on the data files or analysis requirements.
```

---

## Expected Outcomes

### Clinical Data Analysis
- Patient demographics parsed from FHIR-inspired JSON
- CA-125 trend: 1456 → 22 → 389 → 289 U/mL
- Clear platinum resistance pattern identified

### Genomic Analysis
- TP53 R175H hotspot mutation identified
- PIK3CA E545K activating mutation noted
- PTEN loss of heterozygosity detected
- Copy number alterations: MYC, CCNE1, AKT2 amplifications

### TCGA Comparison
- Likely C1 (immunoreactive) or C2 (differentiated) subtype
- BRCA-mutant cohort comparison
- Poor prognosis with platinum resistance
- PI3K/AKT pathway commonly activated

### Multi-Omics Analysis
- Top resistance genes: AKT1, PIK3CA, ABCB1, BCL2L1
- Stouffer's Z-scores highly significant (Z > 3, p < 0.001)
- Directionality preserved (upregulated in resistant)
- FDR < 0.05 for pathway genes

### Spatial Analysis
- 6 distinct regions identified
- Tumor core: high proliferation (MKI67, TOP2A)
- Tumor-stroma interface: EMT markers (VIM, SNAI1)
- Immune regions: CD3+, CD8+ T cells
- Spatial heterogeneity in resistance markers

### Imaging Analysis
- H&E: necrotic regions, tumor architecture
- IF: CD8+ infiltration (quantified)
- Ki67: proliferation index ~60-70%
- Multiplex: cell phenotypes and spatial relationships

### Treatment Recommendations
- PI3K inhibitors (alpelisib + fulvestrant)
- AKT inhibitors (capivasertib + paclitaxel)
- mTOR inhibitors (everolimus)
- Clinical trials for BRCA-mutant resistant HGSOC
- Monitor CA-125, consider serial biopsies

---

## Testing Tips

1. **Run in stages** if needed:
   - First: Clinical + Genomic (Parts 1-2)
   - Then: Multi-omics (Part 3)
   - Finally: Spatial + Imaging (Part 4)
   - Synthesis: Integration (Part 5)

2. **Monitor server usage:**
   - Check Claude Desktop logs to see which servers are invoked
   - Verify all 9 servers are accessible

3. **Validate results:**
   - Compare findings to TESTING_GUIDE.md expected results
   - Check if resistance mechanisms match synthetic data generation

4. **Performance:**
   - Full analysis may take 5-10 minutes
   - Multi-omics Stouffer's analysis is most computationally intensive

---

**Created:** November 12, 2025
**Test Status:** ✅ Ready for execution
**Expected Duration:** 5-10 minutes for complete analysis
