# Prompt Inventory - Precision Medicine MCP

**Discovery Date:** 2026-01-16
**Updated:** 2026-01-21
**Repository:** precision-medicine-mcp
**Total Prompts Found:** 55+

This document inventories all existing prompts found in the repository during Phase 1 discovery.

---

## 1. PatientOne Clinical Workflow Prompts

**Location:** `docs/test-docs/patient-one-scenario/test-prompts/*.md`

### TEST_1: Clinical Data and Genomic Analysis
- **Source File:** `test-prompts/test-1-clinical-genomic.md`
- **MCP Servers:** mcp-mockepic, mcp-fgbio, mcp-tcga
- **Audience:** Researchers, Bioinformaticians, Clinicians
- **Prompt Summary:** Comprehensive clinical and genomic analysis for Stage IV HGSOC patient
- **Key Components:**
  - Retrieve patient demographics and BRCA1 status
  - Analyze CA-125 tumor marker trends
  - Parse somatic variants (TP53, PIK3CA, PTEN)
  - Compare to TCGA ovarian cancer cohort
- **Expected Output:** Patient summary, CA-125 trend analysis, somatic mutations, TCGA subtype
- **Time Estimate:** 5-10 minutes

### TEST_2: Multi-Omics Resistance Analysis (ENHANCED)
- **Source File:** `test-prompts/test-2-multiomics-enhanced.md`
- **MCP Servers:** mcp-multiomics (9 tools)
- **Audience:** Researchers, Bioinformaticians
- **Prompt Summary:** Complete multi-omics preprocessing and analysis pipeline for platinum resistance
- **Key Components:**
  - STEP 0: Data preprocessing (validation, batch correction, QC)
  - STEP 1: Data integration (RNA/Protein/Phospho)
  - STEP 2: Stouffer's meta-analysis
  - STEP 3: Upstream regulator prediction
- **Expected Output:** Preprocessing summary, gene-level results, upstream regulators, pathway interpretation
- **Time Estimate:** 15-25 minutes

### TEST_3: Spatial Transcriptomics Analysis
- **Source File:** `test-prompts/test-3-spatial.md`
- **MCP Servers:** mcp-spatialtools
- **Audience:** Researchers, Bioinformaticians
- **Prompt Summary:** Analyze spatial gene expression patterns using 10x Visium data
- **Key Components:**
  - Load spatial data (900 spots, 6 regions)
  - Gene expression by region
  - Spatial patterns and heterogeneity
  - Generate visualizations (heatmaps, region composition, gene expression matrix)
- **Expected Output:** Spatial structure, gene expression heatmap, spatial findings, visualizations
- **Time Estimate:** 10-15 minutes

### TEST_4: Histology and Imaging Analysis
- **Source File:** `test-prompts/test-4-imaging.md`
- **MCP Servers:** mcp-openimagedata, mcp-deepcell
- **Audience:** Researchers, Pathologists
- **Prompt Summary:** Analyze H&E histology and immunofluorescence images
- **Key Components:**
  - H&E morphology analysis (cellularity, necrosis)
  - CD8 T cell infiltration quantification
  - Ki67 proliferation index
  - Multiplex IF cell phenotyping
  - Generate visualizations (segmentation overlays, heatmaps)
- **Expected Output:** H&E summary, CD8 analysis, Ki67 index, multiplex results, visualizations
- **Time Estimate:** 20-40 minutes

### TEST_5: Integrated Analysis & Clinical Recommendations
- **Source File:** `test-prompts/test-5-integration.md`
- **MCP Servers:** All (synthesis)
- **Audience:** Clinicians, Researchers
- **Prompt Summary:** Synthesize findings from all modalities and generate clinical recommendations
- **Key Components:**
  - Identify primary resistance mechanisms
  - Multi-modal consistency assessment
  - Therapeutic recommendations (targeted, immunotherapy, clinical trials)
  - Biomarkers for monitoring
  - Multi-modal visualization synthesis
- **Expected Output:** Executive summary, resistance mechanisms, treatment recommendations, monitoring strategy, integrated visualizations
- **Time Estimate:** 5-10 minutes

### TEST_6: Clinician-in-the-Loop Review
- **Source File:** `test-prompts/test-6-citl-review.md`
- **MCP Servers:** None (validation workflow)
- **Audience:** Clinicians, Oncologists
- **Prompt Summary:** Formal clinician validation workflow with quality gates and audit trail
- **Key Components:**
  - Generate draft report with automated QC
  - Structured clinician review (findings validation, guideline compliance)
  - Submit review with digital signature
  - Finalize approved report
- **Expected Output:** Draft report, quality checks, signed review, final approved report
- **Time Estimate:** 30-45 minutes (includes 20-30 min manual review)

---

## 2. MCP Server-Specific Example Prompts

### mcp-multiomics Server
**Location:** `servers/mcp-multiomics/README.md`

#### 2.1 Data Preprocessing Prompts

**Validate Multi-Omics Data Quality:**
```
Claude, please validate my multi-omics data quality:
- RNA: /data/pdx_rna_seq.csv
- Protein: /data/pdx_proteomics.csv
- Metadata: /data/sample_metadata.csv (includes Batch column)

Check for batch effects and sample naming issues.
```
- **MCP Servers:** mcp-multiomics
- **Audience:** Researchers, Bioinformaticians
- **Expected Output:** Validation status, batch effects, missing patterns, recommendations

**Preprocess Proteomics Data with Batch Correction:**
```
Claude, preprocess my proteomics data with batch correction:
- RNA: /data/rna_raw.csv
- Protein: /data/protein_raw.csv
- Metadata: /data/metadata.csv (has Batch column)

Apply quantile normalization, KNN imputation, and ComBat batch correction.
Save to /data/preprocessed/
```
- **MCP Servers:** mcp-multiomics
- **Audience:** Researchers, Bioinformaticians
- **Expected Output:** Preprocessed data paths, QC metrics, batch correction results

#### 2.2 Integration and Analysis Prompts

**Integrate Multi-Omics Data:**
```
Claude, please integrate my multi-omics PDX data:
- RNA: /data/pdx_rna_seq.csv
- Protein: /data/pdx_proteomics.csv
- Phospho: /data/pdx_phosphoproteomics.csv
- Metadata: /data/sample_metadata.csv

Apply Z-score normalization and filter features with >50% missing data.
```
- **MCP Servers:** mcp-multiomics
- **Audience:** Researchers, Bioinformaticians
- **Expected Output:** Integrated data, common samples, feature counts, QC metrics

**Combine p-values using Stouffer's Method:**
```
Combine NOMINAL p-values from HAllA using Stouffer's method:
- RNA p-values: [0.001, 0.05, 0.3] (NOMINAL from HAllA)
- Protein p-values: [0.002, 0.04, 0.25] (NOMINAL from HAllA)
- Phospho p-values: [0.01, 0.06, 0.28] (NOMINAL from HAllA)

With log2 fold changes:
- RNA log2FC: [2.5, 1.2, -0.3]
- Protein log2FC: [1.8, 1.5, -0.2]
- Phospho log2FC: [1.2, 0.8, -0.4]

Returns meta_p_values (nominal) and q_values (FDR-corrected).
USE q_values for identifying significant features!
```
- **MCP Servers:** mcp-multiomics
- **Audience:** Researchers, Bioinformaticians
- **Expected Output:** Meta p-values, q-values, significant features

**Predict Upstream Regulators:**
```
From my 50 significant genes (from Stouffer's meta-analysis),
predict upstream regulators:

Input genes with log2FC and p-values:
- AKT1, MTOR, TP53, MYC, FOXO1, etc.

Analyze:
- Kinases (activated/inhibited)
- Transcription factors
- Drug responses

Expected output:
- AKT1/MTOR activated → Alpelisib (PI3K inhibitor) recommendation
- TP53 inhibited → Loss of tumor suppression
```
- **MCP Servers:** mcp-multiomics
- **Audience:** Researchers, Clinicians
- **Expected Output:** Activated kinases, TFs, drug recommendations

#### 2.3 Complete Workflow Prompts

**Complete Clinical PDX Analysis:**
```
Claude, analyze my PDX treatment response data with complete pipeline:

STEP 0: Preprocessing (CRITICAL for real proteomics data)
1. Validate data quality
2. Preprocess with batch correction
3. Visualize QC before/after

STEP 1: Integration
4. Integrate preprocessed data

STEP 2: Association Testing
5. Run HAllA between RNA and Protein
6. Combine NOMINAL p-values with Stouffer's method

STEP 3: Upstream Regulator Analysis
7. Predict regulators from significant genes
8. Generate final report with therapeutic recommendations
```
- **MCP Servers:** mcp-multiomics
- **Audience:** Researchers, Bioinformaticians, Clinicians
- **Expected Output:** Complete analysis with preprocessing, integration, meta-analysis, regulators

---

### mcp-spatialtools Server
**Location:** `servers/mcp-spatialtools/README.md`, `servers/mcp-spatialtools/QUICKSTART.md`

#### 3.1 Quick Analysis Prompts

**Spatial Cell Type Deconvolution:**
```
Analyze the spatial transcriptomics data for Patient-001.
Perform cell type deconvolution and identify key cell populations.
```
- **MCP Servers:** mcp-spatialtools
- **Audience:** Researchers
- **Expected Output:** Cell type composition, summary statistics

**Pathway Enrichment Analysis:**
```
For the upregulated genes [TP53, BRCA1, MYC, KRAS],
perform pathway enrichment analysis using GO_BP database.
```
- **MCP Servers:** mcp-spatialtools
- **Audience:** Researchers, Bioinformaticians
- **Expected Output:** Enriched pathways, p-values, fold enrichment

#### 3.2 Visualization Prompts

**Generate Spatial Heatmaps:**
```
Generate spatial heatmaps for proliferation markers:
- Expression: /data/PAT001-OVC-2025/visium_gene_expression.csv
- Coordinates: /data/PAT001-OVC-2025/visium_spatial_coordinates.csv
- Genes: MKI67, PCNA, TOP2A, AURKA, CDK1, CCNB1
- Use "plasma" colormap
```
- **MCP Servers:** mcp-spatialtools
- **Audience:** Researchers
- **Expected Output:** Multi-panel spatial heatmap figure

**Create Gene-Region Heatmap:**
```
Create a heatmap showing key genes across tumor regions:
- Expression: /data/expression.csv
- Regions: /data/region_annotations.csv
- Genes: MKI67, PCNA, PIK3CA, AKT1, ABCB1, CD3D, CD8A, CD68
```
- **MCP Servers:** mcp-spatialtools
- **Audience:** Researchers
- **Expected Output:** Annotated gene × region heatmap

#### 3.3 Complete Workflow Prompts

**Complete Patient Spatial Analysis:**
```
Claude, analyze spatial data for Patient-001 (ovarian cancer):

STEP 1: Clinical-Spatial Integration
- Get patient clinical data (epic server)
- Retrieve spatial dataset (bridge tool)

STEP 2: Quality Control
- Filter low-quality spots (min 200 UMIs, 100 genes)
- Check spatial patterns with Moran's I for tumor markers

STEP 3: Regional Analysis
- Split into tumor core, margin, and immune zones
- Perform differential expression between regions

STEP 4: Cell Type Deconvolution
- Estimate cell types per spot
- Quantify endothelial cells (bevacizumab targets)
- Assess CD8+ T-cell infiltration (immunotherapy potential)

STEP 5: Pathway Enrichment
- Run enrichment on upregulated genes from tumor margin
- Identify resistance mechanisms

STEP 6: Integrated Report
- Link spatial findings to clinical context
- Provide treatment recommendations
```
- **MCP Servers:** mcp-spatialtools, mcp-epic
- **Audience:** Researchers, Clinicians
- **Expected Output:** Complete spatial analysis report with visualizations and clinical recommendations

**Batch Correction Workflow:**
```
Correct batch effects across 3 sequencing runs:
- Batch 1: /data/run1_expression.csv
- Batch 2: /data/run2_expression.csv
- Batch 3: /data/run3_expression.csv
- Use ComBat method and save to /data/corrected.csv
```
- **MCP Servers:** mcp-spatialtools
- **Audience:** Researchers, Bioinformaticians
- **Expected Output:** Batch-corrected data, variance reduction statistics

---

### mcp-fgbio Server
**Location:** `servers/mcp-fgbio/README.md`

**Fetch Reference Genome:**
```
Can you fetch the hg38 reference genome and save it to /workspace/data/reference?
```
- **MCP Servers:** mcp-fgbio
- **Audience:** Developers, Bioinformaticians
- **Expected Output:** Downloaded genome path, size, MD5 checksum

**Validate FASTQ File:**
```
Please validate the FASTQ file at /data/sample.fastq.gz and check if the average quality is above 25.
```
- **MCP Servers:** mcp-fgbio
- **Audience:** Bioinformaticians
- **Expected Output:** Validation status, quality statistics

**Query Gene Annotations:**
```
What are the genomic coordinates for the BRCA1 gene in hg38?
```
- **MCP Servers:** mcp-fgbio
- **Audience:** Researchers, Bioinformaticians
- **Expected Output:** Gene coordinates, chromosome, strand

---

### mcp-perturbation Server
**Location:** `servers/mcp-perturbation/README.md`, `docs/architecture/perturbation/README.md`

#### 4.1 Treatment Response Prediction

**Predict Checkpoint Inhibitor Response:**
```
Load the GSE184880 ovarian cancer dataset, setup a GEARS model,
train it, and predict how Patient-001's T cells will respond to
checkpoint inhibitor therapy (PD1 + CTLA4 combination).

Show the top 10 upregulated genes.
```
- **MCP Servers:** mcp-perturbation
- **Audience:** Researchers, Clinicians
- **Expected Output:** Predicted gene expression changes, upregulated immune markers (IFNG, TNF), biomarkers
- **Time Estimate:** 8-12 minutes (includes training)

**Screen Multiple Treatment Options:**
```
For Patient-001 with ovarian cancer, predict treatment responses for:
1. Checkpoint inhibitors (PD1 + CTLA4)
2. PARP inhibitors
3. Platinum-based therapy

Rank them by predicted efficacy and identify biomarkers
for the most promising option.
```
- **MCP Servers:** mcp-perturbation
- **Audience:** Clinicians, Researchers
- **Expected Output:** Treatment rankings, predicted gene changes, therapy recommendation
- **Time Estimate:** 10-15 minutes

#### 4.2 Biomarker Discovery

**Identify Resistance Biomarkers:**
```
Using the trained GEARS model, predict which genes are most
differentially expressed when cells are exposed to anti-PD1 therapy.

Perform differential expression analysis and identify the top 20
genes that predict treatment resistance.
```
- **MCP Servers:** mcp-perturbation
- **Audience:** Researchers
- **Expected Output:** Differentially expressed genes, statistical significance, resistance markers
- **Time Estimate:** 3-5 minutes (after model training)

#### 4.3 Complete Workflow Prompts

**Complete Perturbation Prediction Workflow:**
```
Claude, analyze treatment response potential for Patient-001:

STEP 1: Load Data
- Dataset: GSE184880 (ovarian cancer T-cell perturbations)
- Cell type: T-cells
- Extract single-cell RNA-seq data

STEP 2: Model Setup & Training
- Setup GEARS graph neural network model
- Train on perturbation data (20 epochs)
- Report training metrics (Pearson correlation, loss)

STEP 3: Predict Treatment Responses
- Predict response to: PD1 knockout, CTLA4 knockout, combination
- Extract top upregulated and downregulated genes
- Generate UMAP visualization

STEP 4: Differential Expression
- Compare perturbed vs control states
- Statistical testing with FDR correction
- Identify biomarkers (genes with largest changes)

STEP 5: Clinical Interpretation
- Summarize predicted response (likely responder vs non-responder)
- Link findings to immunotherapy literature
- Provide treatment recommendation with supporting evidence
```
- **MCP Servers:** mcp-perturbation
- **Audience:** Researchers, Clinicians
- **Expected Output:** Complete treatment response prediction with biomarkers and clinical recommendation
- **Time Estimate:** 12-18 minutes

---

## 3. UI-Specific Example Prompts

### Streamlit App
**Location:** `ui/streamlit-app/README.md`

**Quick Spatial Analysis:**
```
Analyze the spatial transcriptomics data for Patient-001.
Perform cell type deconvolution and identify key cell populations.
```
- **Interface:** Streamlit Web UI
- **Audience:** Clinicians, Researchers
- **Expected Output:** Cell type composition displayed in chat

**Multi-omics Integration:**
```
Integrate RNA, protein, and phosphorylation data.
Run HAllA association analysis and identify significant correlations.
```
- **Interface:** Streamlit Web UI
- **Audience:** Researchers
- **Expected Output:** Correlation results displayed in chat

**Complete PatientOne Workflow:**
```
For Patient-001 (ovarian cancer):
1. Get clinical data from FHIR
2. Retrieve spatial transcriptomics data
3. Perform cell type deconvolution
4. Run differential expression between tumor core and margin
5. Generate treatment recommendations
```
- **Interface:** Streamlit Web UI
- **Audience:** Clinicians
- **Expected Output:** Comprehensive clinical report in chat

---

### Jupyter Notebook
**Location:** `ui/jupyter-notebook/README.md`

**Quick Test:**
```python
client = MCPClient()
result = client.call_servers(
    prompt="List all available tools from the spatialtools server.",
    servers=["spatialtools"]
)
print(result["response"])
```
- **Interface:** Jupyter Notebook
- **Audience:** Data Scientists, Researchers
- **Expected Output:** List of available tools

**Complete PatientOne Workflow:**
```python
result = client.call_servers(
    prompt="""For Patient-001 (ovarian cancer):
    1. Get clinical data
    2. Retrieve spatial transcriptomics
    3. Perform cell type deconvolution
    4. Run differential expression
    5. Generate treatment recommendations""",
    servers=["spatialtools", "multiomics"]
)
```
- **Interface:** Jupyter Notebook
- **Audience:** Data Scientists, Researchers
- **Expected Output:** Comprehensive analysis results with token usage stats

---

## 4. Quick Start Workflow Prompts

### QUICKSTART.md Examples
**Location:** `servers/mcp-spatialtools/QUICKSTART.md`

**Batch Correction:**
```python
async def correct_batch_effects():
    result = await perform_batch_correction.fn(
        data_files=[
            "/path/to/batch1_expression.csv",
            "/path/to/batch2_expression.csv",
            "/path/to/batch3_expression.csv"
        ],
        batch_labels=["batch1", "batch2", "batch3"],
        parametric=True
    )
```
- **Format:** Python async code
- **Audience:** Developers, Bioinformaticians
- **Expected Output:** Batch-corrected data with variance reduction metrics

**Pathway Enrichment:**
```python
async def analyze_pathways():
    de_genes = ["TP53", "PIK3CA", "AKT1", "PTEN", "BRCA1"]
    result = await perform_pathway_enrichment.fn(
        differential_genes=de_genes,
        all_genes=all_genes,
        fdr_threshold=0.05
    )
```
- **Format:** Python async code
- **Audience:** Developers, Bioinformaticians
- **Expected Output:** Enriched pathways with statistical significance

---

## 5. Summary Statistics

### Prompts by Audience
- **Clinicians:** 11 prompts (TEST_1, TEST_5, TEST_6, UI examples, perturbation prediction)
- **Researchers:** 30+ prompts (all TEST prompts, spatial analysis, visualization, GEARS)
- **Bioinformaticians:** 20+ prompts (multiomics, preprocessing, batch correction)
- **Developers:** 10+ prompts (server testing, API calls, quickstart examples)
- **Funders/Reviewers:** 3 prompts (executive summaries, ROI analysis)
- **Educators/Students:** 5+ prompts (quick tests, tutorial examples)

### Prompts by MCP Server
- **mcp-multiomics:** 12+ prompts (preprocessing, integration, meta-analysis, upstream regulators)
- **mcp-spatialtools:** 15+ prompts (spatial analysis, visualization, cell type deconvolution)
- **mcp-perturbation:** 5 prompts (treatment response prediction, biomarker discovery, GEARS workflows)
- **mcp-fgbio:** 4 prompts (reference genome, FASTQ validation, annotations)
- **mcp-mockepic:** 2 prompts (clinical data retrieval)
- **mcp-tcga:** 2 prompts (TCGA cohort comparison)
- **mcp-openimagedata:** 3 prompts (imaging analysis, H&E)
- **mcp-deepcell:** 3 prompts (cell segmentation, IF analysis)
- **Multiple servers:** 10+ prompts (integrated workflows)

### Prompts by Complexity
- **Quick (<5 min):** 15 prompts (list tools, query gene, validate file)
- **Medium (5-15 min):** 20 prompts (spatial analysis, pathway enrichment, data integration)
- **Comprehensive (15-45 min):** 15 prompts (complete TEST workflows, preprocessing pipelines)

### Prompts by Expected Output Type
- **Data Tables:** 10 prompts (gene expression, statistics)
- **Visualizations:** 12 prompts (heatmaps, spatial plots, segmentation overlays)
- **Clinical Reports:** 8 prompts (patient summaries, treatment recommendations)
- **Technical Results:** 15 prompts (QC metrics, p-values, pathway enrichment)
- **Integrated Multi-Modal:** 5 prompts (synthesis across genomics, spatial, imaging)

---

## 6. Gaps Identified

### Missing Prompt Categories
1. **Educational Prompts:** Need beginner-friendly prompts for classroom/tutorial use
2. **Funder Demo Prompts:** Need 5-7 high-impact prompts for grant reviewers
3. **Developer Testing Prompts:** Need systematic server validation prompts
4. **Troubleshooting Prompts:** Need diagnostic prompts for common errors
5. **Performance Benchmarking Prompts:** Need prompts for cost/time analysis

### Audience Coverage Gaps
- **Funders/Grant Reviewers:** Only 2-3 prompts, need more executive-level demonstrations
- **Educators:** Limited tutorial-style prompts
- **Hospital IT/Admins:** No deployment validation prompts
- **Regulatory/Compliance:** No HIPAA audit prompts

---

## 7. Next Steps

1. **Create funder-demo-prompts.md** with 5-7 high-impact prompts
2. **Create clinical-workflow-prompts.md** with 15 comprehensive prompts
3. **Create developer-test-prompts.md** with 9+ server validation prompts
4. **Create educational-prompts.md** with 8-10 tutorial prompts
5. **Create main README.md** to index and navigate all prompt collections

---

**Document Version:** 1.0
**Date:** 2026-01-16
**Status:** Discovery complete, ready for prompt curation
