# Precision Medicine MCP Servers

AI-Orchestrated Clinical Bioinformatics for Precision Oncology using Model Context Protocol

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![MCP](https://img.shields.io/badge/MCP-2025--06--18-green.svg)](https://modelcontextprotocol.io/)
[![Claude Desktop](https://img.shields.io/badge/Claude-Desktop-orange.svg)](https://claude.ai/download)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

## Statement of Purpose

**Transform clinical decision-making with AI-orchestrated bioinformatics**
- Analyze complete patient profiles—from genomics to spatial transcriptomics—using natural language
- Demonstrate end-to-end precision medicine workflows for Stage IV Ovarian Cancer across 9 specialized MCP servers and 40 tools
- Named 'Patient-One' in memory of a dear friend who passed away in 2025 due to HGSOC

---

## Featured Use Case: PatientOne

**Comprehensive Precision Medicine Workflow for Stage IV Ovarian Cancer**

<kbd><img src="https://github.com/lynnlangit/spatial-mcp/blob/main/architecture/patient-one/patient-one-holistic.png" width=800></kbd>

PatientOne demonstrates how Claude orchestrates ALL 9 MCP servers to analyze a complete patient profile:

**Patient Profile (Synthetic):**
- 58-year-old female with Stage IV High-Grade Serous Ovarian Cancer (HGSOC)
- Platinum-resistant, BRCA1 germline mutation
- Post-surgery, considering experimental therapies

**5 Integrated Data Modalities:**
1. **Clinical** (MockEpic) - Demographics, CA-125 trends, treatment history
2. **Genomic** (FGbio, TCGA) - Somatic mutations (TP53, PIK3CA, PTEN), CNVs, HRD score
3. **Multiomics** (MultiOmics) - RNA/Protein/Phospho from 15 PDX samples (resistant vs sensitive)
4. **Spatial** (SpatialTools) - 10x Visium spatial transcriptomics (900 spots, 31 genes, 6 regions)
5. **Imaging** (OpenImageData, DeepCell) - H&E histology, multiplex IF, cell segmentation

**Key Findings from PatientOne Analysis:**
- PI3K/AKT/mTOR pathway activation in resistant samples
- Immune exclusion phenotype (low CD8+ infiltration)
- High proliferation (Ki67+) in tumor core regions
- Actionable targets: PIK3CA inhibitors, immune checkpoint combinations

**Try PatientOne:** [Quick Start Guide →](manual_testing/PatientOne-OvarianCancer/README.md)  
**Architecture Details:** [PatientOne Documentation →](architecture/patient-one/README.md)  
**Outputs Generated:** [For developers, care teams, and patients →](architecture/patient-one/patient-one-outputs/)  

---

## Component Workflows

PatientOne integrates these specialized workflows:

### 1. Spatial Transcriptomics Analysis
Process 10x Visium spatial RNA-seq: tissue segmentation → alignment → quantification → spatial expression mapping

**Example Prompt:**
```
Process 10x Visium: fetch hg38 → validate FASTQ → extract UMIs → align → quantify → compare TCGA
```

[Spatial Architecture →](architecture/spatial/README.md) | [18 Example Prompts →](docs/spatial/MCP_POC_Example_Prompts.md)

### 2. Multiomics Integration
Integrate RNA, protein, and phosphorylation data using HAllA association testing and Stouffer's meta-analysis

**Example Prompt:**
```
Analyze PDX treatment resistance with RNA, Protein, Phospho data for TP53, MYC, KRAS:
- RNA p-values: [0.001, 0.002, 0.05], log2FC: [2.5, 1.8, 1.2]
- Protein p-values: [0.005, 0.01, 0.03], log2FC: [2.0, 1.6, 1.1]

Combine using Stouffer's method with directionality and FDR correction.
```

[Multiomics Architecture →](architecture/multiomics/README.md) | [Detailed Guide →](docs/multiomics/README.md)

### 3. Clinical-Genomic Analysis
Link patient EHR data with genomic variants, compare to TCGA cohorts, identify molecular subtypes

**Example Prompt:**
```
Query patient PAT001 clinical records → extract genomic variants from VCF → compare to TCGA ovarian cancer cohort → identify molecular subtype
```

---

## Which MCP Servers Are Here?

### Foundational Tools (All Workflows)

**1. FGbio - FASTQ & Genomic Data Processing (4 tools)**
- `fetch_reference_genome` - Download reference genome sequences
- `validate_fastq` - Quality validation of FASTQ files
- `extract_umis` - UMI extraction and processing
- `query_gene_annotations` - Retrieve gene annotation data

**4. Seqera - Nextflow Workflow Management (3 tools)**
- `launch_nextflow_pipeline` - Execute Nextflow workflows via Seqera Platform
- `monitor_workflow_status` - Track pipeline execution status
- `list_available_pipelines` - Query nf-core and custom pipelines

**5. HuggingFace - Genomic Language Models (3 tools)**
- `load_genomic_model` - Load pre-trained genomic language models
- `predict_cell_type` - Cell type classification using foundation models
- `embed_sequences` - Generate embeddings for DNA/RNA sequences

### Spatial Biology Tools

**2. SpatialTools - Spatial Transcriptomics Analysis (8 tools)**
- `filter_quality` - QC filtering of spatial barcodes
- `split_by_region` - Segment data by spatial regions
- `align_spatial_data` - Align reads to reference genome using STAR
- `merge_tiles` - Combine multiple spatial tiles
- `calculate_spatial_autocorrelation` - Calculate spatial gene expression patterns
- `perform_differential_expression` - Differential expression analysis
- `perform_batch_correction` - Batch correction across samples
- `perform_pathway_enrichment` - Pathway enrichment analysis

**3. OpenImageData - Histology Image Processing (3 tools)**
- `fetch_histology_image` - Retrieve tissue histology images
- `register_image_to_spatial` - Align histology images with spatial coordinates
- `extract_image_features` - Extract computer vision features from histology

**6. DeepCell - Cell Segmentation & Analysis (2 tools)**
- `segment_cells` - Deep learning-based cell segmentation
- `classify_cell_states` - Cell state/phenotype classification

### Multiomics & Clinical Tools

**9. MultiOmics - Multi-Omics Integration (9 tools)**
- `validate_multiomics_data` - Quality validation before analysis (batch effects, missing values)
- `preprocess_multiomics_data` - Batch correction, imputation, normalization
- `visualize_data_quality` - QC plots (PCA, correlation, before/after comparison)
- `integrate_omics_data` - Integrate RNA, protein, and phosphorylation data
- `run_halla_analysis` - HAllA hierarchical all-against-all association testing
- `calculate_stouffer_meta` - Combine p-values across omics modalities
- `predict_upstream_regulators` - Identify kinases, TFs, and drug targets
- `create_multiomics_heatmap` - Create integrated heatmap visualization
- `run_multiomics_pca` - PCA on integrated multi-omics data

**8. TCGA - The Cancer Genome Atlas (5 tools)**
- `query_tcga_cohorts` - Search TCGA datasets by cancer type
- `fetch_expression_data` - Download gene expression data
- `compare_to_cohort` - Compare sample expression to TCGA cohort
- `get_survival_data` - Retrieve survival data correlated with expression
- `get_mutation_data` - Retrieve mutation frequencies from cohort

**7. MockEpic - Clinical Data (Mock EHR) (3 tools)**
- `query_patient_records` - Retrieve patient demographics and clinical data
- `link_spatial_to_clinical` - Connect spatial data to clinical outcomes
- `search_diagnoses` - Query ICD-10 diagnosis codes

**Total: 9 servers, 40 tools**

---

## Quick Start

IMPORTANT: In this POC all MCP servers are running locally and are expected to use your local Claude Desktop as their client.

```bash
# Install (5 min)
git clone https://github.com/lynnlangit/spatial-mcp.git
cd spatial-mcp/manual_testing
./install_dependencies.sh

# Configure Claude Desktop to use the custom MCP Servers
cp ../configs/claude_desktop_config.json ~/Library/Application\ Support/Claude/claude_desktop_config.json

# Verify (restart Claude Desktop first)
./verify_servers.sh
```

**Prerequisites:** Python 3.11+, Claude Desktop, 16GB RAM, 50GB disk

---

## Example Client Usage

**Example of MCP servers in action using Claude Desktop:**
<kbd><img src="https://github.com/lynnlangit/spatial-mcp/blob/main/data/images/Claude-client.png" width=800></kbd>

**End-to-end Patient Analysis Prompt:**
```
Full multiomics AND spatial integrated analysis of patient data using all custom mcp servers and associated tools
```
[See full PatientOne workflow prompt →](https://github.com/lynnlangit/spatial-mcp/blob/main/manual_testing/PatientOne-OvarianCancer/implementation/END_TO_END_TEST_PROMPT.md)

---

## Resources

**References:**
- [MCP Specification](https://modelcontextprotocol.io/specification/2025-06-18)
- [FastMCP Docs](https://github.com/modelcontextprotocol/python-sdk)
- [BioinfoMCP Paper](https://arxiv.org/html/2510.02139v1)
- [Spatial Transcriptomics Review](https://academic.oup.com/nar/article/53/12/gkaf536/8174767)

**Acknowledgments:** Model Context Protocol (Anthropic), BioinfoMCP, FGbio, TCGA, Seqera Platform

---


**Last Updated:** December 26, 2025
**Built for the precision medicine community**
