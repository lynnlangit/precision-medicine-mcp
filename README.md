# Spatial + Multiomics MCP(s) Demonstration POC

AI-Orchestrated Spatial Transcriptomics (and Multiomics) Bioinformatics Workflows using Model Context Protocol

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![MCP](https://img.shields.io/badge/MCP-2025--06--18-green.svg)](https://modelcontextprotocol.io/)
[![Claude Desktop](https://img.shields.io/badge/Claude-Desktop-orange.svg)](https://claude.ai/download)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

## What's In It For You?

**Stop writing glue code.**   
Describe your analysis goals in natural language and let Claude orchestrate bioinformatics workflows across 9 specialized MCP servers using up to 36 associated tools.

**Key Benefits:**
- ✅ Replace bash (or workflow) scripts with **conversational requests**
- ✅ Instant access to FASTQ QC, alignment, quantification, multi-omics meta-analysis
- ✅ Reproducible by default - logged, versioned, repeatable
- ✅ Modular architecture - add tools **without** rewriting pipelines  

---

## Which MCP Servers are Here?

### 1. **FGbio** - FASTQ & Genomic Data Processing (4 tools)
- `fetch_reference_genome` - Download reference genome sequences
- `validate_fastq` - Quality validation of FASTQ files
- `extract_umis` - UMI extraction and processing
- `query_gene_annotations` - Retrieve gene annotation data

### 2. **SpatialTools** - Spatial Transcriptomics Analysis (8 tools)
- `filter_quality` - QC filtering of spatial barcodes
- `split_by_region` - Segment data by spatial regions
- `align_spatial_data` - Align reads to reference genome using STAR
- `merge_tiles` - Combine multiple spatial tiles
- `calculate_spatial_autocorrelation` - Calculate spatial gene expression patterns
- `perform_differential_expression` - Differential expression analysis
- `perform_batch_correction` - Batch correction across samples
- `perform_pathway_enrichment` - Pathway enrichment analysis

### 3. **OpenImageData** - Histology Image Processing (3 tools)
- `fetch_histology_image` - Retrieve tissue histology images
- `register_image_to_spatial` - Align histology images with spatial coordinates
- `extract_image_features` - Extract computer vision features from histology

### 4. **Seqera** - Nextflow Workflow Management (3 tools)
- `launch_nextflow_pipeline` - Execute Nextflow workflows via Seqera Platform
- `monitor_workflow_status` - Track pipeline execution status
- `list_available_pipelines` - Query nf-core and custom pipelines

### 5. **HuggingFace** - Genomic Language Models (3 tools)
- `load_genomic_model` - Load pre-trained genomic language models
- `predict_cell_type` - Cell type classification using foundation models
- `embed_sequences` - Generate embeddings for DNA/RNA sequences

### 6. **DeepCell** - Cell Segmentation & Analysis (2 tools)
- `segment_cells` - Deep learning-based cell segmentation
- `classify_cell_states` - Cell state/phenotype classification

### 7. **MockEpic** - Clinical Data (Mock EHR) (3 tools)
- `query_patient_records` - Retrieve patient demographics and clinical data
- `link_spatial_to_clinical` - Connect spatial data to clinical outcomes
- `search_diagnoses` - Query ICD-10 diagnosis codes

### 8. **TCGA** - The Cancer Genome Atlas (5 tools)
- `query_tcga_cohorts` - Search TCGA datasets by cancer type
- `fetch_expression_data` - Download gene expression data
- `compare_to_cohort` - Compare sample expression to TCGA cohort
- `get_survival_data` - Retrieve survival data correlated with expression
- `get_mutation_data` - Retrieve mutation frequencies from cohort

### 9. **MultiOmics** - Multi-Omics Integration (5 tools)
- `integrate_omics_data` - Integrate RNA, protein, and phosphorylation data
- `run_halla_analysis` - HAllA hierarchical all-against-all association testing
- `calculate_stouffer_meta` - Combine p-values across omics modalities
- `create_multiomics_heatmap` - Create integrated heatmap visualization
- `run_multiomics_pca` - PCA on integrated multi-omics data

**Total: 9 servers, 36 tools**

---

## Quick Start

IMPORTANT: In this POC all MCP servers are running locally and are expected to use your local Claude Desktop as their client.  

```bash
# Install (5 min)
git clone https://github.com/your-org/spatial-mcp.git
cd spatial-mcp/manual_testing
./install_dependencies.sh

# Configure Claude Desktop to use the custom MCP Servers
cp ../configs/claude_desktop_config.json ~/Library/Application\ Support/Claude/claude_desktop_config.json

# Verify (restart Claude Desktop first)
./verify_servers.sh
```

**Prerequisites:** Python 3.11+, Claude Desktop, 16GB RAM, 50GB disk

---

## Example Spatial-Multiomics Prompts

**Example of a custom MCP server (mcpFGBio) in action (using Claude Desktop) shown below:** 
<kbd><img src="https://github.com/lynnlangit/spatial-mcp/blob/main/data/images/Claude-client.png" width=800></kbd>

**Multiomics Analysis:**
```
Analyze PDX treatment resistance with RNA, Protein, Phospho data for TP53, MYC, KRAS:
- RNA p-values: [0.001, 0.002, 0.05], log2FC: [2.5, 1.8, 1.2]
- Protein p-values: [0.005, 0.01, 0.03], log2FC: [2.0, 1.6, 1.1]

Combine using Stouffer's method with directionality and FDR correction.
```

**Spatial Analysis:**
```
Process 10x Visium: fetch hg38 → validate FASTQ → extract UMIs → align → quantify → compare TCGA
```

[View 18 example prompts →](docs/spatial/MCP_POC_Example_Prompts.md)

**End-to-end Patient Spatial-Multiomics Analysis (Ovarian Cancer example):**
```
Full multiomics AND spatial integrated analysis of patient data using all custom mcp servers and associated tools  
```
[See full example of full patient workflow prompt →](https://github.com/lynnlangit/spatial-mcp/blob/main/manual_testing/PatientOne-OvarianCancer/Synthetic_sample_data/END_TO_END_TEST_PROMPT.md)

---

## Resources

**References:**
- [MCP Specification](https://modelcontextprotocol.io/specification/2025-06-18)
- [FastMCP Docs](https://github.com/modelcontextprotocol/python-sdk)
- [BioinfoMCP Paper](https://arxiv.org/html/2510.02139v1)
- [Spatial Transcriptomics Review](https://academic.oup.com/nar/article/53/12/gkaf536/8174767)

**Acknowledgments:** Model Context Protocol (Anthropic), BioinfoMCP, FGbio, TCGA, Seqera Platform

---


**Last Updated:** November 12, 2025
**Built with ❤️ for the bioinformatics community**
