# 🏗️ Precision Medicine MCP Server Architectures

Architecture documentation for modality-specific analysis workflows.

---

## 📊 Architecture by Analysis Modality

📋 **[See Individual Server Status →](../servers/README.md#-server-status)** - Detailed tools and documentation for all 10 servers

7 analysis modalities, 10 specialized servers, 55 tools:

| Modality | Servers | Tools | Status | Documentation |
|----------|---------|-------|--------|---------------|
| 🧬 **Clinical Data** | mcp-epic, mcp-mockepic | 7 | ✅ Production/Mock | [clinical/README.md](clinical/README.md) |
| 🧪 **Genomic Cohorts** | mcp-tcga | 5 | ❌ Mocked (GDC-ready) | [genomic/README.md](genomic/README.md) |
| 🖼️ **Imaging** | mcp-openimagedata, mcp-deepcell | 9 | ⚠️ Partial (60%/0%) | [imaging/README.md](imaging/README.md) |
| 🔬 **Multiomics** | mcp-multiomics | 10 | ✅ Production (85%) | [multiomics/README.md](multiomics/README.md) |
| 📍 **Spatial Transcriptomics** | mcp-fgbio, mcp-spatialtools | 18 | ✅ Production (95%) | [spatial-transcriptomics/README.md](spatial-transcriptomics/README.md) |
| 🤖 **AI/ML Inference** | mcp-huggingface | 3 | ❌ Mocked (HF-ready) | [ai-ml/README.md](ai-ml/README.md) |
| ⚙️ **Workflow Orchestration** | mcp-seqera | 3 | ❌ Mocked (Seqera-ready) | [workflow/README.md](workflow/README.md) |

**⚠️ Important:** Not all servers are production-ready. Check [Server Implementation Status](../docs/SERVER_IMPLEMENTATION_STATUS.md) before using.

---

## 🧬 1. Clinical Data Retrieval

**EHR integration for patient clinical context**

**Servers:** mcp-epic (real FHIR, local only) • mcp-mockepic (synthetic, GCP deployed)

**Key Features:**
- FHIR R4 data retrieval (demographics, diagnoses, labs, medications)
- HIPAA Safe Harbor de-identification
- Clinical-spatial outcome linkage

**Workflow:** `Patient EHR → FHIR API → De-identification → Clinical Data`

📖 **[Detailed Architecture →](clinical/README.md)**

---

## 🧪 2. Genomic Cohort Analysis

**TCGA cohort comparison for population-level genomic context**

**Server:** mcp-tcga (33 cancer types, 11,000+ samples)

**Key Features:**
- Gene expression comparison (z-scores, percentiles, p-values)
- Somatic mutation frequency queries
- Survival stratification (Kaplan-Meier, hazard ratios)

**Workflow:** `TCGA Database → Statistical Comparison → Survival Analysis → Integration`

📖 **[Detailed Architecture →](genomic/README.md)**

---

## 🖼️ 3. Imaging Analysis

**Histology and multiplexed immunofluorescence (MxIF) image processing**

**Servers:** mcp-openimagedata (60% real) • mcp-deepcell (mocked)

**Key Workflows:**
- **H&E (Brightfield):** Morphology assessment, necrosis identification
- **MxIF (Fluorescence):** Cell segmentation and quantification

📖 **[Detailed Architecture →](imaging/README.md)**

---

## 🔬 4. Multiomics Integration

**PDX multi-omics data integration with preprocessing and therapeutic target prediction**

**Server:** mcp-multiomics (10 tools, 85% real)

**Key Features:**
- Preprocessing pipeline (batch correction, KNN imputation, QC visualization)
- Association testing (HAllA with chunking)
- Meta-analysis (Stouffer's method)
- Therapeutic targets (kinase/TF/drug prediction)

**Workflow:** `RNA/Protein/Phospho → Validate → Preprocess → Integrate → HAllA → Meta-Analysis → Upstream Regulators`

📖 **[Detailed Architecture →](multiomics/README.md)**

---

## 📍 5. Spatial Transcriptomics

**Spatial gene expression analysis with tissue context**

**Servers:** mcp-fgbio (reference genomes, FASTQ QC) • mcp-spatialtools (spatial analysis, DE, pathway enrichment)

**Key Features:**
- **Analysis Tools (10):** Spatial autocorrelation (Moran's I), differential expression, batch correction, pathway enrichment, cell type deconvolution
- **Visualization Tools (4):** Spatial heatmaps, gene expression heatmaps, region composition charts
- **Bridge Tool:** Integrates with mcp-multiomics

**Workflows:** CSV (current) • FASTQ with STAR alignment (implemented)

📖 **[Detailed Architecture →](spatial-transcriptomics/README.md)**

---

## 🤖 6. AI/ML Model Inference

**Genomic foundation model inference for cell type prediction and sequence embedding**

**Server:** mcp-huggingface (DNABERT-2, Geneformer, Nucleotide-Transformer)

**Key Features:**
- Cell type annotation from single-cell expression data
- DNA/RNA sequence embeddings (768-dimensional vectors)
- Variant effect scoring via embedding distance

**Workflow:** `Single-cell Data / DNA Sequences → Load Model → Generate Predictions → Integration`

📖 **[Detailed Architecture →](ai-ml/README.md)**

---

## ⚙️ 7. Workflow Orchestration

**Nextflow pipeline execution and monitoring via Seqera Platform**

**Server:** mcp-seqera (90+ nf-core workflows)

**Key Features:**
- nf-core pipeline submission (rnaseq, sarek, spatial, etc.)
- Multi-cloud execution (AWS, Azure, GCP, local HPC)
- Real-time progress monitoring and resource tracking

**Workflow:** `Select Pipeline → Configure → Choose Compute → Submit → Monitor → Retrieve Results`

📖 **[Detailed Architecture →](workflow/README.md)**

---

## 🏥 End-to-End Example: PatientOne

**Complete precision medicine workflow combining all 10 MCP servers**

**Use Case:** Stage IV High-Grade Serous Ovarian Cancer (HGSOC), platinum-resistant
**Patient:** PAT001-OVC-2025 (synthetic test case)
**Data Modalities:** Clinical (FHIR) • Genomic (VCF) • Multiomics (RNA/Protein/Phospho) • Spatial (Visium) • Imaging (H&E, MxIF)

**Tests:**
- 🧬 TEST_1: Clinical data retrieval (mcp-epic)
- 🔬 TEST_2: Multiomics integration (mcp-multiomics)
- 📍 TEST_3: Spatial transcriptomics (mcp-spatialtools)
- 🖼️ TEST_4: Imaging analysis (mcp-openimagedata, mcp-deepcell)
- 🎯 TEST_5: Complete end-to-end workflow

📖 **[PatientOne Workflow →](../tests/manual_testing/PatientOne-OvarianCancer/README.md)**
📖 **[PatientOne Architecture →](../tests/manual_testing/PatientOne-OvarianCancer/architecture/README.md)**

---

## 🔗 Related Documentation

### 📋 Operations & Deployment
- [Server Implementation Status](../docs/SERVER_IMPLEMENTATION_STATUS.md) - Production readiness (4/10 ready)
- [GCP Cloud Run Deployment](../docs/deployment/DEPLOYMENT_STATUS.md) - 9 servers deployed ✅
- [Hospital Deployment Guide](../docs/hospital-deployment/) - HIPAA-compliant production setup
- [Cost Analysis](../docs/operations/COST_ANALYSIS.md) - Token costs, ROI analysis

### 🧬 Server Documentation
- [Servers Directory](../servers/) - Individual server READMEs (10 servers)
- [Testing Guide](../tests/README.md) - 167 automated tests ✅

### 📚 Main Docs
- [Main README](../README.md) - Project overview
- [Documentation Index](../docs/README.md) - Complete documentation

---

**Last Updated:** 2026-01-11

**Organization Principle:**
- `architecture/` = High-level design & workflows by modality
- `servers/` = Detailed tool specifications & implementation
- `docs/` = Operational guides & deployment
- `tests/` = End-to-end use cases & validation
