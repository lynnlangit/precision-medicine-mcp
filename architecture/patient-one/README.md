# PatientOne: Comprehensive Precision Medicine Architecture

## Overview

PatientOne demonstrates how Claude orchestrates 9 specialized MCP servers to deliver comprehensive precision medicine analysis for a patient with Stage IV Ovarian Cancer. This end-to-end workflow integrates clinical, genomic, multiomics, spatial, and imaging data to identify resistance mechanisms and actionable treatment targets.

**What makes PatientOne unique:** Unlike traditional bioinformatics pipelines that analyze individual data types in isolation, PatientOne shows how AI can seamlessly integrate across all modalities through natural language—replacing weeks of glue code with conversational requests.

---

## Clinical Context

### Patient Profile (Synthetic Data)

**Patient ID:** PAT001-OVC-2025
**Demographics:** Sarah Anderson (pseudonym), 58-year-old female, Caucasian
**Diagnosis:** Stage IV High-Grade Serous Ovarian Carcinoma (HGSOC)
**FIGO Stage:** IVB (pleural effusion, distant metastases)
**Treatment History:**
- Primary debulking surgery (2023)
- First-line: Carboplatin + Paclitaxel (6 cycles) - Initial response
- Recurrence at 8 months → Platinum-resistant disease
- Second-line: Doxil (pegylated liposomal doxorubicin) - Partial response
- Current status: Stable disease, considering clinical trial enrollment

**Family History:** Mother diagnosed with breast cancer at age 52
**Germline Status:** BRCA1 pathogenic mutation (c.5266dupC, p.Gln1756fs)

**Clinical Markers:**
- CA-125 trajectory:
  - Diagnosis: 1,456 U/mL
  - Post-treatment nadir: 22 U/mL
  - Platinum-resistant recurrence: 389 U/mL
  - Current: 289 U/mL (stable on second-line therapy)

---

## The Precision Medicine Challenge

How do we integrate:
- **Clinical records** from EHR systems (demographics, treatment history, biomarkers)
- **Genomic variants** from whole-exome/genome sequencing (mutations, CNVs)
- **Gene expression** from bulk RNA-seq (transcriptomic profiles)
- **Multi-omics data** from PDX models (RNA/Protein/Phospho resistance signatures)
- **Spatial context** from spatial transcriptomics (tissue microenvironment, immune landscape)
- **Histology imaging** from H&E and IF microscopy (cellular architecture, phenotypes)

...to make **actionable treatment recommendations** that account for:
- Molecular mechanisms of resistance
- Tumor microenvironment composition
- Pathway-level dysregulation
- Available therapeutic options

**Traditional Approach:** Weeks of custom scripts, multiple software tools, manual data wrangling, siloed analysis
**PatientOne with MCP:** Conversational requests that orchestrate 40 tools across 9 servers automatically

### Data Integration Flow

```mermaid
flowchart LR
    subgraph Input["5 Data Modalities"]
        CLIN[Clinical<br/>Demographics<br/>CA-125<br/>Treatment Hx]
        GEN[Genomic<br/>VCF Mutations<br/>CNVs<br/>TCGA Compare]
        OMICS[Multi-Omics<br/>RNA-seq<br/>Proteomics<br/>Phospho]
        SPAT[Spatial<br/>10x Visium<br/>900 spots<br/>6 regions]
        IMG[Imaging<br/>H&E<br/>IF markers<br/>Cell seg]
    end

    subgraph Integration["AI Integration Layer"]
        CLAUDE[Claude Desktop<br/>MCP Orchestration<br/>9 Servers • 40 Tools]
    end

    subgraph Output["Precision Medicine Output"]
        RES[Resistance<br/>Mechanisms]
        TARGETS[Treatment<br/>Targets]
        TRIALS[Clinical<br/>Trials]
    end

    CLIN --> CLAUDE
    GEN --> CLAUDE
    OMICS --> CLAUDE
    SPAT --> CLAUDE
    IMG --> CLAUDE

    CLAUDE --> RES
    CLAUDE --> TARGETS
    CLAUDE --> TRIALS

    style CLAUDE fill:#fff4e1,stroke:#ff9800,stroke-width:3px
    style RES fill:#e8f5e9,stroke:#4caf50,stroke-width:2px
    style TARGETS fill:#e8f5e9,stroke:#4caf50,stroke-width:2px
    style TRIALS fill:#e8f5e9,stroke:#4caf50,stroke-width:2px
```

---

## Architecture Overview

### 5 Integrated Data Modalities

```
┌─────────────────────────────────────────────────────────────┐
│                      PatientOne Workflow                     │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│  1. CLINICAL DATA (MockEpic Server)                         │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│  • Demographics, family history                             │
│  • CA-125 tumor marker trends (4 timepoints)                │
│  • CBC, metabolic panel                                     │
│  • Treatment timeline & response assessment                 │
├─────────────────────────────────────────────────────────────┤
│  2. GENOMIC DATA (FGbio + TCGA Servers)                     │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│  • Somatic mutations: TP53 R175H, PIK3CA E545K, PTEN LOH    │
│  • Copy number: MYC/CCNE1/AKT2/KRAS amplified               │
│  • HRD score: 42 (positive for homologous recombination)    │
│  • TCGA cohort comparison → C1 immunoreactive subtype       │
├─────────────────────────────────────────────────────────────┤
│  3. MULTIOMICS DATA (MultiOmics Server)                     │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│  • 15 PDX samples (7 resistant, 8 sensitive to carboplatin) │
│  • RNA-seq: 1,000 genes × 15 samples                        │
│  • Proteomics: 500 proteins × 15 samples                    │
│  • Phosphoproteomics: 300 sites × 15 samples                │
│  • Analysis: Stouffer's meta-analysis, FDR correction       │
│  • Findings: PI3K/AKT/mTOR pathway activation in resistant  │
├─────────────────────────────────────────────────────────────┤
│  4. SPATIAL DATA (SpatialTools + DeepCell Servers)          │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│  • 10x Visium spatial transcriptomics                       │
│  • 900 spatial spots, 31 genes profiled                     │
│  • 6 tissue regions: tumor_core, proliferative, interface,  │
│    stroma_immune, stroma, necrotic_hypoxic                  │
│  • Spatial heterogeneity in resistance markers              │
│  • Immune exclusion phenotype (CD8+ low/peripheral)         │
├─────────────────────────────────────────────────────────────┤
│  5. IMAGING DATA (OpenImageData + DeepCell Servers)         │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│  • 7 high-resolution TIFF images                            │
│  • H&E histology (20x magnification)                        │
│  • Immunofluorescence: DAPI, CD3, CD8, Ki67, PanCK          │
│  • Multiplex IF: TP53/Ki67/DAPI 3-channel                   │
│  • Cell segmentation, phenotyping, quantification           │
│  • Findings: 70-80% tumor cellularity, Ki67 ~45-55%         │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│               INTEGRATED SYNTHESIS & INSIGHTS                │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│  • Molecular resistance signature: PIK3CA/AKT1/ABCB1 high   │
│  • Microenvironment: Immune exclusion, high stroma          │
│  • Proliferation: Ki67+ tumor core, regional heterogeneity  │
│  • Actionable targets: PI3K inhibitors, immune checkpoint   │
│  • Treatment recommendation: Alpelisib + anti-PD-1 combo    │
└─────────────────────────────────────────────────────────────┘
```

<kbd><img src="https://github.com/lynnlangit/precision-medicine-mcp/blob/main/architecture/patient-one/patient-one-holistic.png"></kbd>

---

## MCP Server Orchestration

### How All 9 Servers Contribute

| Workflow Stage | MCP Servers Engaged | Tools Used | Output |
|---|---|---|---|
| **1. Clinical Retrieval** | MockEpic | `query_patient_records`, `search_diagnoses` | Demographics, CA-125 trends, ICD-10 codes |
| **2. Genomic Analysis** | FGbio, TCGA | `validate_fastq`, `query_gene_annotations`, `compare_to_cohort`, `get_mutation_data` | VCF variants, CNV profile, TCGA subtype |
| **3. Multiomics Integration** | MultiOmics | `integrate_omics_data`, `calculate_stouffer_meta`, `create_multiomics_heatmap` | Resistance gene signatures, pathway activation |
| **4. Spatial Processing** | SpatialTools, DeepCell | `filter_quality`, `split_by_region`, `align_spatial_data`, `segment_cells` | Spatial expression maps, tissue segmentation |
| **5. Histology Analysis** | OpenImageData, DeepCell | `fetch_histology_image`, `register_image_to_spatial`, `classify_cell_states` | Cell counts, phenotype distributions |
| **6. Workflow Orchestration** | Seqera | `launch_nextflow_pipeline`, `monitor_workflow_status` | Reproducible pipeline execution |
| **7. ML Inference** | HuggingFace | `predict_cell_type`, `embed_sequences` | Cell type predictions, sequence embeddings |

**Total Servers:** 9
**Total Tools:** 36
**Integration:** Seamless orchestration through natural language prompts

---

## Data Assets

All synthetic patient data located in: `/data/patient-data/PAT001-OVC-2025/`

### File Inventory (17 files, ~3.2 MB total)

| Modality | Files | Size | Content Description |
|----------|-------|------|---------------------|
| **Clinical** | 2 | 10.7 KB | `clinical_demographics.json`, `ca125_timeline.csv` |
| **Genomic** | 1 | 2.3 KB | `PAT001_somatic_variants.vcf` (12 key variants) |
| **Multiomics** | 4 | 504 KB | `pdx_rna_expression.csv` (1K genes), `pdx_protein.csv` (500), `pdx_phospho.csv` (300), `pdx_metadata.csv` |
| **Spatial** | 3 | 559 KB | `visium_spots.csv` (900 spots), `visium_expression.csv` (31 genes), `visium_regions.csv` (6 regions) |
| **Imaging** | 7 | 2.2 MB | H&E + IF (DAPI, CD3, CD8, Ki67, PanCK) + multiplex |
| **TOTAL** | **17** | **~3.2 MB** | Complete precision medicine dataset |

---

## Key Findings from PatientOne Analysis

### 1. Molecular Resistance Mechanisms

**From Multiomics Integration (MCP-MultiOmics):**
- **PI3K/AKT/mTOR pathway activation** in carboplatin-resistant PDX samples
- Upregulated genes/proteins: `PIK3CA`, `AKT1`, `mTOR`, `RPS6KB1` (Stouffer's combined p < 0.001)
- Drug efflux: `ABCB1` (MDR1) overexpression (log2FC = 2.8, FDR < 0.01)
- Anti-apoptotic: `BCL2L1` upregulation

**From Genomic Analysis (MCP-FGbio + MCP-TCGA):**
- `PIK3CA E545K` activating mutation (allele frequency 38%)
- `TP53 R175H` hotspot mutation (loss of function)
- `PTEN` loss of heterozygosity (tumor suppressor inactivation)
- TCGA subtype: **C1 Immunoreactive** (immune infiltration expected, but...)

### 2. Tumor Microenvironment

**From Spatial Transcriptomics (MCP-SpatialTools):**
- **6 distinct spatial regions** identified
- **Immune exclusion phenotype:** CD8+ T cells enriched at tumor periphery, sparse in core
- **Proliferation gradient:** Ki67/PCNA high in tumor_proliferative region
- **Stroma barrier:** Thick stromal band separating immune cells from tumor

**From Histology Imaging (MCP-OpenImageData + MCP-DeepCell):**
- Tumor cellularity: 70-80%
- Ki67 proliferation index: 45-55% (high)
- CD8+ T cell density: 5-15 cells/mm² (LOW, mostly peripheral)
- CD3+ overall: 30-50 cells/mm² (moderate T cell presence, but not cytotoxic)

### 3. Clinical-Molecular Integration

**From Clinical Data (MCP-MockEpic):**
- **CA-125 response pattern:** Initial deep response (1456 → 22 U/mL) followed by resistance (→ 389 U/mL)
- **BRCA1 germline mutation:** HRD-positive → PARP inhibitor candidate, BUT PIK3CA pathway may confer resistance
- **Platinum-free interval:** 8 months → platinum-resistant category

### 4. Actionable Treatment Recommendations

**Primary Target: PI3K/AKT Pathway**
- Consider: **Alpelisib** (PIK3CA inhibitor) given E545K mutation
- Rationale: Direct target of activating mutation, demonstrated efficacy in PIK3CA-mutant cancers
- Clinical trial: NCT03602859 (alpelisib + paclitaxel in ovarian cancer)

**Secondary Target: Immune Checkpoint**
- Consider: **Anti-PD-1** (pembrolizumab, nivolumab) to overcome immune exclusion
- Rationale: C1 immunoreactive subtype suggests immune-responsive potential
- Combination: Alpelisib + anti-PD-1 may reverse exclusion + activate cytotoxic response

**PARP Inhibitor Re-consideration**
- Given BRCA1 mutation + HRD score 42, PARP inhibitor (olaparib, niraparib) remains option
- Caution: PIK3CA pathway activation may limit efficacy
- Strategy: Sequence after PI3K inhibitor or combination approach

---

## Outputs by Stakeholder

### For Developers (`patient-one-outputs/for-developer/`)
- **MCP_Report_PAT001.pdf:** Technical validation report showing all MCP server calls, data flows, and integration points
- **MCP_Servers_Reference_Guide.pdf:** Complete documentation of 9 servers and 40 tools used
- **Full_Test_Prompt.pdf:** End-to-end prompt that reproduces entire analysis

### For Care Teams (`patient-one-outputs/for-care-team/`)
- **Spatial_Transcriptomics_Analysis.pdf:** Tissue region maps, immune landscape, spatial heterogeneity
- **Histology_Imaging_Analysis.pdf:** Cell segmentation, Ki67 proliferation, CD8 quantification
- **Multiomics_Resistance_Analysis.pdf:** PI3K/AKT pathway activation, resistance gene signatures

### For Patients (`patient-one-outputs/for-patient/`)
- **Medication_Guide.pdf:** Plain-language explanation of recommended therapies
- **Patient_Summary.pdf:** Disease status, test results, next steps
- **Patient_Infographic.pdf:** Visual summary of tumor profile and treatment strategy

---

## Implementation Guide

### Testing Strategy

PatientOne testing is divided into **5 modular tests** to avoid Claude Desktop context limits:

| Test | Focus | Servers | Data Files | Duration | Status |
|------|-------|---------|------------|----------|--------|
| **TEST_1** | Clinical + Genomic | MockEpic, FGbio, TCGA | 3 files | 5-10 min | ✅ |
| **TEST_2** | Multi-Omics | MultiOmics | 4 files | 5-10 min | ✅ |
| **TEST_3** | Spatial Transcriptomics | SpatialTools, DeepCell | 3 files | 5-10 min | ✅ |
| **TEST_4** | Histology & Imaging | OpenImageData, DeepCell | 4 images | 5-10 min | ✅ |
| **TEST_5** | Integration & Recommendations | All servers | Synthesis | 5 min | ✅ |

**Each test is independently runnable and sequentially composable.**

### Quick Start

**Running PatientOne:**

PatientOne supports two modes:
- **DRY_RUN mode** (default): Uses synthetic data for demonstration — works immediately
- **Actual Data mode**: Analyzes your real patient data files — requires configuration

For detailed instructions:
- [PatientOne Quick Start Guide →](../../tests/manual_testing/PatientOne-OvarianCancer/README.md) - How to run the tests
- [Data Modes Configuration Guide →](../../tests/manual_testing/PatientOne-OvarianCancer/DATA_MODES_GUIDE.md) - DRY_RUN vs Actual Data setup

### Test Prompts

All test prompts available in:
- `manual_testing/PatientOne-OvarianCancer/implementation/TEST_1_CLINICAL_GENOMIC.txt`
- `manual_testing/PatientOne-OvarianCancer/implementation/TEST_2_MULTIOMICS.txt`
- `manual_testing/PatientOne-OvarianCancer/implementation/TEST_3_SPATIAL.txt`
- `manual_testing/PatientOne-OvarianCancer/implementation/TEST_4_IMAGING.txt`
- `manual_testing/PatientOne-OvarianCancer/implementation/TEST_5_INTEGRATION.txt`

---

## Why PatientOne Matters

### Paradigm Shift in Precision Medicine

**From:** AI tools for bioinformatics analysis
**To:** AI as orchestrator of complete precision medicine workflows

PatientOne demonstrates that with MCP servers, Claude can seamlessly coordinate across:
- ✅ EHR systems (clinical context)
- ✅ Genomic databases (molecular foundation)
- ✅ Multi-omics platforms (functional landscape)
- ✅ Spatial biology (tissue microenvironment)
- ✅ Medical imaging (cellular morphology)
- ✅ Reference cohorts (comparative context)
- ✅ Workflow engines (reproducibility)
- ✅ ML models (predictive power)

**All in one conversational interface** — replacing weeks of glue code with natural language prompts.

### Real-World Impact

While PatientOne uses synthetic data, the workflow represents a **real clinical decision-making process**:

1. **Tumor Board Preparation:** Integrate all available molecular data before multidisciplinary review
2. **Clinical Trial Matching:** Identify actionable targets and match to available trials
3. **Treatment Selection:** Evidence-based therapy recommendations accounting for resistance mechanisms
4. **Biomarker Monitoring:** Track CA-125 and other markers to assess response

**Future Direction:** Integration with real EHR systems, genomic databases, and clinical workflows for production precision medicine.

---

## Technical Architecture Details

### Data Flow

```
User Prompt → Claude Desktop → MCP Protocol → Server Selection
                                    ↓
                    ┌───────────────────────────────┐
                    │  9 MCP Servers (Local)        │
                    │  Each with specialized tools  │
                    └───────────────────────────────┘
                                    ↓
            ┌───────────────────────────────────────┐
            │  Data Sources (Local Files)           │
            │  /data/patient-data/PAT001-OVC-2025/  │
            └───────────────────────────────────────┘
                                    ↓
            ┌───────────────────────────────────────┐
            │  Tool Execution & Results             │
            │  JSON responses with data             │
            └───────────────────────────────────────┘
                                    ↓
            ┌───────────────────────────────────────┐
            │  Claude Synthesis & Interpretation    │
            │  Integrated insights, recommendations │
            └───────────────────────────────────────┘
                                    ↓
                        User-friendly output
```

### Reproducibility

All PatientOne analyses are:
- **Logged:** Complete MCP server call history
- **Versioned:** Git-tracked data and configurations
- **Repeatable:** Deterministic tool execution (DRY_RUN mode available)
- **Auditable:** Full trace from raw data to recommendations

---

## Related Documentation

- **Main Project:** [Precision Medicine MCP Servers →](../../README.md)
- **Spatial Workflow:** [Architecture →](../spatial/README.md)
- **Multiomics Workflow:** [Architecture →](../multiomics/README.md)
- **Testing Guide:** [PatientOne Implementation →](../../manual_testing/PatientOne-OvarianCancer/README.md)

---

**Last Updated:** December 26, 2025
**Version:** 1.0
**Status:** Demonstration POC with synthetic data
