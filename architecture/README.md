# Spatial Multiomics MCP Servers POC

See each subfolder for architectures. See also main Repo README.md for list of tools for each custom MCP server

## Architecture Workflows

### 1) Multiomics (custom mcp servers + tools)

**Purpose:** PDX multi-omics data integration with Stouffer's meta-analysis

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    MULTIOMICS WORKFLOW ARCHITECTURE                      │
└─────────────────────────────────────────────────────────────────────────┘

                         Claude Desktop (MCP Host)
                                    │
                    ┌───────────────┼───────────────┐
                    │                               │
                    ▼                               ▼
         ┌──────────────────────┐        ┌──────────────────────┐
         │  mcp-multiomics      │        │  mcp-mockepic        │
         │  ─────────────────   │        │  ───────────────     │
         │  • Data Integration  │        │  • Patient Data      │
         │  • HAllA Analysis    │        │  • Clinical Metadata │
         │  • Stouffer's Meta   │        └──────────────────────┘
         │  • Visualizations    │
         └──────────┬───────────┘
                    │
    ┌───────────────┼───────────────────────┐
    │               │                       │
    ▼               ▼                       ▼
┌────────┐    ┌──────────┐         ┌────────────┐
│  RNA   │    │ Protein  │         │  Phospho   │
│  Data  │    │   Data   │         │    Data    │
└────────┘    └──────────┘         └────────────┘
    │               │                       │
    └───────────────┼───────────────────────┘
                    │
            ┌───────▼────────┐
            │  HAllA Analysis│
            │  (per modality)│
            └───────┬────────┘
                    │
         ┌──────────▼──────────┐
         │  Stouffer's Method  │
         │  (meta-analysis)    │
         │  • p-value combine  │
         │  • Directionality   │
         │  • FDR correction   │
         └──────────┬──────────┘
                    │
         ┌──────────▼──────────┐
         │   Visualization     │
         │  • Heatmaps         │
         │  • PCA plots        │
         │  • Pathway results  │
         └─────────────────────┘

Key Features:
  • Integrates RNA, Protein, Phosphorylation data
  • Statistical meta-analysis across modalities
  • Identifies multi-modal resistance signatures
  • Suitable for PDX treatment response studies
```

### 2) Spatial (custom mcp servers + tools)

**Purpose:** Spatial transcriptomics bioinformatics pipeline

```
┌─────────────────────────────────────────────────────────────────────────┐
│               SPATIAL TRANSCRIPTOMICS WORKFLOW ARCHITECTURE              │
└─────────────────────────────────────────────────────────────────────────┘

                         Claude Desktop (MCP Host)
                      AI-Orchestrated Workflow Execution
                                    │
        ┌───────────────────────────┼───────────────────────────┐
        │                           │                           │
┌───────▼──────────┐    ┌──────────▼─────────┐    ┌───────────▼────────┐
│  STAGE 1:        │    │   STAGE 2:         │    │   STAGE 3:         │
│  Data Ingestion  │───▶│   Spatial          │───▶│   Sequence         │
│  & QC            │    │   Segmentation     │    │   Alignment        │
└──────────────────┘    └────────────────────┘    └────────────────────┘
│ • mcp-fgbio      │    │ • mcp-spatialtools │    │ • mcp-fgbio        │
│   - validate_    │    │   - split_by_      │    │   - fetch_ref      │
│     fastq        │    │     region         │    │ • mcp-spatialtools │
│   - extract_     │    │ • mcp-openimagedata│    │   - align_spatial  │
│     umis         │    │   - fetch_         │    │ • mcp-seqera       │
│ • mcp-spatial    │    │     histology      │    │   - launch_nf      │
│   - filter_      │    │   - register_      │    │                    │
│     quality      │    │     image          │    │ Output: BAM files  │
│                  │    │                    │    │         w/ spatial │
│ Input: FASTQ +   │    │ Output: Segmented  │    │         tags       │
│        barcodes  │    │         regions    │    │                    │
└──────────────────┘    └────────────────────┘    └────────────────────┘
        │                           │                           │
        └───────────────────────────┼───────────────────────────┘
                                    │
        ┌───────────────────────────┼───────────────────────────┐
        │                           │                           │
┌───────▼──────────┐    ┌──────────▼─────────┐                │
│  STAGE 4:        │    │   STAGE 5:         │                │
│  Expression      │───▶│   Analysis &       │                │
│  Quantification  │    │   Integration      │                │
└──────────────────┘    └────────────────────┘                │
│ • mcp-spatialtools│   │ • mcp-seqera       │                │
│   - count_umis    │   │   - run_rnaseq     │                │
│ • mcp-deepcell    │   │ • mcp-huggingface  │                │
│   - segment_cells │   │   - predict_cell   │                │
│ • mcp-huggingface │   │     _type          │                │
│   - embed_        │   │ • mcp-mockepic     │                │
│     sequences     │   │   - link_spatial   │                │
│                   │   │ • mcp-tcga         │                │
│ Output: Gene x    │   │   - compare_to     │                │
│         Spot/Cell │   │     _tcga          │                │
│         matrix    │   │                    │                │
│                   │   │ Output: Analysis   │                │
│                   │   │         results,   │                │
│                   │   │         reports    │                │
└───────────────────┘   └────────────────────┘                │
                                    │                          │
                                    ▼                          │
                        ┌────────────────────┐                 │
                        │  Final Deliverable │                 │
                        │  ────────────────  │                 │
                        │  • Spatial maps    │                 │
                        │  • Differential    │                 │
                        │    expression      │                 │
                        │  • Cell types      │                 │
                        │  • Visualizations  │                 │
                        └────────────────────┘                 │
                                                               │
┌──────────────────────────────────────────────────────────────┘
│  MCP Servers Used (8):
│  ├─ mcp-fgbio          (Reference data, FASTQ validation)
│  ├─ mcp-tcga           (Cancer genomics reference)
│  ├─ mcp-spatialtools   (Core spatial processing)
│  ├─ mcp-huggingface    (ML models)
│  ├─ mcp-mockepic       (Clinical data)
│  ├─ mcp-openimagedata  (Histology images)
│  ├─ mcp-seqera         (Workflow orchestration)
│  └─ mcp-deepcell       (Cell segmentation)
└──────────────────────────────────────────────────────────────┘

Key Features:
  • End-to-end FASTQ → Analysis pipeline
  • Spatial coordinate preservation throughout
  • Integration with histology images
  • AI-assisted cell type identification
  • Comparison to TCGA reference cohorts
```

### 3) PatientOne (combine spatial + multiomics mcp servers) - end-to-end use case

**Purpose:** Comprehensive precision medicine analysis (Stage IV Ovarian Cancer)

```
┌─────────────────────────────────────────────────────────────────────────┐
│            PATIENTONE END-TO-END WORKFLOW ARCHITECTURE                   │
│              (Stage IV Ovarian Cancer Use Case)                          │
└─────────────────────────────────────────────────────────────────────────┘

                         Claude Desktop (MCP Host)
                  Complete Precision Medicine Workflow Orchestration
                                    │
        ┌───────────────────────────┼───────────────────────────┐
        │                           │                           │
        │                           │                           │
┌───────▼──────────┐    ┌──────────▼─────────┐    ┌───────────▼────────┐
│  CLINICAL DATA   │    │   GENOMIC DATA     │    │  MULTIOMICS DATA   │
│  ──────────────  │    │   ────────────     │    │  ──────────────    │
│                  │    │                    │    │                    │
│ mcp-mockepic     │    │ mcp-fgbio          │    │ mcp-multiomics     │
│ • Demographics   │    │ • FASTQ validation │    │ • RNA-seq (PDX)    │
│ • CA-125 trends  │    │ • VCF processing   │    │ • Proteomics       │
│ • Treatment Hx   │    │                    │    │ • Phosphoproteomics│
│                  │    │ mcp-tcga           │    │ • Integration      │
│ Output:          │    │ • TCGA comparison  │    │ • Stouffer's meta  │
│ • Patient profile│    │ • Mutation data    │    │                    │
│ • Clinical       │    │                    │    │ Output:            │
│   context        │    │ Output:            │    │ • Resistance genes │
│                  │    │ • Mutations        │    │ • Pathway analysis │
│                  │    │   (TP53, BRCA1,    │    │ • Multi-modal      │
│                  │    │    PIK3CA)         │    │   signatures       │
│                  │    │ • TCGA subtype     │    │                    │
└──────────────────┘    └────────────────────┘    └────────────────────┘
        │                           │                           │
        └───────────────────────────┼───────────────────────────┘
                                    │
        ┌───────────────────────────┼───────────────────────────┐
        │                           │                           │
┌───────▼──────────┐    ┌──────────▼─────────┐    ┌───────────▼────────┐
│  SPATIAL DATA    │    │   IMAGING DATA     │    │  ANALYSIS & Rx     │
│  ────────────    │    │   ────────────     │    │  ─────────────     │
│                  │    │                    │    │                    │
│ mcp-spatialtools │    │ mcp-openimagedata  │    │ Integration of     │
│ • Visium data    │    │ • H&E histology    │    │ ALL data streams   │
│ • 900 spots      │    │ • Immunofluorescence│   │                    │
│ • 31 genes       │    │   (CD3, CD8, KI67) │    │ • Treatment        │
│ • 6 regions      │    │                    │    │   recommendations  │
│                  │    │ mcp-deepcell       │    │ • Pathway targets  │
│ • Spatial        │    │ • Cell segmentation│    │   (PI3K/AKT/mTOR)  │
│   heterogeneity  │    │ • Cell counting    │    │ • Clinical trials  │
│ • Immune         │    │                    │    │ • Monitoring plan  │
│   localization   │    │ Output:            │    │                    │
│                  │    │ • Tissue           │    │ • Synthetic results│
│ Output:          │    │   architecture     │    │   across all       │
│ • Expression     │    │ • Cell phenotypes  │    │   modalities       │
│   maps           │    │ • Immune infiltrate│    │                    │
│ • Region         │    │                    │    │                    │
│   analysis       │    │                    │    │                    │
└──────────────────┘    └────────────────────┘    └────────────────────┘
        │                           │                           │
        └───────────────────────────┼───────────────────────────┘
                                    │
                                    ▼
                    ┌───────────────────────────┐
                    │  PRECISION MEDICINE       │
                    │  RECOMMENDATIONS          │
                    │  ───────────────────────  │
                    │                           │
                    │  Molecular Profile:       │
                    │  • TP53 R175H (hotspot)   │
                    │  • BRCA1 germline mut     │
                    │  • PIK3CA E545K (resist)  │
                    │                           │
                    │  Resistance Mechanisms:   │
                    │  • PI3K/AKT activation    │
                    │  • MDR1 upregulation      │
                    │  • Anti-apoptotic signals │
                    │                           │
                    │  Treatment Recommendations│
                    │  • PI3K inhibitors        │
                    │    (alpelisib)            │
                    │  • AKT inhibitors         │
                    │    (capivasertib)         │
                    │  • mTOR inhibitors        │
                    │  • Clinical trials        │
                    │                           │
                    │  Monitoring Strategy:     │
                    │  • CA-125 every 3 weeks   │
                    │  • Imaging every 6 weeks  │
                    │  • PDX model validation   │
                    └───────────────────────────┘

┌──────────────────────────────────────────────────────────────────────────┐
│  ALL 9 MCP Servers Utilized:                                             │
│  ├─ mcp-mockepic       (EHR/Clinical data)                              │
│  ├─ mcp-fgbio          (Genomic QC & validation)                        │
│  ├─ mcp-tcga           (TCGA cohort comparison)                         │
│  ├─ mcp-multiomics     (PDX multi-omics integration)                    │
│  ├─ mcp-spatialtools   (Spatial transcriptomics)                        │
│  ├─ mcp-openimagedata  (Histology imaging)                              │
│  ├─ mcp-deepcell       (AI cell segmentation)                           │
│  ├─ mcp-seqera         (Workflow orchestration)                         │
│  └─ mcp-huggingface    (ML model inference)                             │
│                                                                           │
│  Synthetic Data (17 files):                                              │
│  • Clinical: 2 files (demographics, labs)                                │
│  • Genomics: 1 file (VCF with mutations)                                 │
│  • Multiomics: 4 files (RNA/Protein/Phospho + metadata)                  │
│  • Spatial: 3 files (coordinates, expression, annotations)               │
│  • Imaging: 7 files (H&E, IF markers, multiplex)                         │
│                                                                           │
│  Test Location: /manual_testing/PatientOne-OvarianCancer/implementation/ │
└──────────────────────────────────────────────────────────────────────────┘

Key Features:
  • Complete precision medicine workflow
  • Synthetic patient: PAT001-OVC-2025 (58yo, Stage IV HGSOC)
  • Integration of 5 data modalities
  • Resistance mechanism identification
  • Treatment recommendations based on molecular profile
  • Demonstrates all MCP servers working together
```

---

**See subfolder READMEs for detailed architecture documentation:**
- `multiomics/README.md` - Multiomics server architecture
- `spatial/README.md` - Spatial transcriptomics pipeline architecture
- `patient-one/README.md` - PatientOne end-to-end use case