# PatientOne Data Modes Configuration Guide

Complete guide to running PatientOne in DRY_RUN mode (synthetic data) vs Actual Data mode (real analysis).

---

## Table of Contents

- [Overview](#overview)
- [Mode Comparison](#mode-comparison)
- [DRY_RUN Mode (Default)](#dry_run-mode-default)
- [Actual Data Mode](#actual-data-mode)
- [Configuration Steps](#configuration-steps)
- [Data File Placement](#data-file-placement)
- [Expected Outputs](#expected-outputs)
- [Troubleshooting](#troubleshooting)
- [Use Case Recommendations](#use-case-recommendations)

---

## Overview

PatientOne supports two distinct operating modes to accommodate different use cases:

**📊 Quick Cost Comparison:**
- **DRY_RUN:** 25-35 min, ~$0.32 total
- **Actual Data:** 2-4 hours, $15-45 total

**[→ See Full Cost Analysis & ROI](../../../COST_ANALYSIS.md)**

---

### DRY_RUN Mode (Default)
**Purpose:** Demonstration, testing, and validation without external dependencies

- ✅ No real data files required
- ✅ No external API calls (TCGA, HuggingFace, etc.)
- ✅ Fast execution (~5-10 min per test)
- ✅ Safe for CI/CD pipelines
- ✅ Returns realistic synthetic responses

### Actual Data Mode
**Purpose:** Production analysis with real patient data

- 🔬 Processes your actual data files
- 🔬 May make external API calls (TCGA, HuggingFace)
- 🔬 Longer execution time (depends on data size)
- 🔬 Requires proper data file placement
- 🔬 Returns real analysis results

---

## Mode Comparison

### Feature Matrix

| Feature | DRY_RUN Mode | Actual Data Mode |
|---------|--------------|------------------|
| **Data Source** | Synthetic responses | Your files in `/data/patient-data/` |
| **External APIs** | None (mocked) | Real calls (TCGA, HuggingFace, etc.) |
| **Execution Time** | Fast (25-35 min total) | Longer (2-4 hours total) |
| **Cost** | **~$0.32 total** | **$15-45 total** |
| **Cost Breakdown** | Claude tokens only (~30K tokens) | Compute ($14-30) + APIs ($0-5) + Claude tokens |
| **Setup Required** | Minimal | Data files + environment config |
| **File I/O** | Minimal (no writes) | Full (reads/writes) |
| **API Keys** | Not needed | May be required (HF_TOKEN, etc.) |
| **Network** | Not needed | Required for external APIs |
| **Reproducibility** | 100% deterministic | Variable (API responses may change) |
| **Value** | Perfect for learning & workflow validation | Production analysis, replaces ~40 hours manual work |
| **Use Cases** | Demo, testing, CI/CD, education | Real analysis, research, clinical decision support |

### Server Behavior Differences

| MCP Server | DRY_RUN Behavior | Actual Data Behavior |
|------------|------------------|----------------------|
| **mcp-fgbio** | Returns mock BAM/VCF validation | Validates real genomic files |
| **mcp-multiomics** | Returns synthetic resistance signatures | Processes real RNA/Protein/Phospho data |
| **mcp-spatialtools** | Returns mock spatial patterns | Analyzes real Visium/Xenium data |
| **mcp-tcga** | Returns mock cohort comparisons | Queries real TCGA database |
| **mcp-huggingface** | Returns mock predictions | Calls HuggingFace API with real models |
| **mcp-deepcell** | Returns mock segmentation | Runs DeepCell segmentation models |
| **mcp-mockepic** | Returns mock EHR data | Queries mock EHR (always synthetic) |
| **mcp-seqera** | Returns mock workflow status | Launches real Nextflow workflows |
| **mcp-openimagedata** | Returns mock image metadata | Processes real histology images |

---

## DRY_RUN Mode (Default)

### When to Use

✅ **Best for:**
- First-time users exploring PatientOne
- Demonstrations and presentations
- Testing MCP server integration
- CI/CD pipeline validation
- Educational workshops
- Verifying Claude Desktop configuration
- Quick sanity checks

❌ **Not suitable for:**
- Real patient analysis
- Research publications
- Clinical decision support
- Production deployments

### Setup (Minimal)

**DRY_RUN mode is the default** — no additional setup required beyond basic installation.

#### 1. Verify Environment Variables

Check your Claude Desktop configuration:

```bash
cat ~/Library/Application\ Support/Claude/claude_desktop_config.json
```

Look for DRY_RUN environment variables (should all be `"true"`):

```json
{
  "mcpServers": {
    "mcp-fgbio": {
      "command": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-fgbio/venv/bin/python",
      "args": ["-m", "mcp_fgbio"],
      "env": {
        "FGBIO_DRY_RUN": "true"
      }
    },
    "mcp-multiomics": {
      "command": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-multiomics/venv/bin/python",
      "args": ["-m", "mcp_multiomics"],
      "env": {
        "MULTIOMICS_DRY_RUN": "true"
      }
    },
    "mcp-spatialtools": {
      "command": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-spatialtools/venv/bin/python",
      "args": ["-m", "mcp_spatialtools"],
      "env": {
        "SPATIAL_DRY_RUN": "true"
      }
    }
    // ... all 9 servers with DRY_RUN=true
  }
}
```

#### 2. Run PatientOne Tests

Simply copy/paste any test prompt from `implementation/TEST_*.txt` into Claude Desktop.

**Example (TEST_1):**
```
I want to run the PatientOne clinical and genomic analysis.

Please analyze patient PAT001-OVC-2025 using:
1. MockEpic for clinical data
2. FGbio for genomic variants
3. TCGA for cohort comparison

This is TEST_1 from the PatientOne workflow.
```

#### 3. Interpret Results

Responses will include `"status": "success (DRY_RUN mode)"` to indicate synthetic data usage.

**Example response snippet:**
```json
{
  "status": "success (DRY_RUN mode)",
  "patient_demographics": {
    "patient_id": "PAT001-OVC-2025",
    "age": 58,
    "diagnosis": "Stage IV HGSOC",
    "ca125_baseline": 1456,
    "ca125_current": 289
  },
  "note": "This is synthetic data for demonstration purposes"
}
```

### Data Files (Optional)

**DRY_RUN mode does NOT require actual data files**, but you can reference them in prompts for context:

```bash
# Example file structure (for reference in prompts)
/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/data/patient-data/PAT001-OVC-2025/
├── clinical/
│   ├── patient_demographics.json
│   └── lab_results.json
├── genomics/
│   └── somatic_variants.vcf
├── multiomics/
│   ├── pdx_rna_seq.csv
│   ├── pdx_proteomics.csv
│   ├── pdx_phosphoproteomics.csv
│   └── sample_metadata.csv
├── spatial/
│   ├── visium_spatial_coordinates.csv
│   ├── visium_gene_expression.csv
│   └── visium_region_annotations.csv
└── imaging/
    ├── PAT001_tumor_HE_20x.tiff
    ├── PAT001_tumor_IF_DAPI.tiff
    ├── PAT001_tumor_IF_CD3.tiff
    ├── PAT001_tumor_IF_CD8.tiff
    ├── PAT001_tumor_IF_KI67.tiff
    ├── PAT001_tumor_IF_PanCK.tiff
    └── PAT001_tumor_multiplex_IF_TP53_KI67_DAPI.tiff
```

**Note:** In DRY_RUN mode, servers will return synthetic responses even if you provide file paths.

---

## Actual Data Mode

### When to Use

✅ **Best for:**
- Real patient data analysis
- Research studies and publications
- Clinical decision support (with proper validation)
- Custom dataset processing
- Production precision medicine workflows
- Validating PatientOne with your own data

⚠️ **Requirements:**
- Real data files properly formatted
- API keys for external services (optional)
- Sufficient compute resources
- Data privacy compliance (HIPAA, GDPR, etc.)

### Setup (Comprehensive)

#### 1. Configure Environment Variables

Edit Claude Desktop configuration to disable DRY_RUN:

```bash
# Open configuration file
nano ~/Library/Application\ Support/Claude/claude_desktop_config.json
```

**Change all DRY_RUN variables from `"true"` to `"false"`:**

```json
{
  "mcpServers": {
    "mcp-fgbio": {
      "command": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-fgbio/venv/bin/python",
      "args": ["-m", "mcp_fgbio"],
      "env": {
        "FGBIO_DRY_RUN": "false",
        "FGBIO_DATA_DIR": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/data",
        "FGBIO_CACHE_DIR": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/data/cache"
      }
    },
    "mcp-multiomics": {
      "command": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-multiomics/venv/bin/python",
      "args": ["-m", "mcp_multiomics"],
      "env": {
        "MULTIOMICS_DRY_RUN": "false",
        "MULTIOMICS_DATA_DIR": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/data/multiomics",
        "MULTIOMICS_CACHE_DIR": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/data/cache/multiomics"
      }
    },
    "mcp-spatialtools": {
      "command": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-spatialtools/venv/bin/python",
      "args": ["-m", "mcp_spatialtools"],
      "env": {
        "SPATIAL_DRY_RUN": "false",
        "SPATIAL_DATA_DIR": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/data",
        "SPATIAL_CACHE_DIR": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/data/cache"
      }
    },
    "mcp-tcga": {
      "command": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-tcga/venv/bin/python",
      "args": ["-m", "mcp_tcga"],
      "env": {
        "TCGA_DRY_RUN": "false"
        // TCGA may require API configuration
      }
    },
    "mcp-huggingface": {
      "command": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-huggingface/venv/bin/python",
      "args": ["-m", "mcp_huggingface"],
      "env": {
        "HF_DRY_RUN": "false",
        "HF_TOKEN": "your_huggingface_token_here"  // Required for actual API calls
      }
    },
    "mcp-deepcell": {
      "command": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-deepcell/venv/bin/python",
      "args": ["-m", "mcp_deepcell"],
      "env": {
        "DEEPCELL_DRY_RUN": "false"
      }
    },
    "mcp-seqera": {
      "command": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-seqera/venv/bin/python",
      "args": ["-m", "mcp_seqera"],
      "env": {
        "SEQERA_DRY_RUN": "false",
        "SEQERA_ACCESS_TOKEN": "your_seqera_token_here"  // Required for Nextflow workflows
      }
    },
    "mcp-mockepic": {
      "command": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-mockepic/venv/bin/python",
      "args": ["-m", "mcp_mockepic"],
      "env": {
        "EPIC_DRY_RUN": "true"  // Note: MockEpic always uses synthetic data
      }
    },
    "mcp-openimagedata": {
      "command": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-openimagedata/venv/bin/python",
      "args": ["-m", "mcp_openimagedata"],
      "env": {
        "IMAGE_DRY_RUN": "false",
        "IMAGE_DATA_DIR": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/data/images",
        "IMAGE_CACHE_DIR": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/data/cache/images"
      }
    }
  }
}
```

**Important Notes:**
- `mcp-mockepic` always uses synthetic EHR data (no real EHR integration)
- `HF_TOKEN` required for HuggingFace model downloads
- `SEQERA_ACCESS_TOKEN` required for Nextflow Tower/Seqera Platform
- Adjust file paths to match your system

#### 2. Obtain API Keys (If Using External Services)

**HuggingFace (for ML models):**
```bash
# Get token from https://huggingface.co/settings/tokens
export HF_TOKEN="hf_your_token_here"

# Or add to claude_desktop_config.json (recommended)
```

**Seqera Platform (for Nextflow workflows):**
```bash
# Get token from https://cloud.seqera.io → Tokens
export SEQERA_ACCESS_TOKEN="your_token_here"
```

**TCGA API:**
- No token required for basic queries
- Rate limits apply: https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/#api-endpoints

#### 3. Restart Claude Desktop

After modifying configuration:

```bash
# macOS: Quit and restart Claude Desktop app
# The servers will reload with new DRY_RUN=false settings
```

Verify servers loaded correctly in Claude Desktop (should see all 9 servers in available tools).

---

## Configuration Steps

### Step-by-Step: Switching from DRY_RUN to Actual Data

#### Step 1: Backup Current Configuration

```bash
cp ~/Library/Application\ Support/Claude/claude_desktop_config.json \
   ~/Library/Application\ Support/Claude/claude_desktop_config.json.backup
```

#### Step 2: Prepare Your Data Files

Place your patient data in the correct directory structure:

```bash
# Create directory for your patient
mkdir -p /Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/data/patient-data/YOUR_PATIENT_ID/

# Create subdirectories for each data modality
mkdir -p /Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/data/patient-data/YOUR_PATIENT_ID/{clinical,genomics,multiomics,spatial,imaging}
```

#### Step 3: Validate Data File Formats

Each data type has specific format requirements:

**Clinical Data (JSON):**
```json
{
  "patient_id": "YOUR_PATIENT_ID",
  "age": 58,
  "sex": "Female",
  "diagnosis": "Stage IV HGSOC",
  "treatment_history": [...],
  "biomarkers": {...}
}
```

**Genomic Data (VCF):**
```vcf
##fileformat=VCFv4.2
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
chr17   7577548 .       G       A       100     PASS    DP=150;AF=0.38
```

**Multiomics Data (CSV):**
```csv
Gene,Sample_01,Sample_02,Sample_03,...
TP53,145.2,98.7,234.5,...
PIK3CA,89.1,203.4,156.8,...
```

**Spatial Data (CSV):**
```csv
barcode,x_coord,y_coord,region
AAACAAGTATCTCCCA-1,100.5,200.3,tumor_core
AAACACCAATAACTGC-1,105.2,198.7,stroma
```

**Imaging Data (TIFF):**
- Supported formats: `.tiff`, `.tif`, `.png`
- Multi-channel TIFFs supported
- Recommended: 16-bit or 8-bit depth

#### Step 4: Update Configuration File

Edit each server's DRY_RUN setting to `"false"` as shown in [Actual Data Mode Setup](#1-configure-environment-variables).

#### Step 5: Test Individual Servers

Before running full PatientOne workflow, test each server:

```bash
# Test FGbio with real VCF
cd /Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-fgbio
FGBIO_DRY_RUN=false venv/bin/python -m mcp_fgbio

# Test MultiOmics with real data
cd /Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-multiomics
MULTIOMICS_DRY_RUN=false venv/bin/python -m mcp_multiomics

# etc. for all servers
```

#### Step 6: Run PatientOne with Real Data

Modify test prompts to reference your patient data:

```
I want to analyze patient [YOUR_PATIENT_ID] using the PatientOne workflow.

Please read the following files from /Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/data/patient-data/[YOUR_PATIENT_ID]/:
- clinical/patient_demographics.json
- genomics/somatic_variants.vcf
- multiomics/rna_expression.csv
- spatial/visium_coordinates.csv
- imaging/tumor_HE_20x.tiff

Then:
1. Use MockEpic to summarize clinical profile
2. Use FGbio to analyze genomic variants
3. Use MultiOmics to identify resistance signatures
4. Use SpatialTools to analyze tissue microenvironment
5. Use OpenImageData to quantify histology features

Synthesize findings and provide treatment recommendations.
```

---

## Data File Placement

### Required Directory Structure

```
/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/data/
├── patient-data/
│   └── [PATIENT_ID]/
│       ├── clinical/
│       │   ├── patient_demographics.json       # Required
│       │   └── lab_results.json               # Optional
│       ├── genomics/
│       │   ├── somatic_variants.vcf           # Required for genomic analysis
│       │   └── copy_number.seg                # Optional
│       ├── multiomics/
│       │   ├── pdx_rna_seq.csv                # Required for multiomics
│       │   ├── pdx_proteomics.csv             # Optional
│       │   ├── pdx_phosphoproteomics.csv      # Optional
│       │   └── sample_metadata.csv            # Required for batch correction
│       ├── spatial/
│       │   ├── visium_spatial_coordinates.csv # Required for spatial
│       │   ├── visium_gene_expression.csv     # Required
│       │   └── visium_region_annotations.csv  # Optional
│       └── imaging/
│           ├── tumor_HE_20x.tiff              # Required for histology
│           ├── tumor_IF_DAPI.tiff             # Optional (immunofluorescence)
│           ├── tumor_IF_CD3.tiff              # Optional
│           ├── tumor_IF_CD8.tiff              # Optional
│           ├── tumor_IF_KI67.tiff             # Optional
│           └── tumor_multiplex_IF.tiff        # Optional (multi-channel)
├── cache/                                      # Auto-created for intermediate files
│   ├── multiomics/
│   ├── images/
│   └── spatial/
└── output/                                     # Analysis results
    └── [PATIENT_ID]/
        ├── reports/
        ├── figures/
        └── processed_data/
```

### File Format Specifications

#### Clinical Data (JSON)

**Minimal example:**
```json
{
  "patient_id": "PAT002-OVC-2025",
  "age": 62,
  "sex": "Female",
  "diagnosis": "Stage IIIC Ovarian Carcinoma",
  "treatment_history": [
    {
      "date": "2024-01-15",
      "treatment": "Carboplatin + Paclitaxel",
      "response": "Complete Response"
    }
  ],
  "biomarkers": {
    "ca125_baseline": 890,
    "ca125_current": 45
  }
}
```

#### Genomic Data (VCF v4.2+)

Must include standard VCF headers and variant calls:
```vcf
##fileformat=VCFv4.2
##source=YourVariantCaller
#CHROM  POS       ID  REF  ALT  QUAL  FILTER  INFO                    FORMAT  TUMOR
chr17   7577548   .   G    A    100   PASS    DP=150;AF=0.38;GENE=TP53  GT:AD  0/1:92,58
chr3    178936091 .   G    A    95    PASS    DP=120;AF=0.45;GENE=PIK3CA GT:AD 0/1:66,54
```

#### Multiomics Data (CSV)

**RNA-seq expression (genes × samples):**
```csv
Gene,Sample_01,Sample_02,Sample_03,Sample_04,Sample_05
TP53,145.2,98.7,234.5,167.3,201.9
PIK3CA,89.1,203.4,156.8,198.2,175.6
AKT1,234.7,289.3,267.1,301.5,278.4
MTOR,156.3,178.9,165.2,187.4,171.8
```

**Sample metadata (required for batch correction):**
```csv
Sample,Treatment,Batch,Response
Sample_01,Sensitive,Batch1,Complete
Sample_02,Resistant,Batch1,Partial
Sample_03,Sensitive,Batch2,Complete
Sample_04,Resistant,Batch2,None
```

#### Spatial Data (CSV)

**Spatial coordinates:**
```csv
barcode,x_coord,y_coord
AAACAAGTATCTCCCA-1,100.5,200.3
AAACACCAATAACTGC-1,105.2,198.7
AAACAGAGCGACTCCT-1,110.8,195.4
```

**Gene expression (spots × genes):**
```csv
barcode,TP53,CD8A,CD3E,KI67,PCNA,COL1A1
AAACAAGTATCTCCCA-1,145,23,45,89,67,234
AAACACCAATAACTGC-1,89,156,178,34,45,456
```

**Region annotations (optional):**
```csv
barcode,region,annotation
AAACAAGTATCTCCCA-1,tumor_core,High proliferation
AAACACCAATAACTGC-1,stroma,Fibroblast enriched
```

#### Imaging Data (TIFF)

- **Format:** TIFF, multi-page TIFF, or PNG
- **Bit depth:** 8-bit or 16-bit
- **Channels:** Single-channel or multi-channel
- **Naming convention:** `[PATIENT_ID]_[tissue]_[stain]_[magnification].tiff`

**Examples:**
- `PAT002_tumor_HE_20x.tiff` - H&E histology at 20x
- `PAT002_tumor_IF_DAPI.tiff` - DAPI nuclear stain
- `PAT002_tumor_multiplex_IF_TP53_KI67_DAPI.tiff` - 3-channel multiplex

---

## Expected Outputs

### DRY_RUN Mode Outputs

**Response Structure:**
```json
{
  "status": "success (DRY_RUN mode)",
  "data": {
    // Synthetic but realistic data
  },
  "note": "This is a demonstration using synthetic data",
  "warnings": [
    "DRY_RUN mode enabled - no real data processing occurred"
  ]
}
```

**Example TEST_1 Output (DRY_RUN):**
```
✅ Clinical Data Retrieved (DRY_RUN)
Patient: PAT001-OVC-2025, 58yo Female, Stage IV HGSOC
CA-125 trajectory: 1456 → 22 → 389 U/mL (platinum-resistant)

✅ Genomic Analysis (DRY_RUN)
Key mutations: TP53 R175H (p.Arg175His), PIK3CA E545K
BRCA1 germline mutation detected
HRD score: 42 (positive)

✅ TCGA Comparison (DRY_RUN)
Molecular subtype: C1 Immunoreactive
Cohort similarity: 87%

⚠️  NOTE: All results are synthetic (DRY_RUN mode)
```

### Actual Data Mode Outputs

**Response Structure:**
```json
{
  "status": "success",
  "data": {
    // Real analysis results from your data
  },
  "files_processed": [
    "/path/to/patient_demographics.json",
    "/path/to/somatic_variants.vcf"
  ],
  "output_files": {
    "report": "/path/to/output/analysis_report.pdf",
    "figures": "/path/to/output/figures/",
    "processed_data": "/path/to/output/processed/"
  }
}
```

**Example TEST_2 Output (Actual Data):**
```
✅ Multiomics Data Loaded
Files: pdx_rna_seq.csv (1000 genes × 15 samples)
       pdx_proteomics.csv (500 proteins × 15 samples)
       sample_metadata.csv (15 samples, 2 batches)

✅ Preprocessing Applied
- Sample alignment: 15 common samples across modalities
- Batch correction: ComBat (PC1-batch correlation 0.78 → 0.15)
- Normalization: Quantile normalization
- Missing values: Median imputation (2.3% missing)

✅ Differential Analysis (Resistant vs Sensitive)
Significant genes (FDR < 0.05): 127
Top pathway: PI3K/AKT/mTOR (p = 1.2e-8)
Upregulated: PIK3CA (log2FC=2.8), AKT1 (log2FC=2.3)

✅ Outputs Saved
- Heatmap: /path/to/output/resistance_heatmap.png
- Gene list: /path/to/output/deg_list.csv
- Full report: /path/to/output/multiomics_report.pdf

✅ Real data analysis complete (MULTIOMICS_DRY_RUN=false)
```

### Performance Comparison

| Test | DRY_RUN Time | Actual Data Time | Data Size |
|------|--------------|------------------|-----------|
| TEST_1 (Clinical + Genomic) | 2-3 min | 5-10 min | ~2 MB |
| TEST_2 (Multiomics) | 3-5 min | 10-20 min | ~500 MB |
| TEST_3 (Spatial) | 3-5 min | 15-30 min | ~600 MB |
| TEST_4 (Imaging) | 2-4 min | 20-40 min | ~2 GB |
| TEST_5 (Integration) | 2-3 min | 5-10 min | N/A (synthesis) |
| **TOTAL** | **15-20 min** | **60-120 min** | **~3 GB** |

---

## Troubleshooting

### DRY_RUN Mode Issues

#### Issue: "Still seeing DRY_RUN warnings after disabling"

**Cause:** Configuration not reloaded or incorrect syntax

**Solution:**
```bash
# 1. Verify configuration
cat ~/Library/Application\ Support/Claude/claude_desktop_config.json | grep DRY_RUN
# Should show "false" for all servers

# 2. Check JSON syntax
python3 -m json.tool ~/Library/Application\ Support/Claude/claude_desktop_config.json
# Should not show errors

# 3. Restart Claude Desktop (fully quit, not just close window)
# macOS: Cmd+Q then reopen
```

#### Issue: "Servers returning mock data even with DRY_RUN=false"

**Cause:** Some servers have additional internal checks

**Solution:**
```bash
# Check server logs
tail -f ~/Library/Logs/Claude/mcp*.log

# Verify environment variables are being passed
# Look for "DRY_RUN: false" in server startup logs
```

### Actual Data Mode Issues

#### Issue: "FileNotFoundError: No such file"

**Cause:** Incorrect file paths or missing data files

**Solution:**
```bash
# 1. Verify absolute paths
ls -la /Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/data/patient-data/[PATIENT_ID]/

# 2. Check file permissions
chmod -R 644 /Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/data/patient-data/

# 3. Update paths in prompt to match your system
pwd  # Check current directory
```

#### Issue: "API Rate Limit Exceeded (TCGA)"

**Cause:** Too many requests to TCGA API

**Solution:**
```bash
# Add delays between TCGA queries
# Or use local TCGA data downloads:
# https://portal.gdc.cancer.gov/
```

#### Issue: "HuggingFace Authentication Failed"

**Cause:** Missing or invalid HF_TOKEN

**Solution:**
```bash
# 1. Get token from https://huggingface.co/settings/tokens
# 2. Add to configuration:
{
  "mcp-huggingface": {
    "env": {
      "HF_TOKEN": "hf_your_actual_token_here",
      "HF_DRY_RUN": "false"
    }
  }
}

# 3. Restart Claude Desktop
```

#### Issue: "Out of Memory during image processing"

**Cause:** Large TIFF images exceeding available RAM

**Solution:**
```bash
# 1. Process images in smaller tiles
# 2. Reduce image resolution
# 3. Increase system RAM allocation
# 4. Use image pyramids for large files
```

#### Issue: "Batch correction failed - singular matrix"

**Cause:** Too few samples per batch or perfect confounding

**Solution:**
```bash
# 1. Check sample distribution across batches
# Batch 1: ≥3 samples recommended
# Batch 2: ≥3 samples recommended

# 2. If <3 samples per batch, skip batch correction:
# In prompt: "skip batch correction due to small sample size"

# 3. Or combine batches if appropriate
```

### Data Validation Issues

#### Issue: "VCF format error"

**Solution:**
```bash
# Validate VCF file
vcf-validator your_variants.vcf

# Common fixes:
# - Ensure ##fileformat=VCFv4.2 header
# - Check for missing columns
# - Verify chromosome naming (chr1 vs 1)
```

#### Issue: "CSV encoding error"

**Solution:**
```bash
# Convert to UTF-8
iconv -f ISO-8859-1 -t UTF-8 input.csv > output.csv

# Remove BOM if present
sed '1s/^\xEF\xBB\xBF//' input.csv > output.csv
```

#### Issue: "TIFF not recognized"

**Solution:**
```bash
# Verify TIFF format
file your_image.tiff
# Should show: TIFF image data

# Convert if needed
convert input.png output.tiff

# For multi-channel TIFFs, ensure proper channel ordering
```

---

## Use Case Recommendations

### Choosing the Right Mode

#### Use DRY_RUN Mode For:

1. **First-time Setup**
   - Learning the PatientOne workflow
   - Verifying MCP server installation
   - Testing Claude Desktop integration

2. **Demonstrations**
   - Conference presentations
   - Educational workshops
   - Stakeholder demos

3. **Development**
   - Testing new features
   - CI/CD pipelines
   - Integration testing
   - Debugging server configurations

4. **Quick Checks**
   - Verifying tool orchestration
   - Testing prompt engineering
   - Exploring available tools

#### Use Actual Data Mode For:

1. **Research Analysis**
   - Processing your own datasets
   - Validating PatientOne methodology
   - Comparative studies
   - Publications requiring real data

2. **Clinical Applications** (with validation)
   - Tumor board preparation
   - Treatment decision support
   - Biomarker discovery
   - Clinical trial matching

3. **Production Workflows**
   - Routine precision medicine analysis
   - High-throughput processing
   - Integration with institutional pipelines
   - Real-time analysis

4. **Validation Studies**
   - Benchmarking against gold standards
   - Cross-platform comparison
   - Method development

### Hybrid Approach

**Recommended workflow for new analyses:**

1. **Start with DRY_RUN** (5-10 min)
   - Test prompts and workflow
   - Verify tool orchestration
   - Identify potential issues

2. **Switch to Single Modality Actual Data** (30-60 min)
   - Test one data type at a time
   - Validate file formats
   - Check outputs

3. **Full Actual Data Analysis** (1-2 hours)
   - Run complete PatientOne workflow
   - Generate publication-quality results
   - Create reports for stakeholders

4. **Return to DRY_RUN for iterations** (5-10 min)
   - Test modifications quickly
   - Experiment with prompts
   - Validate changes

---

## Data Privacy and Compliance

### Important Considerations for Actual Data Mode

⚠️ **WARNING:** When using Actual Data mode with real patient information:

1. **HIPAA Compliance (US)**
   - Ensure proper de-identification (HIPAA Safe Harbor or Expert Determination)
   - Use institutional review board (IRB) approval for research
   - Implement access controls and audit logs
   - Never include Protected Health Information (PHI) in prompts

2. **GDPR Compliance (EU)**
   - Obtain proper consent for data processing
   - Implement data minimization
   - Ensure right to erasure capabilities
   - Document data processing activities

3. **Local Data Processing**
   - All MCP servers run locally (no cloud upload by default)
   - Claude Desktop processes data on your machine
   - External API calls only for TCGA, HuggingFace (can be disabled)

4. **Recommendations:**
   - Use synthetic or de-identified data for demos
   - Store patient data on encrypted drives
   - Implement access controls (file permissions)
   - Maintain audit trails of analyses
   - Follow institutional data governance policies

---

## Additional Resources

### Documentation
- [PatientOne Architecture →](../../../architecture/patient-one/README.md)
- [PatientOne Quick Start →](./README.md)
- [MCP Server Configuration →](../../../configs/README.md)
- [Test Coverage Report →](../../README.md)

### External Links
- [TCGA API Documentation](https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/)
- [HuggingFace API](https://huggingface.co/docs/api-inference/index)
- [Seqera Platform](https://seqera.io/platform/)
- [DeepCell Documentation](https://deepcell.readthedocs.io/)

### Support
- **Issues:** https://github.com/lynnlangit/precision-medicine-mcp/issues
- **Discussions:** https://github.com/lynnlangit/precision-medicine-mcp/discussions

---

**Last Updated:** December 27, 2025
**Version:** 1.0
**Covers:** DRY_RUN and Actual Data modes for all 9 MCP servers
