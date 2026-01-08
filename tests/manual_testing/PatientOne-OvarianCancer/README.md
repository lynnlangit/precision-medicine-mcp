# PatientOne: Quick Start Guide

Comprehensive precision medicine workflow for Stage IV Ovarian Cancer using 9 MCP servers

## Overview

PatientOne demonstrates end-to-end precision medicine analysis integrating:
- **Clinical** data (demographics, CA-125 trends)
- **Genomic** variants (VCF, CNVs, TCGA comparison)
- **Multiomics** (RNA/Protein/Phospho from PDX models)
- **Spatial** transcriptomics (10x Visium, 900 spots)
- **Imaging** (H&E histology, multiplex IF)

**All synthetic data** for demonstration purposes.


---

## ⚠️ IMPORTANT: Research Use Only Disclaimer

```
╔═══════════════════════════════════════════════════════════════════════════╗
║                     PATIENTONE WORKFLOW DISCLAIMER                         ║
╚═══════════════════════════════════════════════════════════════════════════╝

⚠️  RESEARCH USE ONLY - NOT FOR CLINICAL DECISION-MAKING ⚠️

This PatientOne workflow has NOT been clinically validated.
Do NOT use for actual patient care decisions.

CRITICAL LIMITATIONS:

1. SYNTHETIC DATA
   This demo uses synthetic/de-identified data for demonstration.
   Results are for educational purposes only.

2. AI-GENERATED INSIGHTS
   Claude AI orchestrates the analysis but:
   • May misinterpret complex data
   • May hallucinate connections not supported by data
   • Cannot replace expert human interpretation
   • Requires validation by qualified bioinformaticians

3. VALIDATION REQUIRED
   ALL findings must be validated before research use:
   • Statistical findings → Independent cohort
   • Pathway predictions → Mechanistic studies
   • Drug targets → Experimental models
   • Clinical recommendations → Clinical trials

4. NO LIABILITY
   The developers and AI provider assume NO liability for:
   • Treatment decisions based on this analysis
   • Patient outcomes
   • Misinterpretation of results
   • Technical errors or bugs

╔═══════════════════════════════════════════════════════════════════════════╗
║  CONSULT A QUALIFIED ONCOLOGIST BEFORE ANY CLINICAL DECISION               ║
╚═══════════════════════════════════════════════════════════════════════════╝
```

**See [DISCLAIMERS.md](../../../docs/DISCLAIMERS.md) for complete safety guidelines.**

---

## Prerequisites

### System Requirements
- **Python:** 3.11+
- **Claude Desktop:** Latest version ([Download](https://claude.ai/download))
- **RAM:** 16GB recommended
- **Disk:** 50GB free space
- **OS:** macOS, Linux, or Windows with WSL2

### Setup Verification

1. **Check Python version:**
```bash
python3 --version  # Should show 3.11 or higher
```

2. **Verify Claude Desktop configuration:**
```bash
cat ~/Library/Application\ Support/Claude/claude_desktop_config.json
# Should show all 9 MCP servers configured
```

3. **Confirm data files exist:**
```bash
ls -lh ../../data/patient-data/PAT001-OVC-2025/
# Should show 17 files (~3.2 MB total)
```

### First-Time Setup

If you haven't installed the MCP servers yet:

```bash
# Clone repository
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp

# Install dependencies (5-10 min)
cd manual_testing
./install_dependencies.sh

# Configure Claude Desktop
cp ../configs/claude_desktop_config.json ~/Library/Application\ Support/Claude/claude_desktop_config.json

# Restart Claude Desktop
# Verify servers loaded (should see 9 servers in Claude Desktop)

# Test basic server connectivity
./verify_servers.sh
```

---

## Running Modes: DRY_RUN vs Actual Data

PatientOne can run in two modes:

| Mode | Purpose | Data Source | External APIs | Best For |
|------|---------|-------------|---------------|----------|
| **DRY_RUN** (default) | Demo & testing | Synthetic responses | None | Quick demo, CI/CD, learning |
| **Actual Data** | Real analysis | Your files | May connect | Production, research, clinical |

**Quick Mode Selection:**
- **DRY_RUN mode** (default): No setup needed, works immediately with synthetic data
- **Actual Data mode**: Requires data files and configuration — see [📘 Data Modes Guide](./DATA_MODES_GUIDE.md)

💡 **Tip:** Start with DRY_RUN mode to understand the workflow (5 min), then switch to actual data for real analysis.

---

## Try PatientOne in 5 Minutes

### Option 1: Quick Demo (Single Test)

Run **TEST_1** to see clinical + genomic integration:

1. **Open Claude Desktop**

2. **Copy/paste this prompt:**
```
I want to run the PatientOne clinical and genomic analysis (TEST_1).

Please read the following files:
- /Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/clinical/patient_demographics.json
- /Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/clinical/lab_results.json
- /Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/genomics/somatic_variants.vcf

Then:
1. Use mockepic server to summarize patient clinical profile (demographics, BRCA1 status)
2. Use mockepic to analyze CA-125 tumor marker trends (evidence of platinum resistance?)
3. Use fgbio server to parse somatic variants (expect TP53 R175H, PIK3CA E545K, PTEN LOH)
4. Use tcga server to compare mutations to TCGA ovarian cancer cohort and identify molecular subtype
5. Synthesize findings into a clinical summary

This is TEST_1 from the PatientOne workflow.
```

3. **Expected output:**
- Patient demographics (Sarah Anderson, 58yo, Stage IV HGSOC)
- CA-125 trajectory showing initial response then resistance
- Key mutations: TP53 R175H, PIK3CA E545K, PTEN LOH
- TCGA subtype: C1 Immunoreactive
- BRCA1 germline mutation implications

**Duration:** 5-10 minutes

---

### Option 2: Complete Analysis (All 5 Tests)

Run all modular tests sequentially for comprehensive precision medicine analysis:

#### Demonstration Data (Small Synthetic Files)

| Test | Focus | DRY_RUN Time | Real Data Time | DRY_RUN Cost | Real Data Cost (Small Files) |
|------|-------|--------------|----------------|--------------|------------------------------|
| **TEST_1** | Clinical + Genomic | 3-5 min | 10-15 min | ~$1 | ~$1 |
| **TEST_2** | Multi-Omics Resistance | 5-8 min | 15-25 min | ~$1 | $2-4 |
| **TEST_3** | Spatial Transcriptomics | 4-6 min | 45-90 min | ~$1 | $8-17 |
| **TEST_4** | Histology & Imaging | 3-5 min | 20-40 min | ~$1 | $3-7 |
| **TEST_5** | Integration & Recommendations | 5-7 min | 5-10 min | ~$1 | ~$1 |
| **TOTAL** | Complete Analysis | **25-35 min** | **2-4 hours** | **~$1** | **$15-45** |

**Data sizes:** 315 KB spatial, 505 KB multi-omics, 4.1 MB imaging (total ~4.9 MB)

#### Production Data (Realistic Hospital Volumes)

| Test | Focus | Real Data Time (Pre-aligned) | Real Data Time (Raw FASTQ) | Cost (Pre-aligned) | Cost (Raw FASTQ) |
|------|-------|------------------------------|----------------------------|-------------------|------------------|
| **TEST_1** | Clinical + Genomic | 15-30 min | 15-30 min | ~$1-2 | ~$1-2 |
| **TEST_2** | Multi-Omics Resistance | 30-60 min | 30-60 min | $6-20 | $6-20 |
| **TEST_3** | Spatial Transcriptomics | 45-120 min | 90-240 min | $5-30 | $10-40 |
| **TEST_4** | Histology & Imaging | 40-90 min | 40-90 min | $10-36 | $10-36 |
| **TEST_5** | Integration & Recommendations | 10-20 min | 10-20 min | ~$1-2 | ~$1-2 |
| **TOTAL (Compute + API)** | Complete Analysis | **2-4 hours** | **4-8 hours** | **$23-90** | **$28-100** |
| **+ Claude Tokens** | - | - | - | **~$1-2** | **~$1-2** |
| **GRAND TOTAL** | - | **2-4 hours** | **4-8 hours** | **$24-92** | **$29-102** |

**Data sizes:** 100-500 MB spatial (3,000-5,000 spots × 18,000-30,000 genes), 15-20 MB multi-omics processed (or 2.7 GB raw), 500 MB - 2 GB imaging

**Cost Breakdown:**
- **DRY_RUN:** Claude token usage only (~30K tokens) - ~$1
- **Real Data (Demonstration):** Compute ($7-24) + APIs (~$1) + Claude tokens (~$1) = **$8-26 total**
- **Real Data (Production):** Compute ($22-99) + APIs (~$1) + Claude tokens (~$1-2) = **$24-102 total**
  - 3-4× more expensive due to 300-1500× larger data files
  - **Token costs stay low** (~$1-2) because MCP servers return summaries, not raw 3-8 GB files!

📊 **[Full Cost Analysis & ROI →](../../../docs/operations/COST_ANALYSIS.md)**

**Instructions:**
1. Open each `TEST_*.txt` file in `implementation/` directory
2. Copy/paste the prompt into Claude Desktop
3. Review results before proceeding to next test
4. Tests build on each other but are independently runnable

---

## Imaging Modality Reference

Understanding the difference between imaging types is critical for correct analysis:

| Image Type | Microscopy Mode | Staining Method | Analysis Server(s) | Use Case |
|------------|----------------|-----------------|-------------------|----------|
| **H&E** | Brightfield | Chromogenic (Hematoxylin=blue nuclei, Eosin=pink cytoplasm) | OpenImageData | Tissue architecture, morphology, cellularity assessment |
| **IF (single-plex)** | Fluorescence | Single fluorescent antibody | OpenImageData + DeepCell | Protein marker quantification (CD8, Ki67, etc.) |
| **IF (multiplex)** | Fluorescence | Multiple fluorophores (2-7 colors) | OpenImageData + DeepCell | Cell phenotyping, protein co-localization |
| **Spatial RNA-seq** | N/A (sequencing) | Tabular CSV data (no images) | SpatialTools only | Gene expression patterns across tissue |

**Key Differences:**
- **H&E:** Brightfield microscopy with colored (chromogenic) stains - NOT fluorescence
- **IF:** Fluorescence microscopy with fluorescent antibodies - requires different analysis
- **Spatial data:** No images, just CSV files with coordinates and expression values

**When to use DeepCell:**
- ✅ IF images requiring cell segmentation
- ✅ H&E images requiring nuclear/cell segmentation
- ❌ Tabular spatial data (CSV files) - use SpatialTools instead

---

## Test Descriptions

### TEST_1: Clinical + Genomic Analysis
**Servers:** Epic, FGbio, TCGA
**Files:** 3 (patient_demographics.json, lab_results.json, somatic_variants.vcf)

**What it does:**
- Retrieves patient demographics and treatment history
- Analyzes CA-125 tumor marker trajectory
- Identifies somatic mutations and CNVs
- Compares to TCGA ovarian cancer cohort
- Determines molecular subtype

**Key Findings:**
- Platinum-resistant disease (8-month recurrence)
- TP53/PIK3CA/PTEN driver mutations
- C1 immunoreactive subtype
- BRCA1 germline mutation → HRD-positive

---

### TEST_2: Multi-Omics Resistance Analysis
**Servers:** MultiOmics
**Files:** 4 (pdx_rna_seq.csv, pdx_proteomics.csv, pdx_phosphoproteomics.csv, sample_metadata.csv)

**What it does:**
- Integrates RNA/Protein/Phospho data from 15 PDX samples
- Compares resistant vs sensitive samples (7 vs 8)
- Performs Stouffer's meta-analysis with FDR correction
- Identifies dysregulated pathways

**Key Findings:**
- PI3K/AKT/mTOR pathway activation in resistant samples
- PIK3CA, AKT1, mTOR, RPS6KB1 upregulated (p < 0.001)
- Drug efflux: ABCB1 (MDR1) overexpression
- Anti-apoptotic: BCL2L1 upregulation

---

### TEST_3: Spatial Transcriptomics
**Servers:** SpatialTools
**Files:** 3 (visium_gene_expression.csv, visium_spatial_coordinates.csv, visium_region_annotations.csv)

**What it does:**
- Processes 10x Visium spatial RNA-seq **tabular data** (900 spots, 31 genes)
- Identifies 6 tissue regions (tumor_core, proliferative, interface, stroma, etc.)
- Maps spatial expression patterns from CSV files
- Quantifies immune cell distribution
- **Generates visualizations:** Spatial heatmaps, gene expression matrices, autocorrelation plots

**Note:** Uses only tabular CSV data, not images. DeepCell is NOT needed for this test.

**Key Findings:**
- Immune exclusion phenotype (CD8+ low in tumor core)
- High proliferation in tumor_proliferative region (Ki67+, PCNA+)
- Thick stromal barrier separating immune cells from tumor
- Spatial heterogeneity in resistance markers

**Expected Visualizations:**
- Spatial heatmap showing top 6 spatially variable genes across tissue coordinates
- Gene expression heatmap (8 key genes × 6 regions)
- Region composition bar chart
- Spatial autocorrelation plot (Moran's I)

---

### TEST_4: Histology & Imaging
**Servers:** OpenImageData, DeepCell
**Files:** 7 TIFF images (H&E brightfield, IF single-markers, multiplex IF)

**What it does:**
- Analyzes **H&E histology** using brightfield microscopy (tissue architecture, morphology)
- Processes **immunofluorescence (IF) images** using fluorescence microscopy (DAPI, CD3, CD8, Ki67, PanCK)
- Performs cell segmentation on IF images with DeepCell
- Quantifies proliferation and immune infiltration
- **Generates visualizations:** Segmentation overlays, spatial distribution maps, phenotype analyses

**Note:** H&E uses chromogenic stains (brightfield), while IF uses fluorescent antibodies (fluorescence). These require different imaging modalities.

**Key Findings:**
- Tumor cellularity: 70-80%
- Ki67 proliferation index: 45-55% (HIGH)
- CD8+ T cell density: 5-15 cells/mm² (LOW, mostly peripheral)
- CD3+ overall: 30-50 cells/mm² (moderate T cells, but not cytotoxic)

**Expected Visualizations:**
- **CD8 IF:** Segmentation overlay (CD8+ vs CD8-), spatial distribution heatmap
- **Ki67 IF:** Nuclear segmentation overlay (Ki67+ vs Ki67-), proliferation heatmap
- **Multiplex IF:** RGB channel composite, cell phenotype segmentation, scatter plot (TP53 vs KI67)
- **H&E:** Annotated morphology showing necrotic regions and high cellularity areas

---

### TEST_5: Integration & Recommendations
**Servers:** All 9 servers (synthesis)
**Files:** None (builds on previous tests)

**What it does:**
- Synthesizes findings across all 5 modalities
- Integrates molecular, spatial, and clinical insights
- Identifies actionable treatment targets
- Generates precision medicine recommendations

**Key Recommendations:**
- **Primary:** PI3K inhibitor (Alpelisib) targeting PIK3CA E545K mutation
- **Secondary:** Anti-PD-1 immunotherapy to overcome immune exclusion
- **Tertiary:** PARP inhibitor consideration (BRCA1 mutation, but PIK3CA pathway may limit efficacy)
- **Clinical trial:** NCT03602859 (alpelisib + paclitaxel in ovarian cancer)

---

## Expected Outputs

### For Each Test

Claude Desktop will generate:
1. **Data Summary:** Key statistics from loaded files
2. **Tool Execution:** MCP server calls with results
3. **Analysis:** Interpretation and synthesis
4. **Findings:** Bullet-point key discoveries

### Final Integrated Output (After TEST_5)

Comprehensive report including:
- **Executive Summary:** Patient profile and precision medicine strategy
- **Molecular Profile:** Genomic alterations, pathway dysregulation
- **Microenvironment:** Spatial distribution, immune landscape
- **Resistance Mechanisms:** Multi-omics signatures
- **Treatment Plan:** Evidence-based recommendations with rationale

---

## Troubleshooting

### Issue: "MCP servers not found"

**Cause:** Claude Desktop config not loaded or servers not installed

**Fix:**
```bash
# Verify config exists
cat ~/Library/Application\ Support/Claude/claude_desktop_config.json

# If missing, copy template
cp configs/claude_desktop_config.json ~/Library/Application\ Support/Claude/

# Restart Claude Desktop
```

---

### Issue: "Cannot find data files"

**Cause:** Incorrect file paths or data not present

**Fix:**
```bash
# Verify data exists
ls -lh data/patient-data/PAT001-OVC-2025/
# Should show 17 files

# Check absolute path in prompt matches your system
pwd  # Note current directory
# Update file paths in prompts to match your installation
```

---

### Issue: "Server returned error / DRY_RUN warnings"

**Cause:** Servers are in DRY_RUN mode (expected behavior for testing)

**Explanation:**
- All MCP servers are configured with `DRY_RUN=true` by default
- This prevents actual external API calls while demonstrating tool orchestration
- Servers return realistic synthetic responses

**To use your own data:**
- See the [📘 Data Modes Guide](./DATA_MODES_GUIDE.md) for complete instructions on:
  - Switching from DRY_RUN to Actual Data mode
  - Configuring environment variables
  - Setting up data file directories
  - Obtaining required API keys

---

### Issue: "Context limit exceeded"

**Cause:** Trying to run all tests in single prompt

**Fix:**
- Run tests individually (TEST_1 through TEST_5)
- Each test designed to fit within Claude Desktop context limits
- Do NOT combine multiple tests in one prompt
- Clear conversation history between tests if needed

---

### Issue: "Missing Python packages"

**Cause:** Virtual environments not set up correctly

**Fix:**
```bash
cd manual_testing
./install_dependencies.sh

# Verify each server's venv
for server in ../servers/mcp-*/; do
    echo "Checking $server"
    $server/venv/bin/python --version
done
```

---

## File Access Configuration

For detailed guidance on configuring Claude Desktop for file access, see:
- `implementation/CLAUDE_DESKTOP_FILE_ACCESS_GUIDE.md`
- `implementation/CLAUDE_DESKTOP_FIX_SUMMARY.md`

---

## Next Steps

After completing PatientOne:

1. **Explore Individual Workflows:**
   - [Spatial Transcriptomics](../../../architecture/spatial/README.md)
   - [Multiomics Integration](../../../architecture/multiomics/README.md)

2. **Customize for Your Data:**
   - Replace synthetic patient data with your own
   - Adjust MCP server configurations
   - Modify test prompts to match your analysis goals

3. **Production Deployment:**
   - Integrate with real EHR systems
   - Connect to institutional genomic databases
   - Deploy Nextflow pipelines via Seqera Platform

4. **Read Full Architecture:**
   - [PatientOne Architecture](../../../architecture/patient-one/README.md)
   - [Comprehensive Documentation](../../../README.md)

---

## Support

**Issues:** https://github.com/lynnlangit/precision-medicine-mcp/issues
**Documentation:** https://github.com/lynnlangit/precision-medicine-mcp

---

**Last Updated:** December 29, 2025
**Testing Status:** All 5 tests validated with synthetic data
**Data:** 100% synthetic for demonstration purposes
