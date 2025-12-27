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

<kbd><img src="https://github.com/lynnlangit/spatial-mcp/blob/main/architecture/patient-one/patient-one-holistic.png" width=600></kbd>

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

## Try PatientOne in 5 Minutes

### Option 1: Quick Demo (Single Test)

Run **TEST_1** to see clinical + genomic integration:

1. **Open Claude Desktop**

2. **Copy/paste this prompt:**
```
I want to run the PatientOne clinical and genomic analysis.

Please read the following files from /Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/:
- clinical_demographics.json
- ca125_timeline.csv
- PAT001_somatic_variants.vcf

Then:
1. Use MockEpic to summarize patient clinical profile
2. Use FGbio to analyze somatic variants
3. Use TCGA to compare mutations to ovarian cancer cohort
4. Synthesize findings and identify molecular subtype

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

| Test | Focus | Time | Prompt File |
|------|-------|------|-------------|
| **TEST_1** | Clinical + Genomic | 5-10 min | `implementation/TEST_1_CLINICAL_GENOMIC.txt` |
| **TEST_2** | Multi-Omics Resistance | 5-10 min | `implementation/TEST_2_MULTIOMICS.txt` |
| **TEST_3** | Spatial Transcriptomics | 5-10 min | `implementation/TEST_3_SPATIAL.txt` |
| **TEST_4** | Histology & Imaging | 5-10 min | `implementation/TEST_4_IMAGING.txt` |
| **TEST_5** | Integration & Recommendations | 5 min | `implementation/TEST_5_INTEGRATION.txt` |

**Total Duration:** 30-45 minutes

**Instructions:**
1. Open each `TEST_*.txt` file in `implementation/` directory
2. Copy/paste the prompt into Claude Desktop
3. Review results before proceeding to next test
4. Tests build on each other but are independently runnable

---

## Test Descriptions

### TEST_1: Clinical + Genomic Analysis
**Servers:** MockEpic, FGbio, TCGA
**Files:** 3 (clinical_demographics.json, ca125_timeline.csv, PAT001_somatic_variants.vcf)

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
**Files:** 4 (pdx_rna_expression.csv, pdx_protein.csv, pdx_phospho.csv, pdx_metadata.csv)

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
**Servers:** SpatialTools, DeepCell
**Files:** 3 (visium_spots.csv, visium_expression.csv, visium_regions.csv)

**What it does:**
- Processes 10x Visium spatial RNA-seq (900 spots, 31 genes)
- Identifies 6 tissue regions (tumor_core, proliferative, interface, stroma, etc.)
- Maps spatial expression patterns
- Quantifies immune cell distribution

**Key Findings:**
- Immune exclusion phenotype (CD8+ low in tumor core)
- High proliferation in tumor_proliferative region (Ki67+, PCNA+)
- Thick stromal barrier separating immune cells from tumor
- Spatial heterogeneity in resistance markers

---

### TEST_4: Histology & Imaging
**Servers:** OpenImageData, DeepCell
**Files:** 7 TIFF images (H&E, IF single-markers, multiplex)

**What it does:**
- Analyzes H&E histology (tissue architecture)
- Processes immunofluorescence (DAPI, CD3, CD8, Ki67, PanCK)
- Performs cell segmentation and phenotyping
- Quantifies proliferation and immune infiltration

**Key Findings:**
- Tumor cellularity: 70-80%
- Ki67 proliferation index: 45-55% (HIGH)
- CD8+ T cell density: 5-15 cells/mm² (LOW, mostly peripheral)
- CD3+ overall: 30-50 cells/mm² (moderate T cells, but not cytotoxic)

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
- To disable DRY_RUN, edit `configs/claude_desktop_config.json` and change `*_DRY_RUN` to `false`

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
   - [Spatial Transcriptomics](../../architecture/spatial/README.md)
   - [Multiomics Integration](../../architecture/multiomics/README.md)

2. **Customize for Your Data:**
   - Replace synthetic patient data with your own
   - Adjust MCP server configurations
   - Modify test prompts to match your analysis goals

3. **Production Deployment:**
   - Integrate with real EHR systems
   - Connect to institutional genomic databases
   - Deploy Nextflow pipelines via Seqera Platform

4. **Read Full Architecture:**
   - [PatientOne Architecture](../../architecture/patient-one/README.md)
   - [Comprehensive Documentation](../../README.md)

---

## Support

**Issues:** https://github.com/lynnlangit/precision-medicine-mcp/issues
**Documentation:** https://github.com/lynnlangit/precision-medicine-mcp

---

**Last Updated:** December 26, 2025
**Testing Status:** All 5 tests validated with synthetic data
**Data:** 100% synthetic for demonstration purposes
