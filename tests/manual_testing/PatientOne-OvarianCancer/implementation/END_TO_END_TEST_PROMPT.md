# End-to-End Test Prompt: Patient One - Ovarian Cancer

**Purpose:** Comprehensive test of all 9 MCP servers using realistic synthetic patient data

**Patient ID:** PAT001-OVC-2025

---

## ⚠️ IMPORTANT: Use the FIXED Test Prompts

The original comprehensive end-to-end test hit **context limits** in Claude Desktop and had **file access issues** due to containerization.

**✅ Solution:** Use the 5 focused FIXED test prompts instead.

---

## Recommended Testing Approach

### Use These Test Files:

| Test | File | Purpose | Time |
|------|------|---------|------|
| **1** | `TEST_1_CLINICAL_GENOMIC.txt` | Clinical data + genomic analysis | 5-10 min |
| **2** | `TEST_2_MULTIOMICS.txt` | Multi-omics resistance analysis | 5-10 min |
| **3** | `TEST_3_SPATIAL.txt` | Spatial transcriptomics | 5-10 min |
| **4** | `TEST_4_IMAGING.txt` | Histology and imaging | 5-10 min |
| **5** | `TEST_5_INTEGRATION.txt` | Integrated clinical report | 5 min |

**Total time:** 25-45 minutes (run separately in Claude Desktop)

---

## Why FIXED Tests?

### Problem 1: Context Limits ❌
The comprehensive test tried to load all data at once:
- 1000 genes from RNA-seq
- 31 genes from spatial transcriptomics
- 7 imaging files
- Result: **Conversation size limit exceeded**

### Problem 2: File Access ❌
Claude Desktop's MCP servers run in containerized environments and cannot access arbitrary file paths like:
```
/Users/lynnlangit/Documents/GitHub/spatial-mcp/manual_testing/PatientOne-OvarianCancer/...
```

### Solution: FIXED Tests ✅
- **Focused scope:** Each test analyzes 3-8 key features (not hundreds)
- **MCP-accessible paths:** Uses relative paths like `patient-data/PAT001-OVC-2025/`
- **Data copied:** All 17 files copied to `/data/patient-data/PAT001-OVC-2025/`
- **Explicit instructions:** Each test instructs Claude Desktop to use MCP server tools

---

## Quick Start

### Step 1: Read Quick Reference
```bash
cat QUICK_TEST_REFERENCE.md
```

### Step 2: Run Test 1
```bash
cat TEST_1_CLINICAL_GENOMIC.txt
```
Copy the output → Paste into Claude Desktop

### Step 3: Run Tests 2-5
Use a **new conversation** for each test:
```bash
cat TEST_2_MULTIOMICS.txt
cat TEST_3_SPATIAL.txt
cat TEST_4_IMAGING.txt
cat TEST_5_INTEGRATION.txt
```

---

## What Gets Tested?

### TEST 1: Clinical & Genomic ✅
**Servers:** mcp-mockepic, mcp-fgbio, mcp-tcga
- Patient demographics (Sarah Anderson, 58yo, BRCA1 mutation)
- CA-125 tumor marker trends (1456→22→389→289 U/mL)
- Somatic mutations (TP53 R175H, PIK3CA E545K, PTEN LOH)
- TCGA cohort comparison

### TEST 2: Multi-Omics ✅
**Servers:** mcp-multiomics
- 15 PDX samples (7 resistant, 8 sensitive)
- 6 key genes across RNA/Protein/Phospho modalities
- Stouffer's meta-analysis with directionality
- PI3K/AKT pathway activation analysis

### TEST 3: Spatial Transcriptomics ✅
**Servers:** mcp-spatialtools
- 900 spatial spots across 6 tissue regions
- 8 key genes (proliferation, resistance, immune markers)
- Spatial heterogeneity in resistance markers
- Immune cell exclusion patterns

### TEST 4: Imaging ✅
**Servers:** mcp-openimagedata, mcp-deepcell
- H&E histology (tissue architecture, necrosis, cellularity)
- CD8 T cell infiltration quantification
- Ki67 proliferation index
- Multiplex IF cell segmentation (TP53/Ki67/DAPI)

### TEST 5: Integration ✅
**Servers:** None (synthesis only)
- Resistance mechanisms ranked by evidence strength
- Multi-modal consistency analysis
- Top 3 therapeutic recommendations
- Immunotherapy assessment
- Clinical monitoring strategy

---

## Documentation

### For Quick Testing:
- **`QUICK_TEST_REFERENCE.md`** - One-page quick start guide

### For Detailed Testing:
- **`TESTING_STRATEGY.md`** - Comprehensive testing strategy

### For Troubleshooting:
- **`CLAUDE_DESKTOP_FILE_ACCESS_GUIDE.md`** - File access issue explanation
- **`CLAUDE_DESKTOP_FIX_SUMMARY.md`** - Complete solution summary

---

## Verification Checkpoints

After running all 5 tests, you should have confirmed:

**Clinical & Genomic:**
- ✅ Patient: Sarah Anderson, 58yo
- ✅ BRCA1 germline mutation identified
- ✅ CA-125 resistance pattern (1456→22→389→289)
- ✅ TP53 R175H, PIK3CA E545K, PTEN LOH mutations found

**Multi-Omics:**
- ✅ 15 PDX samples analyzed (7R, 8S)
- ✅ PI3K/AKT pathway activated
- ✅ Stouffer's Z-scores >3 for top genes
- ✅ FDR-adjusted p-values <0.05

**Spatial:**
- ✅ 900 spots across 6 regions loaded
- ✅ Resistance markers concentrated in tumor regions
- ✅ Immune cells excluded from tumor
- ✅ Tumor microenvironment: COLD

**Imaging:**
- ✅ High cellularity (~70-80%)
- ✅ Low CD8+ infiltration (~5-15 cells/mm²)
- ✅ High Ki67 proliferation (~45-55%)
- ✅ TP53+ cells are highly proliferative

**Integration:**
- ✅ Top 3 resistance mechanisms identified
- ✅ Therapeutic recommendations provided
- ✅ Immunotherapy assessment (limited efficacy expected)
- ✅ Monitoring strategy defined

---

## Alternative: Natural Language Prompts

If you prefer not to use the prepared test files, you can also use natural language prompts in Claude Desktop:

```
For patient PAT001-OVC-2025:
1. Retrieve patient demographics and CA-125 trends
2. Parse somatic variants (TP53, PIK3CA, PTEN)
3. Analyze multi-omics PDX resistance data
4. Examine spatial gene expression patterns
5. Quantify immune infiltration from imaging

Files are in: patient-data/PAT001-OVC-2025/
```

Claude Desktop will use the appropriate MCP server tools to access the data.

---

## Summary

✅ **Data Ready:** 17 files in `/data/patient-data/PAT001-OVC-2025/`
✅ **Tests Ready:** 5 FIXED test prompts available
✅ **Documentation:** Complete guides for testing and troubleshooting

**Start with:** `TEST_1_CLINICAL_GENOMIC.txt`

**Need help?** See `QUICK_TEST_REFERENCE.md` or `CLAUDE_DESKTOP_FILE_ACCESS_GUIDE.md`

---

**Last Updated:** November 12, 2025
**Status:** Ready for Claude Desktop testing
