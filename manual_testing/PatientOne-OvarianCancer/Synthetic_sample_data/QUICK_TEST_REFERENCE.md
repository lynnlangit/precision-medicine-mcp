# Quick Test Reference Card

## Problem: Context Limit + File Access

The full end-to-end test is too large for Claude Desktop, and Claude Desktop's MCP servers run in containerized environments that can't access arbitrary file paths.

**✅ Solution:** Use the FIXED test prompts that use MCP-accessible data locations.

---

## 5 Focused Tests (Run Separately)

| Test | File | Servers | Time | Context |
|------|------|---------|------|---------|
| **1. Clinical & Genomic** | `TEST_1_CLINICAL_GENOMIC_FIXED.txt` | mockepic, fgbio, tcga | 5-10 min | Low ⭐ |
| **2. Multi-Omics** | `TEST_2_MULTIOMICS_FIXED.txt` | multiomics | 5-10 min | Medium ⭐⭐ |
| **3. Spatial** | `TEST_3_SPATIAL_FIXED.txt` | spatialtools | 5-10 min | Medium ⭐⭐ |
| **4. Imaging** | `TEST_4_IMAGING_FIXED.txt` | openimagedata, deepcell | 5-10 min | Medium-High ⭐⭐⭐ |
| **5. Integration** | `TEST_5_INTEGRATION_FIXED.txt` | none (synthesis) | 5 min | Low ⭐ |

**Total Time:** 25-45 minutes (all 5 tests)

**Note:** The FIXED versions use relative paths (`patient-data/PAT001-OVC-2025/`) that MCP servers can access.

---

## Quick Start

### Step 1: Run Test 1
```bash
cat TEST_1_CLINICAL_GENOMIC_FIXED.txt
```
Copy → Paste into Claude Desktop → Verify results

### Step 2: Run Test 2
```bash
cat TEST_2_MULTIOMICS_FIXED.txt
```
Copy → Paste into NEW conversation → Verify results

### Step 3: Continue with Tests 3, 4, 5
```bash
cat TEST_3_SPATIAL_FIXED.txt
cat TEST_4_IMAGING_FIXED.txt
cat TEST_5_INTEGRATION_FIXED.txt
```
Use a new conversation for each test

**Note:** Test 5 (Integration) should be run AFTER completing Tests 1-4, as it synthesizes all findings.

---

## What Each Test Does

### Test 1: Clinical & Genomic ✅
- Patient demographics (Sarah Anderson, 58yo)
- CA-125 trend (1456→22→389→289)
- Mutations (TP53 R175H, PIK3CA E545K)
- TCGA comparison

### Test 2: Multi-Omics ✅
- 15 PDX samples (7 resistant, 8 sensitive)
- 6 key genes across RNA/Protein/Phospho
- Stouffer's meta-analysis
- PI3K/AKT pathway activation

### Test 3: Spatial ✅
- 6 spatial regions, 900 spots
- 8 key genes across regions
- Resistance marker distribution
- Immune cell localization

### Test 4: Imaging ✅
- H&E tissue architecture
- CD8 T cell quantification
- Ki67 proliferation
- Multiplex IF segmentation

### Test 5: Integration ✅
- Synthesize findings from Tests 1-4
- Resistance mechanisms
- Treatment recommendations
- Monitoring strategy

---

## Verification Checklist

Quick checks for each test:

**Test 1:**
- [ ] Patient is Sarah Anderson, age 58
- [ ] CA-125 values: 1456, 22, 389, 289
- [ ] TP53 R175H found
- [ ] TCGA subtype identified

**Test 2:**
- [ ] 15 PDX samples (7R, 8S)
- [ ] 6 genes analyzed
- [ ] PI3K/AKT pathway activated

**Test 3:**
- [ ] 6 regions, 900 spots
- [ ] Resistance markers spatially heterogeneous
- [ ] Immune exclusion noted

**Test 4:**
- [ ] H&E analyzed
- [ ] CD8 quantified
- [ ] Ki67 index calculated

**Test 5:**
- [ ] Resistance mechanisms identified
- [ ] Treatment recommendations provided

---

## If Still Hitting Context Limits

### Test 2 too large?
Reduce to 3 genes: AKT1, PIK3CA, ABCB1

### Test 3 too large?
Reduce to 4 genes: MKI67, PIK3CA, CD8A, CD68

### Test 4 too large?
Split into 3 conversations:
- 4a: H&E only
- 4b: CD8 IF only
- 4c: Multiplex IF only

---

## Files Location

All test files in:
```
/Users/lynnlangit/Documents/GitHub/spatial-mcp/manual_testing/
PatientOne-OvarianCancer/Synthetic_sample_data/
```

---

**Start Here:** `TEST_1_CLINICAL_GENOMIC_FIXED.txt`

**Read Full Guides:**
- `TESTING_STRATEGY.md` - Detailed testing strategy
- `CLAUDE_DESKTOP_FILE_ACCESS_GUIDE.md` - Why FIXED versions are needed
