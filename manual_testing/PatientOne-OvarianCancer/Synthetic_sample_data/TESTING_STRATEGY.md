# Testing Strategy: Avoiding Context Limits

**Problem:** The comprehensive end-to-end test (COPY_PASTE_PROMPT.txt) hits Claude Desktop's context limit because it tries to load all data files and run all 9 MCP servers at once.

**Solution:** Break testing into 5 smaller, focused tests that can be run independently.

---

## Test Suite Overview

### Test 1: Clinical Data & Genomic Analysis (5-10 minutes)
**File:** `TEST_1_CLINICAL_GENOMIC.txt`
**Servers Used:** mcp-mockepic, mcp-fgbio, mcp-tcga
**Data Files:** 3 small files (2 JSON, 1 VCF)
**Context Usage:** ⭐ Low

**Tests:**
- Patient demographics and family history
- CA-125 trend analysis
- Somatic mutations (TP53, PIK3CA, PTEN)
- TCGA cohort comparison

---

### Test 2: Multi-Omics Resistance Analysis (5-10 minutes)
**File:** `TEST_2_MULTIOMICS.txt`
**Servers Used:** mcp-multiomics
**Data Files:** 4 CSV files (subset of genes)
**Context Usage:** ⭐⭐ Medium

**Tests:**
- PDX sample integration (15 samples)
- Stouffer's meta-analysis for 6 key genes only
  - AKT1, PIK3CA, MTOR, PTEN, ABCB1, BCL2L1
- FDR correction
- Pathway activation (PI3K/AKT)

**Note:** Only analyzes 6 genes instead of all 1000 to reduce context usage

---

### Test 3: Spatial Transcriptomics (5-10 minutes)
**File:** `TEST_3_SPATIAL.txt`
**Servers Used:** mcp-spatialtools
**Data Files:** 3 CSV files (subset of genes)
**Context Usage:** ⭐⭐ Medium

**Tests:**
- Spatial region identification (6 regions)
- Expression of 8 key genes across regions
- Resistance marker spatial heterogeneity
- Immune cell localization
- Immune exclusion analysis

**Note:** Only analyzes 8 genes instead of all 31 to reduce context usage

---

### Test 4: Histology & Imaging (5-10 minutes)
**File:** `TEST_4_IMAGING.txt`
**Servers Used:** mcp-openimagedata, mcp-deepcell
**Data Files:** 4 TIFF images (subset)
**Context Usage:** ⭐⭐⭐ Medium-High

**Tests:**
- H&E tissue architecture
- CD8 T cell quantification
- Ki67 proliferation index
- Multiplex IF cell segmentation

**Note:** Only processes 4 images instead of all 7 to reduce context usage

---

### Test 5: Integration & Recommendations (5 minutes)
**File:** `TEST_5_INTEGRATION.txt`
**Servers Used:** None (synthesis only)
**Data Files:** None (uses results from Tests 1-4)
**Context Usage:** ⭐ Low

**Tests:**
- Synthesize findings across all modalities
- Identify primary resistance mechanisms
- Generate treatment recommendations
- Propose monitoring strategy

**Note:** Run this AFTER completing Tests 1-4

---

## Recommended Testing Approach

### Option 1: Sequential Testing (Recommended)
Run tests in order, one conversation per test:

```
Day 1:
  Conversation 1 → TEST_1_CLINICAL_GENOMIC.txt
  Conversation 2 → TEST_2_MULTIOMICS.txt

Day 2:
  Conversation 3 → TEST_3_SPATIAL.txt
  Conversation 4 → TEST_4_IMAGING.txt

Day 3:
  Conversation 5 → TEST_5_INTEGRATION.txt (synthesizes 1-4)
```

**Advantages:**
- ✅ Never hits context limits
- ✅ Can troubleshoot each MCP server independently
- ✅ Clear verification of each data type
- ✅ Easy to identify which server has issues

**Disadvantages:**
- ⏱️ More conversations to manage
- 📋 Need to track results across conversations

---

### Option 2: Combined Testing (If Context Allows)
Try combining some tests in a single conversation:

```
Conversation 1 → TEST_1 + TEST_2 (Clinical + Multi-omics)
Conversation 2 → TEST_3 + TEST_4 (Spatial + Imaging)
Conversation 3 → TEST_5 (Integration)
```

**Try this first:** Start with Conversation 1. If it hits the limit, fall back to Option 1.

---

### Option 3: Minimal Testing (Quick Validation)
Just verify MCP servers work with existing data:

```
Conversation 1 → TEST_1_CLINICAL_GENOMIC.txt (verify 3 servers)
Conversation 2 → TEST_2_MULTIOMICS.txt (verify multiomics server)
Skip Tests 3-5
```

**Use this for:** Quick verification that MCP servers can read existing files

---

## How to Use Each Test

### Step 1: Open test file
```bash
cat TEST_1_CLINICAL_GENOMIC.txt
```

### Step 2: Copy entire contents

### Step 3: Paste into Claude Desktop

### Step 4: Verify results
Check that:
- ✅ Existing files were read (not new data generated)
- ✅ Expected MCP servers were invoked
- ✅ Data values match synthetic data
- ✅ Analysis completed without errors

### Step 5: Save results
Copy Claude Desktop's response to a file:
```
TEST_1_RESULTS.md
TEST_2_RESULTS.md
...
```

---

## Verification Checklist

### Test 1 (Clinical + Genomic)
- [ ] Patient name: Sarah Elizabeth Anderson
- [ ] Patient age: 58 years
- [ ] CA-125 trend: 1456 → 22 → 389 → 289
- [ ] TP53 mutation: R175H
- [ ] PIK3CA mutation: E545K
- [ ] PTEN LOH identified
- [ ] TCGA subtype identified

### Test 2 (Multi-Omics)
- [ ] 15 PDX samples identified (7R, 8S)
- [ ] 6 genes analyzed across 3 modalities
- [ ] Stouffer's Z-scores calculated
- [ ] FDR correction applied
- [ ] PI3K/AKT pathway identified as activated

### Test 3 (Spatial)
- [ ] 6 spatial regions identified
- [ ] 900 spots total
- [ ] Resistance markers show spatial heterogeneity
- [ ] Immune cells localized
- [ ] Immune exclusion noted

### Test 4 (Imaging)
- [ ] H&E analyzed (necrosis, cellularity)
- [ ] CD8 quantified
- [ ] Ki67 index calculated
- [ ] Multiplex IF segmented

### Test 5 (Integration)
- [ ] Primary resistance mechanisms identified
- [ ] Treatment recommendations provided
- [ ] Monitoring strategy proposed

---

## Troubleshooting

### "Still hitting context limit in Test 2"
**Solution:** Reduce genes further
- Only test: AKT1, PIK3CA, ABCB1 (3 genes instead of 6)

### "Still hitting context limit in Test 3"
**Solution:** Reduce to 4 genes
- Only test: MKI67, PIK3CA, CD8A, CD68

### "Still hitting context limit in Test 4"
**Solution:** Test images one at a time
- Conversation 4a: H&E only
- Conversation 4b: CD8 IF only
- Conversation 4c: Multiplex IF only

### "MCP server not responding"
**Solution:** Check server status
```bash
# In Claude Desktop, ask:
"What MCP servers are available?"

# Verify the specific server works:
"Can you use mcp-multiomics to list available tools?"
```

---

## Expected Total Time

**Sequential Testing (Option 1):** 25-45 minutes total
- Test 1: 5-10 min
- Test 2: 5-10 min
- Test 3: 5-10 min
- Test 4: 5-10 min
- Test 5: 5 min

**Combined Testing (Option 2):** 15-30 minutes total
- Tests 1+2: 10-15 min
- Tests 3+4: 10-15 min
- Test 5: 5 min

**Minimal Testing (Option 3):** 10-20 minutes total
- Test 1: 5-10 min
- Test 2: 5-10 min

---

## Why This Approach Works

### 1. Reduced Data Loading
- Instead of loading ALL genes: Load 6-8 key genes
- Instead of loading ALL images: Load 3-4 key images
- Instead of loading ALL spots: Summarize by region

### 2. Focused Analysis
- Each test has a clear objective
- Easier to verify correctness
- Simpler to debug issues

### 3. Modular Design
- Can skip tests if needed
- Can re-run failed tests independently
- Can combine tests if context allows

### 4. Better Verification
- Clear checkpoints for each data type
- Easy to identify which MCP server has issues
- Results are more manageable

---

## Files Summary

```
Synthetic_sample_data/
├── COPY_PASTE_PROMPT.txt          ⚠️ May hit context limit
├── END_TO_END_TEST_PROMPT.md      ⚠️ May hit context limit
│
├── TESTING_STRATEGY.md            ← This file
│
├── TEST_1_CLINICAL_GENOMIC.txt    ✅ Small, safe to run
├── TEST_2_MULTIOMICS.txt          ✅ Medium, should be safe
├── TEST_3_SPATIAL.txt             ✅ Medium, should be safe
├── TEST_4_IMAGING.txt             ✅ Medium-large, may need adjustment
└── TEST_5_INTEGRATION.txt         ✅ Small, safe to run
```

---

**Recommendation:** Start with **TEST_1_CLINICAL_GENOMIC.txt** to verify the MCP servers can read existing files, then proceed with the other tests based on your testing goals.

**Created:** November 12, 2025
**Purpose:** Provide manageable testing approach that avoids Claude Desktop context limits
