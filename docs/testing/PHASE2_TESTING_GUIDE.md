# Phase 2 Testing Guide - Claude Desktop

**Ready to Test:** Differential Expression + Cell Type Deconvolution

---

## IMPORTANT: This Guide is for LOCAL Testing Only

This guide describes testing MCP servers running **locally** with **Claude Desktop** using **STDIO transport**.

**If you want to test GCP-deployed servers:**
- Claude Desktop **cannot** connect to remote GCP Cloud Run servers
- GCP servers use **SSE transport** over HTTP, not STDIO
- You must use the **Claude API** (Python/TypeScript SDK) for testing deployed servers

**See:** [GCP Testing Guide](deployment/GCP_TESTING_GUIDE.md) for instructions on testing Cloud Run deployments.

---

## Prerequisites

✅ Both servers configured in Claude Desktop:
- `epic` - FHIR clinical data (de-identified)
- `spatialtools` - Spatial transcriptomics analysis (real mode)

✅ Patient-001 data available:
- Clinical: Stage IV HGSOC, platinum-resistant, on Bevacizumab
- Spatial: 900 spots × 31 genes, spatial coordinates, region annotations

---

## Testing Checklist

### ✅ Test 1: Cell Type Deconvolution (Quick Test)

**Prompt:**
```
Perform cell type deconvolution on Patient-001 spatial data. Use the ovarian cancer cell type signatures and tell me what you find.
```

**What to expect:**
- **Steps: 2-3 steps** (efficient with token optimization)
- Analysis of 900 spots across 8 cell types
- Summary statistics showing:
  - Tumor cells: ~20%
  - Endothelial cells: ~17% (bevacizumab targets)
  - Immune cells: ~28% combined
  - Fibroblasts/mesenchymal: ~34%
- Clinical interpretation for treatment planning

**Key findings to look for:**
- Endothelial cell presence validates bevacizumab therapy
- CD8+ T-cell infiltration suggests immunotherapy potential
- Mesenchymal signature explains platinum resistance

**Token Efficiency Note:**
- By default, the tool returns summary statistics only (3KB response)
- Spot-level scores (212KB) are excluded to prevent token bloat
- If you need per-spot details, use `include_spot_scores=True`

---

### ✅ Test 2: Differential Expression (Region Comparison)

**Prompt:**
```
Compare gene expression between tumor core and tumor margin for Patient-001.
Split the data by region first, then run differential expression analysis using the Wilcoxon test.
Show me the top differentially expressed genes.
```

**What to expect:**
- Region segmentation: ~219 core spots, ~218 margin spots
- Statistical testing with FDR correction
- Top genes with log2FC, p-values, q-values
- Clinical interpretation of findings

**Note:** With simulated regions, you may see modest fold changes. Look for:
- Statistical test completion
- FDR correction applied
- Proper handling of multiple genes

---

### ✅ Test 3: Integrated Clinical-Spatial Analysis

**Prompt:**
```
For Patient-001 (ovarian cancer, platinum-resistant, on bevacizumab):

1. Get their clinical data from the epic server
2. Retrieve spatial transcriptomics data using the bridge tool
3. Perform cell type deconvolution
4. Generate a clinical-spatial integration report

Focus on:
- Immune infiltration (CD8+ T-cells)
- Angiogenesis markers (for bevacizumab response)
- EMT signatures (platinum resistance)
```

**What to expect:**
- Epic server retrieves de-identified FHIR data:
  - Demographics (57 years old, female)
  - Condition: Stage IV HGSOC
  - Medications: Bevacizumab (active), Carboplatin/Paclitaxel (completed)
  - Observations: CA-125 elevated (487 U/mL), BRCA negative
- Bridge tool maps patient-001 → PAT001-OVC-2025 spatial data
- Cell type deconvolution quantifies cell populations
- Integrated report with treatment recommendations

---

### ✅ Test 4: Regional Cell Type Heterogeneity

**Prompt:**
```
For Patient-001, split the spatial data into tumor core, tumor margin, and immune infiltrate regions.
Then perform cell type deconvolution on each region separately and compare the results.
Tell me how cell type composition differs across regions.
```

**What to expect:**
- Three region-specific analyses
- Comparison table showing cell type enrichment by region
- Insights on spatial heterogeneity
- Clinical implications (e.g., immune cells at margin, tumor cells in core)

---

### ✅ Test 5: Complete Workflow (Most Comprehensive)

**Prompt:**
```
Generate a comprehensive clinical-spatial analysis report for Patient-001:

1. Clinical Context: Get patient demographics, conditions, medications, and biomarkers from the epic server

2. Spatial Analysis:
   - Perform cell type deconvolution on all spots
   - Split data by regions (tumor core, margin, immune zones)
   - Run differential expression between core and margin

3. Integration:
   - Link cell type findings to clinical context (bevacizumab targets, immune infiltration)
   - Explain how spatial patterns relate to platinum resistance
   - Provide treatment recommendations based on spatial composition

4. Summary:
   - Key findings that impact treatment decisions
   - Spatial heterogeneity insights
   - Monitoring recommendations
```

**What to expect:**
- Full end-to-end workflow using both servers
- Clinical context from FHIR store
- Multiple spatial analyses
- Integrated interpretation
- Actionable treatment recommendations

---

## Expected Outputs by Tool

### Cell Type Deconvolution Output

```json
{
  "status": "success",
  "spots_analyzed": 900,
  "cell_types": ["tumor_cells", "cd8_tcells", "endothelial_cells", ...],
  "summary_statistics": {
    "tumor_cells": {
      "mean": 0.000,
      "median": -0.045,
      "std": 0.999,
      "markers_used": 4
    },
    "cd8_tcells": {
      "mean": 0.000,
      "median": 0.023,
      "std": 0.999,
      "markers_used": 3
    },
    ...
  },
  "dominant_cell_type_distribution": {
    "tumor_cells": 181,      // 20.1% of spots
    "fibroblasts": 164,       // 18.2%
    "endothelial_cells": 151, // 16.8%
    ...
  }
}
```

### Differential Expression Output

```json
{
  "status": "success",
  "test_method": "wilcoxon",
  "total_genes_tested": 31,
  "significant_genes": 0-5,  // Depends on data
  "top_upregulated": [
    {
      "gene": "FAP",
      "log2_fold_change": -0.8086,
      "pvalue": 0.116622,
      "qvalue": 0.901882,
      "significant": false
    }
  ],
  "top_downregulated": [...]
}
```

---

## Troubleshooting

### Issue: "Server not connected"

**Solution:**
1. Restart Claude Desktop
2. Check server logs: `~/Library/Logs/Claude/mcp-server-spatialtools.log`
3. Verify configuration in Claude Desktop settings

### Issue: "DRY_RUN mode warning"

**Solution:**
- Check `claude_desktop_config.json`
- Ensure `SPATIAL_DRY_RUN=false` for spatialtools
- Restart Claude Desktop

### Issue: "No spatial data found for patient"

**Solution:**
- Verify patient ID (use "patient-001")
- Check data exists at: `/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/`

### Issue: Serialization errors

**Solution:**
- Already fixed in latest commits
- If issues persist, check logs and report

---

## Expected Clinical Insights

### From Cell Type Deconvolution:

1. **Tumor Purity**
   - ~20% spots dominated by tumor cells
   - Moderate tumor content for analysis

2. **Immune Infiltration**
   - CD8+ T-cells present (~9% of spots)
   - Suggests immunotherapy potential
   - Consider PD-1/PD-L1 checkpoint inhibitors

3. **Angiogenesis (Bevacizumab Target)**
   - Endothelial cells detected (~17% of spots)
   - Validates bevacizumab therapy
   - Monitor for anti-angiogenic response

4. **EMT/Platinum Resistance**
   - Mesenchymal signature present (~16% of spots)
   - Explains platinum resistance mechanism
   - EMT-targeting therapies may be beneficial

5. **Stromal Component**
   - Fibroblasts (CAFs) detected (~18% of spots)
   - Significant stromal contribution
   - May impact drug penetration

### From Differential Expression:

**Note:** With simulated regions, expect modest differences. In real biological data:
- Core: High proliferation markers (MKI67, PCNA, TOP2A)
- Margin: Immune infiltration (CD8A, CD3D), EMT markers (VIM, SNAI1)
- Significance depends on actual biological variation

---

## Success Criteria

✅ **Cell Type Deconvolution:**
- All 8 cell types analyzed
- Summary statistics generated
- Dominant cell types identified
- Clinical interpretation provided

✅ **Differential Expression:**
- Statistical tests complete (Mann-Whitney U or t-test)
- FDR correction applied (Benjamini-Hochberg)
- Top genes identified
- Results match test script output

✅ **Integration:**
- Epic server retrieves FHIR data successfully
- Bridge tool maps patient → spatial dataset
- Multiple analyses coordinate correctly
- Coherent clinical-spatial report generated

✅ **No Errors:**
- No serialization errors
- No "DRY_RUN" warnings for spatialtools
- Proper error handling if issues occur

---

## What to Look For

### Good Signs ✅
- "status": "success" in all responses
- "mode": "real_analysis" (not "dry_run")
- Reasonable cell type distributions (no cell type at 100%)
- Statistical test p-values calculated
- Clinical interpretations make sense

### Warning Signs ⚠️
- "DRY_RUN mode enabled" messages
- All p-values = 1.0 or NaN
- Serialization errors
- "Server not connected"

---

## Next Steps After Testing

**If all tests pass:**
1. Document any interesting findings
2. Consider adding more patients/datasets
3. Move to Phase 2B (Gene Co-expression) or 2D (Pathway Enrichment)
4. Plan production deployment

**If issues found:**
1. Check server logs
2. Note specific error messages
3. Report back for fixes

---

## Quick Reference

**File Paths:**
- Expression: `/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial/visium_gene_expression.csv`
- Coordinates: `/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial/visium_spatial_coordinates.csv`
- Annotations: `/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial/visium_region_annotations.csv`

**Patient Info:**
- ID: patient-001
- Spatial dataset: PAT001-OVC-2025
- Diagnosis: Stage IV HGSOC
- Treatment: Platinum-resistant, on Bevacizumab
- Biomarkers: CA-125 elevated (487), BRCA negative

**Tools Available:**
- `deconvolve_cell_types` - Cell type estimation
- `perform_differential_expression` - Regional gene comparison
- `calculate_spatial_autocorrelation` - Moran's I spatial patterns
- `filter_quality` - QC filtering
- `split_by_region` - Regional segmentation
- `get_spatial_data_for_patient` - Clinical-spatial bridge

---

**Ready to test! Start with Test 1 (Cell Type Deconvolution) for a quick validation.**
