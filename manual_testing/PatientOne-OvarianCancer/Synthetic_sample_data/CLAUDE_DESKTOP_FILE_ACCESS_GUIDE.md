# Claude Desktop File Access Guide

## Problem: Claude Desktop Cannot See Repository Files

Claude Desktop's MCP servers run in containerized/isolated environments and can only access files in specific configured directories.

## Solution Applied

### ✅ Step 1: Data Copied to MCP-Accessible Location

All 17 synthetic patient data files have been copied to:
```
/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/
```

This directory is within the MCP servers' configured data directories.

### ✅ Step 2: File Structure

```
data/patient-data/PAT001-OVC-2025/
├── clinical/
│   ├── patient_demographics.json
│   └── lab_results.json
├── genomics/
│   └── somatic_variants.vcf
├── multiomics/
│   ├── sample_metadata.csv
│   ├── pdx_rna_seq.csv
│   ├── pdx_proteomics.csv
│   └── pdx_phosphoproteomics.csv
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

---

## How to Use in Claude Desktop

### ❌ DON'T: Use Absolute Paths

Claude Desktop CANNOT access files with paths like:
```
/Users/lynnlangit/Documents/GitHub/spatial-mcp/manual_testing/PatientOne-OvarianCancer/...
```

### ✅ DO: Use MCP Server Tools

Claude Desktop CAN access files through MCP server tools. Here's how:

---

## Updated Test Prompts

**✅ All 5 fixed test prompts are now available:**
- `TEST_1_CLINICAL_GENOMIC_FIXED.txt`
- `TEST_2_MULTIOMICS_FIXED.txt`
- `TEST_3_SPATIAL_FIXED.txt`
- `TEST_4_IMAGING_FIXED.txt`
- `TEST_5_INTEGRATION_FIXED.txt`

### Option 1: Use the Fixed Test Prompts (Recommended)

Copy and paste the content from any of the FIXED test files into Claude Desktop. These prompts:
- Use relative paths that MCP servers can access
- Explicitly instruct to use MCP server tools
- Provide expected results for validation
- Include verification checkpoints

Example:
```bash
cat TEST_1_CLINICAL_GENOMIC_FIXED.txt
# Copy the output and paste into Claude Desktop
```

---

### Option 2: Let Claude Desktop Infer the Location

Just describe what you want and let Claude Desktop use the MCP tools:

```
For patient PAT001-OVC-2025, please:
1. Retrieve patient demographics using mockepic
2. Get CA-125 lab results
3. Parse somatic variants (TP53, PIK3CA, PTEN mutations)
4. Compare to TCGA-OV cohort
```

Claude Desktop will use the appropriate MCP server tools to find and read the files.

---

### Option 3: Provide Relative Paths

Give paths relative to the MCP server's data directory:

```
For patient PAT001-OVC-2025:

1. Clinical data (use mockepic):
   Read patient_demographics.json and lab_results.json
   Files are in: patient-data/PAT001-OVC-2025/clinical/

2. Genomic data (use fgbio):
   Read somatic_variants.vcf
   File is in: patient-data/PAT001-OVC-2025/genomics/

3. Multi-omics data (use multiomics):
   Read files from: patient-data/PAT001-OVC-2025/multiomics/
   - sample_metadata.csv
   - pdx_rna_seq.csv
   - pdx_proteomics.csv
   - pdx_phosphoproteomics.csv

4. Spatial data (use spatialtools):
   Read files from: patient-data/PAT001-OVC-2025/spatial/
   - visium_spatial_coordinates.csv
   - visium_gene_expression.csv
   - visium_region_annotations.csv

5. Imaging data (use openimagedata):
   Read files from: patient-data/PAT001-OVC-2025/imaging/
   - PAT001_tumor_HE_20x.tiff
   - PAT001_tumor_IF_CD8.tiff
   - PAT001_tumor_IF_KI67.tiff
   - PAT001_tumor_multiplex_IF_TP53_KI67_DAPI.tiff
```

---

## MCP Server Configuration Reference

Each MCP server has a configured data directory:

| Server | Config Variable | Path |
|--------|----------------|------|
| spatialtools | SPATIAL_DATA_DIR | `/data/` |
| openimagedata | IMAGE_DATA_DIR | `/data/images/` |
| multiomics | MULTIOMICS_DATA_DIR | `/data/multiomics/` |
| fgbio | FGBIO_REFERENCE_DATA_DIR | `/data/reference/` |

Our patient data is in: `/data/patient-data/PAT001-OVC-2025/`

This is accessible to all servers because it's under `/data/`.

---

## Testing Strategy

### Test 1: Clinical + Genomic

```
For patient PAT001-OVC-2025, use the available MCP servers to:

1. Get patient demographics (name, age, family history, BRCA1 status)
2. Get CA-125 lab results and identify platinum resistance pattern
3. Parse somatic variants VCF (look for TP53 R175H, PIK3CA E545K, PTEN LOH)
4. Compare molecular profile to TCGA-OV cohort

Expected findings:
- Patient: Sarah Anderson, 58yo
- BRCA1 germline mutation
- CA-125: 1456 → 22 → 389 → 289 (shows resistance)
- TP53 R175H, PIK3CA E545K, PTEN LOH present
```

### Test 2: Multi-Omics

```
For patient PAT001-OVC-2025, analyze PDX resistance data:

Use multiomics server to:
1. Load sample metadata (expect 15 samples: 7 resistant, 8 sensitive)
2. Analyze these genes: AKT1, PIK3CA, MTOR, PTEN, ABCB1, BCL2L1
3. Run Stouffer's meta-analysis with directionality
4. Apply FDR correction (α=0.05)
5. Which pathway is activated? (expect PI3K/AKT)

Data location: patient-data/PAT001-OVC-2025/multiomics/
```

### Test 3: Spatial

```
For patient PAT001-OVC-2025, analyze spatial transcriptomics:

Use spatialtools to:
1. Load Visium data (expect 900 spots, 6 regions)
2. Analyze genes: MKI67, PCNA, PIK3CA, AKT1, CD3D, CD8A, CD68, ABCB1
3. Where are resistance markers expressed? (spatial heterogeneity)
4. Where are immune cells located? (immune exclusion?)

Data location: patient-data/PAT001-OVC-2025/spatial/
```

### Test 4: Imaging

```
For patient PAT001-OVC-2025, analyze histology images:

Use openimagedata and deepcell to:
1. Load H&E image and describe tissue architecture
2. Load CD8 IF image and quantify T cell infiltration
3. Load Ki67 IF image and calculate proliferation index
4. Load multiplex IF and segment cells

Data location: patient-data/PAT001-OVC-2025/imaging/
```

### Test 5: Integration

```
For patient PAT001-OVC-2025, synthesize findings from Tests 1-4:

Based on clinical, genomic, multi-omics, spatial, and imaging data:
1. What are the primary resistance mechanisms?
2. Which findings are consistent across all modalities?
3. What are the top 3 therapeutic recommendations?
4. Should immunotherapy be considered? Why or why not?

No data loading required - synthesis only.
```

---

## Verification

### To verify files are accessible:

```
Can you list the available data for patient PAT001-OVC-2025?
Check what files are in the patient-data/PAT001-OVC-2025/ directory.
```

Expected response: Claude Desktop should list all 17 files across the 5 subdirectories.

---

## Troubleshooting

### Issue: "Cannot find file"

**Cause:** File path is not relative to MCP server data directory

**Solution:** Use paths relative to `/data/`, like:
- `patient-data/PAT001-OVC-2025/clinical/patient_demographics.json`

### Issue: "Generating synthetic data"

**Cause:** Prompt doesn't explicitly ask to use existing files

**Solution:** Add:
- "Use the EXISTING patient data for PAT001-OVC-2025"
- "Do NOT generate new data"
- "Read from patient-data/PAT001-OVC-2025/ directory"

### Issue: "Permission denied"

**Cause:** MCP server configuration doesn't allow access

**Solution:** Verify file is in `/data/patient-data/PAT001-OVC-2025/` (should work)

---

## Summary

✅ **Data Location:** `/data/patient-data/PAT001-OVC-2025/`
✅ **17 files copied:** Clinical (2), Genomics (1), Multi-omics (4), Spatial (3), Imaging (7)
✅ **Access Method:** Use MCP server tools (mockepic, fgbio, multiomics, spatialtools, openimagedata)
✅ **Path Format:** Relative to `/data/`, like `patient-data/PAT001-OVC-2025/clinical/`

**Use:** All 5 fixed test prompts (`TEST_1_CLINICAL_GENOMIC_FIXED.txt` through `TEST_5_INTEGRATION_FIXED.txt`)

---

**Created:** November 12, 2025
**Purpose:** Enable Claude Desktop to access synthetic patient data through MCP servers
