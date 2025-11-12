# Test Prompt Update Notes

**Date:** November 12, 2025
**Issue:** Claude Desktop was generating new synthetic data instead of using existing files
**Solution:** Updated prompts to explicitly direct Claude Desktop to use pre-existing data files

---

## Problem

When the original `COPY_PASTE_PROMPT.txt` was used in Claude Desktop, it generated new synthetic patient data instead of reading and analyzing the 17 existing data files that were already created in the repository.

---

## Changes Made

### 1. Updated COPY_PASTE_PROMPT.txt

**Added prominent warning at the beginning:**
```
⚠️ IMPORTANT: All patient data files ALREADY EXIST in this directory.
Please USE THE EXISTING FILES - do NOT generate new synthetic data.
Read the actual files from:
/Users/lynnlangit/Documents/GitHub/spatial-mcp/manual_testing/PatientOne-OvarianCancer/Synthetic_sample_data/
```

**Updated each analysis section to be explicit:**
- Part 1: "READ THE EXISTING FILE: clinical/patient_demographics.json"
- Part 2: "READ THE EXISTING FILE: genomics/somatic_variants.vcf"
- Part 3: "READ THE EXISTING FILES from multiomics/ directory"
- Part 4: "READ THE EXISTING FILES from spatial/ directory"
- Part 4: "READ THE EXISTING IMAGE: imaging/PAT001_tumor_HE_20x.tiff"

**Added reminder at the end:**
```
All data files (JSON, CSV, VCF, TIFF) ALREADY EXIST in:
[directory path]

Do NOT generate new data. Read and analyze the existing files.
There are 17 data files ready for analysis:
- 2 clinical JSON files
- 1 genomics VCF file
- 4 multi-omics CSV files
- 3 spatial CSV files
- 7 imaging TIFF files
```

### 2. Updated END_TO_END_TEST_PROMPT.md

**Added critical warning section:**
```
## ⚠️ CRITICAL: Use Existing Data Files

**All 17 data files ALREADY EXIST** in the directory above. When testing:
- ✅ DO: Read and analyze the existing files
- ❌ DON'T: Generate new synthetic data
- ❌ DON'T: Create new patient data
- ✅ DO: Use the pre-generated files in clinical/, genomics/,
        multiomics/, spatial/, and imaging/ directories
```

**Updated embedded prompt** with same changes as COPY_PASTE_PROMPT.txt

---

## Key Phrases Added

Throughout both files, added emphatic phrases:
- "ALREADY EXISTS"
- "READ THE EXISTING FILE"
- "This file already exists"
- "These files contain real data - read and analyze them"
- "Do NOT generate new data"
- "Use the existing data"

---

## File List Reference

Explicitly listed all 17 existing files:

### Clinical (2 files)
- patient_demographics.json
- lab_results.json

### Genomics (1 file)
- somatic_variants.vcf

### Multi-Omics (4 files)
- sample_metadata.csv
- pdx_rna_seq.csv
- pdx_proteomics.csv
- pdx_phosphoproteomics.csv

### Spatial (3 files)
- visium_spatial_coordinates.csv
- visium_gene_expression.csv
- visium_region_annotations.csv

### Imaging (7 files)
- PAT001_tumor_HE_20x.tiff
- PAT001_tumor_IF_DAPI.tiff
- PAT001_tumor_IF_CD3.tiff
- PAT001_tumor_IF_CD8.tiff
- PAT001_tumor_IF_KI67.tiff
- PAT001_tumor_IF_PanCK.tiff
- PAT001_tumor_multiplex_IF_TP53_KI67_DAPI.tiff

---

## Expected Behavior After Update

When using the updated prompt in Claude Desktop:

1. ✅ Claude Desktop should read patient_demographics.json (not generate demographics)
2. ✅ Claude Desktop should read lab_results.json (not generate lab results)
3. ✅ Claude Desktop should parse somatic_variants.vcf (not generate mutations)
4. ✅ Claude Desktop should load 4 multi-omics CSV files (not generate expression data)
5. ✅ Claude Desktop should load 3 spatial CSV files (not generate spatial data)
6. ✅ Claude Desktop should process 7 existing TIFF images (not generate images)

### MCP Servers Should Be Used:
- mcp-mockepic: To read clinical JSON files
- mcp-fgbio: To parse VCF file
- mcp-tcga: To compare to TCGA-OV cohort
- mcp-multiomics: To integrate CSV files and run Stouffer's analysis
- mcp-spatialtools: To analyze spatial CSV files
- mcp-openimagedata: To process TIFF images
- mcp-deepcell: To segment cells (if available)

---

## Testing Instructions

1. **Open the updated prompt:**
   ```bash
   cat /Users/lynnlangit/Documents/GitHub/spatial-mcp/manual_testing/PatientOne-OvarianCancer/Synthetic_sample_data/COPY_PASTE_PROMPT.txt
   ```

2. **Copy the entire contents**

3. **Paste into Claude Desktop**

4. **Verify in the response:**
   - Check that file paths are mentioned (e.g., "I read patient_demographics.json...")
   - Check that MCP servers are being invoked (look for tool calls in UI)
   - Check that actual data values match the synthetic data (e.g., CA-125: 1456, 22, 389, 289)
   - Check that TP53 R175H mutation is found in VCF
   - Check that 15 PDX samples are identified (7 resistant, 8 sensitive)

5. **Red flags (indicates generating new data):**
   - ❌ "I'll create synthetic patient data..."
   - ❌ "Let me generate some sample data..."
   - ❌ "I'll simulate lab results..."
   - ❌ Generic/different patient details (not Sarah Anderson, age 58)
   - ❌ Different mutation positions or gene names

---

## Verification Checklist

After running the updated prompt in Claude Desktop:

- [ ] Patient name is Sarah Elizabeth Anderson
- [ ] Patient age is 58 years old
- [ ] CA-125 values match: 1456 → 22 → 389 → 289 U/mL
- [ ] TP53 mutation is R175H at chr17:7,578,406
- [ ] PIK3CA mutation is E545K at chr3:178,936,091
- [ ] PTEN LOH at chr10:89,692,940
- [ ] 15 PDX samples identified (7 resistant, 8 sensitive)
- [ ] 900 spatial spots identified
- [ ] 31 genes in spatial expression data
- [ ] 6 spatial regions annotated
- [ ] 7 imaging files processed

If all checkboxes pass ✅ → Test successful, existing data was used

If any checkbox fails ❌ → Claude Desktop may still be generating new data

---

## Troubleshooting

If Claude Desktop still generates new data:

### Option 1: Be even more explicit
Add this to the very beginning of the prompt:
```
STOP - Before proceeding, confirm that you will:
1. Read existing files only
2. NOT generate any synthetic data
3. Use the 17 files that already exist in the directory

Do you understand? If yes, proceed with reading the existing files.
```

### Option 2: Reference specific file contents
Add expected values to the prompt:
```
For example:
- patient_demographics.json should contain: Sarah Elizabeth Anderson, DOB 1966-03-15
- lab_results.json should contain CA-125 value of 1456 at diagnosis
- somatic_variants.vcf should contain TP53 R175H mutation

Please verify you are reading the correct files by confirming these values.
```

### Option 3: Use direct file paths in MCP server calls
If using Claude Code (not Desktop), can directly specify file paths to MCP tools.

---

**Update Status:** ✅ Complete
**Files Updated:** 2 (COPY_PASTE_PROMPT.txt, END_TO_END_TEST_PROMPT.md)
**Files Created:** 1 (PROMPT_UPDATE_NOTES.md - this file)
