# Prompt Library

Reusable prompts for common precision medicine analyses. Copy, paste, and customize for your use case.

---

## How to Use This Library

1. **Browse by category** - Clinical, Genomic, Multi-omics, Spatial, Workflows
2. **Copy the prompt** - Customize patient ID, parameters, thresholds
3. **Paste into Claude Desktop** - Or Streamlit UI, Jupyter Notebook
4. **Review results** - Verify output, iterate if needed

---

## Prompt Categories

### [Clinical & Genomic Prompts](clinical-genomic.md)
Basic queries for clinical data and genomic variants

**Examples:**
- Get patient demographics and clinical history
- Identify pathogenic variants from VCF
- Link clinical phenotype to genotype
- Generate variant interpretation report

**Servers used:** mcp-epic, mcp-mockepic, mcp-fgbio

---

### [Multi-Omics Prompts](multiomics.md)
Integration of RNA, protein, and phosphoproteomics data

**Examples:**
- Load and validate multi-omics datasets
- Run Stouffer meta-analysis across modalities
- Identify upstream regulators
- Pathway enrichment on integrated data

**Servers used:** mcp-multiomics, mcp-fgbio

---

### [Spatial Transcriptomics Prompts](spatial-transcriptomics.md)
Analysis of Visium spatial RNA-seq data

**Examples:**
- Load and QC spatial data
- Spatial differential expression
- Moran's I spatial autocorrelation
- Cell type deconvolution
- Spatial pathway enrichment
- Batch correction with ComBat

**Servers used:** mcp-spatialtools, mcp-openimagedata

---

### [Complete Workflow Prompts](workflows.md)
End-to-end multi-modal analyses

**Examples:**
- PatientOne comprehensive analysis (35 minutes)
- Tumor microenvironment characterization
- Drug resistance mechanism identification
- Biomarker discovery workflow

**Servers used:** All 10 servers

---

## Quick Reference: PatientOne

**Patient ID:** `PAT001-OVC-2025`

**Available Data:**
- Clinical: FHIR resources (mcp-mockepic)
- Genomic: VCF file with TP53, BRCA1 variants
- Multi-omics: RNA, Protein, Phospho (15 samples)
- Spatial: Visium data (900 spots × 31 genes)
- Imaging: H&E slides, multiplex IF

**Data paths (DRY_RUN mode):**
- Clinical: Retrieved via API
- Genomic: `/data/patient-data/PAT001-OVC-2025/genomic/variants.vcf`
- Multi-omics: `/data/patient-data/PAT001-OVC-2025/multiomics/`
- Spatial: `/data/patient-data/PAT001-OVC-2025/spatial/`
- Imaging: `/data/patient-data/PAT001-OVC-2025/imaging/`

---

## Prompt Best Practices

### Be Specific

**❌ Vague:**
```
Analyze PatientOne data
```

**✅ Specific:**
```
Perform spatial pathway enrichment on PatientOne (PAT001-OVC-2025) tumor regions,
focusing on cancer-related KEGG pathways with FDR < 0.05.
```

### Include Parameters

**❌ Missing parameters:**
```
Run differential expression analysis
```

**✅ With parameters:**
```
Run differential expression analysis comparing PatientOne tumor vs. normal samples,
using Mann-Whitney U test with FDR < 0.05 threshold.
```

### Specify Expected Output

**❌ No output guidance:**
```
Find upregulated genes
```

**✅ With output guidance:**
```
Find upregulated genes in PatientOne tumor regions (log2FC > 1, FDR < 0.05).
Return top 20 genes ranked by fold change, with visualization (volcano plot).
```

### Chain Multiple Steps

**❌ Single step:**
```
Load spatial data
```

**✅ Multi-step workflow:**
```
For PatientOne spatial data:
1. Load Visium dataset and summarize (spots, genes, regions)
2. Run spatial differential expression (tumor vs. normal)
3. Perform pathway enrichment on upregulated genes
4. Create spatial visualization showing top pathway
```

---

## Example: Building a Complex Prompt

### Goal
Identify treatment targets for PatientOne using multi-modal data

### Step 1: Break Down the Task
- Load clinical context
- Analyze genomic variants
- Integrate multi-omics data
- Analyze spatial heterogeneity
- Synthesize into treatment recommendations

### Step 2: Write the Prompt
```
Perform comprehensive multi-modal analysis for PatientOne (PAT001-OVC-2025)
to identify top 3 treatment targets:

1. Clinical Context:
   - Get patient demographics, diagnosis, and treatment history from mcp-mockepic
   - Extract relevant biomarkers (CA-125, BRCA status)

2. Genomic Analysis:
   - Load VCF file and identify pathogenic variants (TP53, BRCA1, PIK3CA, PTEN)
   - Interpret clinical significance of each variant

3. Multi-Omics Integration:
   - Load RNA, protein, and phospho data
   - Run Stouffer meta-analysis to identify concordant pathway activations
   - Perform pathway enrichment (KEGG, Hallmark) with FDR < 0.05

4. Spatial Analysis:
   - Load Visium spatial transcriptomics data
   - Perform spatial pathway enrichment focusing on tumor regions
   - Identify spatially variable genes (Moran's I, p < 0.05)

5. Treatment Recommendations:
   - Map activated pathways to FDA-approved drugs
   - Prioritize by strength of evidence (genomic + multi-omics + spatial)
   - Suggest 3 top treatment options with rationale

Return structured report with:
- Clinical summary (1 paragraph)
- Key genomic findings (bullet list)
- Top 3 activated pathways (with p-values)
- Spatial heterogeneity summary
- Treatment recommendations (ranked 1-3 with evidence)
```

### Step 3: Refine Based on Results
If initial results are unclear:
- Add more specific parameters (thresholds, methods)
- Break into smaller sub-prompts
- Request specific visualizations
- Ask for clarification of ambiguous results

---

## Customization Guide

### Replace Patient ID
```
# Template
Analyze PatientOne (PAT001-OVC-2025) spatial data...

# Your patient
Analyze Patient123 (PAT123-BRCA-2025) spatial data...
```

### Adjust Thresholds
```
# Default
...with FDR < 0.05 threshold

# More stringent
...with FDR < 0.01 threshold (Bonferroni correction)

# More lenient (exploratory)
...with p < 0.05 threshold (uncorrected)
```

### Change Methods
```
# Default
...using Mann-Whitney U test

# Parametric alternative
...using Welch's t-test

# For paired samples
...using Wilcoxon signed-rank test
```

### Specify Output Format
```
# Table
Return results as CSV table with columns: gene, log2FC, p_value, FDR

# Visualization
Create volcano plot showing all genes with significant genes highlighted

# Summary statistics
Return summary statistics: total significant genes, top upregulated, top downregulated
```

---

## Common Variables

### Patient IDs
- `PAT001-OVC-2025` - PatientOne (Stage IV ovarian cancer)
- `[YOUR-PATIENT-ID]` - Replace with your patient identifier

### Data Paths (DRY_RUN mode)
- Clinical: Retrieved via API (no path needed)
- Genomic: `/data/patient-data/{PATIENT_ID}/genomic/variants.vcf`
- Multi-omics: `/data/patient-data/{PATIENT_ID}/multiomics/rna_counts.csv`
- Spatial: `/data/patient-data/{PATIENT_ID}/spatial/filtered_feature_matrix/`
- Imaging: `/data/patient-data/{PATIENT_ID}/imaging/`

### Thresholds
- **FDR:** 0.05 (standard), 0.01 (stringent), 0.1 (exploratory)
- **Log2 Fold Change:** 1.0 (2x), 1.5 (3x), 2.0 (4x)
- **Spatial p-value:** 0.05 (standard)

### Methods
- **Differential Expression:** Mann-Whitney U (non-parametric), t-test (parametric)
- **Multiple Testing:** Benjamini-Hochberg (FDR), Bonferroni (strict)
- **Pathway Enrichment:** Fisher's exact test
- **Batch Correction:** ComBat (Empirical Bayes)
- **Meta-Analysis:** Stouffer's Z-score method

---

## Validation & Quality Control

### Check Your Results

**After running a prompt, validate:**

1. **Data loaded correctly**
   - Prompt: "Summarize the data loaded (number of samples, genes, spots)"
   - Expected: Reasonable numbers matching your dataset

2. **Statistical significance**
   - Check p-values and FDR-corrected values
   - Ensure multiple testing correction was applied

3. **Biological plausibility**
   - Do top genes/pathways make sense for the cancer type?
   - Are fold changes in expected direction?

4. **Reproducibility**
   - Re-run same prompt - results should be identical (in non-DRY_RUN mode)
   - Document parameters used

### Common Issues

**"No significant results"**
- Threshold too stringent? Try FDR < 0.1
- Sample size too small? May lack statistical power
- Batch effects? Apply ComBat correction first

**"Too many significant results"**
- Threshold too lenient? Use FDR < 0.01
- Effect size filter? Add log2FC > 1 requirement
- Biologically relevant? Focus on specific pathways

**"Results don't match literature"**
- Different dataset? Synthetic vs. real data
- Different methods? Check statistical test used
- Different thresholds? Verify cutoffs match

---

## Contributing Prompts

Have a useful prompt? Share it!

**Guidelines:**
- Clear objective (1 sentence)
- Copy-paste ready (include all parameters)
- Expected output described
- Tested with PatientOne

**Submit via:**
- GitHub PR to `/docs/prompt-library/`
- Include prompt file (markdown)
- Example output (optional)

---

## Related Resources

- **[PatientOne Guide](../test-docs/patient-one-scenario/README.md)** - Complete walkthrough
- **[Test Prompts](../test-docs/patient-one-scenario/test-prompts/)** - 6 test scenarios
- **[Server Documentation](../../servers/README.md)** - Tool reference
- **[For Researchers](../for-researchers/README.md)** - Analysis workflows

---

**Last Updated:** 2026-01-14
