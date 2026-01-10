# CSV/Tabular Data Workflow

**Status:** ✅ Production - Current Implementation
**Last Updated:** January 10, 2026
**Used In:** PatientOne end-to-end tests (TEST_3_SPATIAL.txt)

---

## Overview

The CSV workflow processes **pre-processed spatial transcriptomics data** in tabular format (CSV files). This is the **current implementation** used in all PatientOne tests and production deployments.

**Why CSV Instead of FASTQ?**
- Faster iteration (no alignment needed)
- Focuses on analysis rather than data processing
- Suitable for pre-processed Visium data
- Enables immediate statistical and visualization insights

---

## Data Format

### Input Files (3 CSV files required)

**Location:** `/data/patient-data/PAT001-OVC-2025/spatial/`

#### 1. `visium_spatial_coordinates.csv`
Spatial location of each spot on the tissue.

**Columns:** `barcode`, `x`, `y`
**Rows:** 900 spots

#### 2. `visium_gene_expression.csv`
Gene expression values for each spot.

**Columns:** `barcode` + 31 gene columns
**Key Genes:**
- Proliferation: MKI67, PCNA, TOP2A
- Resistance: PIK3CA, AKT1, ABCB1, BCL2L1
- Immune: CD3D, CD8A, CD68, PDCD1

#### 3. `visium_region_annotations.csv`
Tissue region assignment for each spot.

**Columns:** `barcode`, `region`
**Regions (6 total):**
1. tumor_core
2. tumor_proliferative
3. tumor_interface
4. stroma_immune
5. stroma
6. necrotic_hypoxic

---

## Workflow Steps

### Step 1: Data Loading
Load 3 CSV files into memory, validate barcode matching and data completeness.

**Result:** 900 spots × 31 genes with (x,y) coordinates and region labels

### Step 2: Differential Expression Analysis
**Tool:** `perform_differential_expression`

Identify genes with statistically significant differences between regions.

**Example:** Compare tumor_core vs stroma_immune using Wilcoxon test

**Output:**
```json
{
  "comparison": "tumor_core vs stroma_immune",
  "method": "wilcoxon",
  "significant_genes": [
    {
      "gene": "MKI67",
      "mean_group1": 420.5,
      "mean_group2": 45.2,
      "fold_change": 9.3,
      "p_value": 0.0001,
      "adjusted_p_value": 0.003
    }
  ]
}
```

**Key Findings (PatientOne):**
- MKI67, PCNA highly expressed in tumor_proliferative
- CD8A LOW in tumor regions (immune exclusion)
- PIK3CA, AKT1 elevated in tumor_core (resistance)

### Step 3: Spatial Autocorrelation
**Tool:** `calculate_spatial_autocorrelation`

Detect spatial clustering patterns using Moran's I statistic.

**Moran's I interpretation:**
- +1: Perfect clustering (similar values near each other)
- 0: Random spatial pattern
- -1: Checkerboard pattern (dissimilar values near each other)

**Example:** Calculate Moran's I for CD8A expression

**Output:**
```json
{
  "gene": "CD8A",
  "morans_i": 0.42,
  "p_value": 0.001,
  "interpretation": "Significant positive spatial autocorrelation"
}
```

**Key Findings (PatientOne):**
- CD8A: Clustered at tumor periphery (Moran's I = 0.42)
- MKI67: Strong clustering in proliferative zones (Moran's I = 0.68)
- Resistance genes: Heterogeneous distribution

### Step 4: Pathway Enrichment
**Tool:** `perform_pathway_enrichment`

Identify enriched biological pathways from gene lists.

**Databases:** GO, KEGG, Reactome, Hallmark

**Example:** Enrich genes upregulated in tumor_core

**Output:**
```json
{
  "pathway": "PI3K-AKT signaling pathway",
  "genes_in_pathway": ["PIK3CA", "AKT1", "mTOR", "RPS6KB1"],
  "p_value": 0.0001,
  "adjusted_p_value": 0.003,
  "enrichment_score": 4.2
}
```

**Key Findings (PatientOne):**
- PI3K/AKT/mTOR pathway activated (drug resistance)
- Cell cycle pathways enriched (proliferation)
- Immune response pathways depleted (exclusion)

### Step 5: Batch Correction (Optional)
**Tool:** `perform_batch_correction`

Remove technical batch effects when analyzing multiple samples.

**Methods:** ComBat, Harmony, Scanorama

**When to use:** Multiple tissue sections, different sequencing runs, cross-patient comparisons

**Note:** PatientOne uses single patient, so batch correction not needed in current tests.

### Step 6: Visualization
**Tools:** 4 visualization tools (see [TOOL_REFERENCE.md](TOOL_REFERENCE.md) for details)

1. **Spatial Heatmap** - Expression overlaid on (x,y) coordinates
2. **Gene Expression Heatmap** - Clustered heatmap (genes × regions)
3. **Region Composition Chart** - Bar chart of spot counts
4. **Spatial Autocorrelation Plot** - Moran's I visualization

**Output:** PNG files with timestamps

### Step 7: Integration with Multi-Omics
**Tool:** `get_spatial_data_for_patient`

Bridge spatial findings to mcp-multiomics for cross-modality integration.

**Output:**
```json
{
  "patient_id": "PAT001-OVC-2025",
  "spatial_metrics": {
    "immune_infiltration": "LOW",
    "cd8_density": 5.2,
    "proliferation_index": 0.52,
    "resistance_score": 0.78,
    "spatial_heterogeneity": "HIGH"
  }
}
```

---

## Expected Results (PatientOne)

### Spatial Patterns
- **Immune Exclusion:** CD8+ cells LOW in tumor core
- **Proliferation:** MKI67, PCNA HIGH in tumor_proliferative region
- **Resistance:** PIK3CA, AKT1 elevated in tumor_core
- **Heterogeneity:** Spatial variation in resistance markers

### Statistical Findings
- 8 genes significantly different between tumor vs stroma
- CD8A clustered (Moran's I = 0.42, p < 0.001)
- PI3K/AKT/mTOR pathway activated (p < 0.001)

### Clinical Implications
- Immune exclusion → immunotherapy may have limited efficacy
- PI3K/AKT activation → consider PI3K inhibitors (alpelisib)
- Spatial heterogeneity → tumor sampling bias considerations

---

## Testing

### PatientOne End-to-End Test
**File:** [TEST_3_SPATIAL.txt](../../tests/manual_testing/PatientOne-OvarianCancer/implementation/TEST_3_SPATIAL.txt)

**Test execution:**
- Duration: 3-5 minutes
- Data: 900 spots, 31 genes, 6 regions
- Validations: Statistical results + 4 visualizations generated

### Unit Tests
**Location:** `/servers/mcp-spatialtools/tests/`
**Coverage:** >80% for production tools

**Running tests:**
```bash
cd servers/mcp-spatialtools
pytest tests/ -v
```

**Key test files:**
- `test_differential_expression.py`
- `test_spatial_autocorrelation.py`
- `test_pathway_enrichment.py`
- `test_batch_correction.py`
- `test_visualizations.py`

---

## See Also
- [TOOL_REFERENCE.md](TOOL_REFERENCE.md) - Detailed tool documentation
- [OVERVIEW.md](OVERVIEW.md) - System architecture
- [PatientOne README](../../tests/manual_testing/PatientOne-OvarianCancer/README.md) - Complete end-to-end workflow
