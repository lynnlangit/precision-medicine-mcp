# Tool Reference - mcp-spatialtools

All 14 tools with parameters and use cases.

---

## Analysis Tools (10)

### 1. filter_quality
Quality filtering of spatial barcodes and reads.
- **Status:** Implemented (FASTQ workflow only)
- **Use:** Filters reads by barcode quality, UMI uniqueness

### 2. split_by_region
Segment spatial data by tissue regions.
- **Status:** ✅ Production
- **Use:** Split spots into user-defined regions

### 3. align_spatial_data
STAR alignment of spatial transcriptomics reads.
- **Status:** Implemented (FASTQ workflow only)
- **Use:** Align reads to reference genome (requires STAR + indices)

### 4. merge_tiles
Combine multiple spatial tiles into single dataset.
- **Status:** ✅ Production
- **Use:** Multi-tile experiments spanning multiple capture areas

### 5. calculate_spatial_autocorrelation
**Status:** ✅ Production (CSV Workflow)

Compute Moran's I or Geary's C spatial autocorrelation statistics.

**Parameters:**
- `gene`: Gene name to analyze
- `method`: "morans_i" or "gearys_c"
- `weight_type`: "distance" or "knn"

**Returns:**
```json
{
  "gene": "CD8A",
  "method": "morans_i",
  "statistic": 0.42,
  "p_value": 0.001,
  "z_score": 3.2
}
```

### 6. perform_differential_expression
**Status:** ✅ Production (CSV Workflow)

Statistical testing between sample groups.

**Parameters:**
- `group1_spots`: List of spot barcodes for group 1
- `group2_spots`: List of spot barcodes for group 2
- `method`: "wilcoxon", "ttest", or "deseq2_style"

**Returns:** List of genes with fold change, p-value, adjusted p-value

### 7. perform_batch_correction
**Status:** ✅ Production (Multi-sample)

Remove batch effects across samples.

**Methods:** ComBat, Harmony, Scanorama
**Use:** Cross-patient or multi-section analysis

### 8. perform_pathway_enrichment
**Status:** ✅ Production (CSV Workflow)

Gene set enrichment analysis.

**Databases:** GO, KEGG, Reactome, Hallmark
**Method:** Fisher's exact test with FDR correction

### 9. deconvolve_cell_types
**Status:** ⚠️ Synthetic (returns mock signatures)

Cell type deconvolution from bulk spatial data.

**Methods:** CIBERSORTx, EPIC, quanTIseq (mocked)

### 10. get_spatial_data_for_patient
**Status:** ✅ Production (Bridge Tool)

Extract spatial metrics for multi-omics integration.

**Returns:** JSON with spatial findings for mcp-multiomics

---

## Visualization Tools (4)

All visualization tools generate publication-quality PNG images (300 DPI, 10×8 inches) with timestamps.

**Output location:** `SPATIAL_OUTPUT_DIR` environment variable
- Local: `/workspace/output/visualizations/`
- GCP: `/app/data/output/visualizations/`

**Filename format:** `{tool_name}_{gene}_{timestamp}.png`

### 11. generate_spatial_heatmap
**Status:** ✅ Production (Added Jan 8, 2026)

Visualize gene expression overlaid on tissue (x,y) coordinates.

**Parameters:**
- `gene` (required): Gene name to visualize
- `colormap` (optional): matplotlib colormap (default: "viridis")
- `output_filename` (optional): Custom filename

**Output:**
```json
{
  "status": "success",
  "output_file": "/path/to/spatial_heatmap_CD8A_20260109_103000.png",
  "gene": "CD8A",
  "spots_plotted": 900
}
```

**Example:** `Generate a spatial heatmap for CD8A expression showing T cell distribution across the tissue.`

### 12. generate_gene_expression_heatmap
**Status:** ✅ Production (Added Jan 8, 2026)

Clustered heatmap showing genes × regions with hierarchical clustering.

**Parameters:**
- `genes` (required): List of gene names (e.g., ["MKI67", "CD8A", "PIK3CA"])
- `cluster_genes` (optional): Cluster genes (default: True)
- `cluster_regions` (optional): Cluster regions (default: True)
- `colormap` (optional): matplotlib colormap (default: "RdBu_r")

**Output:** PNG with clustered heatmap showing mean expression per region

**Example:** `Generate clustered heatmap for 8 key genes (MKI67, PCNA, PIK3CA, AKT1, ABCB1, CD3D, CD8A, CD68) across all 6 tissue regions.`

### 13. generate_region_composition_chart
**Status:** ✅ Production (Added Jan 8, 2026)

Bar chart showing number of spots per tissue region.

**Parameters:**
- `output_filename` (optional): Custom filename

**Output:** Bar chart PNG showing spot counts for each of the 6 regions

**Example:** `Show the distribution of 900 spots across the 6 tissue regions.`

### 14. visualize_spatial_autocorrelation
**Status:** ✅ Production (Added Jan 9, 2026)

Visualize Moran's I spatial autocorrelation results.

**Parameters:**
- `gene` (required): Gene name
- `morans_i` (required): Moran's I statistic value
- `p_value` (required): Statistical significance
- `output_filename` (optional): Custom filename

**Output:** Scatter plot showing spatial pattern with Moran's I annotation

**Example:** `Visualize spatial autocorrelation for CD8A showing Moran's I = 0.42 (p < 0.001), indicating significant clustering of CD8+ T cells.`

---

## Dependencies

**Visualization tools require:**
- matplotlib >= 3.8.0
- seaborn >= 0.13.0 (for heatmaps)
- pandas >= 2.2.0

---

**See Also:**
- [CSV_WORKFLOW.md](CSV_WORKFLOW.md) - Workflow using these tools
- [FASTQ_WORKFLOW.md](FASTQ_WORKFLOW.md) - FASTQ-specific tools
