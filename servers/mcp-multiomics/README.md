# mcp-multiomics: Multi-Omics PDX Data Analysis

MCP server for integrating and analyzing multi-omics data from Patient-Derived Xenograft (PDX) models. Combines RNA, protein, and phosphorylation data using advanced statistical methods including HAllA association testing and Stouffer's meta-analysis.

## Overview

This server enables:
- **Multi-Omics Integration**: Align and normalize RNA, protein, and phosphorylation data
- **Stouffer's Meta-Analysis**: Combine p-values across modalities with directionality
- **HAllA Analysis**: Hierarchical all-against-all association testing (R-based, optional)
- **Visualization**: Multi-omics heatmaps and PCA plots

## Tools

### 1. integrate_omics_data

Integrate multi-omics data from RNA, protein, and phosphorylation datasets.

**Parameters:**
- `rna_path` (required): Path to RNA expression data (CSV/TSV, genes x samples)
- `protein_path` (optional): Path to protein abundance data
- `phospho_path` (optional): Path to phosphorylation data
- `metadata_path` (optional): Path to sample metadata (must include 'Sample' column)
- `normalize` (default: True): Apply Z-score normalization within each modality
- `filter_missing` (default: 0.5): Remove features with >X fraction missing (0.0-1.0)

**Returns:**
- `integrated_data`: Aligned data matrices per modality
- `common_samples`: List of samples present in all modalities
- `feature_counts`: Number of features per modality after filtering
- `metadata`: Sample metadata if provided
- `qc_metrics`: Quality control statistics
- `cache_path`: Path to cached integrated data (for downstream tools)

**Example:**
```
Claude, please integrate my multi-omics PDX data:
- RNA: /data/pdx_rna_fpkm.csv
- Protein: /data/pdx_protein_tmt.csv
- Phospho: /data/pdx_phospho.csv
- Metadata: /data/sample_info.csv

Apply Z-score normalization and filter features with >50% missing data.
```

### 2. calculate_stouffer_meta

Combine p-values across omics modalities using Stouffer's Z-score method.

**Parameters:**
- `p_values_dict` (required): Dict of modality -> list of p-values (one per feature)
- `effect_sizes_dict` (optional): Dict of modality -> list of effect sizes (log2FC or correlation)
- `weights` (optional): Dict of modality -> weight (default: equal weights)
- `use_directionality` (default: True): Incorporate effect size sign into Z-scores

**Returns:**
- `meta_p_values`: Combined p-values per feature
- `meta_z_scores`: Combined Z-scores
- `q_values`: FDR-corrected q-values
- `significant_features`: Features passing FDR threshold with detailed info
- `statistics`: Summary statistics and diagnostics

**Example:**
```
I have differential expression p-values from 3 modalities:
- RNA p-values: [0.001, 0.05, 0.3, 0.8]
- Protein p-values: [0.002, 0.04, 0.25, 0.75]
- Phospho p-values: [0.01, 0.06, 0.28, 0.82]

And corresponding log2 fold changes:
- RNA log2FC: [2.5, 1.2, -0.3, 0.1]
- Protein log2FC: [1.8, 1.5, -0.2, 0.15]
- Phospho log2FC: [1.2, 0.8, -0.4, -0.05]

Please combine these using Stouffer's method with directionality.
```

### 3. run_halla_analysis

Run HAllA hierarchical all-against-all association testing (R-based, requires rpy2).

**Parameters:**
- `data_path` (required): Path to integrated multi-omics data (from integrate_omics_data)
- `modality1` (required): First modality ("rna", "protein", or "phospho")
- `modality2` (required): Second modality ("rna", "protein", or "phospho")
- `fdr_threshold` (default: 0.05): FDR cutoff for significant associations
- `method` (default: "spearman"): Correlation method - "spearman", "pearson", or "mi"

**Returns:**
- `associations`: List of significant feature pairs
- `clusters`: Hierarchical cluster assignments
- `statistics`: p-values, q-values, effect sizes
- `plot_data`: Data for association heatmap

**Example:**
```
Using the integrated data at /workspace/cache/integrated_data.pkl,
please run HAllA to find associations between RNA and Protein modalities.
Use Spearman correlation and FDR threshold of 0.05.
```

### 4. create_multiomics_heatmap

Create integrated heatmap visualization across multiple omics modalities.

**Parameters:**
- `data_path` (required): Path to integrated multi-omics data
- `features` (optional): List of feature names to include (default: top significant features)
- `cluster_rows` (default: True): Apply hierarchical clustering to rows (features)
- `cluster_cols` (default: True): Apply hierarchical clustering to columns (samples)
- `output_path` (optional): Path to save plot (PNG or PDF)

**Returns:**
- `plot_path`: Path to saved visualization
- `cluster_info`: Cluster assignments
- `figure_data`: Plotly JSON for interactive plot

**Example:**
```
Please create a multi-omics heatmap for the integrated data showing
TP53, MYC, EGFR, and KRAS across all modalities. Cluster both rows
and columns and save to /workspace/plots/multiomics_heatmap.png.
```

### 5. run_multiomics_pca

Run Principal Component Analysis on integrated multi-omics data.

**Parameters:**
- `data_path` (required): Path to integrated multi-omics data
- `modalities` (optional): List of modalities to include (default: all available)
- `n_components` (default: 3): Number of principal components to compute
- `scale_features` (default: True): Apply feature scaling before PCA
- `output_path` (optional): Path to save 2D and 3D PCA plots

**Returns:**
- `variance_explained`: Fraction of variance per component
- `loadings`: Feature loadings on each PC
- `sample_coordinates`: PC coordinates for each sample
- `plot_path`: Path to saved visualization

**Example:**
```
Run PCA on the integrated RNA + Protein data to visualize sample clustering.
Compute 3 components and color samples by treatment response.
```

## Resources

### multiomics://config

Get current server configuration and analysis parameters.

**Example:**
```
What are the current multi-omics server settings?
```

### multiomics://example-data

Get information about example datasets and data format requirements.

**Example:**
```
What data formats does the multi-omics server expect?
```

## Installation

```bash
cd servers/mcp-multiomics
python -m venv venv
source venv/bin/activate  # On Windows: venv\\Scripts\\activate
pip install -e ".[dev]"
```

## Configuration

Add to Claude Desktop config (`~/Library/Application Support/Claude/claude_desktop_config.json`):

```json
{
  "mcpServers": {
    "multiomics": {
      "command": "/absolute/path/to/spatial-mcp/servers/mcp-multiomics/venv/bin/python",
      "args": ["-m", "mcp_multiomics"],
      "cwd": "/absolute/path/to/spatial-mcp/servers/mcp-multiomics",
      "env": {
        "PYTHONPATH": "/absolute/path/to/spatial-mcp/servers/mcp-multiomics/src",
        "MULTIOMICS_DATA_DIR": "/absolute/path/to/data/multiomics",
        "MULTIOMICS_CACHE_DIR": "/absolute/path/to/cache/multiomics",
        "MULTIOMICS_DRY_RUN": "true"
      }
    }
  }
}
```

**Important:** Use the full path to the venv Python executable, not just `python`.

For a complete working config with all servers, see `../../configs/claude_desktop_config.json`.

## DRY_RUN Mode

Set `MULTIOMICS_DRY_RUN=true` to use realistic mock data without requiring:
- R installation and rpy2
- HAllA R package
- Large data downloads
- Real file I/O operations

Perfect for development, testing, and demonstrations!

## Data Format Requirements

### RNA Expression Data (CSV/TSV)
```
Gene,Sample_01,Sample_02,Sample_03,...
TP53,1234.5,2341.2,987.3,...
MYC,5432.1,4321.9,3210.4,...
EGFR,876.4,1023.8,1156.2,...
```
- First column: Gene symbols
- Other columns: FPKM, TPM, or normalized counts per sample
- Missing values: Leave empty or use NA/NaN

### Protein Abundance (TMT/Label-Free)
```
Protein,Sample_01,Sample_02,Sample_03,...
TP53,45.2,67.8,34.1,...
MYC,123.4,145.6,98.7,...
```
- First column: Protein names (preferably matching gene symbols)
- Other columns: TMT intensities or LFQ values

### Phosphorylation Data
```
Phosphosite,Sample_01,Sample_02,Sample_03,...
TP53_S15,12.3,23.4,8.7,...
AKT1_S473,89.2,102.3,76.5,...
```
- First column: Phosphosite notation (PROTEIN_RESIDUE_POSITION)
- Other columns: Phosphorylation intensities

### Sample Metadata
```
Sample,Response,Batch,Age,Sex
Sample_01,Resistant,1,65,M
Sample_02,Sensitive,1,58,F
Sample_03,Resistant,2,72,M
```
- Must include `Sample` column matching sample names in omics data
- Recommended columns: Treatment response, batch, clinical variables

## Example Workflows

### Complete Multi-Omics Analysis

```
Claude, I need to analyze PDX treatment response data:

1. Integrate RNA (/data/rna.csv), Protein (/data/protein.csv),
   and Phospho (/data/phospho.csv) with metadata (/data/metadata.csv)

2. For significant features from each modality's differential analysis,
   combine p-values using Stouffer's method with directionality:
   - RNA: p=[0.001, 0.002, 0.05], log2FC=[2.5, 1.8, 1.2]
   - Protein: p=[0.005, 0.01, 0.03], log2FC=[2.0, 1.6, 1.1]
   - Phospho: p=[0.01, 0.015, 0.04], log2FC=[1.8, 1.4, 0.9]

3. Create a heatmap showing the top 50 most significant features

4. Run PCA to visualize resistant vs sensitive samples
```

### HAllA Association Discovery

```
After integrating my multi-omics data, please run HAllA to discover
RNA-Protein associations in treatment-resistant PDX models. I'm interested
in post-transcriptional regulation patterns.
```

### Meta-Analysis Only

```
I have pre-computed differential expression statistics from 3 omics layers.
Please use Stouffer's method to identify features that are consistently
dysregulated across all modalities. Weight RNA 2x higher than other modalities.
```

## Testing

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=src --cov-report=html

# Run specific test suite
pytest tests/test_stouffer.py -v
pytest tests/test_integration.py -v
```

## Integration with Other Servers

**Works with:**
- `mcp-fgbio`: Reference genome data for gene annotations
- `mcp-spatialtools`: Spatial transcriptomics analysis and clustering
- `mcp-tcga`: Compare PDX results to TCGA cohorts
- `mcp-huggingface`: Apply ML models to multi-omics features

## Technical Details

- **Framework**: FastMCP (Python 3.11+)
- **Statistical Methods**:
  - Stouffer's Z-score meta-analysis (scipy)
  - FDR correction (Benjamini-Hochberg)
  - Z-score normalization
- **R Integration**: rpy2 for HAllA (optional)
- **Data Handling**: pandas, numpy
- **Visualization**: matplotlib, seaborn, plotly
- **Performance**:
  - Mock mode: <100ms per tool call
  - Integration: 1-5 seconds for 1000 features
  - Stouffer's: <1 second for 10,000 features

## References

- **Stouffer's Method**: Whitlock MC (2005). Combining probability from independent tests: the weighted Z-method is superior to Fisher's approach. *J Evol Biol* 18(5):1368-73.
- **HAllA**: Rahnavard et al. (2017). High-sensitivity pattern discovery in large, paired multiomic datasets. *Bioinformatics* 33(14):i81-i89.
- **Multi-Omics Integration**: Menyh\u00e1rt et al. (2023). Multi-omics approaches in cancer research with applications in tumor subtyping, prognosis, and diagnosis. *Comput Struct Biotechnol J* 21:1864-1883.

## Troubleshooting

### Server won't start
- Check Python version: Must be 3.11+
- Verify dependencies: Run `pip install -e .`
- Check environment variables: Ensure paths exist

### Integration fails
- Verify file formats: First column must be feature names
- Check sample names: Must match exactly across files
- Enable DRY_RUN: Set `MULTIOMICS_DRY_RUN=true` for testing

### HAllA not available
- Install R and rpy2: `pip install rpy2` (Unix/Mac only)
- Install HAllA R package: See HAllA documentation
- Use DRY_RUN mode: Works without R installation

## License

See the main repository LICENSE file.

## Support

- **Issues**: https://github.com/spatial-mcp/spatial-mcp/issues
- **Documentation**: See main repository docs/
- **MCP Specification**: https://modelcontextprotocol.io/

---

**Built for the Precision Medicine MCP suite** (key component of the [PatientOne](../../architecture/patient-one/README.md) precision medicine workflow)
