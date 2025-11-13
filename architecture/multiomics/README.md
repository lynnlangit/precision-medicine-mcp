# Multiomics MCP Demonstration POC - Architecture Document

**Version:** 1.0  
**Date:** November 11, 2025  
**Status:** Draft for Review

---
## Executive Summary

The **mcp-multiomics** server is the highest-priority addition for PDX multi-omics analysis. It enables integration of RNA, protein, and phosphorylation data through statistical methods like HAllA and Stouffer's meta-analysis.

**Estimated Development Time:** 1-2 weeks  
**Complexity:** Medium-High (R integration required)  
**Dependencies:** Python 3.11+, R 4.3+, rpy2, HAllA R package

---

## Server Architecture

### Server Name
`mcp-multiomics`

### Purpose
Multi-omics data integration and meta-analysis for PDX models

### Key Differentiators
- First MCP server specifically for multi-omics integration
- Implements Stouffer's method for p-value combination
- Integrates with R's HAllA package via rpy2
- Handles missing data across modalities

---

## Core Functionality

### 1. Tools (5 total)

#### Tool 1: `integrate_omics_data`
**Purpose:** Combine RNA, protein, and phospho data into unified format

**Inputs:**
```json
{
  "rna_file": "string (path to RNA FPKM matrix)",
  "protein_file": "string (path to protein TMT matrix)", 
  "phospho_file": "string (optional, path to phospho matrix)",
  "metadata_file": "string (path to sample metadata)",
  "modalities": ["RNA", "PROTEIN", "PHOSPHO"],
  "batch_correct": "boolean (default: true)",
  "normalization": "string (log2, zscore, quantile)"
}
```

**Outputs:**
```json
{
  "integrated_data": {
    "path": "string",
    "features": "integer",
    "samples": "integer",
    "modalities": ["RNA", "PROTEIN", "PHOSPHO"]
  },
  "qc_metrics": {
    "missing_values_pct": "float",
    "features_per_modality": "dict",
    "batch_correction_applied": "boolean"
  }
}
```

**Processing Steps:**
1. Load expression matrices (CSV/TSV format)
2. Load metadata (SampleID, Status, Group, Treatment)
3. Align samples across modalities
4. Handle missing values (imputation or filtering)
5. Apply batch correction (ComBat if needed)
6. Normalize within modality
7. Save integrated dataset

---

#### Tool 2: `run_halla_analysis`
**Purpose:** Hierarchical all-against-all association testing

**Inputs:**
```json
{
  "integrated_data": "string (path to integrated dataset)",
  "phenotype": "string (Status, Group, or continuous variable)",
  "modality": "string (RNA, PROTEIN, PHOSPHO, or ALL)",
  "regression_type": "string (logistic, linear)",
  "fdr_threshold": "float (default: 0.05)",
  "correlation_threshold": "float (default: 0.4)"
}
```

**Outputs:**
```json
{
  "significant_features": [
    {
      "feature_id": "string",
      "feature_name": "string", 
      "data_modality": "string",
      "association": "string",
      "correlation_value": "float",
      "p_value": "float",
      "adjusted_p_value": "float",
      "log2fc": "float"
    }
  ],
  "summary": {
    "total_features_tested": "integer",
    "significant_count": "integer",
    "modality_breakdown": "dict"
  }
}
```

**Processing Steps:**
1. Call R's HAllA package via rpy2
2. Configure regression type based on phenotype
3. Run association testing
4. Adjust p-values (FDR)
5. Calculate log2FC for significant features
6. Return structured results

---

#### Tool 3: `calculate_stouffer_meta`
**Purpose:** Stouffer's inverse normal method for p-value integration

**Inputs:**
```json
{
  "halla_results_rna": "string (path to RNA results)",
  "halla_results_protein": "string (path to protein results)",
  "halla_results_phospho": "string (optional, path to phospho results)",
  "weight_by_sample_size": "boolean (default: false)",
  "min_modalities": "integer (default: 2)"
}
```

**Outputs:**
```json
{
  "meta_analysis_results": [
    {
      "gene_id": "string",
      "gene_name": "string",
      "modalities_detected": ["RNA", "PROTEIN"],
      "meta_z_score": "float",
      "meta_p_value": "float",
      "adjusted_meta_p_value": "float",
      "direction": "string (up, down, mixed)",
      "individual_results": {
        "RNA": {"z": "float", "p": "float", "log2fc": "float"},
        "PROTEIN": {"z": "float", "p": "float", "log2fc": "float"}
      }
    }
  ],
  "summary": {
    "genes_tested": "integer",
    "significant_genes": "integer",
    "multi_modal_genes": "integer"
  }
}
```

**Processing Steps:**
1. Load results from each modality
2. Map features to common gene IDs
3. Filter to genes present in ≥ min_modalities
4. Convert p-values to Z-scores
5. Calculate Stouffer's combined Z: Z_meta = Σ(Z_i) / sqrt(n)
6. Convert back to p-values
7. Adjust for multiple testing
8. Determine directionality

**Mathematical Formula:**
```
Z_meta = (Z_RNA + Z_protein + Z_phospho) / sqrt(3)
p_meta = 2 * Φ(-|Z_meta|)  # where Φ is standard normal CDF
```

---

#### Tool 4: `create_multiomics_heatmap`
**Purpose:** Visualize integrated multi-omics results

**Inputs:**
```json
{
  "integrated_data": "string (path)",
  "meta_results": "string (path to Stouffer results)",
  "features": "string (all, significant, top_n)",
  "top_n": "integer (default: 50)",
  "cluster_rows": "boolean (default: true)",
  "cluster_cols": "boolean (default: true)",
  "annotation_columns": ["Status", "Group"],
  "output_format": "string (png, pdf, html)"
}
```

**Outputs:**
```json
{
  "heatmap_path": "string",
  "feature_order": ["gene1", "gene2", ...],
  "sample_order": ["sample1", "sample2", ...],
  "cluster_assignments": {
    "features": "dict",
    "samples": "dict"
  }
}
```

**Visualization Features:**
- Row annotation: Data modality (RNA/PROTEIN/PHOSPHO)
- Column annotation: Status (Resistant/Sensitive), Group
- Color scale: Z-score normalized expression
- Hierarchical clustering with dendrograms
- Log2FC annotation on side

---

#### Tool 5: `run_multiomics_pca`
**Purpose:** Dimensionality reduction across modalities

**Inputs:**
```json
{
  "integrated_data": "string (path)",
  "features": "string (all, significant, variable)",
  "n_components": "integer (default: 2)",
  "color_by": "string (Status, Group, Treatment)",
  "scale_data": "boolean (default: true)"
}
```

**Outputs:**
```json
{
  "pca_results": {
    "explained_variance": [0.45, 0.23, ...],
    "pc_coordinates": "array",
    "loadings": "array"
  },
  "plot_path": "string",
  "top_contributing_features": [
    {"feature": "gene1", "pc1_loading": 0.23, "pc2_loading": -0.15}
  ]
}
```

---

## 2. Resources (4 total)

### Resource 1: `multiomics://integrated/all_features`
**Type:** Read-only data resource  
**Content:** Complete integrated dataset across all modalities  
**Format:** Parquet or HDF5 for efficient access

### Resource 2: `multiomics://meta_z/stouffer`
**Type:** Read-only results resource  
**Content:** Stouffer's meta-analysis results  
**Format:** JSON with hierarchical structure

### Resource 3: `multiomics://significant/combined`
**Type:** Read-only results resource  
**Content:** Significant features from meta-analysis  
**Format:** JSON array

### Resource 4: `multiomics://metadata/samples`
**Type:** Read-only metadata resource  
**Content:** PDX sample annotations  
**Format:** JSON with sample metadata

---

## 3. Prompts (3 total)

### Prompt 1: `analyze_resistance_signature`
**Template:**
```
"Analyze the multi-omics data comparing {phenotype} samples. 
Load RNA, protein, and phospho data, run HAllA on each modality, 
integrate using Stouffer's method, and identify genes with 
meta p-value < {threshold}. Focus on genes detected in at least 
{min_modalities} modalities."
```

### Prompt 2: `create_integrated_visualization`
**Template:**
```
"Create a comprehensive heatmap showing the top {n} significant 
features from the integrated multi-omics analysis. Cluster samples 
by {phenotype} and annotate with data modality. Include log2FC 
values."
```

### Prompt 3: `explore_multimodal_concordance`
**Template:**
```
"For genes detected in multiple modalities, assess concordance: 
Are changes in RNA expression reflected in protein levels? 
Are phosphorylation changes correlated with protein abundance? 
Show discordant cases."
```

---

## Technical Implementation

### Directory Structure

```
servers/mcp-multiomics/
├── src/
│   ├── mcp_multiomics/
│   │   ├── __init__.py
│   │   ├── server.py              # FastMCP server definition
│   │   ├── tools/
│   │   │   ├── __init__.py
│   │   │   ├── integration.py     # integrate_omics_data
│   │   │   ├── halla.py           # run_halla_analysis
│   │   │   ├── stouffer.py        # calculate_stouffer_meta
│   │   │   ├── visualization.py   # heatmaps and PCA
│   │   │   └── utils.py           # helper functions
│   │   ├── resources/
│   │   │   ├── __init__.py
│   │   │   └── data_resources.py  # Resource providers
│   │   └── r_interface/
│   │       ├── __init__.py
│   │       ├── halla_wrapper.py   # rpy2 wrapper for HAllA
│   │       └── install_packages.R # R package installation
├── tests/
│   ├── test_integration.py
│   ├── test_halla.py
│   ├── test_stouffer.py
│   ├── fixtures/
│   │   ├── sample_rna.csv
│   │   ├── sample_protein.csv
│   │   ├── sample_phospho.csv
│   │   └── sample_metadata.csv
├── pyproject.toml
├── README.md
└── .env.example
```

---

### Key Dependencies

**Python (pyproject.toml):**
```toml
[project]
name = "mcp-multiomics"
version = "0.1.0"
dependencies = [
    "fastmcp>=0.2.0",
    "pydantic>=2.0.0",
    "pandas>=2.0.0",
    "numpy>=1.24.0",
    "scipy>=1.11.0",
    "scikit-learn>=1.3.0",
    "rpy2>=3.5.0",
    "seaborn>=0.12.0",
    "matplotlib>=3.7.0",
    "plotly>=5.17.0",
    "tables>=3.9.0",  # For HDF5
    "pyarrow>=14.0.0",  # For Parquet
    "statsmodels>=0.14.0",
]

[project.optional-dependencies]
dev = [
    "pytest>=7.4.0",
    "pytest-asyncio>=0.21.0",
    "pytest-cov>=4.1.0",
    "black>=23.0.0",
    "ruff>=0.1.0",
    "mypy>=1.7.0",
]
```

**R Packages (install_packages.R):**
```r
# Install required R packages
install.packages("BiocManager")
BiocManager::install(c(
  "HAllA",
  "sva",  # For ComBat batch correction
  "limma",
  "ComplexHeatmap",
  "DESeq2"
))
```

---

### Configuration (config.py)

```python
from pydantic import BaseSettings
from pathlib import Path

class MultiOmicsConfig(BaseSettings):
    # Data directories
    data_dir: Path = Path("/workspace/multiomics_data")
    results_dir: Path = Path("/workspace/multiomics_results")
    cache_dir: Path = Path("/workspace/multiomics_cache")
    
    # R configuration
    r_home: str = "/usr/lib/R"
    r_libs: str = "/usr/local/lib/R/site-library"
    
    # Analysis parameters
    default_fdr_threshold: float = 0.05
    default_correlation_threshold: float = 0.4
    min_modalities_for_meta: int = 2
    
    # Resource limits
    max_features: int = 50000
    max_samples: int = 1000
    timeout_seconds: int = 600
    
    # Logging
    log_level: str = "INFO"
    
    class Config:
        env_prefix = "MULTIOMICS_"
        env_file = ".env"
```

---

### Core Implementation: Stouffer's Method

**File: `src/mcp_multiomics/tools/stouffer.py`**

```python
import numpy as np
import pandas as pd
from scipy import stats
from typing import List, Dict, Optional
from pathlib import Path

class StoufferMetaAnalysis:
    """
    Implements Stouffer's inverse normal method for p-value integration.
    """
    
    def __init__(self, min_modalities: int = 2):
        self.min_modalities = min_modalities
    
    def load_modality_results(
        self, 
        rna_path: Optional[Path] = None,
        protein_path: Optional[Path] = None,
        phospho_path: Optional[Path] = None
    ) -> Dict[str, pd.DataFrame]:
        """Load HAllA results from each modality."""
        results = {}
        
        if rna_path:
            results['RNA'] = pd.read_csv(rna_path)
        if protein_path:
            results['PROTEIN'] = pd.read_csv(protein_path)
        if phospho_path:
            results['PHOSPHO'] = pd.read_csv(phospho_path)
            
        return results
    
    def p_to_z(self, p_values: np.ndarray) -> np.ndarray:
        """Convert p-values to Z-scores (standard normal)."""
        # Handle edge cases
        p_values = np.clip(p_values, 1e-300, 1 - 1e-15)
        
        # Convert to two-tailed Z-scores
        z_scores = stats.norm.ppf(1 - p_values / 2)
        return z_scores
    
    def z_to_p(self, z_scores: np.ndarray) -> np.ndarray:
        """Convert Z-scores to two-tailed p-values."""
        p_values = 2 * (1 - stats.norm.cdf(np.abs(z_scores)))
        return p_values
    
    def combine_z_scores(
        self, 
        z_scores: np.ndarray, 
        weights: Optional[np.ndarray] = None
    ) -> float:
        """
        Stouffer's method: Z_meta = Σ(w_i * Z_i) / sqrt(Σ(w_i^2))
        If no weights provided, uses equal weights.
        """
        if weights is None:
            weights = np.ones(len(z_scores))
        
        numerator = np.sum(weights * z_scores)
        denominator = np.sqrt(np.sum(weights ** 2))
        
        z_meta = numerator / denominator
        return z_meta
    
    def integrate_modalities(
        self, 
        modality_results: Dict[str, pd.DataFrame]
    ) -> pd.DataFrame:
        """
        Integrate results across modalities using Stouffer's method.
        """
        # Create unified gene list
        all_genes = set()
        for modality, df in modality_results.items():
            if 'GeneId' in df.columns:
                all_genes.update(df['GeneId'].unique())
            elif 'gene_id' in df.columns:
                all_genes.update(df['gene_id'].unique())
        
        results = []
        
        for gene in all_genes:
            gene_data = {}
            modalities_present = []
            z_scores = []
            p_values = []
            log2fcs = []
            
            # Collect data from each modality
            for modality, df in modality_results.items():
                gene_col = 'GeneId' if 'GeneId' in df.columns else 'gene_id'
                gene_rows = df[df[gene_col] == gene]
                
                if not gene_rows.empty:
                    row = gene_rows.iloc[0]
                    modalities_present.append(modality)
                    
                    p_val = row['p_value']
                    z_score = self.p_to_z(np.array([p_val]))[0]
                    
                    # Adjust sign based on log2FC direction
                    if 'log2FC' in row:
                        if row['log2FC'] < 0:
                            z_score = -z_score
                    
                    z_scores.append(z_score)
                    p_values.append(p_val)
                    log2fcs.append(row.get('log2FC', np.nan))
            
            # Only process if present in minimum number of modalities
            if len(modalities_present) >= self.min_modalities:
                # Calculate meta Z-score
                z_meta = self.combine_z_scores(np.array(z_scores))
                p_meta = self.z_to_p(np.array([z_meta]))[0]
                
                # Determine direction
                avg_log2fc = np.mean([fc for fc in log2fcs if not np.isnan(fc)])
                direction = "up" if avg_log2fc > 0 else "down"
                
                # Check concordance
                if len(set(np.sign(log2fcs))) > 1:
                    direction = "mixed"
                
                results.append({
                    'gene_id': gene,
                    'modalities_detected': ','.join(modalities_present),
                    'n_modalities': len(modalities_present),
                    'meta_z_score': z_meta,
                    'meta_p_value': p_meta,
                    'direction': direction,
                    'avg_log2fc': avg_log2fc,
                    **{f'{mod}_z': z for mod, z in zip(modalities_present, z_scores)},
                    **{f'{mod}_p': p for mod, p in zip(modalities_present, p_values)},
                    **{f'{mod}_log2fc': fc for mod, fc in zip(modalities_present, log2fcs)}
                })
        
        meta_df = pd.DataFrame(results)
        
        # FDR correction
        from statsmodels.stats.multitest import multipletests
        _, meta_df['adjusted_meta_p_value'], _, _ = multipletests(
            meta_df['meta_p_value'], 
            method='fdr_bh'
        )
        
        # Sort by significance
        meta_df = meta_df.sort_values('meta_p_value')
        
        return meta_df
```

---

### R Interface: HAllA Wrapper

**File: `src/mcp_multiomics/r_interface/halla_wrapper.py`**

```python
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import pandas as pd
from pathlib import Path
from typing import Literal

# Activate automatic pandas conversion
pandas2ri.activate()

class HAllAWrapper:
    """
    Python wrapper for R's HAllA package using rpy2.
    """
    
    def __init__(self):
        # Import R packages
        try:
            self.halla = importr('HAllA')
            self.base = importr('base')
            self.stats = importr('stats')
        except Exception as e:
            raise RuntimeError(
                f"Failed to import R packages. Ensure HAllA is installed: {e}"
            )
    
    def run_halla(
        self,
        expression_data: pd.DataFrame,
        phenotype_data: pd.DataFrame,
        phenotype_name: str,
        regression_type: Literal['logistic', 'linear'] = 'logistic',
        fdr_threshold: float = 0.05
    ) -> pd.DataFrame:
        """
        Run HAllA analysis via R.
        
        Args:
            expression_data: Genes/proteins × samples matrix
            phenotype_data: Sample metadata with phenotype column
            phenotype_name: Column name in phenotype_data
            regression_type: 'logistic' or 'linear'
            fdr_threshold: FDR cutoff for significance
            
        Returns:
            DataFrame with HAllA results
        """
        # Convert pandas DataFrames to R dataframes
        with ro.default_converter + pandas2ri.converter:
            r_expr = ro.conversion.py2rpy(expression_data)
            r_pheno = ro.conversion.py2rpy(phenotype_data)
        
        # Run HAllA
        try:
            # This is pseudocode - actual HAllA API may differ
            r_results = self.halla.halla(
                X=r_expr,
                Y=r_pheno[phenotype_name],
                method=regression_type,
                fdr=fdr_threshold
            )
            
            # Convert results back to pandas
            with ro.default_converter + pandas2ri.converter:
                results_df = ro.conversion.rpy2py(r_results)
            
            return results_df
            
        except Exception as e:
            raise RuntimeError(f"HAllA analysis failed: {e}")
```

---

### FastMCP Server Definition

**File: `src/mcp_multiomics/server.py`**

```python
from fastmcp import FastMCP
from pathlib import Path
import logging
from .tools.integration import integrate_omics_data
from .tools.halla import run_halla_analysis
from .tools.stouffer import calculate_stouffer_meta
from .tools.visualization import create_multiomics_heatmap, run_multiomics_pca
from .config import MultiOmicsConfig

# Initialize config
config = MultiOmicsConfig()

# Setup logging
logging.basicConfig(
    level=config.log_level,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Create MCP server
mcp = FastMCP("multiomics", dependencies=["pandas", "numpy", "scipy", "rpy2"])

@mcp.tool()
async def integrate_omics_data_tool(
    rna_file: str,
    protein_file: str,
    phospho_file: str | None = None,
    metadata_file: str,
    batch_correct: bool = True,
    normalization: str = "log2"
) -> dict:
    """
    Integrate RNA, protein, and phosphorylation data.
    
    This tool combines multiple omics modalities into a unified dataset
    for downstream meta-analysis.
    """
    logger.info(f"Integrating omics data: RNA={rna_file}, Protein={protein_file}")
    
    result = await integrate_omics_data(
        rna_file=Path(rna_file),
        protein_file=Path(protein_file),
        phospho_file=Path(phospho_file) if phospho_file else None,
        metadata_file=Path(metadata_file),
        batch_correct=batch_correct,
        normalization=normalization,
        config=config
    )
    
    return result

@mcp.tool()
async def run_halla_tool(
    integrated_data: str,
    phenotype: str,
    modality: str = "ALL",
    regression_type: str = "logistic",
    fdr_threshold: float = 0.05
) -> dict:
    """
    Run HAllA hierarchical all-against-all association testing.
    
    Identifies features (genes/proteins) significantly associated
    with the phenotype of interest using regression analysis.
    """
    logger.info(f"Running HAllA analysis for {modality} modality")
    
    result = await run_halla_analysis(
        integrated_data=Path(integrated_data),
        phenotype=phenotype,
        modality=modality,
        regression_type=regression_type,
        fdr_threshold=fdr_threshold,
        config=config
    )
    
    return result

@mcp.tool()
async def calculate_stouffer_tool(
    halla_results_rna: str,
    halla_results_protein: str,
    halla_results_phospho: str | None = None,
    min_modalities: int = 2
) -> dict:
    """
    Integrate p-values across modalities using Stouffer's method.
    
    Combines statistical evidence from multiple data modalities to identify
    genes with consistent effects across RNA, protein, and phosphorylation.
    """
    logger.info("Calculating Stouffer's meta-analysis")
    
    result = await calculate_stouffer_meta(
        rna_path=Path(halla_results_rna),
        protein_path=Path(halla_results_protein),
        phospho_path=Path(halla_results_phospho) if halla_results_phospho else None,
        min_modalities=min_modalities,
        config=config
    )
    
    return result

@mcp.tool()
async def create_heatmap_tool(
    integrated_data: str,
    meta_results: str,
    features: str = "significant",
    top_n: int = 50,
    output_format: str = "png"
) -> dict:
    """
    Create multi-omics integrated heatmap.
    
    Visualizes expression patterns across samples and modalities,
    with hierarchical clustering and metadata annotations.
    """
    logger.info(f"Creating heatmap with {features} features")
    
    result = await create_multiomics_heatmap(
        integrated_data=Path(integrated_data),
        meta_results=Path(meta_results),
        features=features,
        top_n=top_n,
        output_format=output_format,
        config=config
    )
    
    return result

@mcp.tool()
async def run_pca_tool(
    integrated_data: str,
    features: str = "significant",
    n_components: int = 2,
    color_by: str = "Status"
) -> dict:
    """
    Perform PCA on integrated multi-omics data.
    
    Reduces dimensionality and visualizes sample clustering
    based on multi-modal molecular profiles.
    """
    logger.info(f"Running PCA with {n_components} components")
    
    result = await run_multiomics_pca(
        integrated_data=Path(integrated_data),
        features=features,
        n_components=n_components,
        color_by=color_by,
        config=config
    )
    
    return result

# Resources
@mcp.resource("multiomics://integrated/all_features")
async def get_integrated_data() -> str:
    """
    Access the complete integrated multi-omics dataset.
    
    Returns the most recently created integrated dataset
    containing RNA, protein, and phosphorylation data.
    """
    # Implementation to return latest integrated dataset
    pass

@mcp.resource("multiomics://meta_z/stouffer")
async def get_stouffer_results() -> str:
    """
    Access Stouffer's meta-analysis results.
    
    Returns genes ranked by meta-Z score with integrated
    p-values across modalities.
    """
    # Implementation to return latest Stouffer results
    pass

# Prompts
@mcp.prompt()
async def analyze_resistance_signature(
    phenotype: str = "Status",
    threshold: float = 0.05,
    min_modalities: int = 2
) -> str:
    """
    Template for analyzing treatment resistance signatures.
    """
    return f"""
    Analyze the multi-omics data comparing {phenotype} samples.
    
    Steps:
    1. Load RNA, protein, and phosphorylation data
    2. Run HAllA on each modality separately
    3. Integrate results using Stouffer's method
    4. Identify genes with meta p-value < {threshold}
    5. Focus on genes detected in at least {min_modalities} modalities
    6. Create visualization showing top significant genes
    
    Generate a summary report with:
    - Number of significant genes per modality
    - Number of multi-modal significant genes
    - Direction of change (up/down/mixed)
    - Heatmap and PCA plot
    """

if __name__ == "__main__":
    mcp.run()
```

---

## Testing Strategy

### Unit Tests

**File: `tests/test_stouffer.py`**

```python
import pytest
import numpy as np
import pandas as pd
from mcp_multiomics.tools.stouffer import StoufferMetaAnalysis

def test_p_to_z_conversion():
    """Test p-value to Z-score conversion."""
    stouffer = StoufferMetaAnalysis()
    
    p_values = np.array([0.05, 0.01, 0.001])
    z_scores = stouffer.p_to_z(p_values)
    
    # Check conversion is correct
    assert z_scores[0] < z_scores[1] < z_scores[2]
    assert z_scores[1] > 2.326  # Approximately for p=0.01

def test_stouffer_combination():
    """Test Stouffer's method for combining Z-scores."""
    stouffer = StoufferMetaAnalysis()
    
    # Three modalities with same effect
    z_scores = np.array([2.0, 2.0, 2.0])
    z_meta = stouffer.combine_z_scores(z_scores)
    
    # Meta Z should be larger due to combination
    assert z_meta > 2.0
    assert np.isclose(z_meta, 2.0 * np.sqrt(3), rtol=0.01)

def test_integration_with_missing_modality():
    """Test integration when gene is missing in one modality."""
    stouffer = StoufferMetaAnalysis(min_modalities=2)
    
    # Mock data: gene present in 2 of 3 modalities
    rna_df = pd.DataFrame({
        'GeneId': ['GENE1', 'GENE2'],
        'p_value': [0.01, 0.05],
        'log2FC': [1.5, -2.0]
    })
    
    protein_df = pd.DataFrame({
        'gene_id': ['GENE1', 'GENE3'],
        'p_value': [0.02, 0.03],
        'log2FC': [1.2, 0.8]
    })
    
    modality_results = {'RNA': rna_df, 'PROTEIN': protein_df}
    integrated = stouffer.integrate_modalities(modality_results)
    
    # GENE1 should be in results (present in 2 modalities)
    assert 'GENE1' in integrated['gene_id'].values
    
    # GENE2 and GENE3 should not (only in 1 modality)
    assert 'GENE2' not in integrated['gene_id'].values
    assert 'GENE3' not in integrated['gene_id'].values
```

### Integration Test

**File: `tests/test_integration.py`**

```python
import pytest
from pathlib import Path
from mcp_multiomics.server import mcp

@pytest.mark.asyncio
async def test_end_to_end_workflow(tmp_path):
    """Test complete workflow from data loading to meta-analysis."""
    
    # Setup test data paths
    test_data_dir = Path(__file__).parent / "fixtures"
    
    # Step 1: Integrate data
    integrate_result = await mcp.call_tool(
        "integrate_omics_data_tool",
        rna_file=str(test_data_dir / "sample_rna.csv"),
        protein_file=str(test_data_dir / "sample_protein.csv"),
        metadata_file=str(test_data_dir / "sample_metadata.csv"),
        batch_correct=False
    )
    
    assert integrate_result["integrated_data"]["samples"] == 15
    
    # Step 2: Run HAllA
    halla_result = await mcp.call_tool(
        "run_halla_tool",
        integrated_data=integrate_result["integrated_data"]["path"],
        phenotype="Status",
        modality="RNA"
    )
    
    assert "significant_features" in halla_result
    assert len(halla_result["significant_features"]) > 0
```

---

## Documentation

### README.md Template

```markdown
# mcp-multiomics

Multi-omics data integration server for Model Context Protocol (MCP).

## Features

- Integration of RNA, protein, and phosphorylation data
- HAllA (Hierarchical All-against-All) association testing
- Stouffer's meta-analysis for p-value combination
- Multi-modal visualization (heatmaps, PCA)
- Support for PDX model cohort analysis

## Installation

```bash
# Clone and install
git clone https://github.com/yourusername/mcp-multiomics
cd mcp-multiomics
pip install -e .

# Install R dependencies
Rscript src/mcp_multiomics/r_interface/install_packages.R
```

## Quick Start

### 1. Configure Claude Desktop

Add to your `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "multiomics": {
      "command": "python",
      "args": ["-m", "mcp_multiomics.server"],
      "env": {
        "MULTIOMICS_DATA_DIR": "/path/to/data",
        "MULTIOMICS_R_HOME": "/usr/lib/R"
      }
    }
  }
}
```

### 2. Example Usage

```python
# In Claude conversation:
"I have RNA FPKM, protein TMT, and phospho data from 15 PDX samples
(7 resistant, 8 sensitive). Load all three modalities, run HAllA to find
resistance-associated features, integrate using Stouffer's method, and
show me genes significant in multiple modalities."
```

## License

MIT
```

---

## Success Metrics

### Phase 1 Completion Criteria

✅ Can load RNA, protein, phospho CSV files  
✅ Handles missing data across modalities  
✅ Implements Stouffer's method correctly  
✅ Passes unit tests (>80% coverage)  
✅ Integration test succeeds end-to-end  
✅ Generates heatmap and PCA visualizations  
✅ Documentation complete (README, docstrings)  
✅ Claude Desktop connects successfully  

---

## Timeline

| Week | Milestone |
|------|-----------|
| 1 | Setup project, implement data integration, write tests |
| 2 | Implement HAllA wrapper, Stouffer's method, visualizations |
| 3 | Integration testing, documentation, demo with test data |

---

## Next Steps After mcp-multiomics

1. **mcp-pathways** - ActivePathways and GSEA integration
2. **mcp-statsintegration** - Advanced statistical methods
3. **mcp-kinaseenrich** - KEA3 kinase enrichment

The multiomics server provides the foundation that these other servers will build upon.
