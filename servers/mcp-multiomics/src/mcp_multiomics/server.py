"""MCP server for multi-omics PDX data analysis.

Provides tools for:
- Multi-omics data integration (RNA, Protein, Phosphorylation)
- HAllA association testing
- Stouffer's meta-analysis for p-value combination
- Multi-omics visualizations
"""

import logging
from typing import Any, Dict, List, Optional

from fastmcp import FastMCP

from .config import config
from .tools.integration import integrate_omics_data_impl
from .tools.stouffer import calculate_stouffer_meta_impl

# Configure logging
logging.basicConfig(
    level=getattr(logging, config.log_level),
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Initialize FastMCP server
mcp = FastMCP("multiomics")


# ============================================================================
# TOOLS
# ============================================================================


@mcp.tool()
def integrate_omics_data(
    rna_path: str,
    protein_path: Optional[str] = None,
    phospho_path: Optional[str] = None,
    metadata_path: Optional[str] = None,
    normalize: bool = True,
    filter_missing: float = 0.5,
) -> Dict[str, Any]:
    """Integrate multi-omics data from RNA, protein, and phosphorylation datasets.

    Aligns samples across modalities, handles missing data, and performs normalization.

    Args:
        rna_path: Path to RNA expression data (CSV or TSV, genes x samples)
        protein_path: Path to protein abundance data (optional)
        phospho_path: Path to phosphorylation data (optional)
        metadata_path: Path to sample metadata (optional, must include 'Sample' column)
        normalize: Apply Z-score normalization within each modality
        filter_missing: Remove features with >X fraction missing (0.0-1.0)

    Returns:
        Dictionary with:
        - integrated_data: Aligned data matrices per modality
        - common_samples: List of samples present in all modalities
        - feature_counts: Number of features per modality
        - metadata: Sample metadata if provided
        - qc_metrics: Quality control statistics

    Example:
        ```
        result = integrate_omics_data(
            rna_path="/data/rna_fpkm.csv",
            protein_path="/data/protein_tmt.csv",
            phospho_path="/data/phospho.csv",
            metadata_path="/data/metadata.csv",
            normalize=True,
            filter_missing=0.5
        )
        # Returns aligned data for 15 samples across 3 modalities
        ```
    """
    logger.info(f"integrate_omics_data called with rna_path={rna_path}")

    if config.dry_run:
        # Mock response for testing
        return {
            "integrated_data": {
                "rna": {"shape": [1000, 15], "path": rna_path},
                "protein": {"shape": [500, 15], "path": protein_path} if protein_path else None,
                "phospho": {"shape": [300, 15], "path": phospho_path} if phospho_path else None,
            },
            "common_samples": [f"Sample_{i:02d}" for i in range(1, 16)],
            "feature_counts": {
                "rna": 1000,
                "protein": 500 if protein_path else 0,
                "phospho": 300 if phospho_path else 0,
            },
            "metadata": {
                "samples": 15,
                "treatment_resistant": 7,
                "treatment_sensitive": 8,
            } if metadata_path else None,
            "qc_metrics": {
                "normalization": "z-score" if normalize else "none",
                "missing_threshold": filter_missing,
                "samples_filtered": 0,
                "features_filtered": {"rna": 50, "protein": 20, "phospho": 15},
            },
            "status": "success (DRY_RUN mode)",
        }

    # Real implementation
    return integrate_omics_data_impl(
        rna_path=rna_path,
        protein_path=protein_path,
        phospho_path=phospho_path,
        metadata_path=metadata_path,
        normalize=normalize,
        filter_missing=filter_missing,
    )


@mcp.tool()
def run_halla_analysis(
    data_path: str,
    modality1: str,
    modality2: str,
    fdr_threshold: float = 0.05,
    method: str = "spearman",
) -> Dict[str, Any]:
    """Run HAllA hierarchical all-against-all association testing.

    Tests associations between features from two omics modalities using
    hierarchical clustering and statistical testing.

    Args:
        data_path: Path to integrated multi-omics data (from integrate_omics_data)
        modality1: First modality ("rna", "protein", or "phospho")
        modality2: Second modality ("rna", "protein", or "phospho")
        fdr_threshold: FDR cutoff for significant associations (default: 0.05)
        method: Correlation method - "spearman", "pearson", or "mi" (mutual information)

    Returns:
        Dictionary with:
        - associations: List of significant feature pairs
        - clusters: Hierarchical cluster assignments
        - statistics: p-values, q-values, effect sizes
        - plot_data: Data for association heatmap

    Example:
        ```
        result = run_halla_analysis(
            data_path="/workspace/cache/integrated_data.pkl",
            modality1="rna",
            modality2="protein",
            fdr_threshold=0.05,
            method="spearman"
        )
        # Returns 47 significant RNA-Protein associations
        ```
    """
    logger.info(f"run_halla_analysis called: {modality1} vs {modality2}")

    if config.dry_run:
        # Mock response
        return {
            "associations": [
                {
                    "feature1": f"{modality1}_gene_{i}",
                    "feature2": f"{modality2}_protein_{i}",
                    "correlation": 0.75 + (i * 0.01),
                    "p_value": 0.001 / (i + 1),
                    "q_value": 0.01 / (i + 1),
                }
                for i in range(1, 48)
            ],
            "clusters": {
                "modality1_clusters": 5,
                "modality2_clusters": 4,
                "total_associations": 47,
            },
            "statistics": {
                "method": method,
                "fdr_threshold": fdr_threshold,
                "significant_pairs": 47,
                "total_tests": 500000,
            },
            "plot_data": {
                "heatmap_ready": True,
                "dendrogram_ready": True,
            },
            "status": "success (DRY_RUN mode)",
        }

    # TODO: Implement R interface with rpy2
    raise NotImplementedError("HAllA analysis not yet implemented")


@mcp.tool()
def calculate_stouffer_meta(
    p_values_dict: Dict[str, List[float]],
    effect_sizes_dict: Optional[Dict[str, List[float]]] = None,
    weights: Optional[Dict[str, float]] = None,
    use_directionality: bool = True,
) -> Dict[str, Any]:
    """Combine p-values across omics modalities using Stouffer's Z-score method.

    Performs meta-analysis by converting p-values to Z-scores, combining them,
    and computing meta p-values. Optionally incorporates effect size directionality.

    Args:
        p_values_dict: Dict of modality -> list of p-values (one per feature)
        effect_sizes_dict: Dict of modality -> list of effect sizes (log2FC or correlation)
        weights: Dict of modality -> weight (default: equal weights or by sample size)
        use_directionality: Incorporate effect size sign into Z-scores

    Returns:
        Dictionary with:
        - meta_p_values: Combined p-values per feature
        - meta_z_scores: Combined Z-scores
        - significant_features: Features passing FDR threshold
        - statistics: Summary statistics and diagnostics

    Example:
        ```
        result = calculate_stouffer_meta(
            p_values_dict={
                "rna": [0.001, 0.05, 0.3],
                "protein": [0.002, 0.04, 0.25],
                "phospho": [0.01, 0.06, 0.28]
            },
            effect_sizes_dict={
                "rna": [2.5, 1.2, -0.3],
                "protein": [1.8, 1.5, -0.2],
                "phospho": [1.2, 0.8, -0.4]
            },
            use_directionality=True
        )
        # Returns meta p-values and identifies 2 significant features
        ```
    """
    logger.info(f"calculate_stouffer_meta called with {len(p_values_dict)} modalities")

    if config.dry_run:
        # Mock response
        n_features = len(next(iter(p_values_dict.values())))
        return {
            "meta_p_values": [0.0001, 0.02, 0.35][:n_features],
            "meta_z_scores": [3.72, 2.33, 0.93][:n_features],
            "significant_features": [
                {
                    "feature_index": 0,
                    "meta_p": 0.0001,
                    "meta_z": 3.72,
                    "q_value": 0.0003,
                    "modality_contributions": {
                        mod: p_values_dict[mod][0] for mod in p_values_dict
                    },
                },
                {
                    "feature_index": 1,
                    "meta_p": 0.02,
                    "meta_z": 2.33,
                    "q_value": 0.04,
                    "modality_contributions": {
                        mod: p_values_dict[mod][1] for mod in p_values_dict
                    },
                },
            ],
            "statistics": {
                "total_features": n_features,
                "significant_features": 2,
                "fdr_threshold": config.fdr_threshold,
                "directionality_used": use_directionality,
                "weights_used": weights is not None,
            },
            "status": "success (DRY_RUN mode)",
        }

    # Real implementation
    return calculate_stouffer_meta_impl(
        p_values_dict=p_values_dict,
        effect_sizes_dict=effect_sizes_dict,
        weights=weights,
        use_directionality=use_directionality,
        fdr_threshold=config.fdr_threshold,
    )


@mcp.tool()
def create_multiomics_heatmap(
    data_path: str,
    features: Optional[List[str]] = None,
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    output_path: Optional[str] = None,
) -> Dict[str, Any]:
    """Create integrated heatmap visualization across multiple omics modalities.

    Generates publication-quality heatmap with hierarchical clustering and
    annotation tracks for treatment groups and data modalities.

    Args:
        data_path: Path to integrated multi-omics data
        features: List of feature names to include (default: top significant features)
        cluster_rows: Apply hierarchical clustering to rows (features)
        cluster_cols: Apply hierarchical clustering to columns (samples)
        output_path: Path to save plot (PNG or PDF)

    Returns:
        Dictionary with:
        - plot_path: Path to saved visualization
        - cluster_info: Cluster assignments
        - figure_data: Plotly JSON for interactive plot

    Example:
        ```
        result = create_multiomics_heatmap(
            data_path="/workspace/cache/integrated_data.pkl",
            features=["TP53", "MYC", "EGFR"],
            cluster_rows=True,
            cluster_cols=True,
            output_path="/workspace/plots/heatmap.png"
        )
        ```
    """
    logger.info(f"create_multiomics_heatmap called: data_path={data_path}")

    if config.dry_run:
        return {
            "plot_path": output_path or "/workspace/plots/multiomics_heatmap.png",
            "cluster_info": {
                "row_clusters": 4,
                "col_clusters": 2,
                "row_linkage": "complete",
                "col_linkage": "complete",
            },
            "figure_data": {
                "type": "heatmap",
                "rows": len(features) if features else 50,
                "cols": 15,
                "format": "plotly_json",
            },
            "status": "success (DRY_RUN mode)",
        }

    # TODO: Implement visualization
    raise NotImplementedError("Heatmap visualization not yet implemented")


@mcp.tool()
def run_multiomics_pca(
    data_path: str,
    modalities: Optional[List[str]] = None,
    n_components: int = 3,
    scale_features: bool = True,
    output_path: Optional[str] = None,
) -> Dict[str, Any]:
    """Run Principal Component Analysis on integrated multi-omics data.

    Performs PCA for dimensionality reduction and sample clustering visualization.
    Can analyze individual modalities or concatenated multi-omics data.

    Args:
        data_path: Path to integrated multi-omics data
        modalities: List of modalities to include (default: all available)
        n_components: Number of principal components to compute (default: 3)
        scale_features: Apply feature scaling before PCA
        output_path: Path to save 2D and 3D PCA plots

    Returns:
        Dictionary with:
        - variance_explained: Fraction of variance per component
        - loadings: Feature loadings on each PC
        - sample_coordinates: PC coordinates for each sample
        - plot_path: Path to saved visualization

    Example:
        ```
        result = run_multiomics_pca(
            data_path="/workspace/cache/integrated_data.pkl",
            modalities=["rna", "protein"],
            n_components=3,
            output_path="/workspace/plots/pca.png"
        )
        # PC1 explains 42% variance, separates Resistant vs Sensitive
        ```
    """
    logger.info(f"run_multiomics_pca called: n_components={n_components}")

    if config.dry_run:
        return {
            "variance_explained": [0.42, 0.23, 0.15][:n_components],
            "loadings": {
                "PC1": {"top_features": ["TP53", "MYC", "EGFR"], "n_features": 1800},
                "PC2": {"top_features": ["BRCA1", "KRAS", "AKT1"], "n_features": 1800},
                "PC3": {"top_features": ["CDK4", "CCND1", "RB1"], "n_features": 1800},
            },
            "sample_coordinates": {
                "samples": 15,
                "dimensions": n_components,
                "clustering_observed": True,
            },
            "plot_path": output_path or "/workspace/plots/multiomics_pca.png",
            "statistics": {
                "total_variance": 0.80,
                "n_features": 1800,
                "modalities_used": modalities or ["rna", "protein", "phospho"],
            },
            "status": "success (DRY_RUN mode)",
        }

    # TODO: Implement PCA analysis
    raise NotImplementedError("PCA analysis not yet implemented")


# ============================================================================
# RESOURCES
# ============================================================================


@mcp.resource("multiomics://config")
def get_config_resource() -> str:
    """Get current server configuration and analysis parameters.

    Returns:
        Formatted configuration settings
    """
    return f"""# Multi-Omics Server Configuration

## Data Paths
- Data Directory: {config.data_dir}
- Cache Directory: {config.cache_dir}

## Analysis Parameters
- Max Features: {config.max_features}
- Min Samples: {config.min_samples}
- FDR Threshold: {config.fdr_threshold}
- Stouffer Weights: {config.stouffer_weights}

## Mode
- DRY_RUN: {config.dry_run}
- R Home: {config.r_home or 'Auto-detect'}

## Performance
- N Jobs: {config.n_jobs}
- Timeout: {config.timeout_seconds}s
"""


@mcp.resource("multiomics://example-data")
def get_example_data_resource() -> str:
    """Get information about example datasets and data format requirements.

    Returns:
        Data format specifications and example file structures
    """
    return """# Multi-Omics Data Format Requirements

## RNA Expression Data (CSV/TSV)
```
Gene,Sample_01,Sample_02,Sample_03,...
TP53,1234.5,2341.2,987.3,...
MYC,5432.1,4321.9,3210.4,...
EGFR,876.4,1023.8,1156.2,...
```

## Protein Abundance (TMT)
```
Protein,Sample_01,Sample_02,Sample_03,...
TP53,45.2,67.8,34.1,...
MYC,123.4,145.6,98.7,...
```

## Phosphorylation Data
```
Phosphosite,Sample_01,Sample_02,Sample_03,...
TP53_S15,12.3,23.4,8.7,...
AKT1_S473,89.2,102.3,76.5,...
```

## Sample Metadata
```
Sample,Treatment,Response,Batch
Sample_01,Drug_A,Resistant,1
Sample_02,Drug_A,Sensitive,1
Sample_03,Drug_B,Resistant,2
```

## Example Dataset Locations
- `/workspace/data/multiomics/example_rna.csv` (1000 genes x 15 samples)
- `/workspace/data/multiomics/example_protein.csv` (500 proteins x 15 samples)
- `/workspace/data/multiomics/example_metadata.csv` (15 samples, 2 treatment groups)
"""


# Entry point for the server
def main():
    """Run the MCP server."""
    logger.info("Starting mcp-multiomics server...")
    logger.info(f"DRY_RUN mode: {config.dry_run}")
    mcp.run()


if __name__ == "__main__":
    main()
