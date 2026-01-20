"""Visualization utilities for perturbation predictions."""

import scanpy as sc
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Optional, List
import logging

logger = logging.getLogger(__name__)


class PerturbationVisualizer:
    """Create visualizations for perturbation analysis."""

    def __init__(self, output_dir: str = "./data/plots"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def plot_baseline_vs_predicted(
        self,
        baseline_adata: ad.AnnData,
        predicted_adata: ad.AnnData,
        plot_type: str = "pca",
        color_by: str = "condition",
        output_path: Optional[str] = None,
        title: Optional[str] = None
    ) -> str:
        """Plot baseline vs predicted cells in reduced dimensions.

        Args:
            baseline_adata: Baseline cell states
            predicted_adata: Predicted cell states
            plot_type: "pca" or "umap"
            color_by: Column to color by
            output_path: Optional path to save figure
            title: Optional plot title

        Returns:
            Path to saved figure
        """
        # Combine datasets
        baseline_adata.obs[color_by] = "baseline"
        predicted_adata.obs[color_by] = "predicted"

        combined = ad.concat([baseline_adata, predicted_adata], label="source")

        # Compute embedding
        if plot_type == "pca":
            sc.tl.pca(combined, n_comps=50)
            embedding_key = "X_pca"
        elif plot_type == "umap":
            # Need neighbors for UMAP
            sc.pp.neighbors(combined, n_neighbors=15, n_pcs=50)
            sc.tl.umap(combined)
            embedding_key = "X_umap"
        else:
            raise ValueError(f"Unknown plot_type: {plot_type}")

        # Create figure
        fig, ax = plt.subplots(figsize=(10, 8))

        # Get embedding coordinates
        coords = combined.obsm[embedding_key]

        # Plot each condition
        for condition in ["baseline", "predicted"]:
            mask = combined.obs[color_by] == condition
            ax.scatter(
                coords[mask, 0],
                coords[mask, 1],
                label=condition,
                alpha=0.6,
                s=20
            )

        ax.set_xlabel(f"{plot_type.upper()} 1")
        ax.set_ylabel(f"{plot_type.upper()} 2")
        ax.legend()

        if title:
            ax.set_title(title)
        else:
            ax.set_title(f"Baseline vs. Predicted ({plot_type.upper()})")

        plt.tight_layout()

        # Save figure
        if output_path is None:
            output_path = self.output_dir / f"baseline_vs_predicted_{plot_type}.png"
        else:
            output_path = Path(output_path)

        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"Saved {plot_type.upper()} plot to {output_path}")
        return str(output_path)

    def plot_latent_space(
        self,
        latent_representation: np.ndarray,
        labels: np.ndarray,
        output_path: Optional[str] = None,
        title: str = "Latent Space Representation"
    ) -> str:
        """Plot latent space using PCA.

        Args:
            latent_representation: Latent vectors (n_cells x n_latent)
            labels: Cell labels for coloring
            output_path: Path to save figure
            title: Plot title

        Returns:
            Path to saved figure
        """
        # Apply PCA to latent space
        from sklearn.decomposition import PCA

        pca = PCA(n_components=2)
        coords_2d = pca.fit_transform(latent_representation)

        # Plot
        fig, ax = plt.subplots(figsize=(10, 8))

        unique_labels = np.unique(labels)
        for label in unique_labels:
            mask = labels == label
            ax.scatter(
                coords_2d[mask, 0],
                coords_2d[mask, 1],
                label=label,
                alpha=0.6,
                s=20
            )

        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
        ax.set_title(title)
        ax.legend()

        plt.tight_layout()

        # Save
        if output_path is None:
            output_path = self.output_dir / "latent_space_pca.png"
        else:
            output_path = Path(output_path)

        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"Saved latent space plot to {output_path}")
        return str(output_path)

    def plot_differential_expression(
        self,
        de_results: dict,
        n_genes: int = 20,
        output_path: Optional[str] = None
    ) -> str:
        """Plot differential expression results as volcano plot.

        Args:
            de_results: DE results from DifferentialExpressionAnalyzer
            n_genes: Number of top genes to label
            output_path: Path to save figure

        Returns:
            Path to saved figure
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

        # Upregulated genes bar plot
        if de_results["upregulated_genes"]:
            up_genes = de_results["upregulated_genes"][:n_genes]
            genes = [g["gene"] for g in up_genes]
            log_fc = [g["log2_fc"] for g in up_genes]

            ax1.barh(genes, log_fc, color="red", alpha=0.7)
            ax1.set_xlabel("Log2 Fold Change")
            ax1.set_title(f"Top {n_genes} Upregulated Genes")
            ax1.invert_yaxis()

        # Downregulated genes bar plot
        if de_results["downregulated_genes"]:
            down_genes = de_results["downregulated_genes"][:n_genes]
            genes = [g["gene"] for g in down_genes]
            log_fc = [g["log2_fc"] for g in down_genes]

            ax2.barh(genes, log_fc, color="blue", alpha=0.7)
            ax2.set_xlabel("Log2 Fold Change")
            ax2.set_title(f"Top {n_genes} Downregulated Genes")
            ax2.invert_yaxis()

        plt.tight_layout()

        # Save
        if output_path is None:
            output_path = self.output_dir / "differential_expression.png"
        else:
            output_path = Path(output_path)

        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"Saved DE plot to {output_path}")
        return str(output_path)

    def plot_gene_expression_comparison(
        self,
        baseline_adata: ad.AnnData,
        predicted_adata: ad.AnnData,
        genes: List[str],
        output_path: Optional[str] = None
    ) -> str:
        """Compare gene expression between baseline and predicted.

        Args:
            baseline_adata: Baseline cells
            predicted_adata: Predicted cells
            genes: List of genes to compare
            output_path: Path to save figure

        Returns:
            Path to saved figure
        """
        n_genes = len(genes)
        fig, axes = plt.subplots(1, n_genes, figsize=(5 * n_genes, 4))

        if n_genes == 1:
            axes = [axes]

        for i, gene in enumerate(genes):
            if gene not in baseline_adata.var_names:
                logger.warning(f"Gene {gene} not found in data")
                continue

            gene_idx = list(baseline_adata.var_names).index(gene)

            # Get expression values
            baseline_expr = baseline_adata.X[:, gene_idx]
            predicted_expr = predicted_adata.X[:, gene_idx]

            if hasattr(baseline_expr, 'A'):
                baseline_expr = baseline_expr.A.flatten()
                predicted_expr = predicted_expr.A.flatten()

            # Create violin plot
            data = {
                "Expression": np.concatenate([baseline_expr, predicted_expr]),
                "Condition": ["Baseline"] * len(baseline_expr) + ["Predicted"] * len(predicted_expr)
            }

            import pandas as pd
            df = pd.DataFrame(data)

            sns.violinplot(data=df, x="Condition", y="Expression", ax=axes[i])
            axes[i].set_title(gene)
            axes[i].set_ylabel("Expression" if i == 0 else "")

        plt.tight_layout()

        # Save
        if output_path is None:
            output_path = self.output_dir / "gene_expression_comparison.png"
        else:
            output_path = Path(output_path)

        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"Saved gene expression comparison to {output_path}")
        return str(output_path)
