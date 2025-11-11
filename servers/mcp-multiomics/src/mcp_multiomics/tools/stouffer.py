"""Stouffer's meta-analysis for combining p-values across omics modalities."""

import logging
from typing import Any, Dict, List, Optional

import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)


class StoufferMetaAnalysis:
    """Stouffer's Z-score method for meta-analysis of p-values.

    Combines p-values from multiple omics modalities using inverse normal method:
    1. Convert p-values to Z-scores: Z = Φ^(-1)(1 - p/2)
    2. Combine Z-scores: Z_meta = Σ(w_i * Z_i) / sqrt(Σ(w_i^2))
    3. Convert back to p-value: p_meta = 2 * (1 - Φ(Z_meta))

    Supports directionality from effect sizes (log2FC, correlations).
    """

    def __init__(self, use_directionality: bool = True):
        """Initialize Stouffer's meta-analysis.

        Args:
            use_directionality: Incorporate effect size sign into Z-scores
        """
        self.use_directionality = use_directionality

    def p_to_z(
        self,
        p_values: np.ndarray,
        effect_sizes: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """Convert p-values to Z-scores (standard normal deviates).

        Args:
            p_values: Array of p-values (0 < p < 1)
            effect_sizes: Optional effect sizes for directionality

        Returns:
            Array of Z-scores
        """
        # Clip p-values to avoid numerical issues
        p_values = np.clip(p_values, 1e-300, 1 - 1e-15)

        # Convert to two-tailed Z-scores
        z_scores = stats.norm.ppf(1 - p_values / 2)

        # Apply directionality from effect sizes
        if self.use_directionality and effect_sizes is not None:
            signs = np.sign(effect_sizes)
            # Negative effect size means Z-score should be negative
            z_scores = z_scores * signs

        return z_scores

    def z_to_p(self, z_scores: np.ndarray) -> np.ndarray:
        """Convert Z-scores back to two-tailed p-values.

        Args:
            z_scores: Array of Z-scores

        Returns:
            Array of two-tailed p-values
        """
        # Two-tailed p-value: p = 2 * (1 - Φ(|Z|))
        p_values = 2 * (1 - stats.norm.cdf(np.abs(z_scores)))
        return p_values

    def combine_z_scores(
        self,
        z_scores: np.ndarray,
        weights: Optional[np.ndarray] = None,
    ) -> float:
        """Combine Z-scores using Stouffer's method.

        Formula: Z_meta = Σ(w_i * Z_i) / sqrt(Σ(w_i^2))

        Args:
            z_scores: Array of Z-scores from different modalities
            weights: Optional weights (e.g., sqrt of sample sizes)

        Returns:
            Combined Z-score
        """
        if weights is None:
            weights = np.ones(len(z_scores))

        # Stouffer's formula with weights
        numerator = np.sum(weights * z_scores)
        denominator = np.sqrt(np.sum(weights ** 2))

        z_meta = numerator / denominator

        return z_meta

    def meta_analyze(
        self,
        p_values_dict: Dict[str, List[float]],
        effect_sizes_dict: Optional[Dict[str, List[float]]] = None,
        weights: Optional[Dict[str, float]] = None,
    ) -> Dict[str, Any]:
        """Perform Stouffer's meta-analysis across omics modalities.

        Args:
            p_values_dict: Dict of modality -> list of p-values
            effect_sizes_dict: Dict of modality -> list of effect sizes
            weights: Dict of modality -> weight (default: equal weights)

        Returns:
            Dictionary with meta-analysis results
        """
        logger.info(f"Starting Stouffer's meta-analysis for {len(p_values_dict)} modalities")

        # Validate input
        modalities = list(p_values_dict.keys())
        n_features = len(p_values_dict[modalities[0]])

        # Check all modalities have same number of features
        for modality in modalities:
            if len(p_values_dict[modality]) != n_features:
                raise ValueError(
                    f"Modality {modality} has {len(p_values_dict[modality])} features, "
                    f"expected {n_features}"
                )

        # Prepare weights array
        if weights is None:
            weights_array = np.ones(len(modalities))
        else:
            weights_array = np.array([weights.get(mod, 1.0) for mod in modalities])

        # Convert to numpy arrays
        p_values_matrix = np.array([p_values_dict[mod] for mod in modalities])

        effect_sizes_matrix = None
        if effect_sizes_dict is not None and self.use_directionality:
            effect_sizes_matrix = np.array([effect_sizes_dict[mod] for mod in modalities])

        # Initialize results
        meta_p_values = []
        meta_z_scores = []

        # Process each feature
        for feat_idx in range(n_features):
            # Get p-values and effect sizes for this feature
            p_vals = p_values_matrix[:, feat_idx]

            eff_sizes = None
            if effect_sizes_matrix is not None:
                eff_sizes = effect_sizes_matrix[:, feat_idx]

            # Convert p-values to Z-scores
            z_scores = self.p_to_z(p_vals, eff_sizes)

            # Combine Z-scores
            z_meta = self.combine_z_scores(z_scores, weights_array)
            meta_z_scores.append(z_meta)

            # Convert back to p-value
            p_meta = self.z_to_p(np.array([z_meta]))[0]
            meta_p_values.append(p_meta)

        # Convert to numpy arrays
        meta_p_values = np.array(meta_p_values)
        meta_z_scores = np.array(meta_z_scores)

        # Apply FDR correction
        reject, q_values, _, _ = multipletests(
            meta_p_values,
            method="fdr_bh",
            alpha=0.05,
        )

        # Identify significant features
        significant_indices = np.where(reject)[0]

        logger.info(f"Found {len(significant_indices)} significant features")

        return {
            "meta_p_values": meta_p_values.tolist(),
            "meta_z_scores": meta_z_scores.tolist(),
            "q_values": q_values.tolist(),
            "significant_features": significant_indices.tolist(),
            "n_significant": len(significant_indices),
        }


def calculate_stouffer_meta_impl(
    p_values_dict: Dict[str, List[float]],
    effect_sizes_dict: Optional[Dict[str, List[float]]] = None,
    weights: Optional[Dict[str, float]] = None,
    use_directionality: bool = True,
    fdr_threshold: float = 0.05,
) -> Dict[str, Any]:
    """Implementation of Stouffer's meta-analysis.

    Args:
        p_values_dict: Dict of modality -> list of p-values
        effect_sizes_dict: Dict of modality -> list of effect sizes
        weights: Dict of modality -> weight
        use_directionality: Incorporate effect size sign
        fdr_threshold: FDR threshold for significance

    Returns:
        Dictionary with meta-analysis results
    """
    # Create analyzer
    analyzer = StoufferMetaAnalysis(use_directionality=use_directionality)

    # Run meta-analysis
    results = analyzer.meta_analyze(
        p_values_dict=p_values_dict,
        effect_sizes_dict=effect_sizes_dict,
        weights=weights,
    )

    # Get number of features
    n_features = len(results["meta_p_values"])

    # Build detailed results for significant features
    significant_features = []
    for idx in results["significant_features"]:
        feature_info = {
            "feature_index": int(idx),
            "meta_p": float(results["meta_p_values"][idx]),
            "meta_z": float(results["meta_z_scores"][idx]),
            "q_value": float(results["q_values"][idx]),
            "modality_contributions": {
                modality: float(p_values_dict[modality][idx])
                for modality in p_values_dict.keys()
            },
        }

        if effect_sizes_dict is not None:
            feature_info["effect_sizes"] = {
                modality: float(effect_sizes_dict[modality][idx])
                for modality in effect_sizes_dict.keys()
            }

        significant_features.append(feature_info)

    # Prepare final result
    result = {
        "meta_p_values": results["meta_p_values"],
        "meta_z_scores": results["meta_z_scores"],
        "q_values": results["q_values"],
        "significant_features": significant_features,
        "statistics": {
            "total_features": n_features,
            "significant_features": results["n_significant"],
            "fdr_threshold": fdr_threshold,
            "directionality_used": use_directionality,
            "weights_used": weights is not None,
            "modalities": list(p_values_dict.keys()),
        },
        "status": "success",
    }

    return result
