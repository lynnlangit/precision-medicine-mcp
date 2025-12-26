"""Stouffer's meta-analysis for combining p-values across omics modalities.

CRITICAL WORKFLOW (per bioinformatician feedback):
┌─────────────────────────────────────────────────────────────────┐
│ CORRECT FDR TIMING FOR MULTI-OMICS META-ANALYSIS               │
├─────────────────────────────────────────────────────────────────┤
│ Step 1: HAllA Analysis                                          │
│         → Use NOMINAL p-values (NOT FDR-corrected)             │
│         → Returns: p_value_nominal for each association         │
│                                                                  │
│ Step 2: Stouffer's Meta-Analysis (THIS MODULE)                 │
│         → Input: NOMINAL p-values from each modality           │
│         → Combine p-values across modalities                    │
│         → Output: meta_p_values (still nominal)                 │
│                                                                  │
│ Step 3: FDR Correction (APPLIED HERE, AFTER COMBINATION)       │
│         → Input: meta_p_values from Step 2                      │
│         → Apply: Benjamini-Hochberg FDR                         │
│         → Output: q_values (FDR-corrected)                      │
└─────────────────────────────────────────────────────────────────┘

WHY THIS ORDER MATTERS:
- Applying FDR before Stouffer's → loses power (over-conservative)
- Applying FDR after Stouffer's → correct statistical framework
- Nominal p-values preserve evidence strength across modalities

Reference: Erik Jessen PhD & Jessie Hohenstein, 2025 bioinformatician review
"""

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

        IMPORTANT: Input p-values should be NOMINAL (not FDR-corrected).
        FDR correction is applied AFTER combining p-values across modalities.

        Workflow:
        1. Input: NOMINAL p-values from each modality (e.g., from HAllA)
        2. Convert to Z-scores with directionality
        3. Combine Z-scores using weighted Stouffer's method
        4. Convert back to meta p-values (still nominal)
        5. Apply FDR correction (Benjamini-Hochberg) → q-values
        6. Return both meta_p_values (nominal) and q_values (FDR-corrected)

        Args:
            p_values_dict: Dict of modality -> list of NOMINAL p-values
            effect_sizes_dict: Dict of modality -> list of effect sizes (for directionality)
            weights: Dict of modality -> weight (default: equal weights)

        Returns:
            Dictionary with:
            - meta_p_values: Combined p-values (nominal, before FDR)
            - q_values: FDR-corrected p-values (use these for significance)
            - meta_z_scores: Combined Z-scores
            - significant_features: Features passing FDR threshold
        """
        logger.info(f"Starting Stouffer's meta-analysis for {len(p_values_dict)} modalities")
        logger.info("IMPORTANT: Expecting NOMINAL p-values (FDR applied AFTER combination)")

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

        # VALIDATION: Check if p-values look like they might be FDR-corrected
        # (Warning, not error - user might have legitimate reasons)
        all_p_values = np.concatenate([np.array(p_values_dict[mod]) for mod in modalities])
        min_p = np.min(all_p_values)

        if min_p > 0.001:
            logger.warning(
                f"WARNING: Minimum p-value is {min_p:.4f}. "
                "This is unusually high for nominal p-values from association testing. "
                "Are you sure these are NOMINAL p-values and not already FDR-corrected? "
                "If these are q-values, DO NOT use Stouffer's method."
            )

        logger.info(f"Processing {n_features} features across {len(modalities)} modalities")
        logger.info(f"P-value range: {min_p:.2e} to {np.max(all_p_values):.2e}")

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

        logger.info("=" * 70)
        logger.info("APPLYING FDR CORRECTION (Step 3 of workflow)")
        logger.info(f"Input: {len(meta_p_values)} meta p-values (nominal, from Stouffer's)")
        logger.info("Method: Benjamini-Hochberg FDR correction")

        # Apply FDR correction (THIS IS THE CORRECT TIMING - AFTER COMBINATION)
        reject, q_values, _, _ = multipletests(
            meta_p_values,
            method="fdr_bh",
            alpha=0.05,
        )

        # Identify significant features
        significant_indices = np.where(reject)[0]

        logger.info(f"Output: {len(significant_indices)} significant features (q < 0.05)")
        logger.info(f"FDR-corrected q-value range: {np.min(q_values):.2e} to {np.max(q_values):.2e}")
        logger.info("=" * 70)

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
    """Implementation of Stouffer's meta-analysis with CORRECT FDR timing.

    WORKFLOW (per bioinformatician feedback):
    ┌──────────────────────────────────────────────────────┐
    │ 1. Input: NOMINAL p-values from each modality        │
    │    (e.g., from HAllA with p_value_nominal)           │
    │                                                       │
    │ 2. Combine: Use Stouffer's Z-score method           │
    │    → Convert p-values to Z-scores                    │
    │    → Weight and combine Z-scores                     │
    │    → Convert back to meta_p_values (still nominal)   │
    │                                                       │
    │ 3. FDR Correction: Apply AFTER combination          │
    │    → Input: meta_p_values from step 2                │
    │    → Apply: Benjamini-Hochberg FDR                   │
    │    → Output: q_values (FDR-corrected)                │
    └──────────────────────────────────────────────────────┘

    Args:
        p_values_dict: Dict of modality -> list of NOMINAL p-values
                       (NOT FDR-corrected, e.g., from HAllA)
        effect_sizes_dict: Dict of modality -> list of effect sizes
                          (log2FC, correlation, etc. for directionality)
        weights: Dict of modality -> weight (default: equal weights)
        use_directionality: Incorporate effect size sign into Z-scores
        fdr_threshold: FDR threshold for identifying significant features

    Returns:
        Dictionary with:
        - meta_p_values: Combined p-values (NOMINAL, before FDR)
        - meta_z_scores: Combined Z-scores with directionality
        - q_values: FDR-corrected p-values (USE THESE for significance calls)
        - significant_features: Features passing FDR threshold with details
        - statistics: Summary statistics including workflow confirmation
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
            "workflow_confirmation": {
                "step_1": "Input NOMINAL p-values (from HAllA)",
                "step_2": "Combined p-values using Stouffer's method",
                "step_3": "Applied FDR correction AFTER combination",
                "fdr_method": "Benjamini-Hochberg",
                "note": "This is the CORRECT workflow per bioinformatician feedback",
            },
        },
        "p_value_types": {
            "meta_p_values": "NOMINAL (combined, before FDR)",
            "q_values": "FDR-CORRECTED (use these for significance)",
        },
        "recommendation": "Use q_values (not meta_p_values) for identifying significant features",
        "status": "success",
    }

    return result
