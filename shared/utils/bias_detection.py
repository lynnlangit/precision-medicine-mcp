"""
Bias Detection Utilities for Precision Medicine Workflows

This module provides tools for detecting and measuring bias in precision medicine
AI/ML workflows, aligned with FDA, AMA, and NIH standards for ethical healthcare AI.

Key capabilities:
- Data representation analysis across diverse populations
- Fairness metrics calculation (demographic parity, equalized odds, calibration)
- Ancestry-aware confidence scoring
- Proxy feature detection
- Output stratification analysis

Usage:
    from shared.utils.bias_detection import (
        check_data_representation,
        calculate_fairness_metrics,
        detect_proxy_features,
        calculate_ancestry_aware_confidence
    )

References:
- FDA AI/ML-Based Software as Medical Device (SaMD) guidance
- AMA Code of Medical Ethics Opinion 2.3.2
- NIH All of Us Research Program diversity requirements
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
from collections import defaultdict
import warnings


# ============================================================================
# REFERENCE DATASET STATISTICS
# ============================================================================

# Diverse reference datasets for representation comparison
# Source: docs/hospital-deployment/ethics/ETHICS_AND_BIAS.md
REFERENCE_DATASETS = {
    "gnomad": {
        "total_genomes": 76156,
        "ancestry_distribution": {
            "european": 0.43,
            "african": 0.21,
            "latino": 0.14,
            "east_asian": 0.09,
            "south_asian": 0.09,
            "other": 0.04
        }
    },
    "all_of_us": {
        "total_participants": 245388,
        "ancestry_distribution": {
            "underrepresented": 0.80,
            "european": 0.20
        }
    },
    "gtex": {
        "total_donors": 948,
        "ancestry_distribution": {
            "european": 0.85,
            "african_american": 0.10,
            "other": 0.05,
            "asian": 0.00  # Known limitation
        }
    },
    "topmed": {
        "total_genomes": 180000,
        "ancestry_distribution": {
            "diverse": 1.0  # Intentionally diverse US population
        }
    }
}


# Risk thresholds from BIAS_AUDIT_CHECKLIST.md
REPRESENTATION_THRESHOLDS = {
    "critical": 0.05,  # <5% = CRITICAL RISK - Do not use
    "high": 0.10,      # <10% = HIGH RISK - Find alternative
    "medium": 0.20,    # <20% = MEDIUM RISK - Document limitation
    "acceptable": 0.20  # >=20% = Acceptable
}

FAIRNESS_DISPARITY_THRESHOLDS = {
    "critical": 0.20,  # >20% = Critical - Do not deploy
    "high": 0.10,      # >10% = Implement mitigation
    "acceptable": 0.10  # <=10% = Acceptable
}


# ============================================================================
# DATA CLASSES
# ============================================================================

@dataclass
class RepresentationAnalysis:
    """Results from data representation analysis."""
    dataset_name: str
    total_samples: int
    ancestry_counts: Dict[str, int]
    ancestry_proportions: Dict[str, float]
    reference_comparison: Dict[str, Dict[str, float]]
    risk_level: str  # "CRITICAL", "HIGH", "MEDIUM", "ACCEPTABLE"
    warnings: List[str]
    recommendations: List[str]


@dataclass
class FairnessMetrics:
    """Results from fairness metrics calculation."""
    metric_name: str
    overall_value: float
    group_values: Dict[str, float]
    max_disparity: float
    disparity_pairs: List[Tuple[str, str, float]]
    risk_level: str
    warnings: List[str]


@dataclass
class ProxyFeatureAnalysis:
    """Results from proxy feature detection."""
    feature_name: str
    importance_score: float
    correlation_with_protected: Dict[str, float]
    is_proxy: bool
    risk_level: str
    recommendation: str


# ============================================================================
# DATA REPRESENTATION ANALYSIS
# ============================================================================

def check_data_representation(
    data: pd.DataFrame,
    ancestry_column: str,
    reference_dataset: str = "gnomad",
    min_acceptable_pct: float = 0.10
) -> RepresentationAnalysis:
    """
    Check if data has adequate representation across diverse ancestries.

    Args:
        data: DataFrame with patient/sample data
        ancestry_column: Column name containing ancestry information
        reference_dataset: Reference dataset for comparison ("gnomad", "all_of_us", etc.)
        min_acceptable_pct: Minimum acceptable representation percentage (default 10%)

    Returns:
        RepresentationAnalysis with risk assessment and recommendations

    Example:
        >>> df = pd.DataFrame({"ancestry": ["european"]*70 + ["african"]*20 + ["asian"]*10})
        >>> analysis = check_data_representation(df, "ancestry")
        >>> print(analysis.risk_level)  # "MEDIUM" (african <20%, asian <20%)
    """
    total_samples = len(data)
    ancestry_counts = data[ancestry_column].value_counts().to_dict()
    ancestry_proportions = {k: v/total_samples for k, v in ancestry_counts.items()}

    # Get reference distribution
    ref_data = REFERENCE_DATASETS.get(reference_dataset, {})
    ref_dist = ref_data.get("ancestry_distribution", {})

    # Compare to reference
    comparison = {}
    for ancestry, ref_prop in ref_dist.items():
        actual_prop = ancestry_proportions.get(ancestry, 0.0)
        comparison[ancestry] = {
            "reference": ref_prop,
            "actual": actual_prop,
            "difference": actual_prop - ref_prop
        }

    # Assess risk level
    warnings_list = []
    recommendations = []
    min_proportion = min(ancestry_proportions.values()) if ancestry_proportions else 0.0

    if min_proportion < REPRESENTATION_THRESHOLDS["critical"]:
        risk_level = "CRITICAL"
        warnings_list.append(
            f"CRITICAL: At least one ancestry group has <5% representation ({min_proportion:.1%})"
        )
        recommendations.append("DO NOT USE this dataset for production deployment")
        recommendations.append(f"Find alternative dataset with better representation (e.g., {reference_dataset})")
    elif min_proportion < REPRESENTATION_THRESHOLDS["high"]:
        risk_level = "HIGH"
        warnings_list.append(
            f"HIGH RISK: At least one ancestry group has <10% representation ({min_proportion:.1%})"
        )
        recommendations.append("Find alternative dataset or combine with additional sources")
        recommendations.append("Document limitation prominently in all outputs")
    elif min_proportion < REPRESENTATION_THRESHOLDS["medium"]:
        risk_level = "MEDIUM"
        warnings_list.append(
            f"MEDIUM RISK: At least one ancestry group has <20% representation ({min_proportion:.1%})"
        )
        recommendations.append("Document limitation and use ancestry-aware confidence scoring")
        recommendations.append("Consider supplementing with additional diverse datasets")
    else:
        risk_level = "ACCEPTABLE"
        recommendations.append("Representation is adequate across ancestry groups")

    # Add specific warnings for underrepresented groups
    for ancestry, prop in ancestry_proportions.items():
        if prop < min_acceptable_pct:
            warnings_list.append(
                f"{ancestry.capitalize()} ancestry: {prop:.1%} (n={ancestry_counts[ancestry]})"
            )

    return RepresentationAnalysis(
        dataset_name=reference_dataset,
        total_samples=total_samples,
        ancestry_counts=ancestry_counts,
        ancestry_proportions=ancestry_proportions,
        reference_comparison=comparison,
        risk_level=risk_level,
        warnings=warnings_list,
        recommendations=recommendations
    )


# ============================================================================
# FAIRNESS METRICS
# ============================================================================

def demographic_parity(
    y_pred: np.ndarray,
    groups: np.ndarray
) -> FairnessMetrics:
    """
    Calculate demographic parity: equal positive prediction rates across groups.

    Metric: P(Ŷ=1 | A=a) should be equal for all groups a

    Args:
        y_pred: Predicted labels (binary: 0 or 1)
        groups: Protected attribute (e.g., ancestry, sex)

    Returns:
        FairnessMetrics with per-group rates and maximum disparity
    """
    unique_groups = np.unique(groups)
    group_rates = {}

    for group in unique_groups:
        mask = groups == group
        positive_rate = (y_pred[mask] == 1).sum() / mask.sum()
        group_rates[str(group)] = positive_rate

    # Calculate disparities between all pairs
    disparities = []
    rates_list = list(group_rates.values())
    groups_list = list(group_rates.keys())

    for i in range(len(rates_list)):
        for j in range(i+1, len(rates_list)):
            disparity = abs(rates_list[i] - rates_list[j])
            disparities.append((groups_list[i], groups_list[j], disparity))

    max_disparity = max([d[2] for d in disparities]) if disparities else 0.0
    overall_rate = (y_pred == 1).sum() / len(y_pred)

    # Assess risk
    warnings_list = []
    if max_disparity > FAIRNESS_DISPARITY_THRESHOLDS["critical"]:
        risk_level = "CRITICAL"
        warnings_list.append(
            f"CRITICAL: Demographic parity disparity {max_disparity:.1%} exceeds 20% threshold"
        )
        warnings_list.append("DO NOT DEPLOY - implement fairness-aware training")
    elif max_disparity > FAIRNESS_DISPARITY_THRESHOLDS["high"]:
        risk_level = "HIGH"
        warnings_list.append(
            f"HIGH RISK: Demographic parity disparity {max_disparity:.1%} exceeds 10% threshold"
        )
        warnings_list.append("Implement mitigation before deployment")
    else:
        risk_level = "ACCEPTABLE"

    return FairnessMetrics(
        metric_name="Demographic Parity",
        overall_value=overall_rate,
        group_values=group_rates,
        max_disparity=max_disparity,
        disparity_pairs=disparities,
        risk_level=risk_level,
        warnings=warnings_list
    )


def equalized_odds(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    groups: np.ndarray
) -> Dict[str, FairnessMetrics]:
    """
    Calculate equalized odds: equal TPR and FPR across groups.

    Metrics:
    - TPR (Sensitivity): P(Ŷ=1 | Y=1, A=a) should be equal for all groups a
    - FPR: P(Ŷ=1 | Y=0, A=a) should be equal for all groups a

    Args:
        y_true: True labels
        y_pred: Predicted labels
        groups: Protected attribute

    Returns:
        Dict with separate FairnessMetrics for TPR and FPR
    """
    unique_groups = np.unique(groups)
    tpr_by_group = {}
    fpr_by_group = {}

    for group in unique_groups:
        mask = groups == group
        y_t = y_true[mask]
        y_p = y_pred[mask]

        # TPR (True Positive Rate / Sensitivity)
        true_positives = ((y_p == 1) & (y_t == 1)).sum()
        actual_positives = (y_t == 1).sum()
        tpr = true_positives / actual_positives if actual_positives > 0 else 0.0

        # FPR (False Positive Rate)
        false_positives = ((y_p == 1) & (y_t == 0)).sum()
        actual_negatives = (y_t == 0).sum()
        fpr = false_positives / actual_negatives if actual_negatives > 0 else 0.0

        tpr_by_group[str(group)] = tpr
        fpr_by_group[str(group)] = fpr

    # Calculate disparities for TPR
    tpr_disparities = []
    tpr_list = list(tpr_by_group.values())
    groups_list = list(tpr_by_group.keys())

    for i in range(len(tpr_list)):
        for j in range(i+1, len(tpr_list)):
            disparity = abs(tpr_list[i] - tpr_list[j])
            tpr_disparities.append((groups_list[i], groups_list[j], disparity))

    max_tpr_disparity = max([d[2] for d in tpr_disparities]) if tpr_disparities else 0.0

    # Calculate disparities for FPR
    fpr_disparities = []
    fpr_list = list(fpr_by_group.values())

    for i in range(len(fpr_list)):
        for j in range(i+1, len(fpr_list)):
            disparity = abs(fpr_list[i] - fpr_list[j])
            fpr_disparities.append((groups_list[i], groups_list[j], disparity))

    max_fpr_disparity = max([d[2] for d in fpr_disparities]) if fpr_disparities else 0.0

    # Assess risk for TPR
    tpr_warnings = []
    if max_tpr_disparity > FAIRNESS_DISPARITY_THRESHOLDS["critical"]:
        tpr_risk_level = "CRITICAL"
        tpr_warnings.append(f"CRITICAL: TPR disparity {max_tpr_disparity:.1%} exceeds 20%")
    elif max_tpr_disparity > FAIRNESS_DISPARITY_THRESHOLDS["high"]:
        tpr_risk_level = "HIGH"
        tpr_warnings.append(f"HIGH RISK: TPR disparity {max_tpr_disparity:.1%} exceeds 10%")
    else:
        tpr_risk_level = "ACCEPTABLE"

    # Assess risk for FPR
    fpr_warnings = []
    if max_fpr_disparity > FAIRNESS_DISPARITY_THRESHOLDS["critical"]:
        fpr_risk_level = "CRITICAL"
        fpr_warnings.append(f"CRITICAL: FPR disparity {max_fpr_disparity:.1%} exceeds 20%")
    elif max_fpr_disparity > FAIRNESS_DISPARITY_THRESHOLDS["high"]:
        fpr_risk_level = "HIGH"
        fpr_warnings.append(f"HIGH RISK: FPR disparity {max_fpr_disparity:.1%} exceeds 10%")
    else:
        fpr_risk_level = "ACCEPTABLE"

    return {
        "tpr": FairnessMetrics(
            metric_name="True Positive Rate (Sensitivity)",
            overall_value=np.mean(list(tpr_by_group.values())),
            group_values=tpr_by_group,
            max_disparity=max_tpr_disparity,
            disparity_pairs=tpr_disparities,
            risk_level=tpr_risk_level,
            warnings=tpr_warnings
        ),
        "fpr": FairnessMetrics(
            metric_name="False Positive Rate",
            overall_value=np.mean(list(fpr_by_group.values())),
            group_values=fpr_by_group,
            max_disparity=max_fpr_disparity,
            disparity_pairs=fpr_disparities,
            risk_level=fpr_risk_level,
            warnings=fpr_warnings
        )
    }


def calibration_by_group(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    groups: np.ndarray,
    n_bins: int = 10
) -> FairnessMetrics:
    """
    Calculate calibration: predicted probabilities match observed frequencies.

    A well-calibrated model should have P(Y=1 | Ŷ=p, A=a) ≈ p for all groups a.

    Args:
        y_true: True labels
        y_prob: Predicted probabilities
        groups: Protected attribute
        n_bins: Number of bins for calibration curve

    Returns:
        FairnessMetrics with calibration error per group
    """
    unique_groups = np.unique(groups)
    calibration_errors = {}

    for group in unique_groups:
        mask = groups == group
        y_t = y_true[mask]
        y_p = y_prob[mask]

        # Create bins
        bins = np.linspace(0, 1, n_bins + 1)
        bin_indices = np.digitize(y_p, bins) - 1
        bin_indices = np.clip(bin_indices, 0, n_bins - 1)

        # Calculate calibration error
        errors = []
        for bin_idx in range(n_bins):
            bin_mask = bin_indices == bin_idx
            if bin_mask.sum() > 0:
                predicted_prob = y_p[bin_mask].mean()
                actual_freq = y_t[bin_mask].mean()
                errors.append(abs(predicted_prob - actual_freq))

        calibration_error = np.mean(errors) if errors else 0.0
        calibration_errors[str(group)] = calibration_error

    # Calculate disparities
    disparities = []
    errors_list = list(calibration_errors.values())
    groups_list = list(calibration_errors.keys())

    for i in range(len(errors_list)):
        for j in range(i+1, len(errors_list)):
            disparity = abs(errors_list[i] - errors_list[j])
            disparities.append((groups_list[i], groups_list[j], disparity))

    max_disparity = max([d[2] for d in disparities]) if disparities else 0.0
    overall_error = np.mean(errors_list)

    # Assess risk
    warnings_list = []
    if max_disparity > FAIRNESS_DISPARITY_THRESHOLDS["critical"]:
        risk_level = "CRITICAL"
        warnings_list.append(f"CRITICAL: Calibration disparity {max_disparity:.3f} exceeds 0.20")
    elif max_disparity > FAIRNESS_DISPARITY_THRESHOLDS["high"]:
        risk_level = "HIGH"
        warnings_list.append(f"HIGH RISK: Calibration disparity {max_disparity:.3f} exceeds 0.10")
    else:
        risk_level = "ACCEPTABLE"

    return FairnessMetrics(
        metric_name="Calibration Error",
        overall_value=overall_error,
        group_values=calibration_errors,
        max_disparity=max_disparity,
        disparity_pairs=disparities,
        risk_level=risk_level,
        warnings=warnings_list
    )


def calculate_fairness_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    y_prob: Optional[np.ndarray],
    groups: np.ndarray,
    metrics: List[str] = ["demographic_parity", "equalized_odds", "calibration"]
) -> Dict[str, Any]:
    """
    Calculate multiple fairness metrics at once.

    Args:
        y_true: True labels
        y_pred: Predicted labels (binary)
        y_prob: Predicted probabilities (optional, required for calibration)
        groups: Protected attribute
        metrics: List of metrics to calculate

    Returns:
        Dict mapping metric names to FairnessMetrics objects
    """
    results = {}

    if "demographic_parity" in metrics:
        results["demographic_parity"] = demographic_parity(y_pred, groups)

    if "equalized_odds" in metrics:
        results["equalized_odds"] = equalized_odds(y_true, y_pred, groups)

    if "calibration" in metrics:
        if y_prob is None:
            warnings.warn("Calibration requires y_prob, skipping")
        else:
            results["calibration"] = calibration_by_group(y_true, y_prob, groups)

    return results


# ============================================================================
# PROXY FEATURE DETECTION
# ============================================================================

def detect_proxy_features(
    X: pd.DataFrame,
    feature_importances: Dict[str, float],
    protected_attributes: pd.DataFrame,
    correlation_threshold: float = 0.5,
    importance_threshold: float = 0.05
) -> List[ProxyFeatureAnalysis]:
    """
    Detect features that may serve as proxies for protected attributes.

    A proxy feature is one that:
    1. Has high correlation with a protected attribute (e.g., zip code → race)
    2. Has non-trivial feature importance in the model (>5%)

    Args:
        X: Feature matrix
        feature_importances: Dict mapping feature names to importance scores
        protected_attributes: DataFrame with protected attributes (ancestry, sex, etc.)
        correlation_threshold: Minimum correlation to flag as proxy (default 0.5)
        importance_threshold: Minimum importance to flag as proxy (default 0.05)

    Returns:
        List of ProxyFeatureAnalysis for flagged features
    """
    proxy_analyses = []

    for feature_name, importance in feature_importances.items():
        if feature_name not in X.columns:
            continue

        feature_values = X[feature_name]
        correlations = {}

        # Calculate correlations with protected attributes
        for protected_col in protected_attributes.columns:
            protected_values = protected_attributes[protected_col]

            # Handle categorical vs numerical
            if pd.api.types.is_numeric_dtype(feature_values) and \
               pd.api.types.is_numeric_dtype(protected_values):
                corr = np.corrcoef(feature_values, protected_values)[0, 1]
            else:
                # Use Cramér's V for categorical associations
                corr = cramers_v(feature_values, protected_values)

            correlations[protected_col] = abs(corr)

        # Check if feature is a proxy
        max_correlation = max(correlations.values())
        is_proxy = (importance > importance_threshold and
                   max_correlation > correlation_threshold)

        if is_proxy:
            # Assess risk level
            if importance > 0.20:
                risk_level = "CRITICAL"
                recommendation = "REMOVE feature immediately and retrain model"
            elif importance > 0.10:
                risk_level = "HIGH"
                recommendation = "Remove feature and retrain, or document justification"
            else:
                risk_level = "MEDIUM"
                recommendation = "Consider removing feature, document if retained"

            proxy_analyses.append(ProxyFeatureAnalysis(
                feature_name=feature_name,
                importance_score=importance,
                correlation_with_protected=correlations,
                is_proxy=is_proxy,
                risk_level=risk_level,
                recommendation=recommendation
            ))

    return proxy_analyses


def cramers_v(x: pd.Series, y: pd.Series) -> float:
    """
    Calculate Cramér's V statistic for categorical association.

    Cramér's V ranges from 0 (no association) to 1 (perfect association).
    """
    contingency_table = pd.crosstab(x, y)
    chi2 = ((contingency_table -
             contingency_table.sum(axis=1).values[:, None] *
             contingency_table.sum(axis=0).values /
             contingency_table.sum().sum()) ** 2 /
            (contingency_table.sum(axis=1).values[:, None] *
             contingency_table.sum(axis=0).values /
             contingency_table.sum().sum())).sum().sum()

    n = contingency_table.sum().sum()
    min_dim = min(contingency_table.shape[0], contingency_table.shape[1]) - 1

    if min_dim == 0:
        return 0.0

    return np.sqrt(chi2 / (n * min_dim))


# ============================================================================
# ANCESTRY-AWARE CONFIDENCE SCORING
# ============================================================================

def calculate_ancestry_aware_confidence(
    base_confidence: float,
    patient_ancestry: str,
    variant_id: str,
    reference_study_counts: Dict[str, int],
    min_studies_thresholds: Dict[str, Tuple[int, float]] = None
) -> Dict[str, Any]:
    """
    Adjust confidence scores based on ancestral representation in reference data.

    This implements the methodology from ETHICS_AND_BIAS.md section 6.3.

    Args:
        base_confidence: Model's base confidence score (0-1)
        patient_ancestry: Patient's ancestry (e.g., "european", "african", "asian")
        variant_id: Variant identifier (for logging)
        reference_study_counts: Number of studies per ancestry for this variant
        min_studies_thresholds: Dict mapping threshold names to (min_studies, penalty)
                                Default: {"high": (20, 0.0), "medium": (5, 0.1), "low": (0, 0.3)}

    Returns:
        Dict with adjusted confidence, warnings, and metadata

    Example:
        >>> confidence_result = calculate_ancestry_aware_confidence(
        ...     base_confidence=0.85,
        ...     patient_ancestry="african",
        ...     variant_id="BRCA1:c.5266dupC",
        ...     reference_study_counts={"african": 2, "european": 50}
        ... )
        >>> print(confidence_result["adjusted_confidence"])  # 0.595 (30% penalty)
        >>> print(confidence_result["warnings"])  # ["Limited data in african ancestry (<5 studies)"]
    """
    if min_studies_thresholds is None:
        min_studies_thresholds = {
            "high": (20, 0.0),    # >=20 studies: no penalty
            "medium": (5, 0.1),   # 5-19 studies: 10% penalty
            "low": (0, 0.3)       # <5 studies: 30% penalty
        }

    ancestral_studies = reference_study_counts.get(patient_ancestry, 0)

    # Determine confidence penalty
    if ancestral_studies >= min_studies_thresholds["high"][0]:
        confidence_penalty = min_studies_thresholds["high"][1]
        confidence_level = "HIGH"
    elif ancestral_studies >= min_studies_thresholds["medium"][0]:
        confidence_penalty = min_studies_thresholds["medium"][1]
        confidence_level = "MEDIUM"
    else:
        confidence_penalty = min_studies_thresholds["low"][1]
        confidence_level = "LOW"

    adjusted_confidence = base_confidence * (1 - confidence_penalty)

    # Generate warnings
    warnings_list = []
    if confidence_level == "LOW":
        warnings_list.append(
            f"Limited data in {patient_ancestry} ancestry (<5 studies) for {variant_id}"
        )
        warnings_list.append(
            f"Confidence reduced by {confidence_penalty*100:.0f}% due to limited ancestral representation"
        )
    elif confidence_level == "MEDIUM":
        warnings_list.append(
            f"Moderate data in {patient_ancestry} ancestry (5-19 studies) for {variant_id}"
        )

    return {
        "variant_id": variant_id,
        "patient_ancestry": patient_ancestry,
        "base_confidence": base_confidence,
        "adjusted_confidence": adjusted_confidence,
        "confidence_level": confidence_level,
        "confidence_penalty": confidence_penalty,
        "ancestral_studies": ancestral_studies,
        "warnings": warnings_list,
        "all_study_counts": reference_study_counts
    }


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def generate_bias_report(
    representation_analysis: RepresentationAnalysis,
    fairness_metrics: Dict[str, Any],
    proxy_features: List[ProxyFeatureAnalysis],
    output_format: str = "text"
) -> str:
    """
    Generate a comprehensive bias audit report.

    Args:
        representation_analysis: Results from check_data_representation()
        fairness_metrics: Results from calculate_fairness_metrics()
        proxy_features: Results from detect_proxy_features()
        output_format: "text", "html", or "json"

    Returns:
        Formatted report string
    """
    if output_format == "text":
        report = []
        report.append("=" * 80)
        report.append("BIAS AUDIT REPORT")
        report.append("=" * 80)
        report.append("")

        # Section 1: Data Representation
        report.append("1. DATA REPRESENTATION ANALYSIS")
        report.append("-" * 80)
        report.append(f"Risk Level: {representation_analysis.risk_level}")
        report.append(f"Total Samples: {representation_analysis.total_samples}")
        report.append("")
        report.append("Ancestry Distribution:")
        for ancestry, count in representation_analysis.ancestry_counts.items():
            prop = representation_analysis.ancestry_proportions[ancestry]
            report.append(f"  - {ancestry}: {count} ({prop:.1%})")
        report.append("")

        if representation_analysis.warnings:
            report.append("Warnings:")
            for warning in representation_analysis.warnings:
                report.append(f"  ⚠ {warning}")
            report.append("")

        if representation_analysis.recommendations:
            report.append("Recommendations:")
            for rec in representation_analysis.recommendations:
                report.append(f"  → {rec}")
            report.append("")

        # Section 2: Fairness Metrics
        report.append("2. FAIRNESS METRICS")
        report.append("-" * 80)
        for metric_name, metric_result in fairness_metrics.items():
            if isinstance(metric_result, dict) and "tpr" in metric_result:
                # Handle equalized_odds (has both TPR and FPR)
                for sub_metric_name, sub_metric in metric_result.items():
                    report.append(f"{sub_metric.metric_name}:")
                    report.append(f"  Risk Level: {sub_metric.risk_level}")
                    report.append(f"  Max Disparity: {sub_metric.max_disparity:.1%}")
                    report.append("  Per-Group Values:")
                    for group, value in sub_metric.group_values.items():
                        report.append(f"    - {group}: {value:.3f}")
                    if sub_metric.warnings:
                        for warning in sub_metric.warnings:
                            report.append(f"  ⚠ {warning}")
                    report.append("")
            else:
                report.append(f"{metric_result.metric_name}:")
                report.append(f"  Risk Level: {metric_result.risk_level}")
                report.append(f"  Max Disparity: {metric_result.max_disparity:.1%}")
                report.append("  Per-Group Values:")
                for group, value in metric_result.group_values.items():
                    report.append(f"    - {group}: {value:.3f}")
                if metric_result.warnings:
                    for warning in metric_result.warnings:
                        report.append(f"  ⚠ {warning}")
                report.append("")

        # Section 3: Proxy Features
        report.append("3. PROXY FEATURE ANALYSIS")
        report.append("-" * 80)
        if proxy_features:
            for proxy in proxy_features:
                report.append(f"Feature: {proxy.feature_name}")
                report.append(f"  Risk Level: {proxy.risk_level}")
                report.append(f"  Importance: {proxy.importance_score:.3f}")
                report.append("  Correlations with Protected Attributes:")
                for attr, corr in proxy.correlation_with_protected.items():
                    report.append(f"    - {attr}: {corr:.3f}")
                report.append(f"  Recommendation: {proxy.recommendation}")
                report.append("")
        else:
            report.append("No proxy features detected.")
            report.append("")

        report.append("=" * 80)
        report.append("END OF REPORT")
        report.append("=" * 80)

        return "\n".join(report)

    elif output_format == "json":
        import json
        return json.dumps({
            "representation_analysis": representation_analysis.__dict__,
            "fairness_metrics": {k: v.__dict__ if hasattr(v, '__dict__') else v
                               for k, v in fairness_metrics.items()},
            "proxy_features": [p.__dict__ for p in proxy_features]
        }, indent=2)

    else:
        raise ValueError(f"Unsupported output format: {output_format}")
