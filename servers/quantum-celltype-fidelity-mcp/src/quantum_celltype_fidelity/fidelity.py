"""Fidelity-based scoring heads for cell type analysis.

Implements various scoring functions based on quantum fidelity:
- Direct fidelity scoring
- Logit-transformed fidelity
- Confidence-weighted scoring
- Immune evasion detection
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass


@dataclass
class FidelityScore:
    """Container for fidelity-based scores.

    Attributes:
        raw_fidelity: Raw quantum fidelity F(|ψ_a⟩, |ψ_b⟩) in [0, 1]
        logit_score: Logit-transformed score for classification
        confidence: Confidence score based on fidelity distribution
        cell_type_a: First cell type
        cell_type_b: Second cell type (if comparing two cells)
    """
    raw_fidelity: float
    logit_score: float
    confidence: float
    cell_type_a: str
    cell_type_b: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None


class FidelityHead:
    """Basic fidelity scoring head.

    Computes raw fidelity scores and provides thresholding
    for binary classification tasks.
    """

    def __init__(self, threshold: float = 0.5):
        """Initialize fidelity head.

        Args:
            threshold: Fidelity threshold for binary decisions
        """
        self.threshold = threshold

    def score(self, fidelity: float) -> FidelityScore:
        """Score a raw fidelity value.

        Args:
            fidelity: Raw fidelity in [0, 1]

        Returns:
            FidelityScore with basic scoring
        """
        # Confidence based on distance from threshold
        confidence = abs(fidelity - self.threshold)

        # Simple logit transformation
        # Avoid log(0) by clipping
        fidelity_clipped = np.clip(fidelity, 1e-10, 1 - 1e-10)
        logit = np.log(fidelity_clipped / (1 - fidelity_clipped))

        return FidelityScore(
            raw_fidelity=fidelity,
            logit_score=float(logit),
            confidence=float(confidence),
            cell_type_a="unknown"
        )

    def predict(self, fidelity: float) -> bool:
        """Binary prediction based on threshold.

        Args:
            fidelity: Raw fidelity score

        Returns:
            True if fidelity >= threshold, False otherwise
        """
        return fidelity >= self.threshold


class LogitFidelityHead:
    """Logit-transformed fidelity scoring head.

    Applies learned temperature scaling and bias to fidelity scores
    for improved calibration in classification tasks.
    """

    def __init__(
        self,
        temperature: float = 1.0,
        bias: float = 0.0,
        threshold: float = 0.0
    ):
        """Initialize logit fidelity head.

        Args:
            temperature: Temperature parameter for scaling (T > 0)
            bias: Bias term for shifting decision boundary
            threshold: Logit threshold for binary decisions
        """
        if temperature <= 0:
            raise ValueError("Temperature must be positive")

        self.temperature = temperature
        self.bias = bias
        self.threshold = threshold

    def compute_logit(self, fidelity: float) -> float:
        """Compute temperature-scaled logit.

        logit(F) = (1/T) * log(F / (1 - F)) + b

        Args:
            fidelity: Raw fidelity in [0, 1]

        Returns:
            Scaled logit score
        """
        # Clip to avoid log(0)
        fidelity_clipped = np.clip(fidelity, 1e-10, 1 - 1e-10)

        # Compute logit
        raw_logit = np.log(fidelity_clipped / (1 - fidelity_clipped))

        # Apply temperature and bias
        scaled_logit = (raw_logit / self.temperature) + self.bias

        return float(scaled_logit)

    def score(
        self,
        fidelity: float,
        cell_type_a: str,
        cell_type_b: Optional[str] = None
    ) -> FidelityScore:
        """Score fidelity with logit transformation.

        Args:
            fidelity: Raw fidelity
            cell_type_a: First cell type
            cell_type_b: Second cell type (optional)

        Returns:
            FidelityScore with logit scoring
        """
        logit = self.compute_logit(fidelity)

        # Confidence based on logit magnitude
        confidence = 1.0 / (1.0 + np.exp(-abs(logit)))

        return FidelityScore(
            raw_fidelity=fidelity,
            logit_score=logit,
            confidence=float(confidence),
            cell_type_a=cell_type_a,
            cell_type_b=cell_type_b
        )

    def predict(self, fidelity: float) -> bool:
        """Binary prediction using logit threshold.

        Args:
            fidelity: Raw fidelity

        Returns:
            True if logit(fidelity) >= threshold
        """
        logit = self.compute_logit(fidelity)
        return logit >= self.threshold

    def calibrate(
        self,
        fidelities: np.ndarray,
        labels: np.ndarray,
        method: str = "temperature_scaling"
    ) -> None:
        """Calibrate head on labeled data.

        Args:
            fidelities: Array of raw fidelity scores
            labels: Binary labels (0/1)
            method: Calibration method ("temperature_scaling")
        """
        if method == "temperature_scaling":
            # Simple grid search for temperature
            from scipy.optimize import minimize_scalar

            def nll(temp):
                """Negative log-likelihood for temperature."""
                logits = [self.compute_logit(f) for f in fidelities]
                logits = np.array(logits) / temp
                probs = 1.0 / (1.0 + np.exp(-logits))
                probs = np.clip(probs, 1e-10, 1 - 1e-10)
                loss = -np.mean(
                    labels * np.log(probs) + (1 - labels) * np.log(1 - probs)
                )
                return loss

            result = minimize_scalar(nll, bounds=(0.1, 10.0), method='bounded')
            self.temperature = result.x
            print(f"Calibrated temperature: {self.temperature:.3f}")

        else:
            raise ValueError(f"Unknown calibration method: {method}")


class ImmuneEvasionDetector:
    """Detector for immune evasion states using fidelity patterns.

    Identifies cells with:
    - Low fidelity to canonical immune cell types
    - High fidelity to exhausted/dysfunctional states
    - Unusual fidelity patterns indicating evasion
    """

    def __init__(
        self,
        immune_cell_types: List[str],
        exhausted_markers: Optional[List[str]] = None,
        evasion_threshold: float = 0.3
    ):
        """Initialize immune evasion detector.

        Args:
            immune_cell_types: List of canonical immune cell types
            exhausted_markers: Cell types/states indicating exhaustion
            evasion_threshold: Fidelity threshold below which evasion is flagged
        """
        self.immune_cell_types = immune_cell_types
        self.exhausted_markers = exhausted_markers or []
        self.evasion_threshold = evasion_threshold

    def compute_evasion_score(
        self,
        fidelity_to_immune: Dict[str, float],
        fidelity_to_exhausted: Optional[Dict[str, float]] = None
    ) -> Tuple[float, Dict[str, Any]]:
        """Compute immune evasion score for a cell.

        Score is high when:
        - Low fidelity to canonical immune types
        - High fidelity to exhausted states (if provided)

        Args:
            fidelity_to_immune: Dict of {immune_cell_type: fidelity}
            fidelity_to_exhausted: Dict of {exhausted_type: fidelity}

        Returns:
            Tuple of (evasion_score, metadata_dict)
        """
        # Average fidelity to immune cells
        if not fidelity_to_immune:
            raise ValueError("Must provide fidelity to at least one immune type")

        avg_immune_fidelity = np.mean(list(fidelity_to_immune.values()))

        # Evasion score: inverse of immune fidelity
        evasion_score = 1.0 - avg_immune_fidelity

        # Boost score if high fidelity to exhausted states
        exhaustion_boost = 0.0
        if fidelity_to_exhausted:
            avg_exhausted_fidelity = np.mean(list(fidelity_to_exhausted.values()))
            exhaustion_boost = avg_exhausted_fidelity * 0.5  # Weight

        final_score = min(evasion_score + exhaustion_boost, 1.0)

        # Metadata
        metadata = {
            "avg_immune_fidelity": float(avg_immune_fidelity),
            "evasion_score": float(evasion_score),
            "exhaustion_boost": float(exhaustion_boost),
            "final_evasion_score": float(final_score),
            "is_evading": final_score >= self.evasion_threshold,
            "immune_fidelities": fidelity_to_immune,
            "exhausted_fidelities": fidelity_to_exhausted or {}
        }

        return float(final_score), metadata

    def detect_batch(
        self,
        fidelity_matrix: np.ndarray,
        cell_type_labels: List[str],
        cell_indices: Optional[List[int]] = None
    ) -> List[Tuple[int, float, Dict[str, Any]]]:
        """Detect immune evasion in a batch of cells.

        Args:
            fidelity_matrix: Pairwise fidelity matrix (n_cells, n_cells)
            cell_type_labels: Cell type label for each cell
            cell_indices: Indices to analyze (if None, analyze all)

        Returns:
            List of (cell_idx, evasion_score, metadata) for each cell
        """
        if cell_indices is None:
            cell_indices = list(range(len(cell_type_labels)))

        results = []

        # Find indices of immune cell types in the matrix
        immune_indices = [
            i for i, label in enumerate(cell_type_labels)
            if label in self.immune_cell_types
        ]

        exhausted_indices = [
            i for i, label in enumerate(cell_type_labels)
            if label in self.exhausted_markers
        ]

        for cell_idx in cell_indices:
            # Get fidelities to immune cells
            fidelity_to_immune = {
                cell_type_labels[i]: fidelity_matrix[cell_idx, i]
                for i in immune_indices
            }

            # Get fidelities to exhausted cells
            fidelity_to_exhausted = None
            if exhausted_indices:
                fidelity_to_exhausted = {
                    cell_type_labels[i]: fidelity_matrix[cell_idx, i]
                    for i in exhausted_indices
                }

            # Compute evasion score
            evasion_score, metadata = self.compute_evasion_score(
                fidelity_to_immune,
                fidelity_to_exhausted
            )

            results.append((cell_idx, evasion_score, metadata))

        return results


def compute_fidelity_summary_stats(
    fidelity_matrix: np.ndarray,
    cell_type_labels: List[str]
) -> Dict[str, Any]:
    """Compute summary statistics for a fidelity matrix.

    Args:
        fidelity_matrix: Pairwise fidelity matrix (n_cells, n_cells)
        cell_type_labels: Cell type label for each cell

    Returns:
        Dictionary with summary statistics
    """
    n_cells = fidelity_matrix.shape[0]

    # Overall statistics
    # Exclude diagonal (self-fidelities)
    mask = ~np.eye(n_cells, dtype=bool)
    off_diagonal_fidelities = fidelity_matrix[mask]

    stats = {
        "n_cells": n_cells,
        "mean_fidelity": float(np.mean(off_diagonal_fidelities)),
        "std_fidelity": float(np.std(off_diagonal_fidelities)),
        "min_fidelity": float(np.min(off_diagonal_fidelities)),
        "max_fidelity": float(np.max(off_diagonal_fidelities)),
        "median_fidelity": float(np.median(off_diagonal_fidelities))
    }

    # Per-cell-type statistics
    unique_types = sorted(set(cell_type_labels))
    per_type_stats = {}

    for cell_type in unique_types:
        # Get indices for this type
        indices = [i for i, label in enumerate(cell_type_labels) if label == cell_type]

        if len(indices) < 2:
            continue

        # Within-type fidelities
        within_fidelities = []
        for i in indices:
            for j in indices:
                if i != j:
                    within_fidelities.append(fidelity_matrix[i, j])

        # Between-type fidelities
        between_fidelities = []
        for i in indices:
            for j in range(n_cells):
                if cell_type_labels[j] != cell_type:
                    between_fidelities.append(fidelity_matrix[i, j])

        per_type_stats[cell_type] = {
            "n_cells": len(indices),
            "within_type_mean": float(np.mean(within_fidelities)),
            "within_type_std": float(np.std(within_fidelities)),
            "between_type_mean": float(np.mean(between_fidelities)),
            "between_type_std": float(np.std(between_fidelities)),
            "separation_score": float(
                np.mean(within_fidelities) - np.mean(between_fidelities)
            )
        }

    stats["per_cell_type"] = per_type_stats

    return stats
