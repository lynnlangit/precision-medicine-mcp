"""Quantum Cell Type Fidelity MCP Server.

Quantum computing-based analysis of cell type fidelity in spatial transcriptomics data.
Uses Qiskit to build parameterized quantum circuits for embedding cell types into
Hilbert space and computing quantum fidelity metrics.
"""

from .circuits import QuCoWECircuit
from .embeddings import QuCoWECellTypeEmbedding
from .fidelity import (
    FidelityScore,
    FidelityHead,
    LogitFidelityHead,
    ImmuneEvasionDetector,
    compute_fidelity_summary_stats
)
from .spatial_context import SpatialContextGenerator, SpatialNeighborhood
from .training import QuCoWETrainer, TrainingConfig, ParameterShiftGradientEstimator

__version__ = "0.1.0"

__all__ = [
    "QuCoWECircuit",
    "QuCoWECellTypeEmbedding",
    "FidelityScore",
    "FidelityHead",
    "LogitFidelityHead",
    "ImmuneEvasionDetector",
    "compute_fidelity_summary_stats",
    "SpatialContextGenerator",
    "SpatialNeighborhood",
    "QuCoWETrainer",
    "TrainingConfig",
    "ParameterShiftGradientEstimator",
]
