"""Basic usage example for quantum cell type fidelity analysis.

This example demonstrates:
1. Creating synthetic spatial transcriptomics data
2. Training quantum embeddings
3. Computing fidelity between cell types
4. Identifying immune evasion states
"""

import numpy as np
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from quantum_celltype_fidelity import (
    QuCoWECircuit,
    QuCoWECellTypeEmbedding,
    QuCoWETrainer,
    TrainingConfig,
    ImmuneEvasionDetector,
    compute_fidelity_summary_stats
)


def create_synthetic_data(n_cells_per_type=50, n_cell_types=5, feature_dim=256):
    """Create synthetic spatial transcriptomics data.

    Args:
        n_cells_per_type: Number of cells per cell type
        n_cell_types: Number of cell types
        feature_dim: Feature dimension

    Returns:
        Dictionary mapping cell_type -> list of feature vectors
    """
    print(f"Creating synthetic data: {n_cell_types} types, {n_cells_per_type} cells each")

    cell_types = [f"CellType_{i}" for i in range(n_cell_types)]
    training_data = {}

    for i, cell_type in enumerate(cell_types):
        # Create features with cell-type-specific mean
        mean = np.random.randn(feature_dim) * 2
        features_list = []

        for _ in range(n_cells_per_type):
            # Add noise around cell type mean
            features = mean + np.random.randn(feature_dim) * 0.5
            features_list.append(features)

        training_data[cell_type] = features_list

    return training_data, cell_types


def main():
    """Run basic usage example."""
    print("=" * 60)
    print("Quantum Cell Type Fidelity Analysis - Basic Example")
    print("=" * 60)

    # Step 1: Create synthetic data
    print("\n[1/5] Creating synthetic data...")
    training_data, cell_types = create_synthetic_data(
        n_cells_per_type=20,
        n_cell_types=4,
        feature_dim=256
    )
    print(f"  Created data for {len(cell_types)} cell types: {cell_types}")

    # Step 2: Initialize quantum embeddings
    print("\n[2/5] Initializing quantum embeddings...")
    embedding = QuCoWECellTypeEmbedding(
        cell_types=cell_types,
        n_qubits=6,  # 2^6 = 64 dimensional Hilbert space
        n_layers=2,
        feature_dim=256,
        backend="cpu"
    )
    print(f"  Circuit: {embedding.circuit.n_qubits} qubits, {embedding.circuit.n_layers} layers")
    print(f"  Hilbert space dimension: {embedding.circuit.hilbert_dim}")
    print(f"  Parameters per cell type: {len(embedding.circuit.theta_params)}")

    # Step 3: Train embeddings
    print("\n[3/5] Training quantum embeddings...")
    config = TrainingConfig(
        n_epochs=10,
        learning_rate=0.05,
        batch_size=8,
        optimizer="adam"
    )
    trainer = QuCoWETrainer(embedding, config)
    training_summary = trainer.train(training_data, verbose=False)

    print(f"  Training complete!")
    print(f"  Final loss: {training_summary['final_loss']:.4f}")
    print(f"  Training time: {training_summary['training_time']:.1f}s")

    # Step 4: Compute fidelity matrix
    print("\n[4/5] Computing fidelity matrix...")
    fidelity_matrix = embedding.compute_fidelity_matrix(training_data)
    print(f"  Matrix shape: {fidelity_matrix.shape}")

    # Compute summary statistics
    cell_labels = []
    for cell_type, features_list in training_data.items():
        cell_labels.extend([cell_type] * len(features_list))

    summary_stats = compute_fidelity_summary_stats(fidelity_matrix, cell_labels)
    print(f"  Mean fidelity: {summary_stats['mean_fidelity']:.3f}")
    print(f"  Std fidelity: {summary_stats['std_fidelity']:.3f}")

    print("\n  Per-cell-type statistics:")
    for cell_type, stats in summary_stats['per_cell_type'].items():
        print(f"    {cell_type}:")
        print(f"      Within-type fidelity: {stats['within_type_mean']:.3f}")
        print(f"      Between-type fidelity: {stats['between_type_mean']:.3f}")
        print(f"      Separation score: {stats['separation_score']:.3f}")

    # Step 5: Identify immune evasion (demonstration)
    print("\n[5/5] Demonstrating immune evasion detection...")

    # For this demo, treat first two cell types as "immune" types
    immune_types = cell_types[:2]
    print(f"  Immune cell types: {immune_types}")

    detector = ImmuneEvasionDetector(
        immune_cell_types=immune_types,
        evasion_threshold=0.3
    )

    # Detect evasion for all cells
    evasion_results = detector.detect_batch(
        fidelity_matrix=fidelity_matrix,
        cell_type_labels=cell_labels
    )

    n_evading = sum(1 for _, _, metadata in evasion_results if metadata['is_evading'])
    print(f"  Found {n_evading} cells in evasion states")

    # Show top 5 evading cells
    sorted_results = sorted(evasion_results, key=lambda x: x[1], reverse=True)
    print("\n  Top 5 evading cells:")
    for i, (cell_idx, evasion_score, metadata) in enumerate(sorted_results[:5]):
        cell_type = cell_labels[cell_idx]
        print(f"    {i+1}. Cell {cell_idx} ({cell_type}): score={evasion_score:.3f}")

    # Summary
    print("\n" + "=" * 60)
    print("Example complete!")
    print(f"  - Trained embeddings for {len(cell_types)} cell types")
    print(f"  - Computed {fidelity_matrix.shape[0]}x{fidelity_matrix.shape[1]} fidelity matrix")
    print(f"  - Identified {n_evading} cells with immune evasion signatures")
    print("=" * 60)


if __name__ == "__main__":
    main()
