"""Generate realistic synthetic multi-omics test fixtures."""

import numpy as np
import pandas as pd
from pathlib import Path

# Set random seed for reproducibility
np.random.seed(42)

# Sample information
N_SAMPLES = 15
N_RESISTANT = 7
N_SENSITIVE = 8

# Feature counts
N_GENES = 1000
N_PROTEINS = 500
N_PHOSPHO = 300

# Sample names
sample_names = [f"Sample_{i:02d}" for i in range(1, N_SAMPLES + 1)]

# Treatment groups (Resistant = 0-6, Sensitive = 7-14)
treatment_groups = (
    ["Resistant"] * N_RESISTANT +
    ["Sensitive"] * N_SENSITIVE
)

# Gene names (using common cancer genes + generic names)
cancer_genes = [
    "TP53", "MYC", "KRAS", "EGFR", "BRCA1", "BRCA2", "PTEN", "AKT1",
    "PIK3CA", "RB1", "CDK4", "CCND1", "ERBB2", "MET", "ALK",
    "BRAF", "NRAS", "HRAS", "FGFR1", "FGFR2", "JAK2", "STAT3"
]
generic_genes = [f"Gene_{i:04d}" for i in range(1, N_GENES - len(cancer_genes) + 1)]
gene_names = cancer_genes + generic_genes

# Protein names (subset of genes)
protein_names = gene_names[:N_PROTEINS]

# Phosphosite names
phospho_names = [f"{gene}_S{np.random.randint(1, 500)}" for gene in gene_names[:N_PHOSPHO]]


def generate_expression_data(
    n_features: int,
    feature_names: list,
    n_resistant: int,
    n_sensitive: int,
    n_differential: int = 50,
    baseline_mean: float = 5.0,
    baseline_std: float = 1.5,
    fold_change_range: tuple = (1.5, 4.0),
) -> pd.DataFrame:
    """Generate log-normal expression data with differential features.

    Args:
        n_features: Number of features (genes/proteins)
        feature_names: List of feature names
        n_resistant: Number of resistant samples
        n_sensitive: Number of sensitive samples
        n_differential: Number of differentially expressed features
        baseline_mean: Mean of log-normal distribution
        baseline_std: Std of log-normal distribution
        fold_change_range: Range of fold changes for differential features

    Returns:
        DataFrame with features as rows, samples as columns
    """
    n_samples = n_resistant + n_sensitive

    # Initialize data matrix
    data = np.zeros((n_features, n_samples))

    # Generate baseline expression (log-normal)
    for i in range(n_features):
        baseline = np.random.lognormal(baseline_mean, baseline_std, n_samples)
        data[i, :] = baseline

    # Add differential expression for subset of features
    differential_indices = np.random.choice(n_features, n_differential, replace=False)

    for idx in differential_indices:
        # Random fold change
        fold_change = np.random.uniform(*fold_change_range)

        # Randomly decide if upregulated or downregulated in resistant
        direction = np.random.choice([-1, 1])

        if direction == 1:
            # Upregulated in resistant
            data[idx, :n_resistant] *= fold_change
        else:
            # Downregulated in resistant (upregulated in sensitive)
            data[idx, n_resistant:] *= fold_change

    # Add small amount of noise
    noise = np.random.normal(0, 0.1, data.shape)
    data = data * (1 + noise)

    # Ensure positive values
    data = np.maximum(data, 0.1)

    # Create DataFrame
    df = pd.DataFrame(
        data,
        index=feature_names[:n_features],
        columns=sample_names,
    )

    return df


def generate_metadata(sample_names: list, treatment_groups: list) -> pd.DataFrame:
    """Generate sample metadata.

    Args:
        sample_names: List of sample names
        treatment_groups: List of treatment response labels

    Returns:
        Metadata DataFrame
    """
    metadata = pd.DataFrame({
        "Sample": sample_names,
        "Response": treatment_groups,
        "Batch": [1, 1, 2, 2, 3, 3, 4, 1, 1, 2, 2, 3, 3, 4, 4][:len(sample_names)],
        "Age": np.random.randint(40, 80, len(sample_names)),
        "Sex": np.random.choice(["M", "F"], len(sample_names)),
    })

    return metadata


def main():
    """Generate all test fixtures."""
    print("Generating synthetic multi-omics test fixtures...")

    # Get fixture directory
    fixture_dir = Path(__file__).parent

    # Generate RNA expression data (FPKM values)
    print(f"Generating RNA data: {N_GENES} genes x {N_SAMPLES} samples")
    rna_data = generate_expression_data(
        n_features=N_GENES,
        feature_names=gene_names,
        n_resistant=N_RESISTANT,
        n_sensitive=N_SENSITIVE,
        n_differential=100,  # 10% differential
        baseline_mean=6.0,
        baseline_std=1.8,
    )
    rna_path = fixture_dir / "sample_rna.csv"
    rna_data.to_csv(rna_path)
    print(f"  Saved to {rna_path}")

    # Generate protein abundance data (TMT values)
    print(f"Generating Protein data: {N_PROTEINS} proteins x {N_SAMPLES} samples")
    protein_data = generate_expression_data(
        n_features=N_PROTEINS,
        feature_names=protein_names,
        n_resistant=N_RESISTANT,
        n_sensitive=N_SENSITIVE,
        n_differential=50,  # 10% differential
        baseline_mean=4.5,
        baseline_std=1.2,
    )
    protein_path = fixture_dir / "sample_protein.csv"
    protein_data.to_csv(protein_path)
    print(f"  Saved to {protein_path}")

    # Generate phosphorylation data
    print(f"Generating Phospho data: {N_PHOSPHO} sites x {N_SAMPLES} samples")
    phospho_data = generate_expression_data(
        n_features=N_PHOSPHO,
        feature_names=phospho_names,
        n_resistant=N_RESISTANT,
        n_sensitive=N_SENSITIVE,
        n_differential=30,  # 10% differential
        baseline_mean=3.5,
        baseline_std=1.5,
        fold_change_range=(2.0, 5.0),  # More dramatic changes
    )
    phospho_path = fixture_dir / "sample_phospho.csv"
    phospho_data.to_csv(phospho_path)
    print(f"  Saved to {phospho_path}")

    # Generate metadata
    print(f"Generating metadata: {N_SAMPLES} samples")
    metadata = generate_metadata(sample_names, treatment_groups)
    metadata_path = fixture_dir / "sample_metadata.csv"
    metadata.to_csv(metadata_path, index=False)
    print(f"  Saved to {metadata_path}")

    # Print summary statistics
    print("\nData Summary:")
    print(f"  RNA range: {rna_data.min().min():.2f} - {rna_data.max().max():.2f}")
    print(f"  Protein range: {protein_data.min().min():.2f} - {protein_data.max().max():.2f}")
    print(f"  Phospho range: {phospho_data.min().min():.2f} - {phospho_data.max().max():.2f}")
    print(f"\nTreatment groups:")
    print(f"  Resistant: {N_RESISTANT} samples")
    print(f"  Sensitive: {N_SENSITIVE} samples")

    print("\nFixtures generated successfully!")


if __name__ == "__main__":
    main()
