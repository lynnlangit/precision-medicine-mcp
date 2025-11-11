"""Utility functions for multi-omics analysis."""

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats


def load_omics_data(file_path: str) -> pd.DataFrame:
    """Load omics data from CSV or TSV file.

    Args:
        file_path: Path to data file

    Returns:
        DataFrame with features as rows, samples as columns
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Data file not found: {file_path}")

    # Detect delimiter
    delimiter = "\t" if file_path.endswith((".tsv", ".txt")) else ","

    # Load data
    df = pd.read_csv(file_path, sep=delimiter, index_col=0)

    return df


def align_samples(
    dataframes: Dict[str, pd.DataFrame]
) -> Tuple[Dict[str, pd.DataFrame], List[str]]:
    """Align samples across multiple omics modalities.

    Args:
        dataframes: Dict of modality -> DataFrame

    Returns:
        Tuple of (aligned_dataframes, common_samples)
    """
    # Find common samples
    sample_sets = [set(df.columns) for df in dataframes.values()]
    common_samples = sorted(list(set.intersection(*sample_sets)))

    if len(common_samples) == 0:
        raise ValueError("No common samples found across modalities")

    # Align dataframes to common samples
    aligned = {
        modality: df[common_samples] for modality, df in dataframes.items()
    }

    return aligned, common_samples


def filter_missing_features(
    df: pd.DataFrame, threshold: float = 0.5
) -> pd.DataFrame:
    """Filter out features with too many missing values.

    Args:
        df: DataFrame with features as rows
        threshold: Maximum fraction of missing values allowed (0.0-1.0)

    Returns:
        Filtered DataFrame
    """
    missing_fraction = df.isna().sum(axis=1) / df.shape[1]
    keep_features = missing_fraction <= threshold
    return df[keep_features]


def normalize_zscore(df: pd.DataFrame) -> pd.DataFrame:
    """Apply Z-score normalization across samples for each feature.

    Args:
        df: DataFrame with features as rows, samples as columns

    Returns:
        Normalized DataFrame
    """
    # Calculate Z-scores for each row (feature)
    # Use ddof=1 for sample standard deviation
    normalized = df.apply(
        lambda row: (row - row.mean()) / row.std(ddof=1) if row.std() > 0 else row,
        axis=1,
    )

    return normalized


def calculate_qc_metrics(
    dataframes: Dict[str, pd.DataFrame],
    original_counts: Dict[str, int],
) -> Dict[str, any]:
    """Calculate quality control metrics for integrated data.

    Args:
        dataframes: Dict of modality -> processed DataFrame
        original_counts: Dict of modality -> original feature count

    Returns:
        Dictionary with QC metrics
    """
    metrics = {
        "samples": list(dataframes.values())[0].shape[1],
        "features_retained": {
            modality: df.shape[0] for modality, df in dataframes.items()
        },
        "features_filtered": {
            modality: original_counts[modality] - df.shape[0]
            for modality, df in dataframes.items()
        },
        "missing_data": {
            modality: float(df.isna().sum().sum() / (df.shape[0] * df.shape[1]))
            for modality, df in dataframes.items()
        },
    }

    return metrics


def save_integrated_data(
    dataframes: Dict[str, pd.DataFrame],
    metadata: Optional[pd.DataFrame],
    output_path: str,
) -> None:
    """Save integrated multi-omics data to pickle file.

    Args:
        dataframes: Dict of modality -> DataFrame
        metadata: Sample metadata DataFrame
        output_path: Output file path
    """
    data_package = {
        "omics_data": dataframes,
        "metadata": metadata,
    }

    import pickle

    with open(output_path, "wb") as f:
        pickle.dump(data_package, f)


def load_integrated_data(file_path: str) -> Tuple[Dict[str, pd.DataFrame], Optional[pd.DataFrame]]:
    """Load integrated multi-omics data from pickle file.

    Args:
        file_path: Path to integrated data file

    Returns:
        Tuple of (omics_dataframes, metadata)
    """
    import pickle

    with open(file_path, "rb") as f:
        data_package = pickle.load(f)

    return data_package["omics_data"], data_package.get("metadata")
