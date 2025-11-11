"""Tests for multi-omics data integration."""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path

from mcp_multiomics.tools.integration import integrate_omics_data_impl
from mcp_multiomics.tools.utils import (
    load_omics_data,
    align_samples,
    filter_missing_features,
    normalize_zscore,
)


class TestDataLoading:
    """Tests for data loading utilities."""

    def test_load_rna_data(self, rna_path):
        """Test loading RNA expression data."""
        df = load_omics_data(rna_path)

        assert isinstance(df, pd.DataFrame)
        assert df.shape[0] == 1000  # 1000 genes
        assert df.shape[1] == 15  # 15 samples
        assert all(col.startswith("Sample_") for col in df.columns)

    def test_load_protein_data(self, protein_path):
        """Test loading protein abundance data."""
        df = load_omics_data(protein_path)

        assert df.shape[0] == 500  # 500 proteins
        assert df.shape[1] == 15  # 15 samples

    def test_load_missing_file(self):
        """Test error handling for missing files."""
        with pytest.raises(FileNotFoundError):
            load_omics_data("/nonexistent/file.csv")


class TestSampleAlignment:
    """Tests for sample alignment across modalities."""

    def test_align_common_samples(self, rna_path, protein_path):
        """Test alignment of common samples."""
        rna_df = load_omics_data(rna_path)
        protein_df = load_omics_data(protein_path)

        dataframes = {"rna": rna_df, "protein": protein_df}
        aligned, common_samples = align_samples(dataframes)

        assert len(common_samples) == 15
        assert all(sample in aligned["rna"].columns for sample in common_samples)
        assert all(sample in aligned["protein"].columns for sample in common_samples)

    def test_align_no_common_samples(self):
        """Test error when no common samples exist."""
        df1 = pd.DataFrame({"S1": [1, 2], "S2": [3, 4]})
        df2 = pd.DataFrame({"S3": [5, 6], "S4": [7, 8]})

        with pytest.raises(ValueError, match="No common samples"):
            align_samples({"df1": df1, "df2": df2})


class TestFeatureFiltering:
    """Tests for missing data filtering."""

    def test_filter_missing_features(self):
        """Test filtering features with missing data."""
        # Create test data with some missing values
        data = pd.DataFrame({
            "S1": [1.0, np.nan, 3.0, np.nan],
            "S2": [2.0, np.nan, 4.0, 5.0],
            "S3": [3.0, 6.0, np.nan, 6.0],
        })

        # Filter with 50% threshold
        filtered = filter_missing_features(data, threshold=0.5)

        # Feature 1 (row index 1): 2/3 missing = 66% > 50%, should be removed
        assert filtered.shape[0] == 3
        assert 1 not in filtered.index

    def test_no_filtering_needed(self, rna_path):
        """Test when no features need filtering."""
        df = load_omics_data(rna_path)
        filtered = filter_missing_features(df, threshold=0.5)

        assert filtered.shape == df.shape


class TestNormalization:
    """Tests for Z-score normalization."""

    def test_zscore_normalization(self):
        """Test Z-score normalization."""
        data = pd.DataFrame({
            "S1": [1.0, 10.0, 100.0],
            "S2": [2.0, 20.0, 200.0],
            "S3": [3.0, 30.0, 300.0],
        })

        normalized = normalize_zscore(data)

        # Check each row is normalized
        for idx in normalized.index:
            row = normalized.loc[idx]
            assert abs(row.mean()) < 1e-10  # Mean should be ~0
            assert abs(row.std() - 1.0) < 1e-10  # Std should be ~1

    def test_zscore_constant_feature(self):
        """Test normalization with constant feature."""
        data = pd.DataFrame({
            "S1": [5.0, 10.0],
            "S2": [5.0, 20.0],
            "S3": [5.0, 30.0],
        })

        normalized = normalize_zscore(data)

        # Constant row should remain constant
        assert all(normalized.iloc[0] == 5.0)


class TestIntegration:
    """Tests for complete integration pipeline."""

    def test_integrate_rna_only(self, rna_path, tmp_path):
        """Test integration with RNA data only."""
        result = integrate_omics_data_impl(
            rna_path=rna_path,
            normalize=True,
            filter_missing=0.5,
        )

        assert result["status"] == "success"
        assert "rna" in result["feature_counts"]
        assert result["feature_counts"]["rna"] > 0
        assert len(result["common_samples"]) == 15
        assert result["qc_metrics"]["normalization"] == "z-score"

    def test_integrate_all_modalities(
        self, rna_path, protein_path, phospho_path, metadata_path
    ):
        """Test integration with all three modalities."""
        result = integrate_omics_data_impl(
            rna_path=rna_path,
            protein_path=protein_path,
            phospho_path=phospho_path,
            metadata_path=metadata_path,
            normalize=True,
            filter_missing=0.5,
        )

        assert result["status"] == "success"

        # Check all modalities present
        assert "rna" in result["feature_counts"]
        assert "protein" in result["feature_counts"]
        assert "phospho" in result["feature_counts"]

        # Check metadata
        assert result["metadata"] is not None
        assert result["metadata"]["samples"] == 15
        assert "treatment_resistant" in result["metadata"]
        assert "treatment_sensitive" in result["metadata"]

        # Check QC metrics
        assert "features_filtered" in result["qc_metrics"]
        assert "missing_data" in result["qc_metrics"]

    def test_integrate_no_normalization(self, rna_path):
        """Test integration without normalization."""
        result = integrate_omics_data_impl(
            rna_path=rna_path,
            normalize=False,
            filter_missing=1.0,  # No filtering
        )

        assert result["qc_metrics"]["normalization"] == "none"
        assert result["qc_metrics"]["missing_threshold"] == 1.0

    def test_integrate_strict_filtering(self, rna_path):
        """Test integration with strict missing data filtering."""
        result = integrate_omics_data_impl(
            rna_path=rna_path,
            normalize=True,
            filter_missing=0.0,  # No missing data allowed
        )

        # With strict filtering, might remove some features
        assert result["feature_counts"]["rna"] <= 1000

    def test_cache_file_created(self, rna_path):
        """Test that integrated data is cached."""
        result = integrate_omics_data_impl(
            rna_path=rna_path,
            normalize=True,
        )

        cache_path = Path(result["cache_path"])
        assert cache_path.exists()
        assert cache_path.suffix == ".pkl"
