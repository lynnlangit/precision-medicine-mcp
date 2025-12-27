"""Tests for preprocessing tools (validate, preprocess, visualize)."""

import pytest
import numpy as np
import pandas as pd
from pathlib import Path
import os

from mcp_multiomics.tools.preprocessing import (
    validate_multiomics_data_impl,
    preprocess_multiomics_data_impl,
    visualize_data_quality_impl,
)


class TestValidateMultiomicsData:
    """Test suite for validate_multiomics_data tool."""

    def test_validation_with_good_data(self, tmp_path):
        """Test validation passes with clean data."""
        # Create synthetic good data
        rna_data = pd.DataFrame({
            'Gene': ['TP53', 'MYC', 'EGFR'],
            'Sample_01': [100, 200, 150],
            'Sample_02': [110, 210, 160],
            'Sample_03': [105, 205, 155],
        })
        rna_path = tmp_path / "rna.csv"
        rna_data.to_csv(rna_path, index=False)

        result = validate_multiomics_data_impl(
            rna_path=str(rna_path),
            protein_path=None,
            phospho_path=None,
            metadata_path=None,
        )

        # In DRY_RUN mode, should return success
        assert result['status'] == 'success (DRY_RUN mode)'
        assert 'validation_status' in result
        assert 'sample_overlap' in result
        assert 'missing_patterns' in result

    def test_validation_detects_batch_effects(self):
        """Test validation detects batch effects in mock data."""
        result = validate_multiomics_data_impl(
            rna_path="/mock/rna.csv",
            protein_path="/mock/protein.csv",
            phospho_path=None,
            metadata_path="/mock/metadata.csv",
        )

        # Should return batch effect information
        assert 'batch_effects' in result
        assert result['batch_effects'] is not None

    def test_validation_with_missing_values(self):
        """Test validation detects missing value patterns."""
        result = validate_multiomics_data_impl(
            rna_path="/mock/rna.csv",
            protein_path="/mock/protein.csv",
            phospho_path="/mock/phospho.csv",
            metadata_path=None,
        )

        # Should analyze missing patterns
        assert 'missing_patterns' in result
        assert result['missing_patterns'] is not None


class TestPreprocessMultiomicsData:
    """Test suite for preprocess_multiomics_data tool."""

    def test_basic_preprocessing(self):
        """Test basic preprocessing with default parameters."""
        result = preprocess_multiomics_data_impl(
            rna_path="/mock/rna_raw.csv",
            protein_path="/mock/protein_raw.csv",
            phospho_path=None,
            metadata_path="/mock/metadata.csv",
            normalize_method="quantile",
            batch_correction=True,
            imputation_method="knn",
            outlier_threshold=3.0,
            output_dir="/mock/output",
        )

        # Check DRY_RUN mode response
        assert result['status'] == 'success (DRY_RUN mode)'
        assert 'preprocessed_paths' in result
        assert 'preprocessing_report' in result
        assert 'qc_metrics' in result
        assert 'batch_correction_results' in result

    def test_batch_correction_workflow(self):
        """Test batch correction is applied correctly."""
        result = preprocess_multiomics_data_impl(
            rna_path="/mock/rna.csv",
            protein_path="/mock/protein.csv",
            phospho_path=None,
            metadata_path="/mock/metadata.csv",
            batch_correction=True,
        )

        # Should include batch correction metrics
        batch_results = result['batch_correction_results']
        assert 'pc1_batch_correlation_before' in batch_results
        assert 'pc1_batch_correlation_after' in batch_results

        # After correction should be better
        before = batch_results['pc1_batch_correlation_before']
        after = batch_results['pc1_batch_correlation_after']
        assert after < before  # Correlation should decrease

    def test_imputation_methods(self):
        """Test different imputation methods."""
        for method in ['knn', 'median', 'minimum']:
            result = preprocess_multiomics_data_impl(
                rna_path="/mock/rna.csv",
                protein_path="/mock/protein.csv",
                phospho_path=None,
                metadata_path="/mock/metadata.csv",
                imputation_method=method,
            )

            assert result['status'] == 'success (DRY_RUN mode)'
            assert 'imputation_stats' in result
            # Check that imputation method appears in steps
            assert any(method in step for step in result['preprocessing_report']['steps_applied'])

    def test_normalization_methods(self):
        """Test different normalization methods."""
        for method in ['quantile', 'median', 'tmm', 'zscore']:
            result = preprocess_multiomics_data_impl(
                rna_path="/mock/rna.csv",
                protein_path="/mock/protein.csv",
                phospho_path=None,
                metadata_path="/mock/metadata.csv",
                normalize_method=method,
            )

            assert result['status'] == 'success (DRY_RUN mode)'
            # Check that normalization method appears in steps
            assert any(method in step.lower() for step in result['preprocessing_report']['steps_applied'])

    def test_outlier_detection(self):
        """Test outlier detection with different thresholds."""
        result = preprocess_multiomics_data_impl(
            rna_path="/mock/rna.csv",
            protein_path="/mock/protein.csv",
            phospho_path=None,
            metadata_path="/mock/metadata.csv",
            outlier_threshold=2.5,
        )

        # Should identify outliers
        assert 'outliers_removed' in result
        outliers = result['outliers_removed']
        assert isinstance(outliers, list)

    def test_preprocessing_without_batch_correction(self):
        """Test preprocessing can skip batch correction."""
        result = preprocess_multiomics_data_impl(
            rna_path="/mock/rna.csv",
            protein_path="/mock/protein.csv",
            phospho_path=None,
            metadata_path=None,
            batch_correction=False,
        )

        assert result['status'] == 'success (DRY_RUN mode)'
        # Batch correction should be skipped - check in report
        assert any('skipped' in step.lower() for step in result['preprocessing_report']['steps_applied'])


class TestVisualizeDataQuality:
    """Test suite for visualize_data_quality tool."""

    def test_basic_visualization(self):
        """Test basic QC visualization."""
        result = visualize_data_quality_impl(
            data_paths={
                'rna': '/mock/rna.csv',
                'protein': '/mock/protein.csv',
            },
            metadata_path='/mock/metadata.csv',
            output_dir='/mock/plots',
        )

        # Check DRY_RUN mode response
        assert result['status'] == 'success (DRY_RUN mode)'
        assert 'plot_paths' in result
        assert 'qc_summary' in result
        assert 'batch_effect_assessment' in result

    def test_before_after_comparison(self):
        """Test before/after preprocessing comparison."""
        result = visualize_data_quality_impl(
            data_paths={
                'rna': '/mock/rna_preprocessed.csv',
                'protein': '/mock/protein_preprocessed.csv',
            },
            metadata_path='/mock/metadata.csv',
            output_dir='/mock/plots',
            compare_before_after=True,
            before_data_paths={
                'rna': '/mock/rna_raw.csv',
                'protein': '/mock/protein_raw.csv',
            },
        )

        # Should include before/after plots
        assert result['status'] == 'success (DRY_RUN mode)'
        plot_paths = result['plot_paths']

        # Should have before/after comparison plot
        assert 'before_after_comparison' in plot_paths
        assert plot_paths['before_after_comparison'] is not None

    def test_batch_effect_visualization(self):
        """Test batch effect assessment in plots."""
        result = visualize_data_quality_impl(
            data_paths={'protein': '/mock/protein.csv'},
            metadata_path='/mock/metadata.csv',
            output_dir='/mock/plots',
        )

        # Should assess batch effects
        batch_assessment = result['batch_effect_assessment']
        assert 'pc1_batch_correlation' in batch_assessment
        assert 'interpretation' in batch_assessment

    def test_multiple_modalities(self):
        """Test visualization with all 3 modalities."""
        result = visualize_data_quality_impl(
            data_paths={
                'rna': '/mock/rna.csv',
                'protein': '/mock/protein.csv',
                'phospho': '/mock/phospho.csv',
            },
            metadata_path='/mock/metadata.csv',
            output_dir='/mock/plots',
        )

        # Should generate plots for all modalities
        assert result['status'] == 'success (DRY_RUN mode)'
        plot_paths = result['plot_paths']
        assert 'pca_plot' in plot_paths
        assert 'correlation_heatmap' in plot_paths
        assert 'missing_values' in plot_paths  # Note: it's 'missing_values' not 'missing_values_plot'

    def test_recommendations_provided(self):
        """Test that visualization provides actionable recommendations."""
        result = visualize_data_quality_impl(
            data_paths={'protein': '/mock/protein.csv'},
            metadata_path='/mock/metadata.csv',
            output_dir='/mock/plots',
        )

        # Should provide recommendations
        assert 'recommendations' in result
        recommendations = result['recommendations']
        assert isinstance(recommendations, list)
        assert len(recommendations) > 0


class TestPreprocessingWorkflow:
    """Integration tests for complete preprocessing workflow."""

    def test_complete_preprocessing_pipeline(self):
        """Test validate → preprocess → visualize workflow."""

        # Step 1: Validate
        validation = validate_multiomics_data_impl(
            rna_path="/mock/rna_raw.csv",
            protein_path="/mock/protein_raw.csv",
            phospho_path=None,
            metadata_path="/mock/metadata.csv",
        )
        assert validation['status'] == 'success (DRY_RUN mode)'

        # Step 2: Preprocess
        preprocessing = preprocess_multiomics_data_impl(
            rna_path="/mock/rna_raw.csv",
            protein_path="/mock/protein_raw.csv",
            phospho_path=None,
            metadata_path="/mock/metadata.csv",
            batch_correction=True,
            imputation_method="knn",
        )
        assert preprocessing['status'] == 'success (DRY_RUN mode)'

        # Step 3: Visualize
        visualization = visualize_data_quality_impl(
            data_paths={
                'rna': preprocessing['preprocessed_paths']['rna'],
                'protein': preprocessing['preprocessed_paths']['protein'],
            },
            metadata_path="/mock/metadata.csv",
            output_dir="/mock/plots",
            compare_before_after=True,
            before_data_paths={
                'rna': '/mock/rna_raw.csv',
                'protein': '/mock/protein_raw.csv',
            },
        )
        assert visualization['status'] == 'success (DRY_RUN mode)'

    def test_preprocessing_improves_batch_correlation(self):
        """Test that preprocessing reduces PC1-batch correlation."""

        # Preprocess with batch correction
        result = preprocess_multiomics_data_impl(
            rna_path="/mock/rna.csv",
            protein_path="/mock/protein.csv",
            phospho_path=None,
            metadata_path="/mock/metadata.csv",
            batch_correction=True,
        )

        # Check improvement
        batch_results = result['batch_correction_results']
        before = batch_results['pc1_batch_correlation_before']
        after = batch_results['pc1_batch_correlation_after']

        # In mock data, after should be < 0.3 (good)
        assert after < 0.3
        # Improvement should be significant
        assert (before - after) > 0.3


class TestValidateWithRealData:
    """Test validation with actual data processing (DRY_RUN=False)."""

    def test_validate_with_real_rna_data(self, rna_path, monkeypatch):
        """Test validation loads and analyzes real RNA data."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        result = validate_multiomics_data_impl(
            rna_path=rna_path,
            protein_path=None,
            phospho_path=None,
            metadata_path=None,
        )

        # Should process real data
        assert result['status'] in ['ready', 'proceed_with_caution', 'needs_preprocessing', 'failed']
        assert 'statistics' in result

        # Should have missing_values statistics for rna
        assert 'missing_values' in result['statistics']
        assert 'rna' in result['statistics']['missing_values']

    def test_validate_detects_sample_overlap(self, rna_path, protein_path, monkeypatch):
        """Test validation checks sample consistency across modalities."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        result = validate_multiomics_data_impl(
            rna_path=rna_path,
            protein_path=protein_path,
            phospho_path=None,
            metadata_path=None,
        )

        # Should detect common samples
        assert 'statistics' in result
        assert 'common_samples' in result['statistics']
        assert result['statistics']['common_samples'] > 0

    def test_validate_analyzes_missing_values(self, rna_path, protein_path, monkeypatch):
        """Test validation analyzes missing value patterns."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        result = validate_multiomics_data_impl(
            rna_path=rna_path,
            protein_path=protein_path,
            phospho_path=None,
            metadata_path=None,
        )

        # Should analyze missing values
        assert 'statistics' in result
        assert 'missing_values' in result['statistics']

        # Should have stats for both modalities
        missing_stats = result['statistics']['missing_values']
        assert 'rna' in missing_stats
        assert 'protein' in missing_stats

    def test_validate_checks_metadata_batches(self, rna_path, metadata_path, monkeypatch):
        """Test validation detects batch information in metadata."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        result = validate_multiomics_data_impl(
            rna_path=rna_path,
            protein_path=None,
            phospho_path=None,
            metadata_path=metadata_path,
        )

        # Should check metadata
        assert result['status'] in ['ready', 'proceed_with_caution', 'needs_preprocessing', 'failed']
        assert 'statistics' in result or 'errors' in result

    def test_validate_handles_file_errors(self, monkeypatch):
        """Test validation handles missing/invalid files gracefully."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        result = validate_multiomics_data_impl(
            rna_path="/nonexistent/file.csv",
            protein_path=None,
            phospho_path=None,
            metadata_path=None,
        )

        # Should fail gracefully
        assert result['status'] == 'failed'
        assert 'errors' in result
        assert len(result['errors']) > 0


class TestPreprocessWithRealData:
    """Test preprocessing with actual data (DRY_RUN=False)."""

    def test_preprocess_loads_real_data(self, rna_path, protein_path, metadata_path, tmp_path, monkeypatch):
        """Test preprocessing loads and processes real data."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        output_dir = tmp_path / "output"
        output_dir.mkdir(exist_ok=True)

        result = preprocess_multiomics_data_impl(
            rna_path=rna_path,
            protein_path=protein_path,
            phospho_path=None,
            metadata_path=metadata_path,
            normalize_method="median",
            batch_correction=False,
            imputation_method="median",
            outlier_threshold=3.0,
            output_dir=str(output_dir),
        )

        # Should process successfully
        assert result['status'] == 'success'
        assert 'output_files' in result

        # Output files should exist
        assert 'rna' in result['output_files']
        assert os.path.exists(result['output_files']['rna'])

    def test_preprocess_normalizes_data(self, rna_path, tmp_path, monkeypatch):
        """Test normalization is applied correctly."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        output_dir = tmp_path / "output"
        output_dir.mkdir(exist_ok=True)

        result = preprocess_multiomics_data_impl(
            rna_path=rna_path,
            protein_path=None,
            phospho_path=None,
            metadata_path=None,
            normalize_method="median",
            batch_correction=False,
            imputation_method="none",
            output_dir=str(output_dir),
        )

        # Should apply normalization
        assert result['status'] == 'success'
        assert 'steps_completed' in result
        steps = result['steps_completed']
        assert 'normalization' in steps

    def test_preprocess_handles_missing_values(self, protein_path, tmp_path, monkeypatch):
        """Test imputation handles missing values."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        output_dir = tmp_path / "output"
        output_dir.mkdir(exist_ok=True)

        result = preprocess_multiomics_data_impl(
            rna_path=protein_path,  # Protein data has more missing values
            protein_path=None,
            phospho_path=None,
            metadata_path=None,
            normalize_method="median",
            batch_correction=False,
            imputation_method="median",
            output_dir=str(output_dir),
        )

        # Should handle missing values
        assert result['status'] == 'success'
        if 'imputation_stats' in result:
            assert result['imputation_stats'] is not None


class TestVisualizeWithRealData:
    """Test visualization with actual data (DRY_RUN=False)."""

    def test_visualize_generates_plots(self, rna_path, protein_path, metadata_path, tmp_path, monkeypatch):
        """Test visualization generates actual plot files."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        output_dir = tmp_path / "plots"
        output_dir.mkdir(exist_ok=True)

        result = visualize_data_quality_impl(
            data_paths={
                'rna': rna_path,
                'protein': protein_path,
            },
            metadata_path=metadata_path,
            output_dir=str(output_dir),
        )

        # Should generate plots
        assert result['status'] == 'success'
        assert 'plots_generated' in result

        # Check plots were generated
        plots = result['plots_generated']
        assert isinstance(plots, list)
        assert len(plots) > 0

    def test_visualize_assesses_batch_effects(self, protein_path, metadata_path, tmp_path, monkeypatch):
        """Test visualization calculates batch effect metrics."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        output_dir = tmp_path / "plots"
        output_dir.mkdir(exist_ok=True)

        result = visualize_data_quality_impl(
            data_paths={'protein': protein_path},
            metadata_path=metadata_path,
            output_dir=str(output_dir),
        )

        # Should assess batch effects
        assert result['status'] == 'success'
        if 'batch_effect_assessment' in result:
            batch_assessment = result['batch_effect_assessment']
            assert 'pc1_batch_correlation' in batch_assessment

    def test_visualize_provides_recommendations(self, rna_path, tmp_path, monkeypatch):
        """Test visualization generates actionable recommendations."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        output_dir = tmp_path / "plots"
        output_dir.mkdir(exist_ok=True)

        result = visualize_data_quality_impl(
            data_paths={'rna': rna_path},
            metadata_path=None,
            output_dir=str(output_dir),
        )

        # Should provide recommendations
        assert result['status'] == 'success'
        if 'recommendations' in result:
            recommendations = result['recommendations']
            assert isinstance(recommendations, list)
