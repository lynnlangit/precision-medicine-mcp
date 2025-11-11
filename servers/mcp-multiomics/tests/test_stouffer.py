"""Tests for Stouffer's meta-analysis."""

import pytest
import numpy as np

from mcp_multiomics.tools.stouffer import (
    StoufferMetaAnalysis,
    calculate_stouffer_meta_impl,
)


class TestStoufferMetaAnalysis:
    """Tests for StoufferMetaAnalysis class."""

    def test_p_to_z_conversion(self):
        """Test p-value to Z-score conversion."""
        analyzer = StoufferMetaAnalysis(use_directionality=False)

        p_values = np.array([0.05, 0.01, 0.001])
        z_scores = analyzer.p_to_z(p_values)

        # All Z-scores should be positive (two-tailed)
        assert all(z_scores > 0)

        # Smaller p-values should give larger Z-scores
        assert z_scores[2] > z_scores[1] > z_scores[0]

    def test_p_to_z_with_directionality(self):
        """Test p-value to Z-score conversion with effect size directionality."""
        analyzer = StoufferMetaAnalysis(use_directionality=True)

        p_values = np.array([0.01, 0.01, 0.01])
        effect_sizes = np.array([2.0, -2.0, 0.5])  # Mixed directions

        z_scores = analyzer.p_to_z(p_values, effect_sizes)

        # Positive effect size -> positive Z-score
        assert z_scores[0] > 0

        # Negative effect size -> negative Z-score
        assert z_scores[1] < 0

        # Magnitudes should be similar for same p-value
        assert abs(abs(z_scores[0]) - abs(z_scores[1])) < 0.01

    def test_z_to_p_conversion(self):
        """Test Z-score to p-value conversion."""
        analyzer = StoufferMetaAnalysis()

        z_scores = np.array([1.96, 2.58, 3.29])
        p_values = analyzer.z_to_p(z_scores)

        # Check approximate p-values for known Z-scores
        assert abs(p_values[0] - 0.05) < 0.01  # Z=1.96 ≈ p=0.05
        assert abs(p_values[1] - 0.01) < 0.005  # Z=2.58 ≈ p=0.01
        assert abs(p_values[2] - 0.001) < 0.001  # Z=3.29 ≈ p=0.001

    def test_combine_z_scores_equal_weights(self):
        """Test combining Z-scores with equal weights."""
        analyzer = StoufferMetaAnalysis()

        z_scores = np.array([2.0, 3.0, 4.0])
        z_meta = analyzer.combine_z_scores(z_scores)

        # With equal weights: Z_meta = sum(Z) / sqrt(n)
        expected = (2.0 + 3.0 + 4.0) / np.sqrt(3)
        assert abs(z_meta - expected) < 0.001

    def test_combine_z_scores_custom_weights(self):
        """Test combining Z-scores with custom weights."""
        analyzer = StoufferMetaAnalysis()

        z_scores = np.array([2.0, 3.0])
        weights = np.array([1.0, 2.0])  # Second modality has more weight

        z_meta = analyzer.combine_z_scores(z_scores, weights)

        # Z_meta = (1*2 + 2*3) / sqrt(1^2 + 2^2) = 8 / sqrt(5)
        expected = 8.0 / np.sqrt(5)
        assert abs(z_meta - expected) < 0.001

    def test_meta_analyze_basic(self, sample_p_values):
        """Test basic meta-analysis."""
        analyzer = StoufferMetaAnalysis(use_directionality=False)

        result = analyzer.meta_analyze(sample_p_values)

        assert len(result["meta_p_values"]) == 4
        assert len(result["meta_z_scores"]) == 4
        assert len(result["q_values"]) == 4
        assert "n_significant" in result

        # First feature has low p-values in all modalities -> should be significant
        assert result["meta_p_values"][0] < 0.05

    def test_meta_analyze_with_directionality(
        self, sample_p_values, sample_effect_sizes
    ):
        """Test meta-analysis with effect size directionality."""
        analyzer = StoufferMetaAnalysis(use_directionality=True)

        result = analyzer.meta_analyze(
            p_values_dict=sample_p_values,
            effect_sizes_dict=sample_effect_sizes,
        )

        assert "meta_z_scores" in result

        # Features with consistent positive effects should have positive Z
        assert result["meta_z_scores"][0] > 0

        # Features with consistent negative effects should have negative Z
        assert result["meta_z_scores"][2] < 0

    def test_meta_analyze_with_weights(self, sample_p_values):
        """Test meta-analysis with custom weights."""
        analyzer = StoufferMetaAnalysis()

        # Weight RNA more heavily
        weights = {"rna": 2.0, "protein": 1.0, "phospho": 1.0}

        result = analyzer.meta_analyze(
            p_values_dict=sample_p_values,
            weights=weights,
        )

        assert len(result["meta_p_values"]) == 4

    def test_meta_analyze_mismatched_features(self):
        """Test error handling for mismatched feature counts."""
        analyzer = StoufferMetaAnalysis()

        p_values_dict = {
            "rna": [0.01, 0.05],
            "protein": [0.02, 0.04, 0.06],  # Different length
        }

        with pytest.raises(ValueError, match="expected 2"):
            analyzer.meta_analyze(p_values_dict)


class TestStoufferImplementation:
    """Tests for calculate_stouffer_meta_impl function."""

    def test_basic_meta_analysis(self, sample_p_values):
        """Test basic meta-analysis implementation."""
        result = calculate_stouffer_meta_impl(
            p_values_dict=sample_p_values,
            use_directionality=False,
        )

        assert result["status"] == "success"
        assert "meta_p_values" in result
        assert "meta_z_scores" in result
        assert "significant_features" in result
        assert "statistics" in result

        # Check statistics
        stats = result["statistics"]
        assert stats["total_features"] == 4
        assert stats["directionality_used"] is False
        assert "modalities" in stats

    def test_meta_analysis_with_directionality(
        self, sample_p_values, sample_effect_sizes
    ):
        """Test meta-analysis with directionality."""
        result = calculate_stouffer_meta_impl(
            p_values_dict=sample_p_values,
            effect_sizes_dict=sample_effect_sizes,
            use_directionality=True,
        )

        assert result["statistics"]["directionality_used"] is True

        # Check that significant features include effect sizes
        if len(result["significant_features"]) > 0:
            feat = result["significant_features"][0]
            assert "effect_sizes" in feat

    def test_meta_analysis_with_weights(self, sample_p_values):
        """Test meta-analysis with custom weights."""
        weights = {"rna": 2.0, "protein": 1.0, "phospho": 1.0}

        result = calculate_stouffer_meta_impl(
            p_values_dict=sample_p_values,
            weights=weights,
        )

        assert result["statistics"]["weights_used"] is True

    def test_meta_analysis_detailed_features(self, sample_p_values):
        """Test detailed information for significant features."""
        result = calculate_stouffer_meta_impl(
            p_values_dict=sample_p_values,
            fdr_threshold=0.1,  # Lenient threshold to get some hits
        )

        if len(result["significant_features"]) > 0:
            feat = result["significant_features"][0]

            # Check feature structure
            assert "feature_index" in feat
            assert "meta_p" in feat
            assert "meta_z" in feat
            assert "q_value" in feat
            assert "modality_contributions" in feat

            # Check modality contributions
            assert "rna" in feat["modality_contributions"]
            assert "protein" in feat["modality_contributions"]
            assert "phospho" in feat["modality_contributions"]

    def test_fdr_correction(self):
        """Test FDR correction is applied."""
        # Create data where some features are clearly significant
        p_values_dict = {
            "rna": [0.001, 0.001, 0.5, 0.6, 0.7],
            "protein": [0.002, 0.002, 0.55, 0.65, 0.75],
        }

        result = calculate_stouffer_meta_impl(
            p_values_dict=p_values_dict,
            fdr_threshold=0.05,
        )

        # Should have q-values
        assert len(result["q_values"]) == 5

        # Q-values should be >= p-values (FDR correction)
        for i in range(5):
            assert result["q_values"][i] >= result["meta_p_values"][i]

    def test_no_significant_features(self):
        """Test when no features are significant."""
        # All high p-values
        p_values_dict = {
            "rna": [0.8, 0.9, 0.95],
            "protein": [0.85, 0.88, 0.92],
        }

        result = calculate_stouffer_meta_impl(
            p_values_dict=p_values_dict,
            fdr_threshold=0.05,
        )

        assert result["statistics"]["significant_features"] == 0
        assert len(result["significant_features"]) == 0
