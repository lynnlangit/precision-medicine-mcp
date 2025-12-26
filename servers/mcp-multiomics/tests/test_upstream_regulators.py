"""Tests for upstream regulator prediction tool."""

import pytest
import numpy as np
from mcp_multiomics.tools.upstream_regulators import predict_upstream_regulators_impl


class TestUpstreamRegulatorPrediction:
    """Test suite for predict_upstream_regulators tool."""

    def test_basic_prediction(self):
        """Test basic upstream regulator prediction."""
        # Create differential gene input (simplified)
        differential_genes = {
            'TP53': {'log2fc': -2.5, 'p_value': 0.001},
            'MYC': {'log2fc': 3.2, 'p_value': 0.0005},
            'AKT1': {'log2fc': 2.1, 'p_value': 0.002},
            'MTOR': {'log2fc': 2.8, 'p_value': 0.001},
            'FOXO1': {'log2fc': -1.8, 'p_value': 0.003},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['kinase', 'transcription_factor', 'drug'],
            fdr_threshold=0.05,
            activation_zscore_threshold=2.0,
        )

        # Check DRY_RUN mode response
        assert result['status'] == 'success (DRY_RUN mode)'
        assert 'kinases' in result
        assert 'transcription_factors' in result
        assert 'drugs' in result
        assert 'statistics' in result

    def test_kinase_prediction(self):
        """Test kinase-specific prediction."""
        differential_genes = {
            'AKT1': {'log2fc': 2.5, 'p_value': 0.001},
            'MTOR': {'log2fc': 2.3, 'p_value': 0.001},
            'GSK3B': {'log2fc': -1.5, 'p_value': 0.002},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['kinase'],
            fdr_threshold=0.05,
        )

        # Should predict kinases
        kinases = result['kinases']
        assert isinstance(kinases, list)
        assert len(kinases) > 0

        # Check kinase structure
        if kinases:
            kinase = kinases[0]
            assert 'name' in kinase  # Uses 'name' not 'regulator_name'
            assert 'activation_state' in kinase
            assert 'p_value' in kinase
            assert 'q_value' in kinase
            assert 'z_score' in kinase  # Uses 'z_score' not 'activation_zscore'

    def test_transcription_factor_prediction(self):
        """Test transcription factor prediction."""
        differential_genes = {
            'TP53': {'log2fc': -2.5, 'p_value': 0.001},
            'MYC': {'log2fc': 3.0, 'p_value': 0.0005},
            'FOXO1': {'log2fc': -1.8, 'p_value': 0.002},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['transcription_factor'],
            fdr_threshold=0.05,
        )

        # Should predict TFs
        tfs = result['transcription_factors']
        assert isinstance(tfs, list)

    def test_drug_response_prediction(self):
        """Test drug response prediction."""
        differential_genes = {
            'AKT1': {'log2fc': 2.5, 'p_value': 0.001},
            'MTOR': {'log2fc': 2.3, 'p_value': 0.001},
            'PIK3CA': {'log2fc': 2.0, 'p_value': 0.002},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['drug'],
            fdr_threshold=0.05,
        )

        # Should predict drugs
        drugs = result['drugs']
        assert isinstance(drugs, list)
        assert len(drugs) > 0

        # Check drug structure
        if drugs:
            drug = drugs[0]
            assert 'name' in drug  # Uses 'name' not 'drug_name'
            assert 'mechanism' in drug
            assert 'clinical_indication' in drug  # Drugs have 'clinical_indication' not 'activation_state'

    def test_activation_state_detection(self):
        """Test activation vs inhibition state detection."""
        # Genes with positive log2FC (activated pathway)
        differential_genes_activated = {
            'AKT1': {'log2fc': 3.0, 'p_value': 0.0001},
            'MTOR': {'log2fc': 2.5, 'p_value': 0.0002},
            'RPS6KB1': {'log2fc': 2.8, 'p_value': 0.0001},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes_activated,
            regulator_types=['kinase'],
            activation_zscore_threshold=2.0,
        )

        # In mock data, should predict activation
        kinases = result['kinases']
        if kinases:
            # At least some kinases should be "Activated"
            activation_states = [k['activation_state'] for k in kinases]
            assert 'Activated' in activation_states

    def test_fdr_threshold_filtering(self):
        """Test FDR threshold filters results correctly."""
        differential_genes = {
            f'GENE{i}': {'log2fc': 2.0, 'p_value': 0.001}
            for i in range(10)
        }

        # Strict threshold
        result_strict = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['kinase'],
            fdr_threshold=0.01,
        )

        # Relaxed threshold
        result_relaxed = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['kinase'],
            fdr_threshold=0.1,
        )

        # Relaxed should have more or equal predictions
        assert len(result_relaxed['kinases']) >= len(result_strict['kinases'])

    def test_activation_zscore_threshold(self):
        """Test activation Z-score threshold filtering."""
        differential_genes = {
            'TP53': {'log2fc': -2.5, 'p_value': 0.001},
            'MYC': {'log2fc': 3.2, 'p_value': 0.0005},
        }

        # High threshold (only strong activators/inhibitors)
        result_high = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['kinase'],
            activation_zscore_threshold=3.0,
        )

        # Low threshold (include weak activators/inhibitors)
        result_low = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['kinase'],
            activation_zscore_threshold=1.0,
        )

        # Lower threshold should include more predictions
        stats_high = result_high['statistics']
        stats_low = result_low['statistics']

        # Both should have counts
        assert 'kinases_tested' in stats_high  # Uses specific type counts not 'total_regulators_tested'
        assert 'kinases_tested' in stats_low

    def test_statistics_summary(self):
        """Test statistics summary is comprehensive."""
        differential_genes = {
            'TP53': {'log2fc': -2.5, 'p_value': 0.001},
            'MYC': {'log2fc': 3.2, 'p_value': 0.0005},
            'AKT1': {'log2fc': 2.1, 'p_value': 0.002},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['kinase', 'transcription_factor', 'drug'],
        )

        stats = result['statistics']
        assert 'total_genes_analyzed' in stats  # Uses 'total_genes_analyzed' not 'input_genes'
        assert 'kinases_tested' in stats  # Has specific counts not 'total_regulators_tested'
        assert 'significant_kinases' in stats
        assert 'significant_tfs' in stats
        assert 'significant_drugs' in stats
        assert 'method' in stats

    def test_method_documentation(self):
        """Test method field documents the approach."""
        differential_genes = {
            'TP53': {'log2fc': -2.5, 'p_value': 0.001},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
        )

        # Should document method used
        assert 'method' in result
        method = result['method']
        assert 'enrichment_test' in method
        assert 'activation_score' in method  # Uses 'activation_score' not 'activation_zscore'

    def test_recommendation_provided(self):
        """Test therapeutic recommendation is provided."""
        differential_genes = {
            'AKT1': {'log2fc': 2.5, 'p_value': 0.001},
            'MTOR': {'log2fc': 2.3, 'p_value': 0.001},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['drug'],
        )

        # Should provide recommendation
        assert 'recommendation' in result
        recommendation = result['recommendation']
        assert isinstance(recommendation, str)
        assert len(recommendation) > 0

    def test_empty_input(self):
        """Test handling of empty differential genes."""
        result = predict_upstream_regulators_impl(
            differential_genes={},
        )

        # Should handle gracefully
        assert result['status'] == 'success (DRY_RUN mode)'
        assert result['statistics']['total_genes_analyzed'] == 0  # Uses 'total_genes_analyzed'

    def test_few_genes_input(self):
        """Test with very few input genes."""
        differential_genes = {
            'TP53': {'log2fc': -2.5, 'p_value': 0.001},
            'MYC': {'log2fc': 3.0, 'p_value': 0.001},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
        )

        # Should still work but may have fewer predictions
        assert result['status'] == 'success (DRY_RUN mode)'
        assert result['statistics']['total_genes_analyzed'] == 2  # Uses 'total_genes_analyzed'

    def test_all_regulator_types(self):
        """Test requesting all regulator types simultaneously."""
        differential_genes = {
            'TP53': {'log2fc': -2.5, 'p_value': 0.001},
            'MYC': {'log2fc': 3.2, 'p_value': 0.0005},
            'AKT1': {'log2fc': 2.1, 'p_value': 0.002},
            'MTOR': {'log2fc': 2.8, 'p_value': 0.001},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['kinase', 'transcription_factor', 'drug'],
        )

        # Should return all three types
        assert 'kinases' in result
        assert 'transcription_factors' in result
        assert 'drugs' in result

        # Statistics should reflect all types
        stats = result['statistics']
        assert stats['significant_kinases'] >= 0
        assert stats['significant_tfs'] >= 0
        assert stats['significant_drugs'] >= 0


class TestUpstreamRegulatorIntegration:
    """Integration tests combining upstream regulators with other tools."""

    def test_stouffer_to_upstream_workflow(self):
        """Test Stouffer's meta-analysis → Upstream regulators workflow."""

        # Simulate significant genes from Stouffer's meta-analysis
        # (genes with q < 0.05)
        significant_genes = {
            'TP53': {'log2fc': -2.5, 'p_value': 0.001},
            'MYC': {'log2fc': 3.2, 'p_value': 0.0005},
            'AKT1': {'log2fc': 2.1, 'p_value': 0.002},
            'MTOR': {'log2fc': 2.8, 'p_value': 0.001},
            'FOXO1': {'log2fc': -1.8, 'p_value': 0.003},
            'PIK3CA': {'log2fc': 2.3, 'p_value': 0.002},
            'EGFR': {'log2fc': 1.9, 'p_value': 0.004},
        }

        # Predict upstream regulators
        result = predict_upstream_regulators_impl(
            differential_genes=significant_genes,
            regulator_types=['kinase', 'transcription_factor', 'drug'],
            fdr_threshold=0.05,
        )

        # Should successfully predict regulators
        assert result['status'] == 'success (DRY_RUN mode)'
        assert result['statistics']['total_genes_analyzed'] == len(significant_genes)  # Uses 'total_genes_analyzed'

        # Should have therapeutic recommendations
        assert 'recommendation' in result
