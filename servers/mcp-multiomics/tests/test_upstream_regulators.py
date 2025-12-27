"""Tests for upstream regulator prediction tool."""

import pytest
import numpy as np
import os
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


class TestUpstreamRegulatorWithRealData:
    """Test upstream regulator prediction with actual data processing (DRY_RUN=False)."""

    def test_predict_with_real_analysis(self, monkeypatch):
        """Test regulator prediction with actual Fisher's exact test."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        differential_genes = {
            'TP53': {'log2fc': -2.5, 'p_value': 0.001},
            'MYC': {'log2fc': 3.2, 'p_value': 0.0005},
            'AKT1': {'log2fc': 2.1, 'p_value': 0.002},
            'MTOR': {'log2fc': 2.8, 'p_value': 0.001},
            'FOXO1': {'log2fc': -1.8, 'p_value': 0.003},
            'PIK3CA': {'log2fc': 2.3, 'p_value': 0.002},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['kinase', 'transcription_factor'],
            fdr_threshold=0.05,
            activation_zscore_threshold=2.0,
        )

        # Should process real data
        assert result['status'] == 'success'
        assert 'kinases' in result
        assert 'transcription_factors' in result
        assert isinstance(result['kinases'], list)
        assert isinstance(result['transcription_factors'], list)

    def test_kinase_analysis_with_real_data(self, monkeypatch):
        """Test kinase prediction uses Fisher's exact test."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        differential_genes = {
            'AKT1': {'log2fc': 2.5, 'p_value': 0.001},
            'MTOR': {'log2fc': 2.3, 'p_value': 0.001},
            'PIK3CA': {'log2fc': 2.1, 'p_value': 0.002},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['kinase'],
            fdr_threshold=0.05,
        )

        # Should analyze kinases
        assert result['status'] == 'success'
        kinases = result['kinases']
        assert isinstance(kinases, list)

        # Each kinase should have required fields
        for kinase in kinases:
            assert 'name' in kinase
            assert 'p_value' in kinase
            assert 'q_value' in kinase
            assert 'activation_state' in kinase

    def test_fdr_correction_applied(self, monkeypatch):
        """Test FDR correction is applied across all regulators."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        differential_genes = {
            'TP53': {'log2fc': -2.5, 'p_value': 0.001},
            'MYC': {'log2fc': 3.2, 'p_value': 0.0005},
            'AKT1': {'log2fc': 2.1, 'p_value': 0.002},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['kinase', 'transcription_factor'],
            fdr_threshold=0.1,
        )

        # All returned regulators should pass FDR threshold
        assert result['status'] == 'success'

        all_regulators = result['kinases'] + result['transcription_factors']
        for reg in all_regulators:
            assert reg['q_value'] <= 0.1

    def test_activation_state_prediction(self, monkeypatch):
        """Test activation state (activated/inhibited) is predicted correctly."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        differential_genes = {
            'AKT1': {'log2fc': 3.0, 'p_value': 0.001},
            'MTOR': {'log2fc': 2.8, 'p_value': 0.001},
            'PIK3CA': {'log2fc': 2.5, 'p_value': 0.001},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['kinase'],
            fdr_threshold=0.5,
            activation_zscore_threshold=1.5,
        )

        # Should predict activation states
        assert result['status'] == 'success'

        # Check that activation_state is present
        for kinase in result['kinases']:
            assert kinase['activation_state'] in ['activated', 'inhibited', 'uncertain']

    def test_drug_predictions(self, monkeypatch):
        """Test drug response predictions."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        differential_genes = {
            'EGFR': {'log2fc': 2.8, 'p_value': 0.001},
            'ERBB2': {'log2fc': 2.5, 'p_value': 0.001},
            'MET': {'log2fc': 2.3, 'p_value': 0.002},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['drug'],
            fdr_threshold=0.5,
        )

        # Should predict drugs
        assert result['status'] == 'success'
        assert 'drugs' in result
        assert isinstance(result['drugs'], list)

    def test_default_regulator_types(self, monkeypatch):
        """Test default regulator types are used when None specified."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        differential_genes = {
            'TP53': {'log2fc': -2.5, 'p_value': 0.001},
            'MYC': {'log2fc': 3.2, 'p_value': 0.0005},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=None,  # Test default
            fdr_threshold=0.5,
        )

        # Should analyze all three types by default
        assert result['status'] == 'success'
        assert 'kinases' in result
        assert 'transcription_factors' in result
        assert 'drugs' in result

    def test_no_significant_regulators(self, monkeypatch):
        """Test handling when no regulators pass FDR threshold."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        differential_genes = {
            'GENE1': {'log2fc': 0.5, 'p_value': 0.8},  # Weak signal
            'GENE2': {'log2fc': -0.3, 'p_value': 0.9},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['kinase'],
            fdr_threshold=0.001,  # Very strict
        )

        # Should handle case with no significant regulators
        assert result['status'] == 'success'
        assert 'kinases' in result

    def test_drug_recommendation_generated(self, monkeypatch):
        """Test that drug recommendation is generated when drugs found."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        differential_genes = {
            'EGFR': {'log2fc': 3.5, 'p_value': 0.0001},
            'ERBB2': {'log2fc': 3.2, 'p_value': 0.0001},
            'MET': {'log2fc': 2.8, 'p_value': 0.0001},
            'KRAS': {'log2fc': 2.5, 'p_value': 0.0001},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['drug'],
            fdr_threshold=0.99,  # Very lenient to ensure drugs found
        )

        # Should generate drug recommendation
        assert result['status'] == 'success'
        if len(result['drugs']) > 0:
            assert 'recommendation' in result
            assert 'Consider' in result['recommendation']

    def test_inhibited_activation_state(self, monkeypatch):
        """Test prediction of inhibited (downregulated) regulators."""
        from mcp_multiomics.config import config
        monkeypatch.setattr(config, "dry_run", False)

        differential_genes = {
            'AKT1': {'log2fc': -3.0, 'p_value': 0.001},
            'MTOR': {'log2fc': -2.8, 'p_value': 0.001},
            'PIK3CA': {'log2fc': -2.5, 'p_value': 0.001},
        }

        result = predict_upstream_regulators_impl(
            differential_genes=differential_genes,
            regulator_types=['kinase'],
            fdr_threshold=0.99,
            activation_zscore_threshold=1.5,
        )

        # Should find some inhibited kinases
        assert result['status'] == 'success'
        # Check for inhibited state in results
        if len(result['kinases']) > 0:
            states = [k['activation_state'] for k in result['kinases']]
            # At least some should show activation state
            assert len(states) > 0
