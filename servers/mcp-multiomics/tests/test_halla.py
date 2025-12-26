"""Tests for HAllA association analysis with chunking strategy."""

import pytest
from mcp_multiomics.tools.halla import run_halla_analysis_impl


class TestHAllAAnalysis:
    """Test suite for HAllA association testing."""

    def test_basic_halla_analysis(self):
        """Test basic HAllA analysis between two modalities."""
        result = run_halla_analysis_impl(
            data_path="/mock/integrated_data.pkl",
            modality1="rna",
            modality2="protein",
            fdr_threshold=0.05,
            method="spearman",
        )

        # Check DRY_RUN mode response
        assert result['status'] == 'success (DRY_RUN mode)'
        assert 'associations' in result
        assert 'statistics' in result
        assert 'nominal_p_values' in result
        assert result['nominal_p_values'] is True

    def test_chunking_strategy(self):
        """Test HAllA uses chunking for large datasets."""
        result = run_halla_analysis_impl(
            data_path="/mock/integrated_data.pkl",
            modality1="rna",
            modality2="protein",
            chunk_size=1000,
        )

        # Should include chunking information
        assert 'chunks_processed' in result
        chunks = result['chunks_processed']
        assert 'total_chunks' in chunks
        assert 'chunk_size' in chunks
        assert chunks['chunk_size'] == 1000

    def test_different_chunk_sizes(self):
        """Test HAllA with different chunk sizes."""
        for chunk_size in [500, 1000, 2000]:
            result = run_halla_analysis_impl(
                data_path="/mock/integrated_data.pkl",
                modality1="rna",
                modality2="protein",
                chunk_size=chunk_size,
            )

            assert result['status'] == 'success (DRY_RUN mode)'
            assert result['chunks_processed']['chunk_size'] == chunk_size

    def test_correlation_methods(self):
        """Test different correlation methods."""
        for method in ['spearman', 'pearson']:
            result = run_halla_analysis_impl(
                data_path="/mock/integrated_data.pkl",
                modality1="rna",
                modality2="protein",
                method=method,
            )

            assert result['status'] == 'success (DRY_RUN mode)'
            assert result['statistics']['method'] == method

    def test_nominal_p_values_returned(self):
        """Test that HAllA returns NOMINAL p-values (not FDR-corrected)."""
        result = run_halla_analysis_impl(
            data_path="/mock/integrated_data.pkl",
            modality1="rna",
            modality2="protein",
        )

        # Should explicitly flag as nominal
        assert result['nominal_p_values'] is True
        assert 'NOMINAL' in result['statistics']['p_value_type']

        # Should recommend FDR after Stouffer's
        assert 'recommendation' in result
        assert 'Stouffer' in result['recommendation']

    def test_associations_structure(self):
        """Test structure of returned associations."""
        result = run_halla_analysis_impl(
            data_path="/mock/integrated_data.pkl",
            modality1="rna",
            modality2="protein",
        )

        associations = result['associations']
        assert isinstance(associations, list)
        assert len(associations) > 0

        # Check first association structure
        assoc = associations[0]
        assert 'feature1' in assoc
        assert 'feature2' in assoc
        assert 'correlation' in assoc
        assert 'p_value_nominal' in assoc  # Explicitly labeled as NOMINAL
        assert 'chunk_id' in assoc

    def test_different_modality_pairs(self):
        """Test HAllA with different modality combinations."""
        modality_pairs = [
            ('rna', 'protein'),
            ('rna', 'phospho'),
            ('protein', 'phospho'),
        ]

        for mod1, mod2 in modality_pairs:
            result = run_halla_analysis_impl(
                data_path="/mock/integrated_data.pkl",
                modality1=mod1,
                modality2=mod2,
            )

            assert result['status'] == 'success (DRY_RUN mode)'

    def test_hierarchical_clustering(self):
        """Test hierarchical clustering of associations."""
        result = run_halla_analysis_impl(
            data_path="/mock/integrated_data.pkl",
            modality1="rna",
            modality2="protein",
        )

        # Should include clustering results
        if 'clusters' in result:
            clusters = result['clusters']
            assert 'modality1_clusters' in clusters or 'clusters_found' in clusters

    def test_chunk_processing_info(self):
        """Test detailed chunk processing information."""
        result = run_halla_analysis_impl(
            data_path="/mock/integrated_data.pkl",
            modality1="rna",
            modality2="protein",
            chunk_size=1000,
        )

        chunks = result['chunks_processed']
        assert 'strategy' in chunks
        assert '5 min' in chunks['strategy']  # Should mention performance benefit
        assert 'total_features_modality1' in chunks
        assert 'total_features_modality2' in chunks

    def test_statistics_summary(self):
        """Test comprehensive statistics summary."""
        result = run_halla_analysis_impl(
            data_path="/mock/integrated_data.pkl",
            modality1="rna",
            modality2="protein",
        )

        stats = result['statistics']
        assert 'method' in stats
        assert 'total_associations_tested' in stats
        assert 'p_value_type' in stats
        assert 'NOMINAL' in stats['p_value_type']


class TestHAllAIntegration:
    """Integration tests for HAllA with other tools."""

    def test_halla_to_stouffer_workflow(self):
        """Test HAllA → Stouffer's meta-analysis workflow."""

        # Run HAllA for RNA-Protein
        halla_rna_protein = run_halla_analysis_impl(
            data_path="/mock/integrated_data.pkl",
            modality1="rna",
            modality2="protein",
            chunk_size=1000,
        )

        # Run HAllA for RNA-Phospho
        halla_rna_phospho = run_halla_analysis_impl(
            data_path="/mock/integrated_data.pkl",
            modality1="rna",
            modality2="phospho",
            chunk_size=1000,
        )

        # Both should return NOMINAL p-values
        assert halla_rna_protein['nominal_p_values'] is True
        assert halla_rna_phospho['nominal_p_values'] is True

        # Extract p-values for Stouffer's
        p_values_protein = [a['p_value_nominal'] for a in halla_rna_protein['associations']]
        p_values_phospho = [a['p_value_nominal'] for a in halla_rna_phospho['associations']]

        # Should have p-values ready for Stouffer's
        assert len(p_values_protein) > 0
        assert len(p_values_phospho) > 0
        assert all(0 <= p <= 1 for p in p_values_protein)
        assert all(0 <= p <= 1 for p in p_values_phospho)

    def test_chunking_performance_estimate(self):
        """Test chunking provides performance estimates."""
        result = run_halla_analysis_impl(
            data_path="/mock/integrated_data.pkl",
            modality1="rna",
            modality2="protein",
            chunk_size=1000,
        )

        chunks = result['chunks_processed']

        # Should give performance context
        assert 'strategy' in chunks
        strategy = chunks['strategy']
        assert 'min' in strategy.lower() or 'chunk' in strategy.lower()
