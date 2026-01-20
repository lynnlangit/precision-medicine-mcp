"""Tests for prediction functionality."""

import pytest
from mcp_perturbation.prediction import (
    DifferentialExpressionAnalyzer,
    PerturbationPredictor,
    get_top_changed_genes
)
from mcp_perturbation.data_loader import DatasetLoader
from mcp_perturbation.scgen_wrapper import ScGenWrapper


@pytest.fixture
def test_adata():
    """Create test dataset."""
    loader = DatasetLoader()
    return loader.load_geo_dataset("GSE184880", normalize=True, n_hvg=1000)


@pytest.fixture
def trained_wrapper(test_adata):
    """Create and train a small model for testing."""
    wrapper = ScGenWrapper()
    wrapper.setup(test_adata, batch_key="condition", labels_key="cell_type")
    wrapper.initialize_model(n_latent=10, n_hidden=64)
    wrapper.train(max_epochs=3, batch_size=32, early_stopping=False)
    return wrapper


class TestDifferentialExpressionAnalyzer:
    """Test differential expression analysis."""

    def test_compute_de(self, test_adata):
        """Test computing differential expression."""
        # Split into baseline and predicted (mock)
        n_half = test_adata.n_obs // 2
        baseline = test_adata[:n_half]
        predicted = test_adata[n_half:]

        analyzer = DifferentialExpressionAnalyzer()
        results = analyzer.compute_differential_expression(
            baseline,
            predicted,
            n_top_genes=20,
            method="wilcoxon"
        )

        assert "upregulated_genes" in results
        assert "downregulated_genes" in results
        assert len(results["upregulated_genes"]) <= 20
        assert results["method"] == "wilcoxon"


class TestPerturbationPredictor:
    """Test perturbation prediction workflows."""

    def test_apply_perturbation(self, test_adata, trained_wrapper, tmp_path):
        """Test applying perturbation to patient data."""
        predictor = PerturbationPredictor(output_dir=str(tmp_path))

        result = predictor.apply_perturbation_to_patient(
            wrapper=trained_wrapper,
            patient_adata=test_adata,
            ctrl_key="control",
            stim_key="tumor",
            celltype_to_predict="T_cells",
            output_name="test_prediction"
        )

        assert result["status"] == "success"
        assert "output_path" in result
        assert result["cell_type"] == "T_cells"
        assert "delta_norm" in result


def test_get_top_changed_genes(test_adata):
    """Test getting top changed genes."""
    # Split data
    n_half = test_adata.n_obs // 2
    baseline = test_adata[:n_half]
    predicted = test_adata[n_half:]

    upregulated, downregulated = get_top_changed_genes(
        baseline,
        predicted,
        n_genes=10
    )

    assert len(upregulated) == 10
    assert len(downregulated) == 10
    assert all(isinstance(g, str) for g in upregulated)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
