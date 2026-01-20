"""Tests for scGen wrapper functionality."""

import pytest
import numpy as np
from mcp_perturbation.scgen_wrapper import ScGenWrapper
from mcp_perturbation.data_loader import DatasetLoader


@pytest.fixture
def test_adata():
    """Create test dataset."""
    loader = DatasetLoader()
    adata = loader.load_geo_dataset("GSE184880", normalize=True, n_hvg=2000)
    return adata


class TestScGenWrapper:
    """Test ScGenWrapper class."""

    def test_setup(self, test_adata):
        """Test model setup."""
        wrapper = ScGenWrapper()
        wrapper.setup(test_adata, batch_key="condition", labels_key="cell_type")

        assert wrapper.adata is not None
        assert wrapper.adata.n_obs == test_adata.n_obs

    def test_initialize_model(self, test_adata):
        """Test model initialization."""
        wrapper = ScGenWrapper()
        wrapper.setup(test_adata, batch_key="condition", labels_key="cell_type")

        config = wrapper.initialize_model(n_latent=10, n_hidden=64)

        assert wrapper.model is not None
        assert config["n_latent"] == 10
        assert config["n_hidden"] == 64

    def test_train_small_epochs(self, test_adata):
        """Test model training with small number of epochs."""
        wrapper = ScGenWrapper()
        wrapper.setup(test_adata, batch_key="condition", labels_key="cell_type")
        wrapper.initialize_model(n_latent=10, n_hidden=64)

        metrics = wrapper.train(max_epochs=5, batch_size=32, early_stopping=False)

        assert "final_loss" in metrics
        assert "epochs_completed" in metrics
        assert metrics["epochs_completed"] <= 5

    def test_save_and_load_model(self, test_adata, tmp_path):
        """Test model saving and loading."""
        wrapper = ScGenWrapper(model_dir=str(tmp_path))
        wrapper.setup(test_adata, batch_key="condition", labels_key="cell_type")
        wrapper.initialize_model(n_latent=10, n_hidden=64)
        wrapper.train(max_epochs=2, batch_size=32, early_stopping=False)

        # Save model
        model_path = wrapper.save("test_model")
        assert "test_model" in model_path

        # Load model
        wrapper2 = ScGenWrapper(model_dir=str(tmp_path))
        wrapper2.load("test_model", test_adata)

        assert wrapper2.model is not None
        assert wrapper2.model_name == "test_model"

    def test_get_latent_representation(self, test_adata):
        """Test extracting latent representations."""
        wrapper = ScGenWrapper()
        wrapper.setup(test_adata, batch_key="condition", labels_key="cell_type")
        wrapper.initialize_model(n_latent=10, n_hidden=64)
        wrapper.train(max_epochs=2, batch_size=32, early_stopping=False)

        latent = wrapper.get_latent_representation(test_adata)

        assert latent.shape[0] == test_adata.n_obs
        assert latent.shape[1] == 10  # n_latent

    def test_compute_delta(self, test_adata):
        """Test computing perturbation vector."""
        wrapper = ScGenWrapper()
        wrapper.setup(test_adata, batch_key="condition", labels_key="cell_type")
        wrapper.initialize_model(n_latent=10, n_hidden=64)
        wrapper.train(max_epochs=2, batch_size=32, early_stopping=False)

        delta_stats = wrapper.compute_delta(
            ctrl_key="control",
            stim_key="tumor",
            cell_type="T_cells"
        )

        assert "delta_norm" in delta_stats
        assert "ctrl_cells" in delta_stats
        assert "stim_cells" in delta_stats
        assert delta_stats["ctrl_cells"] > 0
        assert delta_stats["stim_cells"] > 0

    def test_predict_response(self, test_adata):
        """Test perturbation prediction."""
        wrapper = ScGenWrapper()
        wrapper.setup(test_adata, batch_key="condition", labels_key="cell_type")
        wrapper.initialize_model(n_latent=10, n_hidden=64)
        wrapper.train(max_epochs=3, batch_size=32, early_stopping=False)

        # Predict for T_cells
        predicted_adata, delta = wrapper.predict(
            ctrl_key="control",
            stim_key="tumor",
            celltype_to_predict="T_cells"
        )

        assert predicted_adata is not None
        assert predicted_adata.n_obs > 0
        assert delta is not None
        assert len(delta) == 10  # n_latent


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
