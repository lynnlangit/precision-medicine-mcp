"""Wrapper for scGen model operations."""

import scanpy as sc
import scgen
import anndata as ad
import numpy as np
from pathlib import Path
from typing import Optional, Tuple, Dict
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class PredictionResult:
    """Results from perturbation prediction."""
    predicted_adata: ad.AnnData
    delta: np.ndarray
    output_path: str


class ScGenWrapper:
    """Manages scGen model lifecycle and predictions.

    scGen uses variational autoencoders to learn latent representations
    of cell states, enabling prediction of treatment responses via
    latent space arithmetic:

    Patient_predicted = Patient_baseline + Δ
    where Δ = Mean(Treated_cells) - Mean(Control_cells)
    """

    def __init__(self, model_dir: str = "./data/models"):
        """Initialize wrapper.

        Args:
            model_dir: Directory to save/load models
        """
        self.model_dir = Path(model_dir)
        self.model_dir.mkdir(parents=True, exist_ok=True)
        self.model: Optional[scgen.SCGEN] = None
        self.adata: Optional[ad.AnnData] = None
        self.model_name: Optional[str] = None

    def setup(
        self,
        adata: ad.AnnData,
        batch_key: str = "condition",
        labels_key: str = "cell_type"
    ) -> None:
        """Configure AnnData for scGen.

        Args:
            adata: Annotated data matrix
            batch_key: Column in adata.obs with condition labels (e.g., "control", "treated")
            labels_key: Column in adata.obs with cell type labels
        """
        self.adata = adata.copy()

        # Verify required keys exist
        if batch_key not in adata.obs.columns:
            raise ValueError(f"batch_key '{batch_key}' not found in adata.obs")
        if labels_key not in adata.obs.columns:
            raise ValueError(f"labels_key '{labels_key}' not found in adata.obs")

        # Setup anndata for scGen
        scgen.SCGEN.setup_anndata(
            self.adata,
            batch_key=batch_key,
            labels_key=labels_key
        )
        logger.info(f"Setup AnnData with batch_key={batch_key}, labels_key={labels_key}")

    def initialize_model(
        self,
        n_latent: int = 100,
        n_hidden: int = 800,
        n_layers: int = 2
    ) -> Dict[str, int]:
        """Create scGen model instance.

        Args:
            n_latent: Dimensionality of latent space
            n_hidden: Number of hidden units per layer
            n_layers: Number of hidden layers

        Returns:
            Model configuration summary
        """
        if self.adata is None:
            raise ValueError("Call setup() before initializing model")

        self.model = scgen.SCGEN(
            self.adata,
            n_latent=n_latent,
            n_hidden=n_hidden,
            n_layers=n_layers
        )

        config = {
            "n_latent": n_latent,
            "n_hidden": n_hidden,
            "n_layers": n_layers,
            "n_genes": self.adata.n_vars,
            "n_cells": self.adata.n_obs
        }

        logger.info(f"Initialized scGen model: {config}")
        return config

    def train(
        self,
        max_epochs: int = 100,
        batch_size: int = 32,
        early_stopping: bool = True,
        early_stopping_patience: int = 25,
        lr: float = 0.001
    ) -> Dict[str, float]:
        """Train the scGen model.

        Args:
            max_epochs: Maximum training epochs
            batch_size: Batch size for training
            early_stopping: Enable early stopping
            early_stopping_patience: Patience for early stopping
            lr: Learning rate

        Returns:
            Training metrics
        """
        if self.model is None:
            raise ValueError("Call initialize_model() before training")

        logger.info(f"Training scGen for up to {max_epochs} epochs")

        self.model.train(
            max_epochs=max_epochs,
            batch_size=batch_size,
            early_stopping=early_stopping,
            early_stopping_patience=early_stopping_patience,
            lr=lr
        )

        # Get final metrics
        final_loss = float(self.model.history["elbo_train"][-1])
        epochs_completed = len(self.model.history["elbo_train"])

        metrics = {
            "final_loss": final_loss,
            "epochs_completed": epochs_completed,
            "early_stopped": epochs_completed < max_epochs
        }

        logger.info(f"Training complete: {metrics}")
        return metrics

    def save(self, name: str, overwrite: bool = True) -> str:
        """Save model to disk.

        Args:
            name: Model name (without extension)
            overwrite: Overwrite existing model

        Returns:
            Path to saved model
        """
        if self.model is None:
            raise ValueError("No model to save")

        save_dir = self.model_dir / name
        self.model.save(str(save_dir), overwrite=overwrite)
        self.model_name = name

        logger.info(f"Saved model to {save_dir}")
        return str(save_dir)

    def load(self, name: str, adata: Optional[ad.AnnData] = None) -> None:
        """Load model from disk.

        Args:
            name: Model name (without extension)
            adata: Optional AnnData to use with model
        """
        load_dir = self.model_dir / name

        if not load_dir.exists():
            raise FileNotFoundError(f"Model not found: {load_dir}")

        # Use stored adata if not provided
        if adata is None and self.adata is not None:
            adata = self.adata
        elif adata is None:
            raise ValueError("Must provide adata when loading model")

        self.model = scgen.SCGEN.load(str(load_dir), adata)
        self.adata = adata
        self.model_name = name

        logger.info(f"Loaded model from {load_dir}")

    def predict(
        self,
        ctrl_key: str,
        stim_key: str,
        celltype_to_predict: str,
        restrict_arithmetic_to: str = "all"
    ) -> Tuple[ad.AnnData, np.ndarray]:
        """Generate perturbation prediction using latent space arithmetic.

        Computes: Patient_predicted = Patient_baseline + Δ
        where Δ = Mean(Treated_cells) - Mean(Control_cells)

        Args:
            ctrl_key: Label for control condition in batch_key column
            stim_key: Label for stimulated/treated condition
            celltype_to_predict: Cell type to predict response for
            restrict_arithmetic_to: Cell type to compute Δ from (or "all")

        Returns:
            Tuple of (predicted_adata, delta_vector)
        """
        if self.model is None:
            raise ValueError("No trained model available")

        logger.info(
            f"Predicting {celltype_to_predict} response "
            f"from {ctrl_key} to {stim_key}"
        )

        predicted_adata, delta = self.model.predict(
            ctrl_key=ctrl_key,
            stim_key=stim_key,
            celltype_to_predict=celltype_to_predict,
            restrict_arithmetic_to=restrict_arithmetic_to
        )

        logger.info(f"Prediction complete. Delta norm: {np.linalg.norm(delta):.4f}")

        return predicted_adata, delta

    def get_latent_representation(
        self,
        adata: Optional[ad.AnnData] = None
    ) -> np.ndarray:
        """Get latent space representation of cells.

        Args:
            adata: Data to embed (uses self.adata if None)

        Returns:
            Latent representations (n_cells x n_latent)
        """
        if self.model is None:
            raise ValueError("No trained model available")

        if adata is None:
            adata = self.adata

        if adata is None:
            raise ValueError("No data to embed")

        latent = self.model.get_latent_representation(adata)
        logger.info(f"Extracted latent representation: shape {latent.shape}")

        return latent

    def compute_delta(
        self,
        ctrl_key: str,
        stim_key: str,
        cell_type: Optional[str] = None
    ) -> Dict[str, any]:
        """Compute perturbation vector Δ between conditions.

        Args:
            ctrl_key: Control condition label
            stim_key: Stimulated condition label
            cell_type: Optional cell type to restrict to

        Returns:
            Dictionary with delta statistics
        """
        if self.adata is None or self.model is None:
            raise ValueError("Model not trained or data not loaded")

        # Get latent representations
        latent = self.get_latent_representation(self.adata)

        # Filter by cell type if specified
        if cell_type is not None:
            mask = self.adata.obs["cell_type"] == cell_type
            latent_subset = latent[mask]
            obs_subset = self.adata.obs[mask]
        else:
            latent_subset = latent
            obs_subset = self.adata.obs

        # Get condition labels
        batch_key = self.model.adata_manager.registry["setup_args"]["batch_key"]
        conditions = obs_subset[batch_key]

        # Compute means
        ctrl_mean = latent_subset[conditions == ctrl_key].mean(axis=0)
        stim_mean = latent_subset[conditions == stim_key].mean(axis=0)

        # Compute delta
        delta = stim_mean - ctrl_mean

        result = {
            "delta_norm": float(np.linalg.norm(delta)),
            "delta_mean": float(delta.mean()),
            "delta_std": float(delta.std()),
            "ctrl_cells": int((conditions == ctrl_key).sum()),
            "stim_cells": int((conditions == stim_key).sum()),
        }

        logger.info(f"Computed delta: {result}")
        return result
