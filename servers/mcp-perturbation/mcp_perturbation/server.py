"""MCP Server for perturbation prediction using scGen."""

from mcp.server.fastmcp import FastMCP
from pydantic import BaseModel, Field, ConfigDict
from typing import Optional, Literal, List
import json
import logging
import scanpy as sc

from .data_loader import load_geo_dataset
from .scgen_wrapper import ScGenWrapper
from .prediction import PerturbationPredictor, DifferentialExpressionAnalyzer
from .visualization import PerturbationVisualizer

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize FastMCP server
mcp = FastMCP("perturbation")

# Global state for models and datasets
_models = {}  # name -> ScGenWrapper
_datasets = {}  # dataset_id -> AnnData


# ==================== Input Models ====================

class LoadDatasetInput(BaseModel):
    """Input for loading scRNA-seq dataset."""
    model_config = ConfigDict(str_strip_whitespace=True, extra='forbid')

    dataset_id: str = Field(..., description="GEO accession (e.g., 'GSE184880') or path to .h5ad")
    normalize: bool = Field(default=True, description="Apply normalize_total + log1p")
    n_hvg: int = Field(default=7000, description="Number of highly variable genes", ge=1000, le=20000)
    cell_type_key: str = Field(default="cell_type", description="Column with cell type labels")
    condition_key: str = Field(default="condition", description="Column with treatment condition")


class SetupModelInput(BaseModel):
    """Input for initializing scGen model."""
    model_config = ConfigDict(str_strip_whitespace=True, extra='forbid')

    dataset_id: str = Field(..., description="Dataset ID from load_dataset")
    n_latent: int = Field(default=100, description="Latent space dimensions", ge=10, le=500)
    n_hidden: int = Field(default=800, description="Hidden layer size", ge=128, le=2048)
    n_layers: int = Field(default=2, description="Number of hidden layers", ge=1, le=5)
    model_name: str = Field(default="scgen_model", description="Name for this model")
    batch_key: str = Field(default="condition", description="Column with condition labels")
    labels_key: str = Field(default="cell_type", description="Column with cell type labels")


class TrainModelInput(BaseModel):
    """Input for training scGen model."""
    model_config = ConfigDict(str_strip_whitespace=True, extra='forbid')

    model_name: str = Field(..., description="Model name from setup_model")
    n_epochs: int = Field(default=100, description="Training epochs", ge=10, le=500)
    batch_size: int = Field(default=32, description="Batch size", ge=16, le=256)
    early_stopping: bool = Field(default=True, description="Enable early stopping")
    early_stopping_patience: int = Field(default=25, description="Patience for early stopping", ge=5, le=100)
    learning_rate: float = Field(default=0.001, description="Learning rate", gt=0, le=0.1)


class ComputeDeltaInput(BaseModel):
    """Input for computing perturbation vector."""
    model_config = ConfigDict(str_strip_whitespace=True, extra='forbid')

    model_name: str = Field(..., description="Trained model name")
    source_cell_type: Optional[str] = Field(None, description="Cell type to compute delta from (None = all)")
    control_key: str = Field(default="control", description="Control condition label")
    treatment_key: str = Field(default="tumor", description="Treatment condition label")


class PredictResponseInput(BaseModel):
    """Input for predicting treatment response."""
    model_config = ConfigDict(str_strip_whitespace=True, extra='forbid')

    model_name: str = Field(..., description="Trained model name")
    patient_data_path: str = Field(..., description="Path to patient .h5ad file")
    cell_type_to_predict: str = Field(..., description="Cell type to transform")
    control_key: str = Field(default="control", description="Control condition label")
    treatment_key: str = Field(default="tumor", description="Treatment condition label")
    output_path: Optional[str] = Field(default=None, description="Path to save predictions")


class DEInput(BaseModel):
    """Input for differential expression analysis."""
    model_config = ConfigDict(str_strip_whitespace=True, extra='forbid')

    baseline_path: str = Field(..., description="Baseline .h5ad")
    predicted_path: str = Field(..., description="Predicted .h5ad")
    n_top_genes: int = Field(default=50, ge=10, le=500, description="Number of top genes to return")
    method: Literal["wilcoxon", "t-test"] = Field(default="wilcoxon", description="Statistical test method")


class GetLatentInput(BaseModel):
    """Input for extracting latent representations."""
    model_config = ConfigDict(str_strip_whitespace=True, extra='forbid')

    model_name: str = Field(..., description="Trained model name")
    data_path: str = Field(..., description=".h5ad file to embed")


class VisualizeInput(BaseModel):
    """Input for visualization."""
    model_config = ConfigDict(str_strip_whitespace=True, extra='forbid')

    baseline_path: str = Field(..., description="Baseline .h5ad file path")
    predicted_path: str = Field(..., description="Predicted .h5ad file path")
    plot_type: Literal["pca", "umap"] = Field(default="pca", description="Type of plot")
    color_by: str = Field(default="condition", description="Column to color by")
    output_path: Optional[str] = Field(default=None, description="Path to save figure")


# ==================== Tool Implementations ====================

@mcp.tool()
async def perturbation_load_dataset(params: LoadDatasetInput) -> str:
    """Load scRNA-seq dataset from GEO or local file.

    Downloads and preprocesses single-cell data, applying normalization
    and highly variable gene selection for scGen training.

    Example:
        {"dataset_id": "GSE184880", "normalize": true, "n_hvg": 7000}
    """
    try:
        adata = await load_geo_dataset(
            params.dataset_id,
            normalize=params.normalize,
            n_hvg=params.n_hvg
        )

        # Store dataset
        key = params.dataset_id.replace("/", "_").replace(".h5ad", "")
        _datasets[key] = adata

        # Get unique values for cell types and conditions
        cell_types = list(adata.obs[params.cell_type_key].unique()) if params.cell_type_key in adata.obs.columns else []
        conditions = list(adata.obs[params.condition_key].unique()) if params.condition_key in adata.obs.columns else []

        return json.dumps({
            "status": "success",
            "dataset_id": key,
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "cell_types": cell_types,
            "conditions": conditions
        }, indent=2)

    except Exception as e:
        logger.error(f"Failed to load dataset: {e}")
        return json.dumps({"status": "error", "message": str(e)}, indent=2)


@mcp.tool()
async def perturbation_setup_model(params: SetupModelInput) -> str:
    """Initialize scGen model architecture.

    Sets up the variational autoencoder for learning cell state representations.

    Example:
        {"dataset_id": "GSE184880", "n_latent": 100, "model_name": "my_model"}
    """
    try:
        # Get dataset
        if params.dataset_id not in _datasets:
            return json.dumps({
                "status": "error",
                "message": f"Dataset '{params.dataset_id}' not found. Load it first."
            }, indent=2)

        adata = _datasets[params.dataset_id]

        # Create wrapper and setup
        wrapper = ScGenWrapper()
        wrapper.setup(
            adata=adata,
            batch_key=params.batch_key,
            labels_key=params.labels_key
        )

        # Initialize model
        config = wrapper.initialize_model(
            n_latent=params.n_latent,
            n_hidden=params.n_hidden,
            n_layers=params.n_layers
        )

        # Store model
        _models[params.model_name] = wrapper

        return json.dumps({
            "status": "success",
            "model_name": params.model_name,
            "configuration": config
        }, indent=2)

    except Exception as e:
        logger.error(f"Failed to setup model: {e}")
        return json.dumps({"status": "error", "message": str(e)}, indent=2)


@mcp.tool()
async def perturbation_train_model(params: TrainModelInput) -> str:
    """Train scGen VAE on control/treated data.

    Trains the model to learn latent representations that enable
    perturbation prediction via latent space arithmetic.

    Example:
        {"model_name": "my_model", "n_epochs": 100, "batch_size": 32}
    """
    try:
        if params.model_name not in _models:
            return json.dumps({
                "status": "error",
                "message": f"Model '{params.model_name}' not found. Setup it first."
            }, indent=2)

        wrapper = _models[params.model_name]

        # Train model
        metrics = wrapper.train(
            max_epochs=params.n_epochs,
            batch_size=params.batch_size,
            early_stopping=params.early_stopping,
            early_stopping_patience=params.early_stopping_patience,
            lr=params.learning_rate
        )

        # Save model
        model_path = wrapper.save(params.model_name)

        return json.dumps({
            "status": "success",
            "model_name": params.model_name,
            "model_path": model_path,
            "training_metrics": metrics
        }, indent=2)

    except Exception as e:
        logger.error(f"Failed to train model: {e}")
        return json.dumps({"status": "error", "message": str(e)}, indent=2)


@mcp.tool()
async def perturbation_compute_delta(params: ComputeDeltaInput) -> str:
    """Calculate perturbation vector (Δ) between conditions.

    Computes: Δ = Mean(Treated_cells) - Mean(Control_cells) in latent space.

    Example:
        {"model_name": "my_model", "source_cell_type": "T_cells", "control_key": "control", "treatment_key": "tumor"}
    """
    try:
        if params.model_name not in _models:
            return json.dumps({
                "status": "error",
                "message": f"Model '{params.model_name}' not found."
            }, indent=2)

        wrapper = _models[params.model_name]

        # Compute delta
        delta_stats = wrapper.compute_delta(
            ctrl_key=params.control_key,
            stim_key=params.treatment_key,
            cell_type=params.source_cell_type
        )

        return json.dumps({
            "status": "success",
            "delta_statistics": delta_stats,
            "formula": "Δ = Mean(Treated) - Mean(Control) in latent space"
        }, indent=2)

    except Exception as e:
        logger.error(f"Failed to compute delta: {e}")
        return json.dumps({"status": "error", "message": str(e)}, indent=2)


@mcp.tool()
async def perturbation_predict_response(params: PredictResponseInput) -> str:
    """Predict treatment response using latent space arithmetic.

    Applies the learned perturbation vector to patient cells:
    Patient_predicted = Patient_baseline + Δ

    Example:
        {"model_name": "my_model", "patient_data_path": "./data/patient_001.h5ad", "cell_type_to_predict": "T_cells"}
    """
    try:
        if params.model_name not in _models:
            return json.dumps({
                "status": "error",
                "message": f"Model '{params.model_name}' not found."
            }, indent=2)

        wrapper = _models[params.model_name]

        # Load patient data
        patient_adata = sc.read_h5ad(params.patient_data_path)

        # Create predictor
        predictor = PerturbationPredictor()

        # Make prediction
        result = predictor.apply_perturbation_to_patient(
            wrapper=wrapper,
            patient_adata=patient_adata,
            ctrl_key=params.control_key,
            stim_key=params.treatment_key,
            celltype_to_predict=params.cell_type_to_predict,
            output_name=params.output_path
        )

        return json.dumps(result, indent=2)

    except Exception as e:
        logger.error(f"Failed to predict response: {e}")
        return json.dumps({"status": "error", "message": str(e)}, indent=2)


@mcp.tool()
async def perturbation_differential_expression(params: DEInput) -> str:
    """Compare baseline vs. predicted expression.

    Identifies genes with significant expression changes between
    baseline and predicted cell states.

    Example:
        {"baseline_path": "./data/baseline.h5ad", "predicted_path": "./data/predicted.h5ad", "n_top_genes": 50}
    """
    try:
        # Load data
        baseline_adata = sc.read_h5ad(params.baseline_path)
        predicted_adata = sc.read_h5ad(params.predicted_path)

        # Run DE analysis
        analyzer = DifferentialExpressionAnalyzer()
        de_results = analyzer.compute_differential_expression(
            baseline_adata=baseline_adata,
            predicted_adata=predicted_adata,
            n_top_genes=params.n_top_genes,
            method=params.method
        )

        return json.dumps({
            "status": "success",
            "differential_expression": de_results
        }, indent=2)

    except Exception as e:
        logger.error(f"Failed to compute differential expression: {e}")
        return json.dumps({"status": "error", "message": str(e)}, indent=2)


@mcp.tool()
async def perturbation_get_latent(params: GetLatentInput) -> str:
    """Extract latent representations for visualization.

    Projects cells into the learned latent space for dimensionality
    reduction and visualization.

    Example:
        {"model_name": "my_model", "data_path": "./data/cells.h5ad"}
    """
    try:
        if params.model_name not in _models:
            return json.dumps({
                "status": "error",
                "message": f"Model '{params.model_name}' not found."
            }, indent=2)

        wrapper = _models[params.model_name]

        # Load data
        adata = sc.read_h5ad(params.data_path)

        # Get latent representation
        latent = wrapper.get_latent_representation(adata)

        # Save to obsm
        adata.obsm["X_scgen"] = latent

        # Save updated file
        output_path = params.data_path.replace(".h5ad", "_with_latent.h5ad")
        adata.write_h5ad(output_path)

        return json.dumps({
            "status": "success",
            "output_path": output_path,
            "latent_shape": list(latent.shape),
            "latent_key": "X_scgen"
        }, indent=2)

    except Exception as e:
        logger.error(f"Failed to extract latent: {e}")
        return json.dumps({"status": "error", "message": str(e)}, indent=2)


@mcp.tool()
async def perturbation_visualize(params: VisualizeInput) -> str:
    """Generate PCA/UMAP plots of baseline vs. predicted.

    Creates visualization comparing cell states before and after
    predicted perturbation.

    Example:
        {"baseline_path": "./data/baseline.h5ad", "predicted_path": "./data/predicted.h5ad", "plot_type": "pca"}
    """
    try:
        # Load data
        baseline_adata = sc.read_h5ad(params.baseline_path)
        predicted_adata = sc.read_h5ad(params.predicted_path)

        # Create visualizer
        visualizer = PerturbationVisualizer()

        # Generate plot
        output_path = visualizer.plot_baseline_vs_predicted(
            baseline_adata=baseline_adata,
            predicted_adata=predicted_adata,
            plot_type=params.plot_type,
            color_by=params.color_by,
            output_path=params.output_path
        )

        return json.dumps({
            "status": "success",
            "output_path": output_path,
            "plot_type": params.plot_type
        }, indent=2)

    except Exception as e:
        logger.error(f"Failed to create visualization: {e}")
        return json.dumps({"status": "error", "message": str(e)}, indent=2)


# ==================== Server Entry Point ====================

def main():
    """Run the MCP server."""
    import asyncio
    logger.info("Starting perturbation MCP server...")
    mcp.run()


if __name__ == "__main__":
    main()
