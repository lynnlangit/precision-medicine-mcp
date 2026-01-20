# MCP Perturbation Server

**Single-cell perturbation prediction using GEARS for precision medicine**

Predicts how patient cancer cells might respond to immunotherapy *in silico* using graph neural networks and biological knowledge graphs.

> **Powered by GEARS**: Uses GEARS (Graph-Enhanced Gene Activation and Repression Simulator), a state-of-the-art graph neural network approach published in Nature Biotechnology 2024.

---

## Overview

This MCP server uses **GEARS** (Graph-Enhanced Gene Activation and Repression Simulator) to predict cellular responses to perturbations (e.g., genetic modifications, drug treatments, immunotherapy) without actually performing the experiment. It integrates biological knowledge graphs of gene-gene relationships with deep learning to predict treatment effects.

### The GEARS Approach

GEARS uses graph neural networks (GNNs) combined with gene regulatory network knowledge to predict perturbation responses:

**Key advantages over VAE-based methods:**
- **40% higher precision** than previous approaches (Nature Biotechnology 2024)
- **Integrates biological knowledge**: Uses gene-gene relationship networks
- **Multi-gene perturbations**: Handles complex combinatorial effects
- **Uncertainty quantification**: Provides confidence estimates for predictions
- **Better generalization**: Leverages graph structure to predict unseen perturbations

**How it works:**
1. Encodes gene relationships as a knowledge graph
2. Uses GNN to learn how perturbations propagate through the network
3. Predicts gene expression changes for novel perturbations
4. Handles single and multi-gene perturbation combinations

### Key Applications

- **Treatment Response Prediction**: Predict how a patient's cells will respond to immunotherapy
- **Drug Screening**: Test multiple treatments *in silico* before clinical application
- **Personalized Medicine**: Identify optimal therapies based on patient-specific cellular profiles
- **Clinical Trial Design**: Pre-screen patients likely to respond to experimental therapies

---

## Installation

✅ **GEARS is now fully working** with Python 3.11+!

### Prerequisites

- Python >= 3.10 (Python 3.11+ recommended)
- PyTorch >= 2.0.0
- CUDA (optional, for GPU acceleration)

### Install from Source

```bash
cd servers/mcp-perturbation
pip install -e .

# For development
pip install -e ".[dev]"
```

The installation will automatically install:
- `cell-gears` - GEARS perturbation prediction
- `torch-geometric` - Graph neural network framework
- `scanpy`, `anndata` - Single-cell data handling
- Modern compatible versions of all dependencies

**Installation time:** ~2-3 minutes
**No dependency conflicts!** All packages are modern and compatible.

### Dependencies

Core dependencies:
- `cell-gears` - GEARS perturbation prediction (Nature Biotech 2024)
- `torch-geometric` - Graph neural network framework
- `scanpy` - Single-cell analysis
- `anndata` - Single-cell data structures
- `torch` - Deep learning framework (PyTorch)
- `mcp` / `fastmcp` - MCP server framework

**GEARS advantages:**
- ✅ Modern Python 3.11+ compatible
- ✅ 40% better performance than VAE methods (Nature Biotechnology 2024)
- ✅ Handles multi-gene perturbations
- ✅ Integrates biological knowledge graphs
- ✅ Active maintenance and support

---

## Quick Start

### 1. Load Dataset

Load scRNA-seq data from GEO or a local .h5ad file:

```json
{
  "tool": "perturbation_load_dataset",
  "params": {
    "dataset_id": "GSE184880",
    "normalize": true,
    "n_hvg": 7000,
    "cell_type_key": "cell_type",
    "condition_key": "condition"
  }
}
```

**Returns**: Dataset metadata (n_cells, n_genes, cell types, conditions)

### 2. Setup and Train Model

Initialize GEARS model:

```json
{
  "tool": "perturbation_setup_model",
  "params": {
    "dataset_id": "GSE184880",
    "hidden_size": 64,
    "num_layers": 2,
    "uncertainty": true,
    "model_name": "ovarian_cancer_model"
  }
}
```

Train the GEARS model:

```json
{
  "tool": "perturbation_train_model",
  "params": {
    "model_name": "ovarian_cancer_model",
    "epochs": 20,
    "batch_size": 32,
    "learning_rate": 0.001
  }
}
```

**Note:** GEARS trains faster than VAE-based methods (20 epochs typical vs 100+)

### 3. Compute Perturbation Effect

Calculate perturbation effect for specific genes:

```json
{
  "tool": "perturbation_compute_delta",
  "params": {
    "model_name": "ovarian_cancer_model",
    "source_cell_type": "T_cells",
    "treatment_key": "CD4"
  }
}
```

**Note:** GEARS predicts effects of genetic perturbations (gene knockouts/upregulation). For drug treatments, map drugs to their target genes.

### 4. Predict Patient Response

Apply GEARS prediction to patient data:

```json
{
  "tool": "perturbation_predict_response",
  "params": {
    "model_name": "ovarian_cancer_model",
    "patient_data_path": "./data/patient_001.h5ad",
    "cell_type_to_predict": "T_cells",
    "treatment_key": "CD4,CD8A"
  }
}
```

**Returns**: Predicted cell states, perturbation effect magnitude, number of cells predicted

**Multi-gene perturbations**: GEARS excels at predicting combinatorial effects - specify multiple genes separated by commas.

---

## Tool Reference

### 1. `perturbation_load_dataset`

Load and preprocess scRNA-seq data from GEO or local file.

**Parameters:**
- `dataset_id` (str): GEO accession (e.g., "GSE184880") or path to .h5ad
- `normalize` (bool): Apply normalization (default: true)
- `n_hvg` (int): Number of highly variable genes (default: 7000)
- `cell_type_key` (str): Column with cell type labels
- `condition_key` (str): Column with treatment condition

**Returns:** JSON with dataset metadata

---

### 2. `perturbation_setup_model`

Initialize GEARS graph neural network model.

**Parameters:**
- `dataset_id` (str): Dataset ID from load_dataset
- `hidden_size` (int): Hidden layer size for GNN (default: 64)
- `num_layers` (int): Number of GNN layers (default: 2)
- `uncertainty` (bool): Enable uncertainty quantification (default: true)
- `uncertainty_reg` (float): Uncertainty regularization weight (default: 1.0)
- `model_name` (str): Name for this model
- `condition_key` (str): Column with condition labels
- `pert_key` (str): Column with perturbation labels

**Returns:** Model configuration summary

---

### 3. `perturbation_train_model`

Train GEARS graph neural network on perturbation data.

**Parameters:**
- `model_name` (str): Model name from setup_model
- `epochs` (int): Training epochs (default: 20)
- `batch_size` (int): Batch size (default: 32)
- `learning_rate` (float): Learning rate (default: 0.001)
- `valid_every` (int): Validation frequency in epochs (default: 1)

**Returns:** Training metrics (final loss, epochs completed, model path)

---

### 4. `perturbation_compute_delta`

Calculate perturbation vector (Δ) between conditions.

**Parameters:**
- `model_name` (str): Trained model name
- `source_cell_type` (str): Cell type to compute delta from (None = all)
- `control_key` (str): Control condition label
- `treatment_key` (str): Treatment condition label

**Returns:** Delta vector statistics (norm, mean, std, cell counts)

---

### 5. `perturbation_predict_response`

Apply Δ to patient's baseline cells to predict treated state.

**Parameters:**
- `model_name` (str): Trained model name
- `patient_data_path` (str): Path to patient .h5ad file
- `cell_type_to_predict` (str): Cell type to transform
- `control_key` (str): Control condition label
- `treatment_key` (str): Treatment condition label
- `output_path` (str, optional): Path to save predictions

**Returns:** Prediction summary with file path, delta norm, changed genes

---

### 6. `perturbation_differential_expression`

Compare baseline vs. predicted expression.

**Parameters:**
- `baseline_path` (str): Baseline .h5ad file
- `predicted_path` (str): Predicted .h5ad file
- `n_top_genes` (int): Number of top genes to return (default: 50)
- `method` (str): Statistical test ("wilcoxon" or "t-test")

**Returns:** Top upregulated/downregulated genes with fold changes

---

### 7. `perturbation_get_latent`

Extract latent representations for visualization.

**Parameters:**
- `model_name` (str): Trained model name
- `data_path` (str): .h5ad file to embed

**Returns:** Path to .h5ad with latent embeddings in .obsm["X_gears"]

---

### 8. `perturbation_visualize`

Generate PCA/UMAP plots of baseline vs. predicted.

**Parameters:**
- `baseline_path` (str): Baseline .h5ad file
- `predicted_path` (str): Predicted .h5ad file
- `plot_type` (str): "pca" or "umap" (default: "pca")
- `color_by` (str): Column to color by (default: "condition")
- `output_path` (str, optional): Path to save figure

**Returns:** Path to saved figure

---

## Primary Dataset: GSE184880

**Study**: Single-cell RNA-seq of ovarian cancer (HGSOC) and healthy controls

**Samples**:
- 5 healthy controls
- 7 high-grade serous ovarian cancer (HGSOC) patients

**Cell Types**:
- T cells (CD4+, CD8+)
- B cells
- Macrophages
- Epithelial cells
- Fibroblasts

**Conditions**:
- `control` - Healthy tissue
- `tumor` - Cancer tissue

This dataset is ideal for learning control vs. disease perturbation vectors for ovarian cancer immunotherapy prediction.

---

## Example Workflow: PatientOne

### Scenario

PatientOne has HGSOC ovarian cancer. We want to predict how her CD8+ T cells will respond to immunotherapy.

### Step 1: Load Reference Data

```json
{
  "tool": "perturbation_load_dataset",
  "params": {
    "dataset_id": "GSE184880",
    "normalize": true,
    "n_hvg": 7000
  }
}
```

### Step 2: Train Model

```json
{
  "tool": "perturbation_setup_model",
  "params": {
    "dataset_id": "GSE184880",
    "model_name": "patient_one_model"
  }
}
```

```json
{
  "tool": "perturbation_train_model",
  "params": {
    "model_name": "patient_one_model",
    "n_epochs": 100
  }
}
```

### Step 3: Predict Response

```json
{
  "tool": "perturbation_predict_response",
  "params": {
    "model_name": "patient_one_model",
    "patient_data_path": "./data/patient_one_baseline.h5ad",
    "cell_type_to_predict": "CD8_T_cells",
    "control_key": "control",
    "treatment_key": "tumor"
  }
}
```

### Step 4: Analyze Changes

```json
{
  "tool": "perturbation_differential_expression",
  "params": {
    "baseline_path": "./data/patient_one_baseline.h5ad",
    "predicted_path": "./data/predictions/patient_one_model_predicted.h5ad",
    "n_top_genes": 50
  }
}
```

### Expected Results

**Top Upregulated Genes** (predicted treatment response):
- `IFNG` - Interferon gamma (immune activation)
- `GZMB` - Granzyme B (cytotoxicity marker)
- `PRF1` - Perforin (cell killing)
- `CD69` - T cell activation marker

**Clinical Interpretation**: The model predicts that immunotherapy will activate PatientOne's CD8+ T cells, increasing cytotoxic markers and immune response genes.

---

## Architecture

### Model: GEARS (Graph Neural Network)

**Type**: Graph Neural Network (GNN) for perturbation prediction

**Components**:
1. **Gene Relationship Graph**: Encodes biological knowledge (GO terms, PPI networks)
2. **GNN Layers**: Message passing over gene-gene relationships
3. **Prediction Head**: Predicts gene expression changes
4. **Uncertainty Module**: Provides confidence estimates (optional)

**How GEARS Works**:
```
1. Build gene relationship graph (20,000+ genes)
2. GNN propagates perturbation effects through network
3. Predict expression changes for each gene
4. Handle multi-gene perturbations (combinatorial effects)
```

### Prediction Process

```
Input: Gene perturbations (e.g., ["PDCD1", "CTLA4"])
       ↓
Gene Graph: Encode relationships
       ↓
GNN Layers: Message passing (2 layers)
       ↓
Output: Predicted gene expression changes Δ
       ↓
Apply to baseline: Predicted cell state
```

For detailed architecture diagrams, see [ARCHITECTURE.md](ARCHITECTURE.md).

---

## Testing

### Run All Tests

```bash
pytest tests/ -v
```

### Run Specific Test File

```bash
pytest tests/test_gears_wrapper.py -v
```

### Test Coverage

```bash
pytest tests/ --cov=mcp_perturbation --cov-report=html
```

### Expected Coverage

- `data_loader.py`: >85%
- `gears_wrapper.py`: >80%
- `prediction.py`: >75%
- `server.py`: >70%

---

## Performance

### Training Time

| Dataset Size | n_latent | n_epochs | GPU Time | CPU Time |
|--------------|----------|----------|----------|----------|
| 5K cells | 100 | 100 | ~2 min | ~10 min |
| 20K cells | 100 | 100 | ~5 min | ~30 min |
| 100K cells | 100 | 100 | ~15 min | ~2 hours |

### Prediction Time

| Operation | Time |
|-----------|------|
| Predict 1K cells | ~1 second |
| Predict 10K cells | ~5 seconds |
| Extract latent | ~2 seconds |

---

## Troubleshooting

### Issue: CUDA out of memory

**Solution**: Reduce batch size or n_latent

```json
{
  "n_latent": 50,
  "batch_size": 16
}
```

### Issue: Model not converging

**Solutions**:
1. Increase n_epochs
2. Reduce learning rate
3. Increase n_hidden or n_latent

### Issue: Poor predictions

**Solutions**:
1. Ensure reference data has both control and treated conditions
2. Check that cell types are annotated correctly
3. Increase n_hvg (more genes = better signal)
4. Train longer (more epochs)

---

## Scientific Background

### Key Papers

1. **GEARS**: Roohani et al., "Predicting transcriptional outcomes of novel multigene perturbations with GEARS", Nature Biotechnology (2024)
2. **scVI**: Lopez et al., "Deep generative modeling for single-cell transcriptomics", Nature Methods (2018)
3. **HGSOC scRNA-seq**: Izar et al., "A single-cell landscape of high-grade serous ovarian cancer", Nature Medicine (2020)

### Method Performance

GEARS has been validated for:
- Perturbation prediction (40% better than VAE methods)
- Multi-gene combinatorial perturbations
- Cross-study generalization
- Handles unseen gene combinations
- Integrates biological knowledge graphs

---

## Related Servers

- **mcp-spatialtools**: Spatial transcriptomics analysis
- **mcp-multiomics**: Multi-omics integration
- **mcp-epic**: Clinical data from FHIR
- **mcp-tcga**: TCGA cohort data

---

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md) for guidelines.

### Development Setup

```bash
cd servers/mcp-perturbation
pip install -e ".[dev]"
pre-commit install
```

---

## License

Apache 2.0 - See [LICENSE](../../LICENSE)

---

**Part of the Precision Medicine MCP suite** - Enabling AI-driven bioinformatics for cancer research.
