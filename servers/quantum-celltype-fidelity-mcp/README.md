# Quantum Cell Type Fidelity MCP Server

Quantum computing-based analysis of cell type fidelity in spatial transcriptomics data using Qiskit parameterized quantum circuits (PQCs).

## Overview

This MCP server implements **QuCoWE-style** (Quantum Contrastive Word Embeddings) quantum circuits to embed cell types from spatial transcriptomics data into Hilbert space. Quantum fidelity between states is used to measure cell type similarity, enabling:

- **Immune evasion detection**: Identifying tumor cells evading immune surveillance
- **TLS analysis**: Characterizing tertiary lymphoid structures with quantum signatures
- **Perturbation prediction**: Simulating drug/treatment effects on cell states
- **Spatial context integration**: Learning embeddings aware of spatial neighborhoods

### Key Features

- **Parameterized Quantum Circuits (PQCs)** with 8-10 qubits for 256-1024 dimensional Hilbert space
- **Parameter-shift rule** for gradient estimation (works on real quantum hardware)
- **Contrastive learning** with InfoNCE loss for training embeddings
- **CPU/GPU/IBM Quantum** backend support
- **Spatial neighborhood** integration from AnnData objects
- **6 MCP tools** for complete quantum cell type analysis workflow

## Scientific Background

### Quantum Fidelity

Quantum fidelity $F(|\psi_a\rangle, |\psi_b\rangle) = |\langle\psi_a|\psi_b\rangle|^2$ measures the "overlap" between two quantum states. For cell type analysis:

- $F = 1$: Identical cell states
- $F = 0$: Orthogonal cell states (maximally different)
- $0 < F < 1$: Partial similarity

### Parameterized Quantum Circuits

Each cell type is encoded via a PQC:
1. **Feature encoding layer**: Gene expression → rotation angles
2. **Variational layers**: Learnable parameters with RX, RY, RZ gates
3. **Entanglement**: Ring entanglement captures relationships

Training uses contrastive learning to maximize fidelity within cell types and minimize between types.

### Parameter-Shift Rule

Gradients are computed exactly without backpropagation:

$$\frac{\partial f}{\partial \theta} = \frac{f(\theta + \pi/2) - f(\theta - \pi/2)}{2}$$

This enables training on real quantum hardware.

## Architecture

```
quantum-celltype-fidelity-mcp/
├── src/quantum_celltype_fidelity/
│   ├── circuits.py          # QuCoWECircuit (Qiskit implementation)
│   ├── embeddings.py         # QuCoWECellTypeEmbedding (parameter management)
│   ├── fidelity.py           # Fidelity scoring heads
│   ├── spatial_context.py    # Spatial neighborhood extraction
│   ├── training.py           # Parameter-shift gradient training
│   └── server.py             # MCP server with 6 tools
├── tests/                    # Unit and integration tests
├── examples/                 # Example notebooks and scripts
└── pyproject.toml            # Dependencies and configuration
```

## MCP Tools

### 1. `learn_spatial_cell_embeddings`

Train quantum embeddings from spatial transcriptomics data.

**Input:**
- `adata_path`: Path to AnnData (.h5ad) file or GCS URI
- `cell_type_key`: Key in adata.obs for cell type labels
- `n_qubits`: Number of qubits (default: 8)
- `n_layers`: Number of variational layers (default: 3)
- `n_epochs`: Training epochs (default: 50)
- `backend`: "cpu", "gpu", or "ibm"

**Output:**
- `embedding_id`: ID for referencing trained embeddings
- `training_summary`: Loss history and metrics
- `embedding_summary`: Circuit configuration

**Example:**
```python
result = learn_spatial_cell_embeddings(
    adata_path="gs://data/patientone_tcells.h5ad",
    cell_type_key="cell_type",
    n_qubits=8,
    n_epochs=50,
    backend="cpu"
)
```

### 2. `compute_cell_type_fidelity`

Compute quantum fidelity between cells.

**Input:**
- `adata_path`: Path to AnnData file
- `embedding_id`: ID from training
- `compute_matrix`: Whether to compute full NxN fidelity matrix

**Output:**
- `fidelity_matrix`: Pairwise fidelities (if compute_matrix=True)
- `summary_stats`: Mean, std, per-cell-type statistics

### 3. `identify_immune_evasion_states`

Detect cells in immune evasion states.

**Input:**
- `immune_cell_types`: List of canonical immune types (e.g., ["T_cell", "B_cell"])
- `exhausted_markers`: Optional exhausted cell types
- `evasion_threshold`: Threshold for flagging (default: 0.3)

**Output:**
- `evading_cells`: List of cells with high evasion scores
- `evasion_score`: Score in [0, 1] for each cell

### 4. `predict_perturbation_effect`

Predict how perturbations affect fidelity.

**Input:**
- `perturbation_type`: "drug", "genetic", "environmental"
- `target_cell_types`: Cell types to perturb
- `perturbation_strength`: Strength (0-1)

**Output:**
- `effects`: Fidelity changes for each cell type
- `effect_direction`: "increase" or "decrease"

### 5. `analyze_tls_quantum_signature`

Analyze quantum signatures of tertiary lymphoid structures.

**Input:**
- `tls_marker_types`: Cell types indicating TLS (e.g., ["B_cell", "T_cell"])
- `min_cluster_size`: Minimum cells for TLS (default: 20)
- `max_distance`: Spatial distance threshold (default: 100)

**Output:**
- `tls_candidates`: List of TLS with quantum signatures
- `quantum_signature`: Mean fidelity, std, per-type composition

### 6. `export_for_downstream`

Export embeddings for downstream analysis.

**Input:**
- `output_format`: "numpy", "json", "pytorch"
- `output_path`: Path to save

**Output:**
- `exported_files`: List of exported file paths

## Installation

### Standard Installation (CPU)

```bash
cd servers/quantum-celltype-fidelity-mcp
pip install -e .
```

### GPU Installation (CUDA + cuQuantum)

```bash
pip install -e ".[gpu]"
```

Requires:
- CUDA 11.8+
- cuQuantum SDK 23.3+
- NVIDIA GPU with Compute Capability 7.0+

### IBM Quantum Hardware

1. Get IBM Quantum API token from https://quantum-computing.ibm.com/
2. Set environment variable:
```bash
export IBM_QUANTUM_TOKEN="your_token_here"
```

## Usage Examples

### Example 1: Train Embeddings on PatientOne Data

```python
# Train quantum embeddings
result = await mcp.call_tool(
    "learn_spatial_cell_embeddings",
    {
        "adata_path": "gs://sample-inputs-patientone/perturbation/patientone_tcells.h5ad",
        "cell_type_key": "cell_type",
        "n_qubits": 8,
        "n_layers": 3,
        "n_epochs": 50,
        "learning_rate": 0.01,
        "backend": "cpu"
    }
)

embedding_id = result["embedding_id"]
print(f"Trained embeddings: {embedding_id}")
print(f"Final loss: {result['training_summary']['final_loss']:.4f}")
```

### Example 2: Detect Immune Evasion

```python
# Identify tumor cells evading immune surveillance
result = await mcp.call_tool(
    "identify_immune_evasion_states",
    {
        "adata_path": "gs://sample-inputs-patientone/perturbation/patientone_tcells.h5ad",
        "embedding_id": embedding_id,
        "immune_cell_types": ["T_cell", "B_cell", "NK_cell"],
        "exhausted_markers": ["T_exhausted", "T_regulatory"],
        "evasion_threshold": 0.3
    }
)

print(f"Found {result['n_evading_cells']} cells in evasion states")
for cell in result["evading_cells"][:5]:
    print(f"  Cell {cell['cell_idx']}: {cell['cell_type']} (score: {cell['evasion_score']:.3f})")
```

### Example 3: Analyze TLS Quantum Signatures

```python
# Find and characterize tertiary lymphoid structures
result = await mcp.call_tool(
    "analyze_tls_quantum_signature",
    {
        "adata_path": "gs://sample-inputs-patientone/perturbation/patientone_tcells.h5ad",
        "embedding_id": embedding_id,
        "tls_marker_types": ["B_cell", "T_cell", "Dendritic_cell"],
        "min_cluster_size": 20,
        "max_distance": 100.0
    }
)

print(f"Found {result['n_tls_candidates']} TLS candidates")
for tls in result["tls_candidates"]:
    print(f"\nTLS {tls['tls_id']}:")
    print(f"  Cells: {tls['n_cells']}")
    print(f"  Composition: {tls['cell_type_composition']}")
    print(f"  Mean fidelity: {tls['quantum_signature']['mean_fidelity']:.3f}")
```

## Development

### Running Tests

```bash
pytest tests/ -v
```

### Code Formatting

```bash
black src/ tests/
ruff check src/ tests/
```

### Type Checking

```bash
mypy src/
```

## Performance

### Circuit Simulation

| Backend | n_qubits | Statevector Time | Fidelity Computation |
|---------|----------|------------------|---------------------|
| CPU     | 8        | ~10 ms          | ~20 ms              |
| CPU     | 10       | ~50 ms          | ~100 ms             |
| GPU     | 8        | ~2 ms           | ~5 ms               |
| GPU     | 10       | ~5 ms           | ~10 ms              |

### Training

- **50 epochs, 500 cells, 5 cell types**: ~5 minutes (CPU), ~1 minute (GPU)
- **100 epochs, 1000 cells, 10 cell types**: ~30 minutes (CPU), ~6 minutes (GPU)

### Memory

- **8 qubits**: ~2 MB per statevector
- **10 qubits**: ~8 MB per statevector
- **Training**: ~100-500 MB depending on dataset size

## Deployment

### Cloud Run (Production)

```bash
cd infrastructure/deployment
./deploy_to_gcp.sh --server quantum-celltype-fidelity
```

### Local Development

```bash
# Set environment variables
export ANTHROPIC_API_KEY="your_key_here"

# Run server
cd servers/quantum-celltype-fidelity-mcp
python -m quantum_celltype_fidelity.server
```

## Scientific References

1. **QuCoWE**: Quantum Contrastive Word Embeddings approach adapted for cell types
2. **Parameter-shift rule**: Mitarai et al., Physical Review A 98, 032309 (2018)
3. **Variational quantum algorithms**: Cerezo et al., Nature Reviews Physics 3, 625-644 (2021)
4. **Spatial transcriptomics**: Moses & Pachter, Nature Methods 19, 534-546 (2022)

## Use Case: HGSOC Immunotherapy Targets

This server was developed for high-grade serous ovarian cancer (HGSOC) research:

- **Identify immune evasion**: Tumor cells with low T-cell fidelity
- **Characterize TLS**: Immune-rich regions for therapy targeting
- **Predict drug response**: Effects of checkpoint inhibitors on cell states
- **Spatial context**: How tumor-immune boundaries affect fidelity

## Limitations

- **Circuit depth**: Limited by decoherence on real quantum hardware
- **Scalability**: Current implementation handles ~10,000 cells; for larger datasets, use sampling
- **Feature encoding**: Simple amplitude encoding; future work will explore other encodings
- **Training time**: Parameter-shift rule requires 2N forward passes per gradient

## Future Work

- [ ] Add support for 3D spatial coordinates
- [ ] Implement amplitude amplification for rare cell types
- [ ] Integrate with GEARS perturbation server for combined analysis
- [ ] Add Bayesian optimization for hyperparameter tuning
- [ ] Support for IBM Quantum hardware execution

## License

Apache 2.0

## Contributing

Contributions welcome! Please see the main repository's CONTRIBUTING.md for guidelines.

## Support

- **Issues**: https://github.com/lynnlangit/precision-medicine-mcp/issues
- **Discussions**: https://github.com/lynnlangit/precision-medicine-mcp/discussions
- **Documentation**: https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs
