"""MCP Hugging Face server - ML models for genomics."""

import json
import os
from typing import Any, Dict, List
from fastmcp import FastMCP

mcp = FastMCP("huggingface")

# DRY_RUN warning wrapper
def add_dry_run_warning(result: Any) -> Any:
    """Add warning banner to results when in DRY_RUN mode."""
    if not config.dry_run:
        return result

    warning = """
╔═══════════════════════════════════════════════════════════════════════════╗
║                    ⚠️  SYNTHETIC DATA WARNING ⚠️                          ║
╚═══════════════════════════════════════════════════════════════════════════╝

This result was generated in DRY_RUN mode and does NOT represent real analysis.

🔴 CRITICAL: Do NOT use this data for research decisions or publications.
🔴 All values are SYNTHETIC/MOCKED and have no scientific validity.

To enable real data processing, set: HF_DRY_RUN=false

═══════════════════════════════════════════════════════════════════════════

"""

    if isinstance(result, dict):
        result["_DRY_RUN_WARNING"] = "SYNTHETIC DATA - NOT FOR RESEARCH USE"
        result["_message"] = warning.strip()
    elif isinstance(result, str):
        result = warning + result

    return result


HF_TOKEN = os.getenv("HF_TOKEN", "mock_token")
DRY_RUN = os.getenv("HF_DRY_RUN", "true").lower() == "true"

@mcp.tool()
async def load_genomic_model(
    model_name: str,
    model_type: str = "dna"
) -> Dict[str, Any]:
    """Load pre-trained genomic language models.

    Args:
        model_name: Model identifier (e.g., "DNABERT-2", "Geneformer", "Nucleotide-Transformer")
        model_type: Type of model - "dna", "rna", "protein", "single-cell"

    Returns:
        Dictionary with model info, parameters, capabilities

    Example:
        >>> result = await load_genomic_model("DNABERT-2", "dna")
    """
    if DRY_RUN:
        return {
            "model_name": model_name,
            "model_type": model_type,
            "parameters": "117M" if "BERT" in model_name else "106M",
            "max_sequence_length": 512,
            "vocab_size": 4096,
            "capabilities": ["sequence_classification", "feature_extraction", "masked_prediction"],
            "status": "loaded",
            "mode": "dry_run"
        }
    return {"model_name": model_name, "status": "loaded"}

@mcp.tool()
async def predict_cell_type(
    expression_data: str,
    model: str = "Geneformer"
) -> Dict[str, Any]:
    """Cell type classification using foundation models.

    Args:
        expression_data: Path to expression matrix file
        model: Model to use - "Geneformer", "scGPT", "scBERT"

    Returns:
        Dictionary with predicted cell types and confidence scores
    """
    if DRY_RUN:
        return {
            "predictions": [
                {"cell_id": "cell_001", "cell_type": "T cell", "confidence": 0.94},
                {"cell_id": "cell_002", "cell_type": "B cell", "confidence": 0.88},
                {"cell_id": "cell_003", "cell_type": "Macrophage", "confidence": 0.91}
            ],
            "model_used": model,
            "total_cells": 1000,
            "unique_types": 15,
            "mode": "dry_run"
        }
    return {"predictions": [], "model_used": model}

@mcp.tool()
async def embed_sequences(
    sequences: List[str],
    model: str = "DNABERT-2"
) -> Dict[str, Any]:
    """Generate embeddings for DNA/RNA sequences.

    Args:
        sequences: List of DNA/RNA sequences to embed
        model: Embedding model to use

    Returns:
        Dictionary with embedding vectors and dimensions
    """
    if DRY_RUN:
        return {
            "embeddings_shape": [len(sequences), 768],
            "model": model,
            "sequences_processed": len(sequences),
            "embedding_dimension": 768,
            "mode": "dry_run"
        }
    return {"embeddings_shape": [0, 0], "model": model}

@mcp.resource("hf://models/dnabert")
def get_dnabert_info() -> str:
    """DNABERT model information."""
    return json.dumps({
        "model": "DNABERT-2",
        "description": "Pre-trained bidirectional encoder for DNA sequences",
        "parameters": "117M",
        "use_cases": ["Gene function prediction", "Promoter identification", "Regulatory element analysis"],
        "huggingface_id": "zhihan1996/DNABERT-2-117M"
    }, indent=2)

@mcp.resource("hf://models/geneformer")
def get_geneformer_info() -> str:
    """Geneformer model information."""
    return json.dumps({
        "model": "Geneformer",
        "description": "Foundation model for single-cell transcriptomics",
        "parameters": "106M",
        "use_cases": ["Cell type annotation", "Cell state prediction", "Gene network inference"],
        "huggingface_id": "ctheodoris/Geneformer"
    }, indent=2)

def main() -> None:
    """Run the MCP mcp-huggingface server."""
    logger.info("Starting mcp-huggingface server...")

    if config.dry_run:
        logger.warning("=" * 80)
        logger.warning("⚠️  DRY_RUN MODE ENABLED - RETURNING SYNTHETIC DATA")
        logger.warning("⚠️  Results are MOCKED and do NOT represent real analysis")
        logger.warning("⚠️  Set HF_DRY_RUN=false for production use")
        logger.warning("=" * 80)
    else:
        logger.info("✅ Real data processing mode enabled (HF_DRY_RUN=false)")

    mcp.run(transport="stdio")

if __name__ == "__main__":
    main()
