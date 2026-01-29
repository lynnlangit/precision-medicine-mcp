"""MCP Server Configuration

Contains URLs and metadata for all 12 deployed GCP Cloud Run MCP servers.
"""

from typing import Dict, List, TypedDict


class MCPServerConfig(TypedDict):
    """MCP server configuration."""
    name: str
    url: str
    description: str
    status: str
    tools_count: int


# MCP Server URLs - Set to local for testing, Cloud Run for production
# LOCAL: http://localhost:PORT/sse (requires ngrok/tunnel for Claude API to reach)
# CLOUD RUN: https://mcp-SERVER-ondu7mwjpa-uc.a.run.app/sse
USE_LOCAL_SERVERS = False  # Must be False when using Claude API beta MCP client

MCP_SERVERS: Dict[str, MCPServerConfig] = {
    "fgbio": {
        "name": "fgbio",
        "url": "http://localhost:8001/sse" if USE_LOCAL_SERVERS else "https://mcp-fgbio-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Genomic reference data and FASTQ validation",
        "status": "production",
        "tools_count": 4
    },
    "multiomics": {
        "name": "multiomics",
        "url": "http://localhost:8002/sse" if USE_LOCAL_SERVERS else "https://mcp-multiomics-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Multi-omics integration (RNA/Protein/Phospho)",
        "status": "production",
        "tools_count": 9
    },
    "spatialtools": {
        "name": "spatialtools",
        "url": "http://localhost:8000/sse" if USE_LOCAL_SERVERS else "https://mcp-spatialtools-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Spatial transcriptomics analysis",
        "status": "production",
        "tools_count": 10
    },
    "tcga": {
        "name": "tcga",
        "url": "https://mcp-tcga-ondu7mwjpa-uc.a.run.app/sse",
        "description": "TCGA cancer genomics data",
        "status": "mock",
        "tools_count": 5
    },
    "openimagedata": {
        "name": "openimagedata",
        "url": "https://mcp-openimagedata-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Medical imaging datasets",
        "status": "mock",
        "tools_count": 4
    },
    "seqera": {
        "name": "seqera",
        "url": "https://mcp-seqera-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Nextflow workflow management",
        "status": "mock",
        "tools_count": 5
    },
    "huggingface": {
        "name": "huggingface",
        "url": "https://mcp-huggingface-ondu7mwjpa-uc.a.run.app/sse",
        "description": "AI/ML models for genomics",
        "status": "mock",
        "tools_count": 4
    },
    "deepcell": {
        "name": "deepcell",
        "url": "https://mcp-deepcell-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Cell segmentation",
        "status": "mock",
        "tools_count": 3
    },
    "mockepic": {
        "name": "mockepic",
        "url": "https://mcp-mockepic-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Mock EHR/FHIR data",
        "status": "mock",
        "tools_count": 5
    },
    "perturbation": {
        "name": "perturbation",
        "url": "http://localhost:8080/sse" if USE_LOCAL_SERVERS else "https://mcp-perturbation-ondu7mwjpa-uc.a.run.app/sse",
        "description": "GEARS perturbation prediction for treatment response",
        "status": "production",
        "tools_count": 8
    },
    "quantum-celltype-fidelity": {
        "name": "quantum-celltype-fidelity",
        "url": "http://localhost:3010/sse" if USE_LOCAL_SERVERS else "https://mcp-quantum-celltype-fidelity-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Quantum computing for cell type validation and immune evasion detection",
        "status": "production",
        "tools_count": 6
    }
}


def get_server_config(server_names: List[str]) -> List[Dict]:
    """Get MCP server configs for Claude API.

    Args:
        server_names: List of server names to include

    Returns:
        List of MCP server configurations for Claude API
    """
    configs = []
    for name in server_names:
        if name in MCP_SERVERS:
            configs.append({
                "type": "url",
                "url": MCP_SERVERS[name]["url"],
                "name": name
            })
    return configs


def get_tools_config(server_names: List[str]) -> List[Dict]:
    """Get tools configuration for Claude API.

    Args:
        server_names: List of server names to include

    Returns:
        List of tool configurations for Claude API
    """
    return [
        {"type": "mcp_toolset", "mcp_server_name": name}
        for name in server_names
        if name in MCP_SERVERS
    ]


def get_server_categories() -> Dict[str, List[str]]:
    """Get servers organized by category.

    Returns:
        Dict of category -> list of server names
    """
    return {
        "Production Servers (Real Analysis)": [
            "fgbio",
            "multiomics",
            "spatialtools",
            "perturbation",
            "quantum-celltype-fidelity"
        ],
        "Mock Servers (Workflow Demo)": [
            "tcga",
            "openimagedata",
            "seqera",
            "huggingface",
            "deepcell",
            "mockepic"
        ]
    }


# Example prompts for different use cases
EXAMPLE_PROMPTS = {
    "Spatial Analysis": "Analyze the spatial transcriptomics data for Patient-001. Perform cell type deconvolution and identify key cell populations.",

    "Multi-omics Integration": "Integrate RNA, protein, and phosphorylation data. Run HAllA association analysis and identify significant correlations.",

    "Genomic QC": "Validate the FASTQ files and check quality metrics. What is the average quality score and read length?",

    "Pathway Enrichment": "For the upregulated genes [TP53, BRCA1, MYC, KRAS], perform pathway enrichment analysis using GO_BP database.",

    "Complete PatientOne Workflow": """For Patient-001 (ovarian cancer):
1. Get clinical data from FHIR
2. Retrieve spatial transcriptomics data
3. Perform cell type deconvolution
4. Run differential expression between tumor core and margin
5. Generate treatment recommendations""",

    "Batch Correction": "I have 3 batches of proteomics data with batch effects. Apply ComBat batch correction and verify PC1 no longer correlates with batch.",

    "Predict Treatment Response": "Load the GSE184880 ovarian cancer dataset, setup a GEARS model, train it, and predict how Patient-001's T cells will respond to checkpoint inhibitor therapy.",

    "Immunotherapy Prediction": "Predict the response of T cells to anti-PD1 and anti-CTLA4 combination therapy. Which genes are most upregulated?",

    "Drug Screening": "For Patient-001 with ovarian cancer, test responses to: 1) Checkpoint inhibitors (PD1/CTLA4), 2) PARP inhibitors, 3) Platinum therapy. Which shows the best predicted response?",

    "Quantum Cell Type Fidelity": "Train quantum embeddings on the T-cell spatial transcriptomics data. Compute fidelity scores and identify cells with immune evasion states near the tumor boundary.",

    "Immune Evasion Detection": "Using quantum fidelity analysis, identify T-cells that are evading immune surveillance. What is the evasion score for cells near the tumor margin?",

    "TLS Analysis": "Analyze the quantum signatures of tertiary lymphoid structures in the spatial data. Which TLS candidates show the highest quantum coherence?",

    "Quantum + GEARS Validation": "First predict T-cell response to checkpoint inhibitors using GEARS. Then encode the predicted gene expression changes into quantum states and compute fidelity changes. Do the quantum and GEARS predictions agree?",
}
