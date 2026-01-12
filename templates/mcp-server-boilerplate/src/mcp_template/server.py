"""
MCP Server for {{MODALITY_DESCRIPTION}}

Provides tools for [describe your modality] including:
- [Tool category 1]
- [Tool category 2]
- [Tool category 3]
"""

import os
from pathlib import Path
from typing import Optional
import pandas as pd
import numpy as np
from fastmcp import FastMCP

# Initialize server
mcp = FastMCP("{{SERVER_NAME}}")

# Configuration from environment
DRY_RUN = os.getenv("{{SERVER_NAME_UPPER}}_DRY_RUN", "true").lower() == "true"
DATA_DIR = os.getenv("{{SERVER_NAME_UPPER}}_DATA_DIR", "/data/{{SERVER_NAME}}")
CACHE_DIR = os.getenv("{{SERVER_NAME_UPPER}}_CACHE_DIR", "/data/cache/{{SERVER_NAME}}")


@mcp.tool()
async def {{TOOL_EXAMPLE}}(
    data_file: str,
    param1: str,
    param2: Optional[int] = None
) -> dict:
    """
    [Brief one-line description of what this tool does]

    [Detailed description explaining the purpose, methodology, and use cases]

    Args:
        data_file: Path to input data file (CSV, TSV, etc.)
        param1: Description of param1
        param2: Optional description of param2 (default: None)

    Returns:
        Dictionary with:
        - status: "success" or "DRY_RUN"
        - result_key1: Description of result
        - result_key2: Description of result

    Raises:
        IOError: If input file not found
        ValueError: If invalid parameters

    Example:
        >>> result = await {{TOOL_EXAMPLE}}(
        ...     data_file="/data/sample.csv",
        ...     param1="value1",
        ...     param2=10
        ... )
        >>> print(result["status"])
        success
    """
    if DRY_RUN:
        # Return synthetic/mocked data for demonstration
        return {
            "status": "DRY_RUN",
            "message": "Simulated [tool name] execution",
            "result_key1": "synthetic_value",
            "result_key2": 42,
        }

    # Real implementation
    try:
        # 1. Load and validate data
        data = pd.read_csv(data_file)

        # 2. Perform analysis/computation
        result = _perform_analysis(data, param1, param2)

        # 3. Return structured results
        return {
            "status": "success",
            "data_file": str(data_file),
            "result_key1": result["key1"],
            "result_key2": result["key2"],
        }

    except FileNotFoundError:
        return {
            "status": "error",
            "error": f"Data file not found: {data_file}",
            "suggestion": "Check file path and ensure data directory is mounted"
        }
    except Exception as e:
        return {
            "status": "error",
            "error": str(e),
            "suggestion": "Check input file format and parameters"
        }


def _perform_analysis(data: pd.DataFrame, param1: str, param2: Optional[int]) -> dict:
    """
    Helper function for core analysis logic.

    Separating implementation from FastMCP tool makes testing easier.

    Args:
        data: Input dataframe
        param1: Analysis parameter
        param2: Optional parameter

    Returns:
        Dictionary with analysis results
    """
    # Implement your analysis logic here
    result = {
        "key1": "computed_value",
        "key2": len(data),
    }
    return result


@mcp.tool()
async def list_available_data(
    data_dir: Optional[str] = None
) -> dict:
    """
    List available data files for this modality.

    Utility tool to help users discover available datasets.

    Args:
        data_dir: Optional custom data directory (default: from env)

    Returns:
        Dictionary with:
        - data_directory: Path to data directory
        - files: List of available data files
        - file_count: Number of files found

    Example:
        >>> result = await list_available_data()
        >>> print(result["file_count"])
        12
    """
    if DRY_RUN:
        return {
            "status": "DRY_RUN",
            "data_directory": DATA_DIR,
            "files": ["sample1.csv", "sample2.csv", "sample3.csv"],
            "file_count": 3,
        }

    # Real implementation
    target_dir = data_dir or DATA_DIR
    data_path = Path(target_dir)

    if not data_path.exists():
        return {
            "status": "error",
            "error": f"Data directory not found: {target_dir}",
            "suggestion": f"Set {{SERVER_NAME_UPPER}}_DATA_DIR environment variable"
        }

    # Find relevant data files
    files = [
        str(f.relative_to(data_path))
        for f in data_path.glob("**/*")
        if f.is_file() and f.suffix in [".csv", ".tsv", ".txt"]
    ]

    return {
        "status": "success",
        "data_directory": str(data_path),
        "files": sorted(files)[:50],  # Limit to first 50
        "file_count": len(files),
    }


# Add more tools here following the same pattern...
# Each tool should:
# 1. Have clear docstring with Args, Returns, Example
# 2. Handle DRY_RUN mode with synthetic data
# 3. Implement real logic with error handling
# 4. Return structured dictionary results


if __name__ == "__main__":
    mcp.run()
