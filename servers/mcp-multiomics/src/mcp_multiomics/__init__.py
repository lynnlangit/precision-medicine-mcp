"""MCP server for multi-omics PDX data analysis."""

__version__ = "0.1.0"

from .config import config
from .server import mcp

__all__ = ["config", "mcp"]
