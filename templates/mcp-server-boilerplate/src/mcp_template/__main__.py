"""Entry point for mcp-{{SERVER_NAME}} server."""

from .server import mcp

if __name__ == "__main__":
    mcp.run()
