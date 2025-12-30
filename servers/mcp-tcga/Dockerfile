FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Create directory structure
RUN mkdir -p /app/servers/mcp-tcga /app/shared/utils

# Copy shared utilities (from repo root)
COPY shared/utils/ /app/shared/utils/

# Copy server code (from repo root)
COPY servers/mcp-tcga/ /app/servers/mcp-tcga/

# Install Python dependencies
WORKDIR /app/servers/mcp-tcga
RUN pip install --no-cache-dir -e .

# Set environment variables for SSE transport
ENV MCP_TRANSPORT=sse
ENV MCP_PORT=3003
ENV TCGA_DRY_RUN=false
ENV TCGA_DATA_DIR=/app/data/tcga
ENV TCGA_CACHE_DIR=/app/data/cache/tcga

# Create data directories
RUN mkdir -p /app/data/tcga /app/data/cache/tcga

# Expose port
EXPOSE 3003

# Run the server
CMD ["python", "-m", "mcp_tcga"]
