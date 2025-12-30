FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Create directory structure
RUN mkdir -p /app/servers/mcp-deepcell /app/shared/utils

# Copy shared utilities (from repo root)
COPY shared/utils/ /app/shared/utils/

# Copy server code (from repo root)
COPY servers/mcp-deepcell/ /app/servers/mcp-deepcell/

# Install Python dependencies
WORKDIR /app/servers/mcp-deepcell
RUN pip install --no-cache-dir -e .

# Set environment variables for SSE transport
ENV MCP_TRANSPORT=sse
ENV MCP_PORT=3007
ENV DEEPCELL_DRY_RUN=false
ENV DEEPCELL_DATA_DIR=/app/data/deepcell
ENV DEEPCELL_CACHE_DIR=/app/data/cache/deepcell

# Create data directories
RUN mkdir -p /app/data/deepcell /app/data/cache/deepcell

# Expose port
EXPOSE 3007

# Run the server
CMD ["python", "-m", "mcp_deepcell"]
