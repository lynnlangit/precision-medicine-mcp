FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Create directory structure
RUN mkdir -p /app/servers/mcp-spatialtools /app/shared/utils

# Copy shared utilities (from repo root)
COPY shared/utils/ /app/shared/utils/

# Copy server code (from repo root)
COPY servers/mcp-spatialtools/ /app/servers/mcp-spatialtools/

# Install Python dependencies
WORKDIR /app/servers/mcp-spatialtools
RUN pip install --no-cache-dir -e .

# Set environment variables for SSE transport
ENV MCP_TRANSPORT=sse
ENV MCP_PORT=3002
ENV SPATIAL_DRY_RUN=false
ENV SPATIAL_DATA_DIR=/app/data/spatial
ENV SPATIAL_CACHE_DIR=/app/data/cache/spatial

# Create data directories
RUN mkdir -p /app/data/spatial /app/data/cache/spatial

# Expose port
EXPOSE 3002

# Run the server
CMD ["python", "-m", "mcp_spatialtools"]
