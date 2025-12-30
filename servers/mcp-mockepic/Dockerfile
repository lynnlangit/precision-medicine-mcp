FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Create directory structure
RUN mkdir -p /app/servers/mcp-mockepic /app/shared/utils

# Copy shared utilities (from repo root)
COPY shared/utils/ /app/shared/utils/

# Copy server code (from repo root)
COPY servers/mcp-mockepic/ /app/servers/mcp-mockepic/

# Install Python dependencies
WORKDIR /app/servers/mcp-mockepic
RUN pip install --no-cache-dir -e .

# Set environment variables for SSE transport
ENV MCP_TRANSPORT=sse
ENV MCP_PORT=3008
ENV EPIC_DRY_RUN=false
ENV EPIC_DATA_DIR=/app/data/epic
ENV EPIC_CACHE_DIR=/app/data/cache/epic

# Create data directories
RUN mkdir -p /app/data/epic /app/data/cache/epic

# Expose port
EXPOSE 3008

# Run the server
CMD ["python", "-m", "mcp_mockepic"]
