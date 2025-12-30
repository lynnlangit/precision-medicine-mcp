FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Create directory structure
RUN mkdir -p /app/servers/mcp-seqera /app/shared/utils

# Copy temporarily staged shared utilities
COPY _shared_temp/utils/ /app/shared/utils/

# Copy server code
COPY . /app/servers/mcp-seqera/

# Install Python dependencies
WORKDIR /app/servers/mcp-seqera
RUN pip install --no-cache-dir -e .

# Set environment variables for SSE transport
ENV MCP_TRANSPORT=sse
ENV MCP_PORT=3005
ENV PYTHONPATH=/app/shared/utils:${PYTHONPATH}
ENV SEQERA_DRY_RUN=false
ENV SEQERA_DATA_DIR=/app/data/seqera
ENV SEQERA_CACHE_DIR=/app/data/cache/seqera

# Create data directories
RUN mkdir -p /app/data/seqera /app/data/cache/seqera

# Expose port
EXPOSE 3005

# Run the server
CMD ["python", "-m", "mcp_seqera"]
