FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    curl \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Create directory structure
RUN mkdir -p /app/servers/mcp-quantum-celltype-fidelity /app/shared/utils

# Copy temporarily staged shared utilities (if needed)
COPY _shared_temp/utils/ /app/shared/utils/ 2>/dev/null || true

# Copy server code
COPY . /app/servers/mcp-quantum-celltype-fidelity/

# Install Python dependencies
WORKDIR /app/servers/mcp-quantum-celltype-fidelity

# Upgrade pip, setuptools, and wheel
RUN pip install --no-cache-dir --upgrade "pip>=24.0" "setuptools>=69.0" "wheel>=0.42.0"

# Install core dependencies first to avoid conflicts
# Note: qiskit-aer removed to avoid build issues; CPU simulation via qiskit works fine
RUN pip install --no-cache-dir \
    "fastmcp>=0.1.0" \
    "qiskit>=1.0.0" \
    "numpy>=1.24.0,<2.0" \
    "scipy>=1.11.0" \
    "pandas>=2.0.0" \
    "scikit-learn>=1.3.0" \
    "anndata>=0.10.0" \
    "fsspec>=2023.1.0" \
    "gcsfs>=2023.1.0"

# Install the package in editable mode (for any additional dependencies)
RUN pip install --no-cache-dir -e .

# Set environment variables for SSE transport
ENV MCP_TRANSPORT=sse
ENV QUANTUM_BACKEND=cpu
ENV QUANTUM_DATA_DIR=/app/data/quantum
ENV QUANTUM_CACHE_DIR=/app/data/cache/quantum

# Create data directories
RUN mkdir -p /app/data/quantum /app/data/cache/quantum

# Expose port (Cloud Run will set PORT environment variable)
EXPOSE 8080

# Run the server
CMD ["python", "-m", "quantum_celltype_fidelity.server"]
