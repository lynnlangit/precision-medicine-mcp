# =============================================================================
# Precision Medicine - TensorFlow Base Image
# =============================================================================
# Base container for TensorFlow-based machine learning tools
# =============================================================================
# OPTIMIZED FOR SIZE - Reduced from ~18GB to ~8-10GB
#
# Switch between GPU and CPU versions by commenting/uncommenting the FROM line:
# - GPU version: ~8-10GB (includes CUDA support)
# - CPU version: ~4-5GB (CPU only, much smaller)
# =============================================================================

# GPU version (CUDA support) - ACTIVE
FROM tensorflow/tensorflow:2.15.0

# CPU version (smaller, no GPU) - Uncomment to use instead of GPU version
# FROM tensorflow/tensorflow:2.15.0-cpu

LABEL maintainer="Precision Medicine Platform"
LABEL description="TensorFlow base image for ML/DL tools (Optimized)"
LABEL version="2.0.0"

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# =============================================================================
# System Dependencies & Python Packages
# =============================================================================
# Combine steps to reduce layers and enable cleanup
# =============================================================================
WORKDIR /app

# Install system dependencies and Python packages in one layer for size optimization
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Runtime dependencies (kept)
    libpng16-16 \
    libjpeg62-turbo \
    libtiff6 \
    libhdf5-103-1 \
    libopenblas0 \
    libgomp1 \
    curl \
    ca-certificates \
    # Build dependencies (removed after use)
    build-essential \
    gcc \
    g++ \
    libpng-dev \
    libjpeg-dev \
    libtiff-dev \
    libhdf5-dev \
    libopenblas-dev \
    && \
    # Upgrade pip
    pip install --no-cache-dir --upgrade pip setuptools wheel && \
    # Install all Python packages in one layer
    pip install --no-cache-dir \
        # Scientific computing (many come with TensorFlow base, use compatible versions)
        numpy>=1.24.0,<2.0.0 \
        scipy>=1.11.0 \
        pandas>=2.0.0 \
        scikit-learn>=1.3.0 \
        scikit-image>=0.21.0 \
        # Deep learning
        keras>=2.15.0 \
        tensorflow-hub>=0.15.0 \
        # Computer vision
        opencv-python-headless>=4.8.0 \
        Pillow>=10.0.0 \
        # Imaging formats
        tifffile>=2023.7.0 \
        imageio>=2.31.0 \
        # Visualization
        matplotlib>=3.7.0 \
        seaborn>=0.12.0 && \
    # Cleanup to reduce image size
    apt-get remove -y build-essential gcc g++ libpng-dev libjpeg-dev libtiff-dev libhdf5-dev libopenblas-dev && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /root/.cache/pip && \
    find /usr/local/lib/python*/dist-packages -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true

# =============================================================================
# REMOVED PACKAGES (for size optimization):
# =============================================================================
# The following packages were removed to reduce image size from ~18GB to ~8-10GB:
#
# - transformers>=4.36.0       (~2-3GB with dependencies)
# - datasets>=2.16.0           (~1GB with caching)
# - cellpose>=3.0.0            (~500MB-1GB)
# - git, wget                  (build tools, not needed at runtime)
#
# To add them back, uncomment this RUN statement:
# RUN pip install --no-cache-dir transformers>=4.36.0 datasets>=2.16.0 cellpose>=3.0.0
# =============================================================================

# =============================================================================
# Runtime
# =============================================================================
ENV TF_CPP_MIN_LOG_LEVEL=2

RUN groupadd -r biouser && useradd -r -g biouser biouser
RUN chown -R biouser:biouser /app

USER biouser
WORKDIR /app

CMD ["python"]
