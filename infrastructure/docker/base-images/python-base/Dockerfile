# =============================================================================
# Precision Medicine - Python Base Image
# =============================================================================
# Base container for Python-based bioinformatics tools
# =============================================================================
FROM python:3.12-slim-bookworm

LABEL maintainer="Precision Medicine Platform"
LABEL description="Python base image for bioinformatics tools"
LABEL version="1.0.0"

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1
ENV PIP_NO_CACHE_DIR=1

# =============================================================================
# System Dependencies
# =============================================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    gfortran \
    libopenblas-dev \
    liblapack-dev \
    libhdf5-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    ca-certificates \
    curl \
    wget \
    git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# =============================================================================
# Python Packages
# =============================================================================
WORKDIR /app

RUN pip install --upgrade pip setuptools wheel

# Scientific computing
RUN pip install \
    numpy>=1.24.0 \
    scipy>=1.11.0 \
    pandas>=2.0.0 \
    scikit-learn>=1.3.0 \
    statsmodels>=0.14.0

# Bioinformatics
RUN pip install \
    biopython>=1.81 \
    pysam>=0.22.0 \
    pyvcf3>=1.0.0 \
    h5py>=3.9.0

# Visualization
RUN pip install \
    matplotlib>=3.7.0 \
    seaborn>=0.12.0

# =============================================================================
# Runtime
# =============================================================================
RUN groupadd -r biouser && useradd -r -g biouser biouser
RUN chown -R biouser:biouser /app

USER biouser
WORKDIR /app

CMD ["python"]
