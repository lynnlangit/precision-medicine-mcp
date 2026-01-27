# =============================================================================
# Precision Medicine - R Base Image (Optimized)  
# =============================================================================

FROM rocker/r-ver:4.3.2

LABEL maintainer="Precision Medicine Platform"
LABEL description="R base image for bioinformatics tools"
LABEL version="1.0.0"

ENV DEBIAN_FRONTEND=noninteractive

# =============================================================================
# System Dependencies (unchanged)
# =============================================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    gfortran \
    libopenblas-dev \
    liblapack-dev \
    libhdf5-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    git \
    curl \
    wget \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# =============================================================================
# R Configuration - Enable parallel compilation
# =============================================================================
WORKDIR /app

RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"))' >> /usr/local/lib/R/etc/Rprofile.site \
    && echo 'options(Ncpus = parallel::detectCores())' >> /usr/local/lib/R/etc/Rprofile.site

# =============================================================================
# R Packages - Split into layers, minimal dependencies
# =============================================================================

# Layer 1: Core packages (cached separately)
RUN R -e "install.packages(c('Rcpp', 'Matrix', 'data.table'), dependencies = FALSE)"

# Layer 2: BiocManager only
RUN R -e "install.packages('BiocManager', dependencies = FALSE)" \
    && R -e "BiocManager::install(version = '3.18', ask = FALSE, update = FALSE)"

# Layer 3: Bioconductor core (no suggested deps)
RUN R -e "BiocManager::install(c( \
    'Biobase', \
    'GenomicRanges', \
    'SummarizedExperiment' \
), ask = FALSE, update = FALSE, dependencies = c('Depends', 'Imports', 'LinkingTo'))"

# Layer 4: Analysis packages
RUN R -e "BiocManager::install(c( \
    'DESeq2', \
    'edgeR', \
    'limma' \
), ask = FALSE, update = FALSE, dependencies = c('Depends', 'Imports', 'LinkingTo'))"

# Layer 5: clusterProfiler (heavy, separate layer)
RUN R -e "BiocManager::install('clusterProfiler', \
    ask = FALSE, update = FALSE, dependencies = c('Depends', 'Imports', 'LinkingTo'))"

# Layer 6: Tidyverse (optional - consider removing if not needed)
# This alone adds 10-15 min build time
RUN R -e "install.packages('tidyverse', dependencies = FALSE)"

# =============================================================================
# Runtime
# =============================================================================
RUN groupadd -r biouser && useradd -r -g biouser biouser \
    && chown -R biouser:biouser /app

USER biouser
WORKDIR /app

CMD ["R"]