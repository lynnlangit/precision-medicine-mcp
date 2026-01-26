# Change Log - Base Docker Images

## 2026-01-25 - Initial Release

### Overview
Created optimized Docker base images for bioinformatics and machine learning workflows.

### New Images

#### python-base
- Python 3.12 slim base image
- Bioinformatics libraries: biopython, pysam, pyvcf3
- Data science stack: numpy, scipy, pandas, scikit-learn
- Visualization: matplotlib, seaborn
- Non-root user execution (biouser)

#### r-base
- R 4.3.2 with Bioconductor packages
- Genomics: DESeq2, edgeR, limma, GenomicRanges
- Pathway analysis: clusterProfiler
- Data manipulation: tidyverse, data.table
- Optimized layer compilation

#### tensorflow-base
- TensorFlow 2.15.0 with Keras
- **Size optimized**: Reduced from ~18GB to ~8-10GB (GPU) or ~4-5GB (CPU)
- Includes: scikit-learn, scikit-image, opencv, Pillow
- GPU/CPU version switching via Dockerfile comments
- Removed Hugging Face transformers/datasets (~3GB savings)
- Removed cellpose (~500MB-1GB savings)

### Security
- All images run as non-root user (biouser)
- Minimal attack surface with slim base images
- Regular rebuild recommended to get security patches

### Notes
- Images are for LOCAL USE ONLY
- Do not push to remote registries
