# FASTQ Alignment Workflow

**Status:** ⚠️ Implemented but not tested in PatientOne
**Last Updated:** January 9, 2026

---

## Overview

The FASTQ workflow processes raw spatial transcriptomics sequencing data from FASTQ files through alignment to generate expression matrices.

**Status:** Implemented in mcp-spatialtools but not used in current PatientOne tests.

**Why not used?**
- Requires STAR aligner installation (large binary ~3GB)
- Requires genome reference indices (~30GB for hg38)
- Longer processing time (~10 min per sample)
- PatientOne focuses on analysis workflows, not data processing

---

## Workflow Steps

### 1. Quality Filtering
**Tool:** `filter_quality`

Filters reads based on:
- Barcode quality scores
- UMI uniqueness
- Spatial barcode whitelist matching

### 2. Alignment
**Tool:** `align_spatial_data`

Uses STAR aligner to:
- Align reads to reference genome (hg38, mm10, etc.)
- Generate BAM files with spatial tags
- Produce alignment statistics

**Requirements:**
- STAR aligner v2.7.11+
- Genome reference indices
- 30+ GB RAM for human genome

### 3. Expression Quantification
Counts UMIs per spot per gene to generate expression matrix.

**Output:** CSV file compatible with CSV workflow

---

## Installation Requirements

### STAR Aligner
```bash
# Via conda (recommended)
conda install -c bioconda star

# Via homebrew (macOS)
brew install star

# Verify
STAR --version
```

### Genome Indices
Download and build STAR indices for your reference genome.
See [mcp-spatialtools/INSTALL_STAR.md](../../servers/mcp-spatialtools/INSTALL_STAR.md)

---

## When to Use FASTQ Workflow

**Use FASTQ workflow when:**
- Processing raw sequencing data from sequencer
- Need custom alignment parameters
- Performing quality control from scratch
- Building new reference genomes

**Use CSV workflow when:**
- Working with pre-processed Visium data
- Focus on analysis, not data processing
- Faster iteration needed
- Limited computational resources

---

**See Also:**
- [CSV_WORKFLOW.md](CSV_WORKFLOW.md) - Current production workflow
- [mcp-spatialtools INSTALL_STAR.md](../../servers/mcp-spatialtools/INSTALL_STAR.md)
