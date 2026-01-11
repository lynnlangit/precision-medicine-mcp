# Spatial Transcriptomics Servers

**Last Updated:** January 10, 2026

This document provides quick links to the 3 core servers used in spatial transcriptomics workflows.

---

## Server Overview

📋 **[See Complete Server Status →](../../servers/README.md#-server-status)** - All 10 MCP servers with tools, status, and documentation

**Spatial Transcriptomics Servers:**
- **mcp-fgbio** - Reference data and FASTQ validation (4 tools, 95% real) - [README](../../servers/mcp-fgbio/README.md)
- **mcp-spatialtools** - Spatial analysis and visualization (14 tools, 95% real) - [README](../../servers/mcp-spatialtools/README.md)
- **mcp-openimagedata** - Image data access (5 tools, 60% real) - [README](../../servers/mcp-openimagedata/README.md)

---

## Server Roles

### 1. mcp-fgbio (Phase 1: Foundation)
**Purpose:** Reference data and FASTQ validation

**Key Tools:**
- Fetch reference genomes and annotations
- Validate FASTQ files
- Extract UMIs
- Provide gene annotations

**Documentation:** [/servers/mcp-fgbio/README.md](../../servers/mcp-fgbio/README.md)

---

### 2. mcp-spatialtools (Phase 2: Spatial Processing)
**Purpose:** Spatial transcriptomics data analysis

**Key Capabilities:**
- Quality filtering and region segmentation
- Spatial autocorrelation (Moran's I, Geary's C)
- Differential expression analysis
- Batch correction
- Pathway enrichment
- Cell type deconvolution
- **NEW:** 4 visualization tools (heatmaps, region charts, autocorrelation plots)

**Documentation:** [/servers/mcp-spatialtools/README.md](../../servers/mcp-spatialtools/README.md)

---

### 3. mcp-openimagedata (Phase 2: Histology Integration)
**Purpose:** Histology image processing and spatial registration

**Key Capabilities:**
- H&E and MxIF image loading
- Spatial registration with transcriptomics
- Feature extraction
- **NEW:** MxIF RGB composites and H&E annotation overlays

**Documentation:** [/servers/mcp-openimagedata/README.md](../../servers/mcp-openimagedata/README.md)

---

## Complete Server Status

For deployment status of ALL 10 MCP servers (including multiomics, TCGA, seqera, etc.):

**See:** [/docs/SERVER_IMPLEMENTATION_STATUS.md](../../docs/SERVER_IMPLEMENTATION_STATUS.md)

---

## Related Documentation

- **Spatial Workflows:** [CSV_WORKFLOW.md](CSV_WORKFLOW.md), [FASTQ_WORKFLOW.md](FASTQ_WORKFLOW.md)
- **PatientOne Tests:** [/tests/manual_testing/PatientOne-OvarianCancer/](../../tests/manual_testing/PatientOne-OvarianCancer/)
- **Main Architecture:** [../README.md](../README.md)
