# Spatial Transcriptomics Documentation

This directory contains documentation for the core spatial transcriptomics workflow and the initial 3-server implementation (Phase 1-2).

---

## Documents in This Directory

### Implementation & Testing

| File | Description | Phase |
|------|-------------|-------|
| **FINAL_IMPLEMENTATION_SUMMARY.md** | Complete implementation summary for Phase 1 (mcp-fgbio) | Phase 1 |
| **FINAL_TEST_RESULTS.md** | Comprehensive test results for Phase 1 | Phase 1 |
| **PHASE_2_SUMMARY.md** | Phase 2 implementation (mcp-spatialtools, mcp-openimagedata) | Phase 2 |

### User Guides

| File | Description |
|------|-------------|
| **MCP_POC_Example_Prompts.md** | 18 ready-to-use test prompts for Claude Desktop |
| **setup_guide.md** | Detailed installation and configuration guide |

---

## Quick Reference

### What is Spatial Transcriptomics?

Spatial transcriptomics combines gene expression profiling with spatial information, allowing researchers to:
- Map gene expression to specific tissue locations
- Identify spatially distinct cell populations
- Understand tissue architecture and microenvironments
- Correlate molecular signatures with histology

### Phase 1-2 Servers (Core Workflow)

**Phase 1: Foundation**
- **mcp-fgbio** (4 tools) - Reference data, FASTQ validation, UMI extraction, gene annotations

**Phase 2: Spatial Processing**
- **mcp-spatialtools** (8 tools) - Spatial data processing and advanced analysis
- **mcp-openimagedata** (3 tools) - Histology image processing and spatial registration

### Example Workflows

For complete workflow examples, see **MCP_POC_Example_Prompts.md** which includes:
- FASTQ validation and QC
- UMI extraction
- Spatial barcode processing
- Gene expression quantification
- Differential expression analysis
- Spatial clustering
- Visualization generation

---

## Related Documentation

### Current Project
- **[Multi-Omics Documentation](../../multiomics/README.md)** - Phase 3 multi-omics integration (mcp-multiomics server)
- **[Main README](../../../README.md)** - Project overview and quick start
- **[Manual Testing](../../../tests/manual_testing/README.md)** - Testing and verification guides
- **[Configuration](../../../desktop-configs/README.md)** - Claude Desktop configuration

### Phase 3+ Servers (Not in this directory)
Phase 3 introduced 6 additional servers for advanced analysis:
- mcp-seqera (Nextflow workflows)
- mcp-huggingface (ML genomics models)
- mcp-deepcell (Cell segmentation)
- mcp-mockepic (EHR integration)
- mcp-tcga (TCGA data)
- mcp-multiomics (Multi-omics integration)

Documentation for these servers is in their respective directories.

---

## Key Technologies (Spatial Workflow)

- **10x Visium** - Spatial gene expression platform
- **FGbio** - Genomic reference data and FASTQ processing
- **STAR** - RNA-seq alignment
- **Scanpy** - Single-cell/spatial analysis in Python
- **Squidpy** - Spatial molecular data analysis

---

## Getting Started

1. **Install all servers:**
   ```bash
   cd manual_testing
   ./install_dependencies.sh
   ```

2. **Configure Claude Desktop:**
   ```bash
   cp configs/claude_desktop_config.json ~/Library/Application\ Support/Claude/claude_desktop_config.json
   ```

3. **Try example prompts:**
   - Open `MCP_POC_Example_Prompts.md`
   - Copy a prompt to Claude Desktop
   - Watch the multi-server orchestration

---

## Document Timeline

| Document | Created | Phase | Status |
|----------|---------|-------|--------|
| setup_guide.md | Oct 24, 2025 | Phase 1 | ✅ Current |
| FINAL_IMPLEMENTATION_SUMMARY.md | Oct 24, 2025 | Phase 1 | ✅ Complete |
| FINAL_TEST_RESULTS.md | Oct 24, 2025 | Phase 1 | ✅ Complete |
| PHASE_2_SUMMARY.md | Oct 24, 2025 | Phase 2 | ✅ Complete |
| MCP_POC_Example_Prompts.md | Oct 24, 2025 | All Phases | ✅ Current (updated for 9 servers) |

---

**Directory Created:** November 11, 2025
**Contents:** Phase 1-2 spatial transcriptomics documentation
**Status:** ✅ Complete - Historical reference for initial implementation
