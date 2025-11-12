# Spatial MCP Demonstration POC

AI-Orchestrated Spatial Transcriptomics (and Multiomics) Bioinformatics Workflows using Model Context Protocol

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![MCP](https://img.shields.io/badge/MCP-2025--06--18-green.svg)](https://modelcontextprotocol.io/)
[![Claude Desktop](https://img.shields.io/badge/Claude-Desktop-orange.svg)](https://claude.ai/download)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

## What's In It For You?

**Stop writing glue code.** Describe your analysis goals in natural language and let Claude orchestrate bioinformatics pipelines across 9 specialized MCP servers with 36 tools.

**Key Benefits:**
- ✅ Replace bash scripts with conversational requests
- ✅ Instant access to FASTQ QC, alignment, quantification, multi-omics meta-analysis
- ✅ Reproducible by default - logged, versioned, repeatable
- ✅ Modular architecture - add tools without rewriting pipelines

**Example:** [18 ready-to-use prompts](docs/spatial/MCP_POC_Example_Prompts.md)

```
Claude, I have 10x Visium spatial transcriptomics data. Please:
1. Validate FASTQ files
2. Extract UMIs and spatial barcodes
3. Filter spots with <200 genes
4. Run differential expression
5. Perform pathway enrichment
6. Compare to TCGA cohorts
```

<img src="https://github.com/lynnlangit/spatial-mcp/blob/main/data/images/Claude-client.png" width=800>

*Single prompt → 6 tools across 4 servers → Complete analysis (no pipeline code)*

---

## Overview

**9 MCP Servers | 36 Tools | 58 Tests (100% Pass) | 80%+ Coverage**

![Architecture](https://github.com/lynnlangit/spatial-mcp/blob/main/architecture/spatial-mcp-arch.png)

| Server | Tools | Purpose |
|--------|-------|---------|
| mcp-fgbio | 4 | Reference data & FASTQ processing |
| mcp-spatialtools | 8 | Spatial processing & analysis |
| mcp-openimagedata | 3 | Image processing & registration |
| mcp-seqera | 3 | Nextflow orchestration |
| mcp-huggingface | 3 | ML genomics models |
| mcp-deepcell | 2 | Cell segmentation |
| mcp-mockepic | 3 | EHR integration |
| mcp-tcga | 5 | TCGA data integration |
| mcp-multiomics | 5 | Multi-omics meta-analysis |

---

## Quick Start

```bash
# Install (5 min)
git clone https://github.com/your-org/spatial-mcp.git
cd spatial-mcp/manual_testing
./install_dependencies.sh

# Configure Claude Desktop
cp ../configs/claude_desktop_config.json ~/Library/Application\ Support/Claude/claude_desktop_config.json

# Verify (restart Claude Desktop first)
./verify_servers.sh
```

**Prerequisites:** Python 3.11+, Claude Desktop, 16GB RAM, 50GB disk

---

## Example Workflows

**Multi-omics analysis:**
```
Analyze PDX treatment resistance with RNA, Protein, Phospho data for TP53, MYC, KRAS:
- RNA p-values: [0.001, 0.002, 0.05], log2FC: [2.5, 1.8, 1.2]
- Protein p-values: [0.005, 0.01, 0.03], log2FC: [2.0, 1.6, 1.1]

Combine using Stouffer's method with directionality and FDR correction.
```

**Spatial pipeline:**
```
Process 10x Visium: fetch hg38 → validate FASTQ → extract UMIs → align → quantify → compare TCGA
```

[View all 18 prompts →](docs/spatial/MCP_POC_Example_Prompts.md)

---

## Status

✅ **Production Ready** | 58 tests passing | 80.5% coverage

- Phase 1: Foundation (mcp-fgbio)  
- Phase 2: Core Processing (mcp-spatialtools, mcp-openimagedata)  
- Phase 3: Advanced Analysis (6 servers: seqera, huggingface, deepcell, mockepic, tcga, multiomics)  

---

## Technology

**Core:** FastMCP, Claude Desktop, Python 3.11+
**Bioinformatics:** FGbio, STAR, samtools, bedtools
**ML:** Hugging Face, PyTorch, DeepCell
**Workflows:** Nextflow (Seqera Platform)
**Statistics:** Stouffer's meta-analysis, scipy, statsmodels

---

## License & Resources

**License:** Apache 2.0 ([LICENSE](LICENSE))

**References:**
- [MCP Specification](https://modelcontextprotocol.io/specification/2025-06-18)
- [FastMCP Docs](https://github.com/modelcontextprotocol/python-sdk)
- [BioinfoMCP Paper](https://arxiv.org/html/2510.02139v1)
- [Spatial Transcriptomics Review](https://academic.oup.com/nar/article/53/12/gkaf536/8174767)

**Acknowledgments:** Model Context Protocol (Anthropic), BioinfoMCP, FGbio, TCGA, Seqera Platform

---

**Built with ❤️ for the bioinformatics community**
