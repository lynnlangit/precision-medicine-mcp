# Spatial MCP Demonstration POC

AI-Orchestrated Spatial Transcriptomics Bioinformatics Pipeline using Model Context Protocol

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![MCP](https://img.shields.io/badge/MCP-2025--06--18-green.svg)](https://modelcontextprotocol.io/)
[![Claude Desktop](https://img.shields.io/badge/Claude-Desktop-orange.svg)](https://claude.ai/download)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

## What's In It For You?

**Are you tired of writing glue code to connect bioinformatics tools?** This POC demonstrates a fundamentally different approach to spatial transcriptomics analysis: instead of manually chaining together FASTQ validators, aligners, expression quantifiers, and statistical tools through brittle shell scripts, you **describe your analysis goals in natural language** and let Claude orchestrate the entire pipeline.

**Why this matters for bioinformaticians:**
- ✅ **No more bash scripts from hell** - Replace complex pipeline code with conversational analysis requests
- ✅ **Instant access to 36 genomics tools** - From FASTQ QC to multi-omics meta-analysis, all via natural language
- ✅ **Reproducible by default** - Every analysis is logged, versioned, and repeatable
- ✅ **Modular & extensible** - Add new tools as MCP servers without rewriting your pipeline
- ✅ **Tested with real workflows** - [18 complete example prompts](docs/MCP_POC_Example_Prompts.md) from QC to publication-ready analysis

**Try it yourself:**
```
Claude, I have 10x Visium spatial transcriptomics data. Please:
1. Validate my FASTQ files (sample_R1.fastq.gz, sample_R2.fastq.gz)
2. Extract UMIs and spatial barcodes
3. Filter spots with <200 genes detected
4. Perform differential expression between tumor and normal regions
5. Run pathway enrichment on upregulated genes
6. Compare key marker genes to TCGA breast cancer cohorts
```
<img src="https://github.com/lynnlangit/spatial-mcp/blob/main/data/images/Claude-client.png" width=800>

This single prompt orchestrates **6 different tools across 4 MCP servers**, generates statistical results, biological interpretations, and TCGA comparisons—all without writing a single line of pipeline code.

---

## Overview

This project demonstrates the power of **Model Context Protocol (MCP)** in orchestrating complex bioinformatics workflows for spatial transcriptomics and multi-omics analysis. Using Claude Desktop as the AI orchestrator, the system coordinates **9 specialized MCP servers** with **36 tools** to process spatial genomics data through a complete analysis pipeline.

### Key Features

- 🤖 **AI-Driven Orchestration** - Claude coordinates multi-server workflows via natural language
- 🧬 **Modular Architecture** - 9 specialized servers, each with single responsibility
- 🔒 **Production-Ready** - Comprehensive testing (58 tests, 100% pass rate, 80%+ coverage)
- 🚀 **Scalable Design** - Containerized, cloud-ready deployment
- 📊 **Industry Tools** - FGbio, TCGA, Hugging Face, Seqera Platform, multi-omics integration

### Architecture

![Spatial MCP Architecture](https://github.com/lynnlangit/spatial-mcp/blob/main/architecture/spatial-mcp-arch.png)

**Pipeline Flow:** Ingest & QC → Segment Spatial → Align Reads → Quantify Expression → Analyze Results

### MCP Servers

| Server | Tools | Purpose |
|--------|-------|---------|
| **mcp-fgbio** | 4 | Genomic reference data & FASTQ processing |
| **mcp-spatialtools** | 8 | Core spatial processing + advanced analysis |
| **mcp-openimagedata** | 3 | Histology image processing & spatial registration |
| **mcp-seqera** | 3 | Nextflow workflow orchestration via Seqera Platform |
| **mcp-huggingface** | 3 | ML models for genomics (DNABERT, Geneformer, scGPT) |
| **mcp-deepcell** | 2 | Deep learning cell segmentation and phenotyping |
| **mcp-mockepic** | 3 | Mock Epic EHR integration with synthetic patient data |
| **mcp-tcga** | 5 | TCGA cancer genomics data integration |
| **mcp-multiomics** | 5 | Multi-omics PDX data integration & Stouffer's meta-analysis |
| **TOTAL** | **36** | **All servers operational** ✅ |

---

## Quick Start

### Installation (5 minutes)

```bash
# Clone repository
git clone https://github.com/your-org/spatial-mcp.git
cd spatial-mcp

# Install all 9 MCP servers
cd manual_testing
./install_dependencies.sh

# Verify installation
./verify_servers.sh
# Expected: "Servers working: 9/9, Total tools: 36"
```

**Prerequisites:** Python 3.11+, Claude Desktop, 16GB+ RAM, 50GB disk space

### Configure Claude Desktop

```bash
# Copy configuration file
cp configs/claude_desktop_config.json ~/Library/Application\ Support/Claude/claude_desktop_config.json

# Restart Claude Desktop
```

### Verify Setup

In Claude Desktop:
```
What MCP servers are available?
```

Expected: 9 servers listed with 36 total tools

---

## Documentation

### Getting Started
- 🚀 **[Manual Testing Guide](manual_testing/README.md)** - Installation, verification, and testing
- 📖 **[Setup Guide](docs/setup_guide.md)** - Detailed installation and configuration
- ⚙️ **[Configuration Guide](configs/README.md)** - Claude Desktop setup and environment variables

### Using the System
- 📋 **[Example Prompts](docs/MCP_POC_Example_Prompts.md)** - 18 ready-to-use test prompts
- 🧬 **[Multi-Omics Guide](docs/multiomics/README.md)** - PDX analysis and Stouffer's meta-analysis
- 🎯 **[Test Results](manual_testing/TEST_RESULTS_ALL_SERVERS.md)** - Complete test suite results

### Technical Documentation
- 🏗️ **[Architecture Document](architecture/Spatial_MCP_POC_Architecture.md)** - Full technical architecture
- 🎨 **[Architecture Diagram](architecture/Spatial_MCP_Architecture_Diagram.html)** - Visual overview
- 🔧 **[Server READMEs](servers/)** - Individual server documentation

---

## Example Usage

### Simple Queries

**Fetch reference genome:**
```
Claude, download the hg38 reference genome and tell me about its characteristics.
```

**Validate FASTQ data:**
```
I have a FASTQ file at /data/sample.fastq.gz. Validate it and check quality.
```

### Complete Workflows

**Multi-omics treatment resistance analysis:**
```
I'm analyzing PDX treatment resistance with RNA, Protein, and Phospho data.
For genes TP53, MYC, KRAS:
- RNA p-values: [0.001, 0.002, 0.05], log2FC: [2.5, 1.8, 1.2]
- Protein p-values: [0.005, 0.01, 0.03], log2FC: [2.0, 1.6, 1.1]

Combine using Stouffer's method with directionality and FDR correction.
```

**Spatial transcriptomics pipeline:**
```
Process 10x Visium data:
1. Fetch hg38 reference
2. Validate FASTQ files
3. Extract UMIs
4. Align reads
5. Quantify expression
6. Compare to TCGA cohorts
```

See [18 complete example prompts](docs/MCP_POC_Example_Prompts.md) for more workflows.

---

## Project Status

**Current:** ✅ **Production Ready - All 9 Servers Operational**

- ✅ Phase 1: Foundation (mcp-fgbio, testing framework)
- ✅ Phase 2: Core Processing (mcp-spatialtools, mcp-openimagedata)
- ✅ Phase 3: Advanced Analysis (6 additional servers: seqera, huggingface, deepcell, mockepic, tcga, multiomics)

**Test Results:** 58 tests, 100% pass rate, 80.5% average coverage

**Next:** Optional Phase 4 enhancements (HAllA implementation, advanced visualizations, CI/CD pipeline)

---

## Technology Stack

- **MCP Framework:** FastMCP (Python), Claude Desktop
- **Bioinformatics:** FGbio, STAR, samtools, bedtools
- **Machine Learning:** Hugging Face Transformers, PyTorch, DeepCell
- **Workflows:** Nextflow via Seqera Platform
- **Statistics:** Stouffer's meta-analysis, scipy, statsmodels
- **Multi-Omics:** pandas, numpy, scikit-learn
- **Testing:** pytest (58 tests, 100% pass rate)

---

## Contributing

Contributions welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass (`pytest`)
5. Submit a pull request

See [Manual Testing Guide](manual_testing/README.md) for development setup.

---

## License

Apache License 2.0 - see [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- **Model Context Protocol** - Anthropic's open standard for AI-data integration
- **BioinfoMCP** - Inspiration for bioinformatics tool integration
- **FGbio** - Fulcrum Genomics bioinformatics toolkit
- **TCGA** - The Cancer Genome Atlas program
- **Seqera Platform** - Nextflow workflow orchestration

---

## Resources

- [MCP Specification](https://modelcontextprotocol.io/specification/2025-06-18)
- [FastMCP Documentation](https://github.com/modelcontextprotocol/python-sdk)
- [BioinfoMCP Paper](https://arxiv.org/html/2510.02139v1)
- [Spatial Transcriptomics Review](https://academic.oup.com/nar/article/53/12/gkaf536/8174767)

**Issues & Discussions:** [GitHub Issues](https://github.com/your-org/spatial-mcp/issues) | [GitHub Discussions](https://github.com/your-org/spatial-mcp/discussions)

---

**Built with ❤️ for the bioinformatics community**
