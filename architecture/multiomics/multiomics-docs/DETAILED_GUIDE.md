# mcp-multiomics Documentation

This directory contains all documentation related to the **mcp-multiomics** server, the 9th MCP server in the Precision Medicine MCP suite.

---

## Overview

**mcp-multiomics** is a Model Context Protocol (MCP) server for multi-omics Patient-Derived Xenograft (PDX) data analysis. It provides tools for integrating and analyzing RNA, Protein, and Phosphorylation data using statistical meta-analysis methods.

**Part of the Precision Medicine MCP suite:** This server is a key component of comprehensive precision medicine workflows, demonstrated in the [PatientOne](../../architecture/patient-one/README.md) example where it identifies treatment resistance mechanisms in Stage IV Ovarian Cancer.

### Key Features

- **Multi-Modality Integration**: Load and align RNA, Protein, and Phosphorylation data
- **Stouffer's Meta-Analysis**: Combine p-values across modalities with directionality
- **Statistical Rigor**: FDR correction, weighted combinations, effect size preservation
- **Ready for Extension**: Mock implementations for HAllA, PCA, and visualization tools

### Server Location

```
/servers/mcp-multiomics/
```

### Documentation in This Directory

| File | Description |
|------|-------------|
| **README.md** | This file - Documentation index and quick reference |
| **Implementation_Plan_mcp_multiomics.md** | Original design document and implementation plan |
| **MULTIOMICS_INTEGRATION_SUMMARY.md** | Integration summary with test results and deployment guide |
| **PYDANTIC_BOOLEAN_FIX.md** | Technical details of Pydantic V2 boolean parsing fix |

---

## Quick Reference

### Tools Provided (5 tools)

1. **integrate_omics_data** ✅ Full Implementation
   - Load and align multi-omics data files
   - Z-score normalization across modalities
   - Missing data filtering with configurable thresholds
   - Quality control metrics and caching

2. **calculate_stouffer_meta** ✅ Full Implementation
   - Stouffer's Z-score meta-analysis
   - P-value to Z-score conversion (inverse normal)
   - Directionality from effect sizes (log2 fold changes)
   - FDR correction (Benjamini-Hochberg)
   - Weighted or unweighted combinations

3. **run_halla_analysis** 🔨 Mock (Ready for R/rpy2)
   - HAllA hierarchical association testing
   - Multi-omics relationship discovery

4. **create_multiomics_heatmap** 🔨 Mock (Ready for Implementation)
   - Multi-omics visualization
   - Heatmap generation across modalities

5. **run_multiomics_pca** 🔨 Mock (Ready for Implementation)
   - Principal Component Analysis
   - Dimensionality reduction across modalities

### Resources Provided (2 resources)

1. **multiomics://cache/integrated** - Cached integrated data
2. **multiomics://cache/meta_analysis** - Cached meta-analysis results

---

## Documentation Guide

### 1. Implementation Plan

**File:** `Implementation_Plan_mcp_multiomics.md`

**Read this if you want to:**
- Understand the original design goals
- Learn about the statistical methods (Stouffer's, HAllA, Fisher's)
- See the planned tool specifications
- Review the data format requirements
- Understand the caching strategy

**Contents:**
- Project rationale and scientific background
- Detailed tool specifications
- Data format examples (CSV with samples as columns)
- R/rpy2 integration plan for HAllA
- Statistical method comparisons
- Success metrics and validation plan

### 2. Integration Summary

**File:** `MULTIOMICS_INTEGRATION_SUMMARY.md`

**Read this if you want to:**
- Get an overview of what was built
- See test coverage and results
- Understand how to deploy the server
- Learn example usage patterns
- Review project totals and metrics

**Contents:**
- Complete server details (5 tools, 2 resources, 84% coverage)
- Tool implementation status
- Test suite summary (29 passing tests)
- Documentation update checklist
- Deployment instructions
- Example workflows
- Next steps for enhancement

### 3. Pydantic Boolean Fix

**File:** `PYDANTIC_BOOLEAN_FIX.md`

**Read this if you want to:**
- Understand the DRY_RUN mode bug
- Learn how Pydantic V2 parses boolean environment variables
- See the technical solution with code examples
- Troubleshoot configuration issues
- Understand the model_validator approach

**Contents:**
- Problem description (string "false" → boolean True)
- Root cause analysis
- Complete solution with code
- Verification tests and results
- Claude Desktop configuration
- Technical deep dive on Pydantic V2 Settings
- Alternative solutions considered

---

## Server Status

**Version:** 1.0.0 (Initial Release)
**Status:** ✅ **Production Ready**
**Test Coverage:** 84% (29/29 tests passing)
**Claude Desktop:** Configured and operational

### Statistics

| Metric | Value |
|--------|-------|
| **Total Lines of Code** | ~2,500 |
| **Tools Implemented** | 5 (2 full, 3 mock) |
| **Resources** | 2 |
| **Tests** | 29 (all passing) |
| **Coverage** | 84% |
| **Integration Tests** | 14 |
| **Stouffer Tests** | 15 |
| **Mock Tools** | 3 (ready for R/rpy2) |

---

## Configuration

### Environment Variables

The mcp-multiomics server uses these environment variables:

| Variable | Purpose | Default | Claude Desktop |
|----------|---------|---------|----------------|
| `MULTIOMICS_DATA_DIR` | Multi-omics data files | `/workspace/data/multiomics` | `/Users/.../precision-medicine-mcp/data/multiomics` |
| `MULTIOMICS_CACHE_DIR` | Cached results | `/workspace/cache/multiomics` | `/Users/.../precision-medicine-mcp/data/cache/multiomics` |
| `MULTIOMICS_DRY_RUN` | Mock execution mode | `true` | `false` |
| `MULTIOMICS_R_HOME` | R installation path | Auto-detect | Auto-detect |
| `MULTIOMICS_LOG_LEVEL` | Logging level | `INFO` | `INFO` |
| `MULTIOMICS_MAX_FEATURES` | Max features per modality | `5000` | `5000` |
| `MULTIOMICS_MIN_SAMPLES` | Min samples required | `3` | `3` |

### Claude Desktop Config

**Location:** `~/Library/Application Support/Claude/claude_desktop_config.json`

**Reference Config:** `/configs/claude_desktop_config.json`

```json
"multiomics": {
  "command": "/path/to/precision-medicine-mcp/servers/mcp-multiomics/venv/bin/python",
  "args": ["-m", "mcp_multiomics"],
  "cwd": "/path/to/precision-medicine-mcp/servers/mcp-multiomics",
  "env": {
    "PYTHONPATH": "/path/to/precision-medicine-mcp/servers/mcp-multiomics/src",
    "MULTIOMICS_DATA_DIR": "/path/to/precision-medicine-mcp/data/multiomics",
    "MULTIOMICS_CACHE_DIR": "/path/to/precision-medicine-mcp/data/cache/multiomics",
    "MULTIOMICS_DRY_RUN": "false"
  }
}
```

---

## Quick Start

### Installation

The server is installed as part of the Precision Medicine MCP project setup:

```bash
cd /path/to/precision-medicine-mcp/manual_testing
./install_dependencies.sh
```

Or manually:

```bash
cd /path/to/precision-medicine-mcp/servers/mcp-multiomics
python3.11 -m venv venv
source venv/bin/activate
pip install -e ".[dev]"
```

### Running Tests

```bash
cd /path/to/precision-medicine-mcp/servers/mcp-multiomics

# All tests
MULTIOMICS_DRY_RUN="false" venv/bin/python -m pytest tests/ -v

# Integration tests only
MULTIOMICS_DRY_RUN="false" venv/bin/python -m pytest tests/test_integration.py -v

# Stouffer tests only
MULTIOMICS_DRY_RUN="false" venv/bin/python -m pytest tests/test_stouffer.py -v

# With coverage report
MULTIOMICS_DRY_RUN="false" venv/bin/python -m pytest tests/ -v --cov
```

### Example Usage in Claude Desktop

**Example 1: Multi-Omics Integration**

```
I have PDX multi-omics data:
- RNA: /data/multiomics/pdx_rna.csv
- Protein: /data/multiomics/pdx_protein.csv
- Phospho: /data/multiomics/pdx_phospho.csv
- Metadata: /data/multiomics/sample_metadata.csv

Please integrate these datasets with Z-score normalization.
```

**Example 2: Stouffer's Meta-Analysis**

```
I have differential analysis results for genes TP53, MYC, KRAS
across RNA and Protein modalities:

RNA:
- p-values: [0.001, 0.002, 0.05]
- log2FC: [2.5, 1.8, 1.2]

Protein:
- p-values: [0.005, 0.01, 0.03]
- log2FC: [2.0, 1.6, 1.1]

Combine these using Stouffer's method with directionality at FDR < 0.05.
```

**Example 3: Full Workflow**

```
Analyze treatment resistance in PDX models:

1. Integrate RNA, Protein, and Phospho data from /data/multiomics/
2. Run differential analysis (using external tools)
3. Combine p-values using Stouffer's weighted method
4. Identify genes consistently dysregulated across all modalities (FDR < 0.05)
5. Generate heatmap of significant features
```

---

## Data Formats

### Input Data Format

**CSV files with samples as columns:**

```csv
feature_id,sample1,sample2,sample3,...
TP53,8.5,9.2,7.8,...
MYC,10.1,9.8,10.5,...
KRAS,7.2,8.1,7.9,...
```

### Metadata Format

```csv
sample_id,treatment,timepoint,batch
sample1,Resistant,Day7,Batch1
sample2,Sensitive,Day7,Batch1
sample3,Resistant,Day14,Batch2
```

### Output Format

**Stouffer's Meta-Analysis Results:**

```json
{
  "features": ["TP53", "MYC", "KRAS"],
  "meta_p_values": [1.2e-6, 3.5e-5, 0.0008],
  "meta_z_scores": [4.85, 4.23, 3.35],
  "q_values": [3.6e-6, 5.2e-5, 0.0008],
  "significant_features": ["TP53", "MYC"],
  "n_significant": 2,
  "directionality_used": true,
  "weighted": true
}
```

---

## Testing

### Test Structure

```
tests/
├── conftest.py              # Pytest configuration and shared fixtures
├── test_integration.py      # Data integration tests (14 tests)
├── test_stouffer.py         # Meta-analysis tests (15 tests)
└── fixtures/
    ├── generate_fixtures.py # Synthetic data generator
    ├── sample_rna.csv       # 1000 genes × 15 samples
    ├── sample_protein.csv   # 500 proteins × 15 samples
    ├── sample_phospho.csv   # 300 sites × 15 samples
    └── sample_metadata.csv  # 15 samples (7 Resistant, 8 Sensitive)
```

### Test Coverage

- **Integration Tests (14)**: Data loading, alignment, normalization, filtering, caching
- **Stouffer Tests (15)**: P-value conversions, Z-score combinations, FDR correction, directionality

### Test Data

Realistic synthetic data with:
- 1000 RNA features
- 500 Protein features
- 300 Phosphorylation sites
- 15 samples (7 Resistant, 8 Sensitive)
- Known differential patterns for validation

---

## Scientific Background

### Stouffer's Method

**Purpose:** Combine p-values from independent tests while preserving directionality of effects.

**Method:**
1. Convert p-values to Z-scores using inverse normal transformation
2. Optionally incorporate effect size directionality (signs from log2FC)
3. Combine Z-scores using weighted or unweighted mean
4. Convert combined Z-score back to p-value
5. Apply FDR correction (Benjamini-Hochberg)

**Advantages:**
- Preserves effect directionality
- Handles both concordant and discordant effects
- More powerful than Fisher's method for large effect sizes
- Supports sample-size weighting

**Reference:** Whitlock MC (2005). Combining probability from independent tests: the weighted Z-method is superior to Fisher's approach. *J Evol Biol* 18(5):1368-73.

### HAllA (Planned)

**Purpose:** Discover hierarchical associations between multi-omics layers.

**Method:**
- Hierarchical clustering of features within each modality
- All-against-all association testing between clusters
- FDR correction for multiple testing
- Association network visualization

**Reference:** Rahnavard et al. (2017). High-sensitivity pattern discovery in large, paired multiomic datasets. *Bioinformatics* 33(14):i81-i89.

### Multi-Omics Integration

**Purpose:** Integrate complementary data types to understand complex biological systems.

**Modalities Supported:**
- **RNA (Transcriptomics)**: Gene expression levels
- **Protein (Proteomics)**: Protein abundance
- **Phospho (Phosphoproteomics)**: Protein phosphorylation sites

**Reference:** Menyhárt et al. (2023). Multi-omics approaches in cancer research with applications in tumor subtyping, prognosis, and diagnosis. *Comput Struct Biotechnol J* 21:1864-1883.

---

## Development

### Project Structure

```
servers/mcp-multiomics/
├── src/mcp_multiomics/
│   ├── __init__.py
│   ├── __main__.py
│   ├── server.py           # FastMCP server with tool implementations
│   ├── config.py           # Pydantic Settings configuration
│   ├── tools/
│   │   ├── __init__.py
│   │   ├── integration.py  # Data integration tools
│   │   ├── stouffer.py     # Stouffer's meta-analysis
│   │   └── utils.py        # Shared utilities
│   ├── resources/          # MCP resources
│   └── r_interface/        # R/rpy2 integration (future)
├── tests/                  # Test suite (29 tests, 84% coverage)
├── pyproject.toml          # Package configuration
├── .env.example            # Environment variable template
└── README.md               # Server-specific documentation
```

### Key Technologies

- **FastMCP**: MCP server framework
- **Pydantic V2**: Configuration and validation
- **NumPy**: Numerical computations
- **SciPy**: Statistical functions
- **Pandas**: Data manipulation
- **Pytest**: Testing framework

### Adding New Tools

1. Create tool function in `src/mcp_multiomics/server.py`
2. Add implementation in `src/mcp_multiomics/tools/`
3. Write tests in `tests/test_<toolname>.py`
4. Update documentation
5. Update test fixtures if needed

### Code Style

- **Type Hints**: All functions fully typed
- **Docstrings**: Google style with examples
- **Testing**: Pytest with fixtures
- **Coverage**: Maintain > 80%
- **Formatting**: Follow existing patterns

---

## Troubleshooting

### Issue: DRY_RUN mode always enabled

**Solution:** See `PYDANTIC_BOOLEAN_FIX.md` for complete details.

```bash
# Verify config
MULTIOMICS_DRY_RUN="false" python -c "from mcp_multiomics.config import config; print(f'DRY_RUN: {config.dry_run}')"
# Should print: DRY_RUN: False
```

### Issue: Environment variables not loading

**Solution:** Check `env_prefix` in config.py and variable names:

```python
# config.py should have:
model_config = SettingsConfigDict(
    env_prefix="MULTIOMICS_",  # Required!
)
```

### Issue: Tests failing

**Solution:** Run with DRY_RUN disabled:

```bash
MULTIOMICS_DRY_RUN="false" pytest tests/ -v
```

### Issue: Claude Desktop not loading server

**Solutions:**
1. Verify venv path in config: `ls /path/to/venv/bin/python`
2. Check JSON syntax: `python -m json.tool ~/Library/Application\ Support/Claude/claude_desktop_config.json`
3. Restart Claude Desktop completely (Quit and relaunch)

---

## Future Enhancements

### Phase 4 (Planned)

1. **HAllA Implementation**
   - R/rpy2 integration
   - Association network generation
   - Visualization tools

2. **Advanced Visualization**
   - Multi-omics heatmaps
   - PCA plots with modality colors
   - Association networks

3. **Additional Meta-Analysis Methods**
   - Fisher's combined probability test
   - Minimum p-value method
   - Rank aggregation methods

4. **Network Analysis**
   - Multi-layer network integration
   - Pathway enrichment
   - Network topology analysis

5. **Machine Learning**
   - Feature selection across modalities
   - Predictive modeling
   - Integration with scikit-learn

---

## References

### Scientific Publications

1. **Stouffer's Method**
   Whitlock MC (2005). Combining probability from independent tests: the weighted Z-method is superior to Fisher's approach. *Journal of Evolutionary Biology* 18(5):1368-1373.

2. **HAllA**
   Rahnavard A, Mann B, Giri A, Chatterjee R, Crandall KA (2017). High-sensitivity pattern discovery in large, paired multiomic datasets. *Bioinformatics* 33(14):i81-i89.

3. **Multi-Omics Integration**
   Menyhárt O, Weltz B, Győrffy B (2023). Multi-omics approaches in cancer research with applications in tumor subtyling, prognosis, and diagnosis. *Computational and Structural Biotechnology Journal* 21:1864-1883.

4. **FDR Correction**
   Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B* 57(1):289-300.

### Technical Documentation

- **FastMCP**: https://github.com/jlowin/fastmcp
- **Pydantic Settings**: https://docs.pydantic.dev/latest/concepts/pydantic_settings/
- **MCP Protocol**: https://modelcontextprotocol.io/

---

## Changelog

### Version 1.0.0 (November 11, 2025)

**Added:**
- Initial release of mcp-multiomics server
- `integrate_omics_data` tool (full implementation)
- `calculate_stouffer_meta` tool (full implementation with directionality)
- `run_halla_analysis` tool (mock, ready for R integration)
- `create_multiomics_heatmap` tool (mock)
- `run_multiomics_pca` tool (mock)
- Comprehensive test suite (29 tests, 84% coverage)
- Pydantic V2 configuration with environment variable support
- Claude Desktop integration

**Fixed:**
- Pydantic V2 boolean parsing issue (string "false" → boolean False)
- Environment variable loading with env_prefix
- Directory creation moved out of __init__

**Documentation:**
- Server README
- Implementation plan
- Integration summary
- Pydantic boolean fix guide

---

## Support

For issues, questions, or contributions:

1. **Server Issues**: Create issue in precision-medicine-mcp repository
2. **Documentation**: Suggest edits via pull request
3. **Scientific Questions**: Refer to publications in References section

---

**Last Updated:** November 11, 2025
**Status:** ✅ Production Ready
**Author:** Claude (Sonnet 4.5)
**Project:** Precision Medicine MCP - Multi-Omics Integration Server
