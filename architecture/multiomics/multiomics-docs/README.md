# Multi-Omics Integration Documentation

This directory contains documentation for the **mcp-multiomics** server (Phase 3), which provides multi-omics PDX data integration and Stouffer's meta-analysis capabilities.

---

## Documents in This Directory

### Active Documentation

| File | Description | Lines |
|------|-------------|-------|
| **README.md** | This file - Quick reference and navigation | 241 |
| **DETAILED_GUIDE.md** | Comprehensive server guide (tools, testing, config, examples) | 584 |
| **MULTIOMICS_INTEGRATION_SUMMARY.md** | Original integration summary (see README for current info) | 217 |

### Archived Documentation

| File | Description |
|------|-------------|
| **archive/Implementation_Plan_mcp_multiomics.md** | Pre-implementation design document |

### Related Files

| File | Location | Description |
|------|----------|-------------|
| **TROUBLESHOOTING.md** | `/servers/mcp-multiomics/` | Technical fixes (Pydantic boolean parsing, etc.) |
| **README.md** | `/servers/mcp-multiomics/` | Server-specific documentation |

---

## Quick Reference

### What is mcp-multiomics?

A Model Context Protocol server for **Patient-Derived Xenograft (PDX) multi-omics analysis**:
- Integrate RNA, Protein, and Phosphorylation data
- Combine p-values across modalities with Stouffer's meta-analysis
- Preserve effect size directionality
- Apply FDR correction (Benjamini-Hochberg)
- Identify genes consistently dysregulated across all modalities

### Tools Provided (5)

| Tool | Status | Description |
|------|--------|-------------|
| **integrate_omics_data** | ✅ Full | Load and align multi-omics data with normalization |
| **calculate_stouffer_meta** | ✅ Full | Stouffer's Z-score meta-analysis with directionality |
| **run_halla_analysis** | 🔨 Mock | HAllA hierarchical association testing (R/rpy2 ready) |
| **create_multiomics_heatmap** | 🔨 Mock | Multi-omics visualization |
| **run_multiomics_pca** | 🔨 Mock | PCA across modalities |

### Test Coverage

- **Tests:** 29 (all passing)
- **Coverage:** 84%
- **Integration Tests:** 14 (data loading, alignment, normalization)
- **Stouffer Tests:** 15 (conversions, combinations, FDR, directionality)

---

## Example Usage

### Simple Meta-Analysis

```
I have p-values from RNA and Protein analyses for genes TP53, MYC, KRAS:

RNA p-values: [0.001, 0.002, 0.05]
RNA effect sizes (log2FC): [2.5, 1.8, 1.2]

Protein p-values: [0.005, 0.01, 0.03]
Protein effect sizes (log2FC): [2.0, 1.6, 1.1]

Please combine these using Stouffer's method with directionality.
```

### PDX Treatment Resistance Workflow

```
I'm analyzing PDX treatment resistance with RNA, Protein, and Phospho data.
For candidate genes TP53, MYC, KRAS, EGFR, PIK3CA:

RNA-seq results:
- p-values: [0.0001, 0.0005, 0.002, 0.015, 0.045]
- log2FC: [3.2, 2.8, 2.1, -1.8, 1.5]
- samples: 15

Proteomics results:
- p-values: [0.0008, 0.001, 0.008, 0.020, 0.060]
- log2FC: [2.5, 2.3, 1.9, -1.5, 1.2]
- samples: 15

Phosphoproteomics results:
- p-values: [0.002, 0.003, 0.015, 0.035, 0.080]
- log2FC: [2.0, 1.8, 1.6, -1.2, 0.9]
- samples: 12

Use Stouffer's weighted meta-analysis with directionality and FDR < 0.05.
Which genes are significantly dysregulated across all modalities?
```

---

## Configuration

### Environment Variables

| Variable | Purpose | Default | Production |
|----------|---------|---------|------------|
| `MULTIOMICS_DATA_DIR` | Data files | `/workspace/data/multiomics` | Custom path |
| `MULTIOMICS_CACHE_DIR` | Cached results | `/workspace/cache/multiomics` | Custom path |
| `MULTIOMICS_DRY_RUN` | Mock mode | `true` | `false` |

### Claude Desktop Setup

```json
"multiomics": {
  "command": "/path/to/precision-medicine-mcp/servers/mcp-multiomics/venv/bin/python",
  "args": ["-m", "mcp_multiomics"],
  "env": {
    "MULTIOMICS_DRY_RUN": "false"
  }
}
```

**Important:** Set `MULTIOMICS_DRY_RUN="false"` to enable directionality features.

---

## Scientific Background

### Stouffer's Method

**Purpose:** Combine p-values from independent tests while preserving effect directionality.

**Steps:**
1. Convert p-values to Z-scores (inverse normal transformation)
2. Apply effect size directionality (signs from log2FC)
3. Combine Z-scores (weighted or unweighted)
4. Convert combined Z to p-value
5. Apply FDR correction

**Advantages:**
- Preserves effect directionality (upregulated vs downregulated)
- Handles concordant and discordant effects
- More powerful than Fisher's method for large effects
- Supports sample-size weighting

**Reference:** Whitlock MC (2005). *J Evol Biol* 18(5):1368-73.

---

## Getting Started

### For Users

1. **Read:** `MULTIOMICS_INTEGRATION_SUMMARY.md` - Overview and deployment
2. **Configure:** Set `MULTIOMICS_DRY_RUN="false"` in Claude Desktop config
3. **Test:** Use example prompts above in Claude Desktop

### For Developers

1. **Read:** `DETAILED_GUIDE.md` - Complete technical reference (584 lines)
2. **Design:** `Implementation_Plan_mcp_multiomics.md` - Original specification
3. **Testing:** Run `pytest tests/` in server directory

### For Troubleshooting

1. **Common Issues:** `PYDANTIC_BOOLEAN_FIX.md` - Boolean parsing, env vars
2. **Test Results:** Check `../../manual_testing/TEST_RESULTS_ALL_SERVERS.md`

---

## Key Features

### Statistical Rigor ✅
- Proper inverse normal transformation
- Effect size directionality preservation
- Benjamini-Hochberg FDR correction
- Sample-size weighted combinations

### Data Handling ✅
- Multi-modality support (RNA, Protein, Phospho)
- Automatic sample alignment
- Missing data filtering
- Z-score normalization
- Integrated data caching

### Production Ready ✅
- 29/29 tests passing (100%)
- 84% code coverage
- Full type hints (mypy clean)
- DRY_RUN mode for testing
- Comprehensive error handling

---

## Related Documentation

### Server Documentation
- **[Server README](../../servers/mcp-multiomics/README.md)** - Server-specific documentation
- **[Test Results](../../manual_testing/Solution-Testing/TEST_RESULTS_ALL_SERVERS.md)** - All server test results

### Project Documentation
- **[Main README](../../README.md)** - Project overview
- **[Spatial Docs](../spatial/README.md)** - Phase 1-2 spatial transcriptomics docs
- **[Configuration](../../configs/README.md)** - Claude Desktop setup

### Technical Resources
- **[Architecture](../../architecture/spatial/README.md)** - Spatial MCP architecture
- **[Manual Testing](../../manual_testing/README.md)** - Testing and verification

---

## Status

| Metric | Value |
|--------|-------|
| **Version** | 1.0.0 |
| **Status** | ✅ Production Ready |
| **Tests** | 29/29 passing |
| **Coverage** | 84% |
| **Claude Desktop** | ✅ Configured |
| **Last Updated** | November 11, 2025 |

---

## Quick Links

📖 **[DETAILED_GUIDE.md](DETAILED_GUIDE.md)** - Complete 584-line technical reference
🔧 **[Troubleshooting](../../servers/mcp-multiomics/TROUBLESHOOTING.md)** - Technical fixes and debugging
📋 **[Integration Summary](MULTIOMICS_INTEGRATION_SUMMARY.md)** - Original integration summary (archived)
📐 **[Implementation Plan](archive/Implementation_Plan_mcp_multiomics.md)** - Original design document (archived)

---

**Created:** November 11, 2025
**Purpose:** Multi-omics PDX analysis and Stouffer's meta-analysis
**Phase:** 3 (Advanced Analysis)
**Server:** mcp-multiomics
