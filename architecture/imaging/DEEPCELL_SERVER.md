# mcp-deepcell: Deep Learning Cell Segmentation for MxIF

MCP server for deep learning-based cell segmentation in multiplexed immunofluorescence (MxIF) and other fluorescence microscopy images, with visualization tools for validating segmentation quality and phenotype analysis.

## Overview

`mcp-deepcell` provides AI-accessible cell segmentation and classification tools through the Model Context Protocol (MCP). This server uses the **open-source DeepCell-TF library** (https://github.com/vanvalenlab/deepcell-tf) to enable automated cell detection, phenotype classification, and visualization generation for MxIF workflows.

### Primary Use Case: MxIF Analysis

**MxIF (Multiplexed Immunofluorescence)** enables imaging of multiple protein markers (2-7+) on a single tissue section. DeepCell segments cells from fluorescence channels (e.g., DAPI nuclear stain) and enables quantitative phenotyping:
- Single-cell marker quantification (CD8+, Ki67+, TP53+, etc.)
- Multi-marker co-expression analysis (e.g., TP53+/Ki67+ double-positive cells)
- Spatial distribution of cell phenotypes
- Cell counting and density measurements

### Key Features

- 🔬 **MxIF Cell Segmentation** - Deep learning-based nuclear, membrane, and cytoplasm segmentation for fluorescence images
- 🎨 **Segmentation Overlays** - Visualize cell boundaries overlaid on original fluorescence images
- 📊 **Phenotype Visualization** - Color cells by marker expression (positive/negative, co-expression patterns)
- 🧬 **Cell Classification** - Classify cell states based on fluorescence intensity (proliferating, quiescent, apoptotic, etc.)
- ⚡ **DRY_RUN Mode** - Test workflows with synthetic data

## Installation

### Prerequisites

- Python 3.11 or higher
- pip package manager

### Local Setup

1. **Create a virtual environment:**

```bash
cd servers/mcp-deepcell
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

2. **Install dependencies:**

```bash
pip install -e ".[dev]"
```

3. **Set up environment variables:**

Create a `.env` file in the server directory:

```bash
# Output directory for visualizations
DEEPCELL_OUTPUT_DIR=/workspace/output

# Development mode (uses mocks instead of real segmentation)
DEEPCELL_DRY_RUN=true
```

## Usage

### Running the Server

**Standalone mode (stdio):**

```bash
python -m mcp_deepcell
```

**With Claude Desktop:**

Add to your Claude Desktop configuration (`~/Library/Application Support/Claude/claude_desktop_config.json` on macOS):

```json
{
  "mcpServers": {
    "deepcell": {
      "command": "/path/to/spatial-mcp/servers/mcp-deepcell/venv/bin/python",
      "args": ["-m", "mcp_deepcell"],
      "cwd": "/path/to/spatial-mcp/servers/mcp-deepcell",
      "env": {
        "PYTHONPATH": "/path/to/spatial-mcp/servers/mcp-deepcell/src",
        "DEEPCELL_OUTPUT_DIR": "/workspace/output",
        "DEEPCELL_DRY_RUN": "false"
      }
    }
  }
}
```

**Important:**
- Use the full path to the venv Python executable
- Set `DEEPCELL_DRY_RUN=false` for real segmentation
- Requires deep learning model weights (loaded automatically)

## Available Tools

### 1. segment_cells

Deep learning-based cell segmentation for microscopy images.

**Parameters:**
- `image_path` (string): Path to microscopy image (TIFF, PNG, etc.)
- `model_type` (string, optional): Model to use - "membrane", "nuclear", or "cytoplasm" (default: "membrane")
- `min_cell_size` (integer, optional): Minimum cell size in pixels (default: 100)

**Returns:**
```json
{
  "segmentation_mask": "/path/to/image.seg.tif",
  "cells_detected": 1247,
  "mean_cell_area": 850,
  "model_used": "membrane",
  "processing_time_seconds": 15.3,
  "quality_metrics": {
    "mean_confidence": 0.89,
    "cells_filtered": 23
  }
}
```

**Example usage with Claude:**
```
Segment cells in the Ki67 IF image:
- Image: /data/IF_Ki67.tiff
- Use nuclear segmentation model
- Filter out cells smaller than 50 pixels
```

**Output:** Label image where each cell has a unique ID (0 = background, 1 = cell 1, 2 = cell 2, etc.)

---

### 2. classify_cell_states

Classify cell states/phenotypes based on segmentation and expression data.

**Parameters:**
- `segmentation_mask` (string): Path to segmentation mask (output from segment_cells)
- `expression_data` (string): Path to expression matrix or intensity measurements

**Returns:**
```json
{
  "classifications": [
    {"cell_id": 1, "state": "proliferating", "confidence": 0.92},
    {"cell_id": 2, "state": "quiescent", "confidence": 0.87},
    {"cell_id": 3, "state": "apoptotic", "confidence": 0.95}
  ],
  "total_cells": 1247,
  "states_identified": ["proliferating", "quiescent", "apoptotic", "senescent"]
}
```

**Example usage with Claude:**
```
Classify cell states based on Ki67/TP53 expression:
- Segmentation: /data/IF_Ki67.tiff.seg.tif
- Expression: /data/cell_intensities.csv
```

**Available states:** proliferating, quiescent, apoptotic, senescent, necrotic

---

## Visualization Tools

The following visualization tools generate publication-quality PNG images:

### 3. generate_segmentation_overlay

Generate segmentation overlay visualization showing cell boundaries on original image.

**Parameters:**
- `original_image_path` (string): Path to original microscopy image
- `segmentation_mask_path` (string): Path to segmentation mask (label image)
- `output_filename` (string, optional): Custom output filename (default: auto-generated with timestamp)
- `overlay_color` (string, optional): Color for cell boundaries - "yellow", "green", "red", "cyan", "magenta" (default: "yellow")
- `overlay_alpha` (float, optional): Transparency of overlay (0=transparent, 1=opaque, default: 0.4)

**Returns:**
```json
{
  "status": "success",
  "output_file": "/workspace/output/visualizations/segmentation_overlay_20250108_143500.png",
  "cells_visualized": 1247,
  "description": "Segmentation overlay showing 1247 segmented cells with yellow boundaries overlaid on original image.",
  "visualization_type": "segmentation_overlay",
  "overlay_color": "yellow"
}
```

**Example usage with Claude:**
```
Generate segmentation overlay for CD8 IF image:
- Original: /data/IF_CD8.tiff
- Segmentation: /data/IF_CD8.tiff.seg.tif
- Use green boundaries with 50% transparency
```

**Output:** Side-by-side figure showing:
1. Original microscopy image
2. Image with colored cell boundaries overlaid

**Use cases:**
- Validate segmentation quality
- Identify over-segmentation or under-segmentation
- Visualize cell density and morphology

---

### 4. generate_phenotype_visualization

Generate cell phenotype visualization coloring cells by marker expression.

**Parameters:**
- `original_image_path` (string): Path to original microscopy image
- `segmentation_mask_path` (string): Path to segmentation mask
- `marker_positive_cells` (list): List of cell IDs that are marker-positive (e.g., [1, 5, 12, 18, ...])
- `output_filename` (string, optional): Custom output filename (default: auto-generated)
- `positive_color` (string, optional): Color for marker-positive cells (default: "green")
- `negative_color` (string, optional): Color for marker-negative cells (default: "red")

**Returns:**
```json
{
  "status": "success",
  "output_file": "/workspace/output/visualizations/phenotype_viz_20250108_143530.png",
  "positive_cells": 623,
  "negative_cells": 624,
  "total_cells": 1247,
  "percent_positive": 49.96,
  "description": "Phenotype visualization showing 623 marker-positive cells (green) and 624 marker-negative cells (red). 49.96% of cells are positive.",
  "visualization_type": "phenotype_coloring",
  "positive_color": "green",
  "negative_color": "red"
}
```

**Example usage with Claude:**
```
Visualize Ki67+ proliferating cells:
- Original: /data/IF_Ki67.tiff
- Segmentation: /data/IF_Ki67.tiff.seg.tif
- Positive cells: [1, 5, 12, 18, 25, ...] (from classification)
- Color positive cells yellow, negative cells blue
```

**Output:** Side-by-side figure showing:
1. Original microscopy image (grayscale)
2. Cells colored by phenotype with legend showing counts

**Use cases:**
- Visualize marker-positive cell distribution (CD8+, Ki67+, etc.)
- Quantify spatial patterns of phenotypes
- Generate figures for publications
- Validate classification results

---

## Example Workflows

### Complete Cell Segmentation and Visualization Workflow

**End-to-end analysis:**

```
Claude, analyze the CD8 immunofluorescence image:

STEP 1: Segment cells
- Image: /data/IF_CD8.tiff
- Use membrane segmentation model
- Filter cells < 100 pixels

STEP 2: Generate segmentation overlay
- Show green boundaries on original image
- Validate segmentation quality

STEP 3: Classify CD8+ cells
- Use intensity threshold to identify positive cells
- Cells with mean intensity > 50 are CD8+

STEP 4: Generate phenotype visualization
- Color CD8+ cells green, CD8- cells red
- Calculate percent positive

STEP 5: Report results
- Total cells detected
- CD8+ count and percentage
- Spatial distribution patterns
```

### Multiplex IF Phenotype Analysis

**Multi-marker phenotyping:**

```
Claude, analyze multiplex IF (DAPI + Ki67 + TP53):

STEP 1: Segment nuclei
- Image: /data/Multiplex_DAPI.tiff (nuclear channel)
- Use nuclear segmentation model

STEP 2: Measure marker intensities
- For each segmented cell:
  - Measure mean Ki67 intensity in nucleus
  - Measure mean TP53 intensity in nucleus

STEP 3: Classify phenotypes
- TP53+/Ki67+: Proliferating mutant cells
- TP53+/Ki67-: Quiescent mutant cells
- TP53-/Ki67+: Proliferating wild-type cells
- TP53-/Ki67-: Quiescent wild-type cells

STEP 4: Visualize each phenotype
- Generate 4 separate phenotype visualizations
- Or use multi-color overlay

STEP 5: Quantify percentages
- Report distribution of 4 phenotypes
```

## Available Resources

### model://deepcell/membrane

Get information about the membrane segmentation model.

**Example:**
```
What model is used for membrane segmentation?
```

Returns model architecture, training data, accuracy, and use cases.

## Data Format Requirements

### Input Images

- **Formats:** TIFF, PNG, JPEG
- **Bit depth:** 8-bit or 16-bit grayscale, or RGB
- **Size:** Up to 4096 × 4096 pixels (larger images may be tiled)
- **Channels:** Single-channel (grayscale) or multi-channel (RGB)

### Segmentation Masks

- **Format:** TIFF or PNG
- **Type:** Label image (integer values)
- **Values:**
  - 0 = background
  - 1 = cell 1
  - 2 = cell 2
  - ...
  - N = cell N

## Development

### Running Tests

**Run all tests:**
```bash
pytest
```

**Run with coverage:**
```bash
pytest --cov=src/mcp_deepcell --cov-report=html
```

**Run specific test:**
```bash
pytest tests/test_segmentation.py -v
```

### Code Quality

**Format code:**
```bash
black src/ tests/
```

## Architecture

### Implementation Status: 30% Real (Mocked for Demonstration)

- **segment_cells**: 30% real (basic segmentation implemented, uses simplified model)
- **classify_cell_states**: 30% real (intensity-based classification)
- **generate_segmentation_overlay**: 100% real (scikit-image boundary detection)
- **generate_phenotype_visualization**: 100% real (matplotlib coloring)

**Note:** Full DeepCell integration requires TensorFlow/PyTorch models. Current implementation uses simplified segmentation for demonstration. Visualization tools are fully functional.

### Design Principles

1. **Visual Validation** - Always generate visualizations to validate results
2. **Flexible Coloring** - Customizable colors for different use cases
3. **Publication Quality** - 300 DPI PNG outputs suitable for papers
4. **Error Handling** - Graceful degradation with informative error messages

### Directory Structure

```
mcp-deepcell/
├── src/
│   └── mcp_deepcell/
│       ├── __init__.py
│       ├── __main__.py
│       └── server.py          # Main server implementation
├── tests/
│   ├── conftest.py            # Pytest fixtures
│   └── test_server.py         # Tool tests
├── pyproject.toml             # Project configuration
└── README.md
```

## Configuration Reference

| Environment Variable | Default | Description |
|---------------------|---------|-------------|
| `DEEPCELL_OUTPUT_DIR` | `/workspace/output` | Directory for output files |
| `DEEPCELL_DRY_RUN` | `false` | Enable mock mode (no real segmentation) |
| `DEEPCELL_LOG_LEVEL` | `INFO` | Logging level |

## Troubleshooting

### Server won't start

- **Check Python version:** Must be 3.11+
- **Verify dependencies:** Run `pip install -e .`
- **Check environment variables:** Ensure output directory exists or is writable

### Segmentation fails

- **Enable dry-run mode:** Set `DEEPCELL_DRY_RUN=true` for testing
- **Check image format:** Ensure image is readable by PIL
- **Check file paths:** Ensure absolute paths are used
- **Review logs:** Check stderr output for details

### Visualization errors

- **Check segmentation mask:** Must be label image (integer values)
- **Verify image dimensions:** Original and mask must have same dimensions
- **Check output directory:** Must be writable
- **Check color names:** Use valid color names (yellow, green, red, etc.)

## License

See the main repository LICENSE file.

## Support

- **Issues:** https://github.com/lynnlangit/precision-medicine-mcp/issues
- **Documentation:** See main repository docs/
- **MCP Specification:** https://modelcontextprotocol.io/

## Related Servers

Part of the Precision Medicine MCP suite:

- **mcp-spatialtools** - Spatial transcriptomics analysis and visualization
- **mcp-openimagedata** - Histology image processing and multiplex IF
- **mcp-epic** - Clinical EHR data (FHIR)
- **mcp-multiomics** - Multi-omics integration
- **mcp-fgbio** - Genomic reference data

---

**Built for the Precision Medicine MCP suite** - Enables AI-driven cell segmentation and phenotype analysis integrated with spatial and clinical data.
