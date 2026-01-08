# mcp-openimagedata: Histology Image Processing

MCP server for histology and immunofluorescence image processing, spatial registration, and visualization generation for microscopy analysis.

## Overview

`mcp-openimagedata` provides AI-accessible histology image processing tools through the Model Context Protocol (MCP). This server enables automated image registration, quality assessment, multiplex channel compositing, and morphology annotation for H&E and immunofluorescence microscopy.

### Key Features

- 🖼️ **Image Registration** - Align H&E with IF images for spatial correlation
- ✅ **Quality Assessment** - Evaluate image quality and detect artifacts
- 🎨 **Multiplex Compositing** - Combine 2-7 fluorescence channels into RGB composites
- 📍 **Morphology Annotation** - Annotate necrotic regions and high cellularity areas on H&E
- 🔬 **Format Conversion** - Convert between imaging formats
- ⚡ **DRY_RUN Mode** - Test workflows with synthetic data

## Installation

### Prerequisites

- Python 3.11 or higher
- pip package manager

### Local Setup

1. **Create a virtual environment:**

```bash
cd servers/mcp-openimagedata
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
# Data directories
IMAGE_DATA_DIR=/workspace/data/images
IMAGE_CACHE_DIR=/workspace/cache/images

# Output directory for visualizations
IMAGE_OUTPUT_DIR=/workspace/output

# Development mode (uses mocks instead of real processing)
IMAGE_DRY_RUN=true
```

## Usage

### Running the Server

**Standalone mode (stdio):**

```bash
python -m mcp_openimagedata
```

**With Claude Desktop:**

Add to your Claude Desktop configuration (`~/Library/Application Support/Claude/claude_desktop_config.json` on macOS):

```json
{
  "mcpServers": {
    "openimagedata": {
      "command": "/path/to/spatial-mcp/servers/mcp-openimagedata/venv/bin/python",
      "args": ["-m", "mcp_openimagedata"],
      "cwd": "/path/to/spatial-mcp/servers/mcp-openimagedata",
      "env": {
        "PYTHONPATH": "/path/to/spatial-mcp/servers/mcp-openimagedata/src",
        "IMAGE_DATA_DIR": "/workspace/data/images",
        "IMAGE_OUTPUT_DIR": "/workspace/output",
        "IMAGE_DRY_RUN": "false"
      }
    }
  }
}
```

**Important:**
- Use the full path to the venv Python executable
- Set `IMAGE_DRY_RUN=false` for real image processing
- Ensure sufficient disk space for large TIFF files

## Available Tools

### 1. register_images

Spatially register H&E and IF images for correlation analysis.

**Parameters:**
- `fixed_image_path` (string): Path to reference image (e.g., H&E)
- `moving_image_path` (string): Path to image to be aligned (e.g., IF)
- `registration_method` (string, optional): Method - "affine", "rigid", or "deformable" (default: "affine")
- `output_path` (string, optional): Path for registered output

**Returns:**
```json
{
  "registered_image": "/output/IF_registered.tiff",
  "transformation_matrix": [[1.02, 0.01, 15], [0.01, 1.01, -8], [0, 0, 1]],
  "registration_quality": {
    "correlation": 0.87,
    "mutual_information": 1.24
  },
  "method_used": "affine"
}
```

**Example usage with Claude:**
```
Register the Ki67 IF image to the H&E reference:
- Fixed (H&E): /data/HE_slide.tif
- Moving (IF): /data/IF_Ki67.tiff
- Use affine registration
```

**Output:** Registered image aligned to reference coordinate system.

---

### 2. assess_image_quality

Evaluate image quality and detect common artifacts.

**Parameters:**
- `image_path` (string): Path to microscopy image
- `modality` (string, optional): Imaging modality - "brightfield" or "fluorescence" (default: "fluorescence")

**Returns:**
```json
{
  "quality_score": 0.85,
  "metrics": {
    "sharpness": 0.89,
    "signal_to_noise": 23.4,
    "dynamic_range": 0.78
  },
  "artifacts_detected": ["edge_effects", "out_of_focus_regions"],
  "recommendation": "Good quality - minor edge effects in upper-left corner"
}
```

**Example usage with Claude:**
```
Assess quality of the CD8 IF image:
- Image: /data/IF_CD8.tiff
- Modality: fluorescence
```

**Detectable artifacts:** out_of_focus, photobleaching, saturation, edge_effects, tissue_folding

---

### 3. convert_format

Convert between imaging file formats.

**Parameters:**
- `input_path` (string): Path to input image
- `output_format` (string): Target format - "tiff", "png", "ome-tiff"
- `compression` (string, optional): Compression - "none", "lzw", "jpeg" (default: "lzw")

**Returns:**
```json
{
  "output_file": "/output/converted_image.tiff",
  "input_format": "png",
  "output_format": "tiff",
  "compression": "lzw",
  "file_size_mb": 45.2
}
```

**Example usage with Claude:**
```
Convert H&E slide from PNG to OME-TIFF:
- Input: /data/HE_slide.png
- Format: ome-tiff
- Compression: lzw
```

---

## Visualization Tools

The following visualization tools generate publication-quality PNG images:

### 4. generate_multiplex_composite

Generate RGB composite from multiplex immunofluorescence channels.

**Parameters:**
- `channel_paths` (list): List of paths to channel images (1-7 channels, grayscale TIFF/PNG)
- `channel_names` (list): List of marker names matching channel order (e.g., ["DAPI", "Ki67", "TP53"])
- `channel_colors` (list, optional): List of colors for each channel - "blue", "green", "red", "cyan", "magenta", "yellow", "white" (default: ["blue", "green", "red", ...])
- `output_filename` (string, optional): Custom output filename (default: auto-generated with timestamp)
- `normalize` (boolean, optional): Normalize each channel to 0-1 range (default: True)

**Returns:**
```json
{
  "status": "success",
  "output_file": "/workspace/output/visualizations/multiplex_composite_20250108_144200.png",
  "channels_combined": 3,
  "channel_info": [
    {"name": "DAPI", "color": "blue"},
    {"name": "Ki67", "color": "green"},
    {"name": "TP53", "color": "red"}
  ],
  "image_dimensions": {"width": 2048, "height": 2048},
  "description": "RGB composite combining 3 channels: DAPI (blue), Ki67 (green), TP53 (red).",
  "visualization_type": "multiplex_composite"
}
```

**Example usage with Claude:**
```
Create RGB composite from multiplex IF channels:
- Channels:
  * /data/Multiplex_DAPI.tiff (nuclear stain)
  * /data/Multiplex_Ki67.tiff (proliferation marker)
  * /data/Multiplex_TP53.tiff (tumor suppressor)
- Channel colors: blue, green, red
- Normalize intensities for each channel
```

**Output:** Multi-panel figure showing:
1. Individual grayscale channels (one per marker)
2. Final RGB composite with all channels combined

**Use cases:**
- Visualize co-localization of multiple markers
- Identify double-positive cells (e.g., TP53+/Ki67+ = yellow in red+green composite)
- Generate publication-quality multiplex IF figures
- Validate channel alignment and crosstalk

**Color combinations:**
- **2 channels:** Blue + Green, Green + Red, Blue + Red
- **3 channels:** Blue (DAPI) + Green (marker 1) + Red (marker 2)
- **4+ channels:** Use cyan, magenta, yellow for additional markers

**Technical details:**
- Each channel is normalized to 0-1 range (if `normalize=True`)
- RGB composite is created by weighted addition: `RGB = Σ(channel_i × color_i)`
- Output is clipped to [0, 255] and saved as 8-bit RGB PNG
- Supports up to 7 channels (though 2-4 is typical)

---

### 5. generate_he_annotation

Generate annotated H&E morphology visualization with region overlays.

**Parameters:**
- `he_image_path` (string): Path to H&E histology image
- `necrotic_regions` (list, optional): List of bounding boxes for necrotic regions (default: [])
  - Each region: `{"x": int, "y": int, "width": int, "height": int}`
- `high_cellularity_regions` (list, optional): List of bounding boxes for high cellularity regions (default: [])
  - Each region: `{"x": int, "y": int, "width": int, "height": int}`
- `output_filename` (string, optional): Custom output filename (default: auto-generated)
- `necrotic_color` (string, optional): Color for necrotic region boxes (default: "red")
- `cellularity_color` (string, optional): Color for cellularity region boxes (default: "green")

**Returns:**
```json
{
  "status": "success",
  "output_file": "/workspace/output/visualizations/he_annotation_20250108_144230.png",
  "necrotic_regions_count": 2,
  "cellularity_regions_count": 3,
  "description": "H&E annotation showing 2 necrotic regions and 3 high cellularity regions.",
  "visualization_type": "he_morphology_annotation"
}
```

**Example usage with Claude:**
```
Annotate H&E slide with pathology findings:
- Image: /data/HE_slide.tif
- Necrotic regions (red dashed boxes):
  * {"x": 500, "y": 300, "width": 200, "height": 150}
  * {"x": 1200, "y": 800, "width": 180, "height": 140}
- High cellularity regions (green solid boxes):
  * {"x": 200, "y": 200, "width": 250, "height": 200}
  * {"x": 900, "y": 400, "width": 300, "height": 250}
  * {"x": 1500, "y": 600, "width": 220, "height": 180}
```

**Output:** Side-by-side figure showing:
1. Original H&E image (unmodified)
2. Annotated H&E with colored bounding boxes and legend

**Use cases:**
- Visualize pathologist annotations on H&E
- Highlight regions of interest for downstream analysis
- Document necrotic areas and tumor cellularity
- Generate annotated figures for pathology reports

**Region types:**
- **Necrotic regions** (default: red dashed boxes):
  - Cell death zones
  - Hypoxic areas
  - Therapy-induced necrosis

- **High cellularity regions** (default: green solid boxes):
  - Dense tumor areas
  - Lymphocyte infiltration
  - Proliferative zones

**Technical details:**
- Bounding boxes drawn as matplotlib rectangles
- Dashed lines for necrotic (pathology), solid for cellularity (biology)
- Legend shows counts for each region type
- Coordinates are in pixel units (x, y from top-left origin)

---

## Example Workflows

### Complete Multiplex IF Analysis Workflow

**End-to-end multiplex imaging analysis:**

```
Claude, analyze the multiplex IF slide (DAPI + Ki67 + TP53):

STEP 1: Quality assessment
- Check all 3 channels for artifacts
- Verify signal-to-noise ratio

STEP 2: Generate RGB composite
- Combine channels: DAPI (blue) + Ki67 (green) + TP53 (red)
- Normalize intensities
- Look for yellow cells (Ki67+/TP53+ double-positive)

STEP 3: Cell segmentation
- Use DAPI channel for nuclear segmentation (deepcell server)
- Generate segmentation mask

STEP 4: Measure marker intensities
- For each segmented cell:
  - Measure mean Ki67 intensity
  - Measure mean TP53 intensity

STEP 5: Phenotype classification and visualization
- Classify cells: TP53+/Ki67+, TP53+/Ki67-, TP53-/Ki67+, TP53-/Ki67-
- Generate phenotype visualization (deepcell server)

STEP 6: Report results
- Distribution of 4 phenotypes
- Spatial patterns (clustered vs. dispersed)
- Clinical interpretation
```

### H&E Registration and Annotation Workflow

**Correlate H&E morphology with IF markers:**

```
Claude, correlate H&E and CD8 IF for Patient-001:

STEP 1: Register images
- Fixed: /data/HE_slide.tif (reference)
- Moving: /data/IF_CD8.tiff
- Use affine registration

STEP 2: Assess H&E quality
- Check for tissue folding or artifacts

STEP 3: Identify regions on H&E
- Necrotic regions: Low cellularity, eosinophilic debris
- High cellularity regions: Dense tumor nests

STEP 4: Annotate H&E
- Draw bounding boxes around:
  * Necrotic areas (red dashed)
  * High cellularity areas (green solid)

STEP 5: Correlate with CD8 IF
- Overlay CD8+ cell locations on annotated H&E
- Quantify CD8+ density in each region type

STEP 6: Report findings
- CD8+ infiltration in high vs. low cellularity regions
- Necrotic areas exclude immune cells
- Spatial correlation between morphology and immune landscape
```

### Multi-Channel Composite Workflow

**Generate composite for 5-plex IF:**

```
Claude, create composite from 5-plex IF panel:

Channels:
1. DAPI (nuclear) → blue
2. CD3 (T-cells) → green
3. CD8 (cytotoxic T-cells) → red
4. PD-L1 (checkpoint) → cyan
5. Ki67 (proliferation) → magenta

Generate multiplex composite:
- Normalize all channels
- Save individual channels and composite
- Look for:
  * CD3+/CD8+ cells (green + red = yellow)
  * PD-L1+/Ki67+ cells (cyan + magenta)
  * Spatial distribution patterns

Expected output:
- 6-panel figure (5 individual + 1 composite)
- Multi-colored cells indicating co-expression
```

## Available Resources

### image://formats

Get information about supported imaging formats and requirements.

**Example:**
```
What image formats are supported?
```

### image://registration-methods

Get details about available registration algorithms.

**Example:**
```
What registration methods are available?
```

## Data Format Requirements

### Input Images

**Fluorescence (IF):**
- **Formats:** TIFF (16-bit preferred), PNG (8-bit)
- **Channels:** Single-channel grayscale per marker
- **Size:** Up to 4096 × 4096 pixels per channel
- **Bit depth:** 8-bit or 16-bit

**Brightfield (H&E):**
- **Formats:** TIFF, PNG, JPEG
- **Channels:** RGB color image
- **Size:** Variable (whole slide images supported via tiling)
- **Resolution:** Typically 0.25-0.5 microns/pixel

### Bounding Box Format

For region annotations:
```json
{
  "x": 500,        // X coordinate of top-left corner (pixels)
  "y": 300,        // Y coordinate of top-left corner (pixels)
  "width": 200,    // Width in pixels
  "height": 150    // Height in pixels
}
```

Origin is top-left corner (0, 0). Positive X goes right, positive Y goes down.

## Development

### Running Tests

**Run all tests:**
```bash
pytest
```

**Run with coverage:**
```bash
pytest --cov=src/mcp_openimagedata --cov-report=html
```

**Run specific test:**
```bash
pytest tests/test_registration.py -v
```

### Code Quality

**Format code:**
```bash
black src/ tests/
```

**Lint code:**
```bash
ruff check src/ tests/
```

**Type checking:**
```bash
mypy src/
```

## Architecture

### Implementation Status: 30% Real (Mocked for Demonstration)

- **register_images**: 30% real (basic registration implemented, uses simplified algorithm)
- **assess_image_quality**: 30% real (basic metrics implemented)
- **convert_format**: 30% real (PIL-based conversion)
- **generate_multiplex_composite**: 100% real (full implementation with color blending)
- **generate_he_annotation**: 100% real (matplotlib bounding box overlay)

**Note:** Full image registration requires OpenCV or ITK. Current implementation uses simplified methods for demonstration. Visualization tools are fully functional.

### Design Principles

1. **Visual Validation** - Always generate visualizations for quality control
2. **Flexible Compositing** - Support 1-7 channels with customizable colors
3. **Publication Quality** - 300 DPI PNG outputs suitable for papers
4. **Modular Design** - Tools work independently or as part of workflows

### Directory Structure

```
mcp-openimagedata/
├── src/
│   └── mcp_openimagedata/
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
| `IMAGE_DATA_DIR` | `/workspace/data/images` | Directory for image datasets |
| `IMAGE_CACHE_DIR` | `/workspace/cache/images` | Directory for cached files |
| `IMAGE_OUTPUT_DIR` | `/workspace/output` | Directory for output files |
| `IMAGE_DRY_RUN` | `false` | Enable mock mode (no real processing) |
| `IMAGE_LOG_LEVEL` | `INFO` | Logging level |

## Troubleshooting

### Server won't start

- **Check Python version:** Must be 3.11+
- **Verify dependencies:** Run `pip install -e .`
- **Check environment variables:** Ensure directories exist or are writable

### Image registration fails

- **Enable dry-run mode:** Set `IMAGE_DRY_RUN=true` for testing
- **Check image formats:** Ensure both images are readable
- **Verify dimensions:** Images should be roughly same size
- **Review logs:** Check stderr output for details

### Multiplex composite looks wrong

- **Check channel order:** Ensure channels match marker names
- **Verify normalization:** Try `normalize=true` if intensities differ greatly
- **Check channel alignment:** Ensure channels are spatially aligned
- **Adjust colors:** Try different color combinations for clarity

### H&E annotation not visible

- **Check bounding box coordinates:** Ensure boxes are within image bounds
- **Verify box sizes:** Very small boxes may not be visible
- **Check colors:** Use contrasting colors (red, green work well on H&E)
- **Increase linewidth:** Modify code to use thicker lines if needed

## License

See the main repository LICENSE file.

## Support

- **Issues:** https://github.com/lynnlangit/precision-medicine-mcp/issues
- **Documentation:** See main repository docs/
- **MCP Specification:** https://modelcontextprotocol.io/

## Related Servers

Part of the Precision Medicine MCP suite:

- **mcp-deepcell** - Cell segmentation and phenotype visualization
- **mcp-spatialtools** - Spatial transcriptomics analysis
- **mcp-epic** - Clinical EHR data (FHIR)
- **mcp-multiomics** - Multi-omics integration
- **mcp-fgbio** - Genomic reference data

---

**Built for the Precision Medicine MCP suite** - Enables AI-driven histology image processing integrated with spatial, cellular, and clinical data.
