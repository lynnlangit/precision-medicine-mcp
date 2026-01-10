# MxIF (Multiplexed Immunofluorescence) Workflow

**Status:** ⚠️ Partially Implemented (openimagedata: 60% real, deepcell: 0% mocked)
**Last Updated:** January 10, 2026
**Used In:** PatientOne TEST_4_IMAGING

---

## Overview

Multiplexed immunofluorescence (MxIF) workflow for quantitative cell phenotyping using **fluorescence microscopy** with multiple antibody markers.

**Pipeline:** mcp-openimagedata (load/composite) → mcp-deepcell (segment/quantify)

**Purpose:** Quantitative single-cell analysis with spatial context

---

## What is MxIF?

**Multiplexed Immunofluorescence:**
- Images **multiple protein markers** (2-7+) on a single tissue section
- Uses **fluorescent antibodies** (not chromogenic stains)
- Repeated rounds of: stain → image → inactivate dye → restain
- Preserves **spatial relationships** across all markers
- Enables **multi-marker co-expression** analysis (e.g., TP53+/Ki67+ double-positive cells)

**Key difference from H&E:**
- MxIF = Fluorescence microscopy → quantitative protein expression
- H&E = Brightfield microscopy → morphology visualization

**Technology:** Uses DeepCell-TF library for AI-based cell segmentation
**Reference:** https://github.com/vanvalenlab/deepcell-tf

---

## Data Formats

### Single-Channel IF

**File:** Grayscale TIFF (e.g., `PAT001_tumor_IF_CD8.tiff`)

**Specifications:**
- **Format:** Grayscale TIFF (1 channel)
- **Microscopy:** Fluorescence (single marker)
- **Marker:** One protein (CD8, Ki67, CD3, etc.)
- **Dimensions:** 512×512 pixels (demo), 1000-2000 pixels (production)
- **File size:** ~268 KB per channel

**Examples:**
- `PAT001_tumor_IF_CD8.tiff` - CD8+ T cells
- `PAT001_tumor_IF_KI67.tiff` - Ki67+ proliferating cells

### Multiplex IF (MxIF)

**File:** RGB TIFF with 3 channels (e.g., `PAT001_tumor_multiplex_IF_TP53_KI67_DAPI.tiff`)

**Specifications:**
- **Format:** RGB TIFF (3 channels)
- **Microscopy:** Fluorescence (multiplex)
- **Markers:** 3 proteins in separate channels
  - Channel 1 (Red): TP53
  - Channel 2 (Green): Ki67
  - Channel 3 (Blue): DAPI (nuclear stain)
- **Dimensions:** 512×512 pixels
- **File size:** ~63 KB (3 channels)

---

## Workflow Steps

### Step 1: Load IF Image
**Tool:** `fetch_histology_image` (mcp-openimagedata)

Load single-channel or multiplex fluorescence TIFF image.

**Parameters:**
- `image_path`: Path to IF TIFF file
- `image_type`: "IF" or "MXIF"

**Output:**
```json
{
  "status": "success",
  "image_path": "/path/to/PAT001_tumor_IF_CD8.tiff",
  "format": "Grayscale TIFF",
  "dimensions": [512, 512],
  "channels": 1,
  "marker": "CD8"
}
```

### Step 2: Generate Composite (Multiplex Only)
**Tool:** `generate_multiplex_composite` (mcp-openimagedata)

For multiplex images, create RGB composite visualization.

**Parameters:**
- `image_path`: Path to multiplex IF TIFF
- `channels`: List of channel names (e.g., ["TP53", "Ki67", "DAPI"])
- `colors`: RGB color mapping (e.g., ["red", "green", "blue"])
- `output_filename`: Custom filename (optional)

**Output:** RGB PNG composite showing all markers overlaid

**File format:** `multiplex_composite_{markers}_{timestamp}.png`

**Example Output:**
- Red: TP53 expression
- Green: Ki67 expression
- Blue: DAPI nuclei
- Yellow: TP53+/Ki67+ double-positive cells (co-expression)

### Step 3: Cell Segmentation
**Tool:** `segment_cells` (mcp-deepcell)

**Status:** ⚠️ Currently mocked - returns synthetic segmentation masks

Segment individual cells from fluorescence images using DeepCell-TF.

**Parameters:**
- `image_path`: Path to IF image (grayscale or RGB)
- `nuclear_channel`: Channel with nuclear stain (e.g., DAPI)
- `segmentation_model`: "mesmer" (default) or "nuclear"

**Output:**
```json
{
  "status": "success",
  "cells_detected": 245,
  "segmentation_mask": "/path/to/segmentation_mask.npy",
  "cell_boundaries": "numpy array with cell IDs"
}
```

**Future:** Will use real DeepCell-TF implementation for production

### Step 4: Quantify Marker Expression
**Tool:** `classify_cell_states` (mcp-deepcell)

**Status:** ⚠️ Currently mocked - returns random classifications

Quantify marker expression per cell and classify phenotypes.

**Parameters:**
- `segmentation_mask`: Output from segment_cells
- `marker_channels`: Channels to quantify (e.g., ["CD8", "Ki67"])
- `threshold`: Intensity threshold for positive calls

**Output:**
```json
{
  "total_cells": 245,
  "cd8_positive": 12,
  "ki67_positive": 112,
  "cd8_ki67_double_positive": 3,
  "percent_cd8": 4.9,
  "percent_ki67": 45.7
}
```

### Step 5: Generate Visualization
**Tool:** `generate_segmentation_overlay` (mcp-deepcell)

**Status:** ⚠️ Currently mocked - returns placeholder overlays

Create overlay showing segmentation boundaries and phenotypes.

**Parameters:**
- `original_image`: Original IF image
- `segmentation_mask`: Cell boundaries
- `phenotype_labels`: Cell classifications

**Output:** PNG with colored cell boundaries showing:
- Green: CD8+ T cells
- Red: Ki67+ proliferating cells
- Yellow: Double-positive cells

---

## Expected Results (PatientOne)

### Single-Channel IF: CD8

**File:** `PAT001_tumor_IF_CD8.tiff`

**Findings:**
- **CD8+ cells detected:** ~12 cells (5-15 cells/mm²)
- **Distribution:** Periphery-dominant (immune exclusion phenotype)
- **Infiltration:** LOW - minimal T cell infiltration into tumor core

**Clinical relevance:** Immune exclusion suggests limited efficacy for checkpoint inhibitors

### Single-Channel IF: Ki67

**File:** `PAT001_tumor_IF_KI67.tiff`

**Findings:**
- **Ki67+ cells detected:** ~112 cells (45-55% proliferation index)
- **Distribution:** Heterogeneous across tissue
- **Proliferation:** HIGH - active tumor growth

**Clinical relevance:** High proliferation supports aggressive biology

### Multiplex IF: TP53/Ki67/DAPI

**File:** `PAT001_tumor_multiplex_IF_TP53_KI67_DAPI.tiff`

**Findings:**
- **Total cells:** ~245 (DAPI+ nuclei)
- **TP53+:** ~180 cells (73%) - mutant TP53 overexpression
- **Ki67+:** ~112 cells (46%) - proliferation
- **TP53+/Ki67+ double-positive:** ~85 cells (35%) - proliferating cells with TP53 mutation

**Clinical relevance:**
- Confirms TP53 mutation at protein level
- Co-expression with Ki67 shows actively dividing mutant cells
- Spatial heterogeneity of resistance markers

---

## Single-Cell Spatial Context

MxIF preserves spatial relationships:
- **Which cells are neighbors?** (tumor-immune interactions)
- **Are markers co-expressed?** (TP53+/Ki67+ double-positive)
- **Spatial gradients?** (marker expression decreasing with distance)
- **Tumor microenvironment?** (immune cells at periphery vs center)

---

## Current Implementation Status

| Component | Status | Implementation |
|-----------|--------|---------------|
| Image loading | ✅ Real | openimagedata fetch_histology_image |
| Composite generation | ✅ Real | openimagedata generate_multiplex_composite |
| Cell segmentation | ❌ Mock | deepcell segment_cells (synthetic masks) |
| Phenotype classification | ❌ Mock | deepcell classify_cell_states (random) |
| Visualization | ❌ Mock | deepcell generate_segmentation_overlay |

**Future:** Real DeepCell-TF implementation will enable production-grade cell segmentation and quantification.

---

## MxIF vs Other Modalities

| Modality | Purpose | Markers | Quantitative? | Server(s) |
|----------|---------|---------|---------------|-----------|
| **H&E** | Morphology | None (chromogenic stains) | No (visual) | openimagedata |
| **IF (single)** | Single protein | 1 marker | Yes | openimagedata + deepcell |
| **MxIF** | Multi-protein | 2-7+ markers | Yes | openimagedata + deepcell |
| **Spatial Transcriptomics** | Gene expression | Whole transcriptome | Yes | spatialtools |

---

## See Also

- [HE_WORKFLOW.md](HE_WORKFLOW.md) - H&E morphology workflow (no segmentation)
- [OPENIMAGEDATA_SERVER.md](OPENIMAGEDATA_SERVER.md) - Image loading and visualization
- [DEEPCELL_SERVER.md](DEEPCELL_SERVER.md) - Cell segmentation and phenotyping
- [GLOSSARY.md](GLOSSARY.md) - Imaging terminology
- [PatientOne TEST_4](../../tests/manual_testing/PatientOne-OvarianCancer/implementation/TEST_4_IMAGING.txt) - Complete imaging test
