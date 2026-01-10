# Imaging Glossary

Imaging-specific terminology for histology and immunofluorescence analysis.

---

## B

**Brightfield Microscopy** - Standard light microscopy using transmitted white light. Used for H&E staining. Produces RGB color images from chromogenic dyes.

---

## C

**Cellularity** - Percentage of viable tumor cells in a tissue sample. Typically assessed visually from H&E slides. Important for determining tumor burden.

**Chromogenic Stain** - Colored dye that absorbs light (not fluorescent). Examples: Hematoxylin (blue), Eosin (pink). Used in H&E staining.

**Co-expression** - When a single cell expresses multiple markers simultaneously. Example: TP53+/Ki67+ double-positive cells show both mutant TP53 and proliferation.

---

## D

**DAPI** - 4',6-diamidino-2-phenylindole. Blue fluorescent nuclear stain that binds DNA. Used as nuclear counterstain in IF/MxIF imaging to identify all cells.

**DeepCell-TF** - Open-source deep learning library for cell segmentation in microscopy images. Based on TensorFlow. Repository: https://github.com/vanvalenlab/deepcell-tf

---

## E

**Eosin** - Pink/red chromogenic dye that stains cytoplasm, collagen, and muscle fibers. Second component of H&E staining.

---

## F

**Fluorescence Microscopy** - Microscopy using fluorescent dyes/proteins that emit light when excited. Used for IF and MxIF imaging. Produces grayscale (single channel) or multi-channel images.

---

## H

**H&E (Hematoxylin & Eosin)** - Standard histology stain using two chromogenic dyes:
- **Hematoxylin:** Stains nuclei blue/purple
- **Eosin:** Stains cytoplasm pink/red
- **Microscopy:** Brightfield (NOT fluorescence)
- **Purpose:** Visual morphology assessment

**Hematoxylin** - Blue/purple chromogenic dye that stains cell nuclei and other basophilic structures. First component of H&E staining.

**HGSOC** - High-Grade Serous Ovarian Carcinoma. Aggressive ovarian cancer subtype characterized by TP53 mutations, papillary architecture, and high nuclear grade.

**Histology** - Microscopic study of tissue architecture and morphology. Typically uses H&E staining on formalin-fixed paraffin-embedded (FFPE) tissue sections.

---

## I

**IF (Immunofluorescence)** - Microscopy technique using fluorescently-labeled antibodies to detect specific proteins in tissue. Single-channel IF visualizes one marker per image.

**Immune Exclusion** - Phenotype where immune cells (e.g., CD8+ T cells) are present at tumor periphery but cannot penetrate tumor core. Associated with poor immunotherapy response.

---

## K

**Ki67** - Nuclear protein marker of cell proliferation. Ki67+ cells are actively dividing. Proliferation index = % Ki67+ cells. High Ki67 (>40%) indicates aggressive tumor growth.

---

## M

**Magnification** - Optical enlargement of microscopy images. Common levels:
- **20×:** Wide-field tissue overview (typical for demos)
- **40×:** High-resolution cellular detail
- **60-100×:** Oil immersion for subcellular structures

**Mesmer** - DeepCell segmentation model for whole-cell segmentation in multiplexed imaging. Uses both nuclear and membrane markers to identify cell boundaries.

**Morphology** - Study of tissue and cellular structure/shape. Assessed visually from H&E images by pathologists.

**MxIF (Multiplexed Immunofluorescence)** - Fluorescence imaging of multiple protein markers (2-7+) on a single tissue section. Enables:
- Multi-marker co-expression analysis
- Spatial protein relationships
- Quantitative phenotyping
- **Technology:** Repeated rounds of stain → image → inactivate → restain

---

## N

**Necrosis** - Cell death visible in tissue as pale/acellular regions with loss of nuclear detail. In H&E images, appears as bright pink (eosinophilic) areas with no blue nuclei.

---

## P

**Phenotype** - Observable characteristics of a cell defined by marker expression. Examples:
- CD8+ phenotype = Cytotoxic T cells
- Ki67+ phenotype = Proliferating cells
- TP53+/Ki67+ phenotype = Proliferating cells with TP53 mutation

**Proliferation Index** - Percentage of cells positive for Ki67 proliferation marker. Indicates tumor growth rate.

---

## S

**Segmentation** - Computational process of identifying individual cell boundaries in microscopy images. Required for quantitative single-cell analysis.

**Segmentation Mask** - Image where each pixel is labeled with a cell ID (0 = background, 1 = cell 1, 2 = cell 2, etc.). Output of cell segmentation algorithms.

**Staining** - Process of adding dyes/antibodies to tissue sections to visualize specific structures or proteins:
- **Chromogenic staining:** H&E (brightfield)
- **Fluorescent staining:** IF/MxIF (fluorescence)

---

## T

**TIFF** - Tagged Image File Format. Standard format for microscopy images. Supports:
- **Grayscale:** Single-channel IF (1 channel)
- **RGB:** H&E (3 channels), multiplex IF (3+ channels)
- **Multi-page:** Z-stacks, time series

**TP53** - Tumor protein p53, tumor suppressor. Mutations in TP53 lead to protein accumulation detectable by IF. Present in >96% of HGSOC cases.

---

## Imaging Modality Comparison

| Term | Microscopy | Staining | Format | Quantitative | Purpose |
|------|-----------|----------|--------|--------------|---------|
| **H&E** | Brightfield | Chromogenic | RGB | No (visual) | Morphology |
| **IF (single)** | Fluorescence | Fluorescent antibody | Grayscale | Yes | Single protein |
| **MxIF** | Fluorescence | Multiple antibodies | Multi-channel | Yes | Multiple proteins |

---

**See Also:**
- [Spatial Transcriptomics Glossary](../spatial-transcriptomics/GLOSSARY.md) - Spatial analysis terms
- [HE_WORKFLOW.md](HE_WORKFLOW.md) - H&E workflow guide
- [MXIF_WORKFLOW.md](MXIF_WORKFLOW.md) - MxIF workflow guide
