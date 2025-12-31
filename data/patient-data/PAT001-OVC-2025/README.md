# PAT001-OVC-2025: Synthetic Patient Data

## Overview

This directory contains **100% synthetic patient data** created for demonstration and testing purposes. The data represents a realistic but completely fabricated Stage IV High-Grade Serous Ovarian Carcinoma (HGSOC) patient case.

**Patient:** Sarah Anderson (synthetic)
**Diagnosis:** Stage IV HGSOC, Platinum-Resistant
**Purpose:** Demonstrate precision medicine MCP server workflows

---

## Data Origin & Creation

### ⚠️ CRITICAL: ALL DATA IS SYNTHETIC

**This is NOT real patient data.** All files were artificially generated using:

1. **Statistical models** - Random number generators with clinically plausible distributions
2. **Expert knowledge** - Curated to reflect realistic ovarian cancer biology
3. **Published literature** - Based on TCGA ovarian cancer studies and clinical trial data
4. **No PHI/PII** - Contains zero protected health information or personally identifiable information

**Source:** Created by the Precision Medicine MCP development team for demonstration purposes.

---

## Data Modalities (5 types, 19 files)

### 1. Clinical Data (`clinical/`)

| File | Size | Description | Source |
|------|------|-------------|--------|
| `clinical_demographics.json` | 532 B | Patient demographics, stage, diagnosis | Synthetic |
| `ca125_timeline.csv` | 365 B | CA-125 tumor marker trajectory | Modeled after clinical patterns |

**Generation method:**
- Demographics: Randomly generated name, age, dates
- CA-125 values: Gaussian noise around typical platinum-resistant trajectory

---

### 2. Genomic Data (`genomics/`)

| File | Size | Description | Source |
|------|------|-------------|--------|
| `PAT001_somatic_variants.vcf` | 2.3 KB | Somatic mutations (TP53, PIK3CA, PTEN) | Based on TCGA-OV data |

**Generation method:**
- Selected common ovarian cancer driver mutations from TCGA
- VCF format with realistic variant quality scores
- Germline BRCA1 mutation added to simulate HRD-positive phenotype

**Key mutations (realistic for HGSOC):**
- TP53 R175H (missense, chr17:7577538)
- PIK3CA E545K (missense, chr3:178936091)
- PTEN LOH (loss of heterozygosity, chr10)
- BRCA1 germline (c.5266dupC, pathogenic)

---

### 3. Multi-Omics Data (`multiomics/`)

| File | Size | Description | Source |
|------|------|-------------|--------|
| `pdx_rna_expression.csv` | 15 KB | RNA-seq for 15 PDX samples | Synthetic |
| `pdx_protein.csv` | 12 KB | Proteomics for 15 PDX samples | Synthetic |
| `pdx_phospho.csv` | 11 KB | Phosphoproteomics for 15 PDX samples | Synthetic |
| `pdx_metadata.csv` | 426 B | Sample annotations (resistant vs sensitive) | Synthetic |

**Generation method:**
- Random sampling from log-normal distributions
- Resistant samples: Elevated PI3K/AKT/mTOR pathway genes
- Sensitive samples: Lower PI3K/AKT activation
- Noise added to simulate biological variability

**Simulated biology:**
- 7 resistant samples (high PIK3CA, AKT1, MTOR, ABCB1)
- 8 sensitive samples (lower pathway activation)
- Fold-changes calibrated to published PDX resistance studies

---

### 4. Spatial Transcriptomics (`spatial/`)

| File | Size | Description | Source |
|------|------|-------------|--------|
| `visium_spots.csv` | 34 KB | Spatial coordinates for 900 spots | Synthetic |
| `visium_expression.csv` | 117 KB | Gene expression (31 genes × 900 spots) | Synthetic |
| `visium_regions.csv` | 18 KB | Tissue region annotations | Synthetic |
| `visium_spatial_data.csv` | 146 KB | Combined coordinates + expression | Synthetic |

**Generation method:**
- Hexagonal Visium array layout (900 spots)
- 31 curated cancer-related genes
- Spatial patterns using Gaussian random fields:
  - Tumor core: High proliferation (Ki67, PCNA, TOP2A)
  - Tumor-stroma interface: Immune markers (CD3D, CD8A)
  - Necrotic regions: Hypoxia markers (HIF1A, CA9, VEGFA)
  - Stroma: Fibroblast markers (COL1A1, FAP, ACTA2)

**Tissue regions simulated (6):**
1. `tumor_core` - Central tumor, high proliferation
2. `tumor_proliferative` - Actively dividing tumor cells
3. `tumor_interface` - Invasion front
4. `stroma` - Fibroblast-rich connective tissue
5. `stroma_immune` - Immune infiltration
6. `necrotic_hypoxic` - Oxygen-deprived regions

**Genes profiled (31):**
- Proliferation: Ki67, PCNA, TOP2A, MYC
- Tumor suppressors: TP53, PTEN
- Oncogenes: PIK3CA, AKT1, MTOR
- Resistance: ABCB1 (MDR1), BCL2L1
- Hypoxia: HIF1A, CA9, VEGFA
- Immune: CD3D, CD8A, CD4, FOXP3, CD68
- Stromal: COL1A1, COL3A1, ACTA2, FAP
- Epithelial: EPCAM, KRT8, KRT18

---

### 5. Imaging Data (`imaging/`)

| File | Size | Description | Source |
|------|------|-------------|--------|
| `HE_tissue_section.tif` | 1.2 MB | H&E histology (placeholder) | Synthetic placeholder |
| `IF_DAPI.tif` | 512 KB | Nuclear staining (placeholder) | Synthetic placeholder |
| `IF_CD3.tif` | 512 KB | T cell marker (placeholder) | Synthetic placeholder |
| `IF_CD8.tif` | 512 KB | Cytotoxic T cell marker (placeholder) | Synthetic placeholder |
| `IF_Ki67.tif` | 512 KB | Proliferation marker (placeholder) | Synthetic placeholder |
| `IF_PanCK.tif` | 512 KB | Epithelial marker (placeholder) | Synthetic placeholder |
| `IF_multiplex.tif` | 2.1 MB | 6-plex IF image (placeholder) | Synthetic placeholder |

**Note:** Image files are placeholder TIFFs. In DRY_RUN mode, MCP servers return synthetic segmentation/analysis results without processing actual images.

---

## Data Quality & Realism

### Realistic Features ✅

1. **Clinically plausible mutations** - Common HGSOC drivers (TP53, PIK3CA, BRCA1)
2. **Biologically coherent** - Resistant samples show expected pathway activation
3. **Spatial organization** - Tumor heterogeneity and microenvironment structure
4. **Temporal dynamics** - CA-125 trajectory matches platinum-resistant pattern

### Synthetic Limitations ❌

1. **No batch effects** - Real data would have technical variation
2. **Simplified biology** - Real tumors have >20,000 genes, not 31
3. **Idealized spatial patterns** - Real tissue is more heterogeneous
4. **No sequencing errors** - Real VCF would have quality issues
5. **Perfect correlations** - Real multi-omics has more noise

---

## Intended Use Cases

### ✅ Appropriate Uses

- **Workflow development** - Test MCP server integration
- **CI/CD testing** - Automated pipeline validation
- **Demonstration** - Show capabilities to stakeholders
- **Education** - Teach precision medicine concepts
- **Software testing** - Verify analysis algorithms

### ❌ Inappropriate Uses

- **Clinical decision-making** - NOT real patient data
- **Publications** - Cannot be cited as real data
- **Regulatory submissions** - Not suitable for FDA/EMA
- **Benchmarking** - Use established datasets instead
- **Training ML models** - Will not generalize to real data

---

## Data License & Attribution

**License:** Apache 2.0 (same as repository)

**Attribution:** If using this synthetic dataset in presentations or educational materials, please cite:

```
Precision Medicine MCP Synthetic Patient Dataset (PAT001-OVC-2025)
Generated for demonstration purposes
https://github.com/lynnlangit/precision-medicine-mcp
2025-12
```

---

## Real Data Alternatives

For **actual** ovarian cancer datasets, use these public resources:

### Genomics
- **TCGA-OV** - The Cancer Genome Atlas Ovarian Cancer
  - https://portal.gdc.cancer.gov/projects/TCGA-OV
  - 489 patients, WGS/RNA-seq/clinical
  - Free, public access via GDC API

### Spatial Transcriptomics
- **10x Genomics Public Datasets**
  - https://www.10xgenomics.com/datasets
  - Human ovarian cancer samples
  - Visium platform, download FASTQ/BAM

### Multi-Omics
- **CPTAC Ovarian Cancer**
  - https://proteomics.cancer.gov/programs/cptac
  - Proteomics + genomics + phosphoproteomics
  - 174 HGSOC samples

### Imaging
- **The Cancer Imaging Archive (TCIA)**
  - https://www.cancerimagingarchive.net/
  - H&E, CT, MRI for ovarian cancer
  - Free download with registration

---

## Technical Details

### File Formats

- **JSON** - Clinical demographics (RFC 8259 compliant)
- **CSV** - Tabular data (comma-delimited, UTF-8)
- **VCF** - Variant Call Format v4.2 (SAMtools spec)
- **TIFF** - Tagged Image File Format (16-bit grayscale)

### Data Size

**Total:** ~4.5 MB uncompressed
- Clinical: <1 KB
- Genomics: 2.3 KB
- Multi-omics: 38 KB
- Spatial: 315 KB
- Imaging: 4.1 MB (placeholders)

### Checksums

For data integrity verification:
```bash
cd /path/to/spatial-mcp/data/patient-data/PAT001-OVC-2025
find . -type f -exec sha256sum {} \; > checksums.txt
```

---

## Related Documentation

- **PatientOne Workflow:** `tests/manual_testing/PatientOne-OvarianCancer/README.md`
- **Data Modes Guide:** `tests/manual_testing/PatientOne-OvarianCancer/DATA_MODES_GUIDE.md`
- **Architecture:** `architecture/patient-one/README.md`
- **Main README:** `README.md`

---

## Support

**Issues:** https://github.com/lynnlangit/precision-medicine-mcp/issues
**Questions:** See repository documentation

---

**Last Updated:** 2025-12-29
**Dataset Version:** 1.0
**Status:** 100% Synthetic - Research/Demo Use Only
