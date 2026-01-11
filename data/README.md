# Data Directory

This directory contains sample data for testing and demonstrating the Precision Medicine MCP servers.

## 📂 Sample Patient Data Directory Structure

- **`patient-data/PAT001-OVC-2025/`** - Complete synthetic patient datasets for PatientOne workflow
  - **`clinical/`** - clinical patient data (JSON)
  - **`genomics/`** - variants - vcf data (VCF)
  - **`imaging/`** - histology and IF imaging data (TIFF)
  - **`multiomics/`** - multi-omics analysis data - RNA/Protein/Phospho (CSV)
  - **`imaging/`** - Visium output (CSV)

## 🏥 PatientOne Synthetic Dataset

For complete details on the PatientOne demonstration dataset (Stage IV Ovarian Cancer):

**📖 [Patient Data Documentation →](patient-data/PAT001-OVC-2025/README.md)**

This includes:
- Clinical data (demographics, CA-125 timeline)
- Genomic variants (VCF with key mutations)
- Multi-omics data (RNA-seq, proteomics, phosphoproteomics)
- Spatial transcriptomics (10x Visium, 900 spots, 31 genes)
- Imaging data (H&E, immunofluorescence)

**⚠️ All data is 100% synthetic** - Created for demonstration and testing purposes only.

---

**See also:** [PatientOne Architecture →](../architecture/patient-one/README.md)
