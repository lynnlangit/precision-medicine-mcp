# Automated Patient Report Generator

**One-command clinical analysis reports for spatial transcriptomics patients**

Generate comprehensive molecular analysis reports integrating FHIR clinical data with spatial transcriptomics analysis in seconds.

## Quick Start

```bash
# Using spatialtools virtual environment (recommended)
/Users/lynnlangit/Documents/GitHub/spatial-mcp/servers/mcp-spatialtools/venv/bin/python3 \
  scripts/generate_patient_report.py --patient-id patient-001 --output-dir ./results

# Or create an alias for convenience
alias analyze_patient='/Users/lynnlangit/Documents/GitHub/spatial-mcp/servers/mcp-spatialtools/venv/bin/python3 /Users/lynnlangit/Documents/GitHub/spatial-mcp/scripts/generate_patient_report.py'

# Then use it like this
analyze_patient --patient-id patient-001 --output-dir ./results
```

## What It Does

The automated report generator:

1. ✅ **Fetches clinical data** from GCP Healthcare API (FHIR)
   - Patient demographics
   - Diagnoses/conditions
   - Laboratory observations
   - Medications

2. ✅ **Loads spatial transcriptomics data**
   - Gene expression matrix
   - Tissue region annotations
   - Spatial coordinates

3. ✅ **Runs Phase 2 molecular analyses**
   - Differential expression (tumor vs stroma)
   - Spatial autocorrelation (Moran's I)
   - Cell type deconvolution
   - Pathway analysis

4. ✅ **Generates clinical summary**
   - Molecular findings
   - Drug resistance markers
   - Hypoxic zones
   - Treatment recommendations

5. ✅ **Saves structured results**
   - CSV files for each analysis
   - Clinical summary (text)
   - Metadata (JSON)

## Output Files

After running the script, you'll get:

```
results/
└── patient-001/
    ├── differential_expression.csv           # DEGs with log2FC and FDR
    ├── spatial_autocorrelation.csv           # Moran's I for all genes
    ├── cell_deconvolution.csv                # Cell type scores by region
    ├── clinical_summary.txt                  # Human-readable report
    ├── metadata.json                         # Analysis metadata
    ├── volcano_plot.png                      # Differential expression visualization
    ├── spatial_heatmap.png                   # Spatial gene expression patterns (top 6 genes)
    ├── cell_composition_heatmap.png          # Cell type enrichment by region
    ├── spatial_autocorrelation_plot.png      # Moran's I bar plot (top 15 genes)
    └── summary_figure.png                    # Multi-panel summary (publication-ready)
```

**Visualizations (300 DPI, PNG format):**
- ✅ **Volcano Plot** - Differential expression with labeled top genes
- ✅ **Spatial Heatmaps** - Expression overlaid on tissue coordinates
- ✅ **Cell Composition** - Heatmap of cell type signatures by region
- ✅ **Spatial Autocorrelation** - Bar chart of Moran's I values
- ✅ **Summary Figure** - 6-panel overview (DEGs, spatial patterns, cell types, stats)

## Example Output

### Clinical Summary (clinical_summary.txt)

```
================================================================================
CLINICAL ANALYSIS REPORT
================================================================================

Patient: Jane TestPatient
Gender: female
DOB: 1968-03-15
Report Date: 2025-12-29 11:21:47

DIAGNOSES:
  - Stage IV High-Grade Serous Ovarian Carcinoma (HGSOC), Platinum-Resistant

CURRENT MEDICATIONS:
  - Bevacizumab (Avastin) (active)
  - Paclitaxel (completed)
  - Carboplatin (completed)

MOLECULAR ANALYSIS SUMMARY:
--------------------------------------------------------------------------------

1. DIFFERENTIAL EXPRESSION (tumor_core vs stroma):
   Total DEGs: 17
   Upregulated: 13
   Downregulated: 4

   Top 5 upregulated genes in tumor:
      TP53: log2FC = 4.654, FDR = 5.04e-20
      KRT8: log2FC = 4.345, FDR = 2.68e-18
      ABCB1: log2FC = 4.285, FDR = 4.20e-18
      BCL2L1: log2FC = 3.683, FDR = 2.95e-19
      MKI67: log2FC = 3.564, FDR = 1.22e-14

2. SPATIAL ORGANIZATION:
   Spatially variable genes: 31

   Top 5 spatially clustered genes:
      HIF1A: Moran's I = 0.1411, p = 0.00e+00
      BCL2L1: Moran's I = 0.1269, p = 0.00e+00
      CD3D: Moran's I = 0.1103, p = 0.00e+00

3. TUMOR MICROENVIRONMENT:
   Cell type enrichment by region:

   tumor_core:
      fibroblasts: 67.1
      immune_cells: 81.7
      hypoxic: 60.6
      resistant: 716.8

   stroma:
      fibroblasts: 748.6
      immune_cells: 76.3
      hypoxic: 70.2
      resistant: 61.2

CLINICAL IMPLICATIONS:
--------------------------------------------------------------------------------

Drug Resistance Markers Detected:
  - PIK3CA: 3.11× fold change
  - AKT1: 3.24× fold change
  - MTOR: 3.17× fold change
  - ABCB1: 4.29× fold change

  ⚠️  Consider:
     - PI3K/AKT pathway inhibitors (Alpelisib, Capivasertib)
     - MDR reversal agents

Hypoxic Regions Identified:
  - HIF1A: Moran's I = 0.1411 (spatially clustered)
  - VEGFA: Moran's I = 0.1053 (spatially clustered)
  - CA9: Moran's I = 0.1026 (spatially clustered)

  ⚠️  Consider:
     - Hypoxia-targeting agents (Evofosfamide, TH-302)
     - VEGF inhibitors (Bevacizumab)

================================================================================
DISCLAIMER:
This analysis is for RESEARCH PURPOSES ONLY.
NOT validated for clinical decision-making.
All treatment decisions must be made by qualified oncologists.
================================================================================
```

### Differential Expression (differential_expression.csv)

| gene | mean_tumor | mean_stroma | log2_fold_change | p_value | fdr |
|------|------------|-------------|------------------|---------|-----|
| TP53 | 432.56 | 17.45 | 4.654 | 5.04e-20 | 5.04e-20 |
| KRT8 | 389.23 | 19.12 | 4.345 | 2.68e-18 | 2.68e-18 |
| ABCB1 | 378.91 | 19.45 | 4.285 | 4.20e-18 | 4.20e-18 |

### Spatial Autocorrelation (spatial_autocorrelation.csv)

| gene | morans_i | z_score | p_value |
|------|----------|---------|---------|
| HIF1A | 0.1411 | 85.11 | 0.00 |
| BCL2L1 | 0.1269 | 76.58 | 0.00 |
| CD3D | 0.1103 | 66.68 | 0.00 |

## Usage Examples

### Example 1: Analyze Patient-001

```bash
analyze_patient --patient-id patient-001 --output-dir ./results
```

**What happens:**
1. Fetches Jane TestPatient data from GCP Healthcare API
2. Loads spatial data from `/data/patient-data/PAT001-OVC-2025/spatial/`
3. Runs complete Phase 2 analysis
4. Generates report in `./results/patient-001/`

**Runtime:** ~10 seconds

### Example 2: Batch Process Multiple Patients

```bash
for patient in patient-001 patient-002 patient-003; do
  analyze_patient --patient-id $patient --output-dir ./results
done
```

### Example 3: Custom Output Location

```bash
analyze_patient \
  --patient-id patient-001 \
  --output-dir /path/to/clinical/reports/$(date +%Y%m%d)
```

## Requirements

**Dependencies:**
- Python 3.10+
- pandas
- numpy
- scipy
- requests
- Google Cloud SDK (for FHIR data)

**Data Requirements:**
- Patient FHIR data in GCP Healthcare API
- Spatial transcriptomics data:
  - `visium_gene_expression.csv`
  - `visium_region_annotations.csv`
  - `visium_spatial_coordinates.csv`

## Configuration

### GCP Healthcare API

The script expects:
- **Project:** `precision-medicine-poc`
- **Region:** `us-central1`
- **Dataset:** `precision-medicine-dataset`
- **FHIR Store:** `identified-fhir-store`

Configure with:
```bash
export GCP_PROJECT_ID="precision-medicine-poc"
export GCP_REGION="us-central1"
export GOOGLE_APPLICATION_CREDENTIALS="/path/to/service-account-key.json"
```

### Data Directory Structure

```
/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/
└── PAT001-OVC-2025/
    └── spatial/
        ├── visium_gene_expression.csv
        ├── visium_region_annotations.csv
        └── visium_spatial_coordinates.csv
```

## Command-Line Options

```bash
python scripts/generate_patient_report.py [OPTIONS]

Required:
  --patient-id TEXT      Patient ID (e.g., patient-001)

Optional:
  --output-dir PATH      Output directory (default: ./results)
  --help                 Show help message
```

## Integration with Workflows

### Use in Scripts

```python
import subprocess

def generate_report(patient_id, output_dir):
    """Generate patient report."""
    cmd = [
        "/path/to/venv/bin/python3",
        "scripts/generate_patient_report.py",
        "--patient-id", patient_id,
        "--output-dir", output_dir
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        print(f"✅ Report generated for {patient_id}")
        return f"{output_dir}/{patient_id}/clinical_summary.txt"
    else:
        print(f"❌ Error: {result.stderr}")
        return None
```

### Integrate with Claude Desktop

In Claude Desktop, you can ask:

```
Generate a clinical report for patient-001 using the automated report generator
```

Claude can run the script and summarize the results for you.

## Performance

**Patient-001 Analysis:**
- Data loading: 0.5 seconds
- Differential expression: 1 second
- Spatial autocorrelation: 2 seconds
- Cell deconvolution: 0.5 seconds
- Report generation: 0.5 seconds

**Total: ~5 seconds for complete analysis**

## Troubleshooting

### Issue: "ModuleNotFoundError: No module named 'pandas'"

**Solution:** Use the spatialtools virtual environment:
```bash
/Users/lynnlangit/Documents/GitHub/spatial-mcp/servers/mcp-spatialtools/venv/bin/python3 \
  scripts/generate_patient_report.py --patient-id patient-001 --output-dir ./results
```

### Issue: "Could not find spatial data for patient-XXX"

**Solution:** Check data directory structure:
```bash
ls -la /Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/
```

Ensure patient data exists in correct format.

### Issue: "FHIR data fetch failed"

**Solution:**
1. Check GCP credentials:
   ```bash
   gcloud auth list
   gcloud config get-value project
   ```

2. Verify FHIR store:
   ```bash
   gcloud healthcare fhir-stores list \
     --dataset=precision-medicine-dataset \
     --location=us-central1
   ```

3. Script will continue with spatial analysis only if FHIR fetch fails

### Issue: "Permission denied" for output directory

**Solution:** Create output directory or use different path:
```bash
mkdir -p ./results
# or
analyze_patient --patient-id patient-001 --output-dir ~/reports
```

## Next Steps

**Enhance the report generator:**

1. **Add PDF generation**
   - Install: `pip install reportlab matplotlib`
   - Generate publication-quality PDFs

2. **Add visualizations**
   - Spatial heatmaps
   - PCA plots
   - Pathway enrichment dot plots

3. **Upload to GCS**
   - Save results to `gs://precision-medicine-poc-results/`
   - Enable team access

4. **Email reports**
   - Automatically send to clinicians
   - Include summary + attachments

5. **Compare multiple patients**
   - Cohort analysis
   - Survival correlations
   - Subgroup identification

## Examples in the Wild

### Research Use Case

```bash
# Analyze all patients in a clinical trial
for patient in $(cat patient_list.txt); do
  echo "Analyzing $patient..."
  analyze_patient --patient-id $patient --output-dir ./trial_results/
done

# Aggregate results
python aggregate_trial_results.py --input-dir ./trial_results/
```

### Clinical Decision Support

```bash
# Generate report for tumor board
analyze_patient --patient-id PAT042 --output-dir ./tumor_board/$(date +%Y%m%d)/

# Review summary
cat ./tumor_board/$(date +%Y%m%d)/PAT042/clinical_summary.txt
```

## Visualizations Guide

The automated report includes 5 publication-quality visualizations (300 DPI PNG files):

### 1. Volcano Plot (volcano_plot.png)

**Purpose:** Visualize differential expression results

**Content:**
- X-axis: log2(Fold Change) - tumor vs stroma
- Y-axis: -log10(FDR) - statistical significance
- Red points: Upregulated in tumor (log2FC > 1, FDR < 0.05)
- Blue points: Downregulated in tumor (log2FC < -1, FDR < 0.05)
- Gray points: Not significant
- Top 10 genes labeled

**Threshold lines:**
- Vertical: |log2FC| = 1
- Horizontal: FDR = 0.05

**Clinical Use:** Quickly identify most significant molecular changes

### 2. Spatial Heatmap (spatial_heatmap.png)

**Purpose:** Show spatial expression patterns for top genes

**Content:**
- 6 panels (2×3 grid)
- Top 6 spatially variable genes (highest Moran's I)
- Each panel shows:
  - Gene expression overlaid on tissue coordinates
  - Viridis colormap (purple=low, yellow=high expression)
  - Moran's I value in title
  - Array row/column coordinates

**Clinical Use:** Identify spatially organized gene expression (hypoxic zones, tumor boundaries, etc.)

### 3. Cell Composition Heatmap (cell_composition_heatmap.png)

**Purpose:** Visualize cell type enrichment across tissue regions

**Content:**
- Heatmap: Cell types (rows) × Tissue regions (columns)
- Color intensity: Signature score (yellow-orange-red scale)
- Annotated with numeric values
- Cell types: fibroblasts, immune_cells, hypoxic, resistant

**Clinical Use:** Understand tumor microenvironment composition

### 4. Spatial Autocorrelation Plot (spatial_autocorrelation_plot.png)

**Purpose:** Rank genes by spatial clustering strength

**Content:**
- Horizontal bar chart
- Top 15 genes with highest Moran's I
- Steel blue bars
- Moran's I values labeled on bars
- Genes sorted by spatial autocorrelation strength

**Clinical Use:** Identify genes with strongest spatial organization

### 5. Summary Figure (summary_figure.png)

**Purpose:** Multi-panel overview for presentations/publications

**Content (6 panels):**

**Panel 1 (Top Left, 2/3 width):** Top 10 DEGs
- Horizontal bar chart
- Red = upregulated, Blue = downregulated
- log2(Fold Change) values

**Panel 2 (Top Right):** Cell Type Enrichment
- Heatmap of cell type scores by region
- Compact version of cell composition heatmap

**Panel 3 & 4 (Middle Left/Center):** Top 2 Spatial Genes
- Spatial expression patterns
- HIF1A and BCL2L1 typically shown
- Expression overlaid on coordinates

**Panel 5 (Middle Right):** Top 8 Spatial Genes
- Bar chart of Moran's I values
- Compact version of autocorrelation plot

**Panel 6 (Bottom Left):** Region Distribution
- Pie chart showing proportion of spots in each tissue region

**Panel 7 (Bottom Right):** Analysis Summary
- Text box with key statistics:
  - Number of genes/spots/regions
  - Number of DEGs and SVGs
  - Drug resistance markers detected
  - Patient ID and date

**Clinical Use:** Single-page summary for tumor boards, presentations, or publications

## Visualization Specifications

**Format:** PNG (Portable Network Graphics)
**Resolution:** 300 DPI (publication quality)
**Color Schemes:**
- Differential expression: Red/Blue (diverging)
- Spatial expression: Viridis (sequential, colorblind-friendly)
- Cell composition: YlOrRd (yellow-orange-red, sequential)
- Bar charts: Steel blue or red/blue

**Size Estimates:**
- Volcano plot: ~150 KB
- Spatial heatmap: ~1.8 MB (6 panels, high resolution)
- Cell composition: ~220 KB
- Spatial autocorrelation: ~140 KB
- Summary figure: ~1.0 MB (multi-panel)

**Total:** ~3.3 MB for all visualizations

## API Reference

### PatientReportGenerator Class

```python
class PatientReportGenerator:
    def __init__(self, patient_id: str, output_dir: str):
        """Initialize report generator."""

    def fetch_fhir_data(self):
        """Fetch patient clinical data from GCP Healthcare API."""

    def load_spatial_data(self):
        """Load spatial transcriptomics data."""

    def perform_differential_expression(self):
        """Run differential expression analysis."""

    def calculate_spatial_autocorrelation(self):
        """Calculate spatial autocorrelation for all genes."""

    def perform_cell_deconvolution(self):
        """Perform cell type deconvolution."""

    def generate_clinical_summary(self):
        """Generate clinical interpretation summary."""

    def generate_report(self):
        """Run complete analysis pipeline."""
```

## License

Research use only. Not for clinical decision-making.

## Support

- **Documentation:** See README files in `servers/mcp-*` directories
- **Issues:** Report at GitHub repository
- **Questions:** Ask in Claude Desktop or Claude Code

---

**🎉 One command. Complete clinical analysis. Seconds, not hours.**

```bash
analyze_patient --patient-id patient-001 --output-dir ./results
```
