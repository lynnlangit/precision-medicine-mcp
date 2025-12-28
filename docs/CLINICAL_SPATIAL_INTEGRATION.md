# Clinical-Spatial Integration Guide

**Phase 1C Complete** - Connecting FHIR clinical data with spatial transcriptomics analysis

---

## Overview

The clinical-spatial integration pipeline enables seamless analysis connecting:
- **Clinical Data**: Patient demographics, diagnoses, medications, lab results (via mcp-epic)
- **Spatial Transcriptomics**: Gene expression with tissue coordinates (via mcp-spatialtools)

This integration provides **clinically-relevant spatial analysis** by automatically identifying genes of interest based on patient conditions and treatments.

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    Claude Desktop Interface                      │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ├─────────────────────┐
                              │                     │
                              ▼                     ▼
                    ┌──────────────────┐  ┌──────────────────┐
                    │   mcp-epic       │  │ mcp-spatialtools │
                    │  (FHIR Data)     │  │  (Spatial Data)  │
                    └──────────────────┘  └──────────────────┘
                              │                     │
                              │    Bridge Tool      │
                              │  ┌──────────────┐   │
                              └─>│ Patient ID   │<──┘
                                 │   Mapping    │
                                 └──────────────┘
                                        │
                                        ▼
                              ┌──────────────────┐
                              │ Genes of Interest│
                              │   + Analyses     │
                              └──────────────────┘
                                        │
                                        ▼
                              ┌──────────────────┐
                              │ Spatial Analysis │
                              │   (Moran's I)    │
                              └──────────────────┘
                                        │
                                        ▼
                              ┌──────────────────┐
                              │ Integrated Report│
                              └──────────────────┘
```

---

## Key Components

### 1. Clinical-Spatial Bridge Tool

**Function**: `get_spatial_data_for_patient`
**Purpose**: Link patient clinical data to spatial transcriptomics datasets

**Inputs**:
- `patient_id`: FHIR patient identifier
- `conditions`: List of diagnoses (e.g., "ovarian cancer", "HGSOC")
- `medications`: Current/past treatments
- `biomarkers`: Lab results and molecular markers
- `tissue_type`: Tissue sample type ("tumor", "normal", "margin")

**Outputs**:
- `spatial_dataset`: Mapped spatial dataset ID
- `files`: Paths to expression, coordinates, and annotation files
- `genes_of_interest`: Genes relevant to patient's clinical context
- `suggested_analyses`: Recommended spatial analyses

**Mapping Logic**:

```python
# Condition → Genes
"ovarian cancer" → [Ki67, TP53, BRCA1, BRCA2, EPCAM, CA125, PAX8]
"HGSOC" → [TP53, Ki67, FOXM1, MYC, CCNE1]
"platinum-resistant" → [ABCB1, ERCC1, GSTP1, BRCA1]

# Treatment → Biomarkers
"Bevacizumab" → [VEGFA, CD31, HIF1A, KDR]  # Anti-angiogenic
"Carboplatin" → [ERCC1, XPA, BRCA1, BRCA2] # DNA repair

# Biomarker → Genes
"CA-125" → [MUC16, MSLN, HE4]
"BRCA" → [BRCA1, BRCA2, RAD51, PARP1]
```

### 2. Spatial Autocorrelation Analysis

**Function**: `calculate_spatial_autocorrelation`
**Purpose**: Quantify spatial clustering patterns of gene expression

**Implementation**: Real Moran's I calculation using scipy

**Algorithm**:
1. Build spatial weights matrix based on distance threshold
2. Row-standardize weights (each row sums to 1)
3. Calculate Moran's I statistic
4. Compute z-score and p-value for significance testing

**Moran's I Formula**:
```
I = (n / W) × (Σᵢ Σⱼ wᵢⱼ(xᵢ - x̄)(xⱼ - x̄)) / (Σᵢ(xᵢ - x̄)²)

where:
- n = number of spots
- W = sum of all weights
- wᵢⱼ = spatial weight between spots i and j
- xᵢ = expression value at spot i
- x̄ = mean expression across all spots
```

**Interpretation**:
- `I > 0`: Clustered pattern (nearby spots have similar expression)
- `I ≈ 0`: Random pattern (no spatial structure)
- `I < 0`: Dispersed pattern (nearby spots have dissimilar expression)

**Parameters**:
- `distance_threshold`: Maximum distance for neighbors (default: 1500 pixels for Visium)
- `genes`: List of genes to analyze
- `method`: "morans_i" (other methods can be added)

---

## Usage Examples

### Example 1: Query Clinical Data and Retrieve Spatial Dataset

**In Claude Desktop:**
```
Get clinical data for patient-001 using the epic server, then retrieve
the corresponding spatial transcriptomics data using the bridge tool.
```

**What happens**:
1. Epic server queries GCP Healthcare API FHIR store
2. Returns de-identified patient data (demographics, conditions, medications, observations)
3. Bridge tool maps patient-001 → PAT001-OVC-2025 spatial dataset
4. Identifies 13 relevant genes based on clinical context
5. Suggests analyses (e.g., "Evaluate angiogenesis markers for treatment response")

### Example 2: Spatial Autocorrelation for Clinical Genes

**In Claude Desktop:**
```
Calculate spatial autocorrelation for genes relevant to patient-001's
ovarian cancer diagnosis and current bevacizumab treatment.
```

**What happens**:
1. Bridge tool identifies genes: MKI67 (proliferation), CD8A (immune), VIM (EMT), etc.
2. Loads spatial expression data from PAT001-OVC-2025
3. Calculates Moran's I for each gene with statistical significance
4. Returns clustering patterns with p-values

**Example Output**:
```
Gene    | Moran's I | Z-score | p-value | Significance | Pattern
--------|-----------|---------|---------|--------------|------------------
MKI67   |    0.0451 |  27.111 |  0.0000 |    ***       | Weakly patterned
CD8A    |    0.0664 |  39.635 |  0.0000 |    ***       | Weakly patterned
VIM     |    0.0574 |  34.369 |  0.0000 |    ***       | Weakly patterned
VEGFA   |    0.0823 |  48.912 |  0.0000 |    ***       | Weakly patterned
```

### Example 3: Integrated Clinical-Spatial Report

**In Claude Desktop:**
```
Generate a clinical-spatial analysis report for patient-001 integrating
FHIR clinical data with spatial transcriptomics findings.
```

**Generated Report**:
```
PATIENT 001 - OVARIAN CANCER SPATIAL ANALYSIS

Clinical Summary:
  • Diagnosis: Stage IV High-Grade Serous Ovarian Carcinoma
  • Treatment: Platinum-resistant, currently on Bevacizumab
  • Biomarkers: CA-125 elevated (487 U/mL), BRCA negative

Spatial Findings:
  MKI67 (Proliferation):
    Moran's I: 0.0451 (p < 0.001)
    Clinical Relevance: Distributed proliferation pattern consistent
                        with heterogeneous tumor biology

  CD8A (T-cell Infiltration):
    Moran's I: 0.0664 (p < 0.001)
    Clinical Relevance: Diffuse immune infiltration may indicate
                        potential immunotherapy response

  VIM (Mesenchymal Marker):
    Moran's I: 0.0574 (p < 0.001)
    Clinical Relevance: Distributed EMT features correlate with
                        platinum resistance mechanisms

Integrated Conclusions:
  • Spatial patterns reveal tumor heterogeneity
  • Immune infiltration supports bevacizumab continuation
  • EMT signatures align with platinum-resistant phenotype
```

---

## Testing

### Standalone Tests

**Test Bridge Tool**:
```bash
cd servers/mcp-spatialtools
SPATIAL_DATA_DIR="../../data" SPATIAL_DRY_RUN="false" \
  venv/bin/python test_bridge_tool.py
```

**Test Spatial Autocorrelation**:
```bash
cd servers/mcp-spatialtools
SPATIAL_DATA_DIR="../../data" SPATIAL_DRY_RUN="false" \
  venv/bin/python test_spatial_autocorrelation.py
```

**Test Integrated Workflow**:
```bash
cd servers/mcp-spatialtools
SPATIAL_DATA_DIR="../../data" SPATIAL_DRY_RUN="false" \
  venv/bin/python test_integrated_workflow.py
```

### Claude Desktop Testing

1. **Ensure servers are configured** in `claude_desktop_config.json`:
   ```json
   {
     "mcpServers": {
       "epic": {
         "env": {
           "GCP_PROJECT_ID": "precision-medicine-poc",
           "DEIDENTIFY_ENABLED": "true"
         }
       },
       "spatialtools": {
         "env": {
           "SPATIAL_DATA_DIR": "/path/to/data",
           "SPATIAL_DRY_RUN": "false"
         }
       }
     }
   }
   ```

2. **Restart Claude Desktop**

3. **Test with prompts**:
   - "Get clinical data for patient-001"
   - "Get spatial data for patient-001 with clinical context"
   - "Calculate spatial autocorrelation for MKI67, CD8A, VIM"
   - "Generate a clinical-spatial report"

---

## Data Requirements

### Clinical Data (FHIR Resources)

Required in GCP Healthcare API FHIR Store:
- **Patient**: Demographics (de-identified)
- **Condition**: Diagnoses (ICD-10 codes, clinical descriptions)
- **MedicationStatement**: Current and past treatments
- **Observation**: Lab results, biomarkers

### Spatial Data

Required files in `data/patient-data/{PATIENT_ID}/spatial/`:
- `visium_gene_expression.csv`: Expression matrix (spots × genes)
- `visium_spatial_coordinates.csv`: Spot coordinates (x, y)
- `visium_region_annotations.csv`: Tissue region labels (optional)

**File Format**:
```
visium_gene_expression.csv:
    barcode,MKI67,CD8A,VIM,EPCAM,TP53,...
    AAACAAGTATCTCCCA-1,5.2,3.1,7.8,2.4,6.5,...
    AAACACCAATAACTGC-1,4.8,2.9,8.2,3.1,5.9,...

visium_spatial_coordinates.csv:
    barcode,x,y,in_tissue
    AAACAAGTATCTCCCA-1,1234.5,2345.6,1
    AAACACCAATAACTGC-1,1256.7,2367.8,1
```

---

## Patient ID Mapping

The bridge tool maps FHIR patient IDs to spatial dataset IDs:

```python
PATIENT_SPATIAL_MAP = {
    "patient-001": "PAT001-OVC-2025",  # Ovarian cancer patient
    # Add more mappings as datasets are added
}
```

To add a new patient:
1. Upload FHIR data to GCP Healthcare API
2. Place spatial data in `data/patient-data/{SPATIAL_ID}/spatial/`
3. Add mapping to `PATIENT_SPATIAL_MAP` in `server.py`

---

## Clinical Gene Mappings

### Condition-Based Genes

| Condition | Genes of Interest | Rationale |
|-----------|-------------------|-----------|
| Ovarian Cancer | Ki67, TP53, BRCA1/2, EPCAM, CA125, PAX8 | Proliferation, tumor suppressor, epithelial markers |
| HGSOC | TP53, Ki67, FOXM1, MYC, CCNE1 | Cell cycle, oncogenes common in HGSOC |
| Platinum-Resistant | ABCB1, ERCC1, GSTP1, BRCA1 | Drug efflux, DNA repair mechanisms |

### Treatment-Based Biomarkers

| Treatment | Biomarkers | Mechanism |
|-----------|------------|-----------|
| Bevacizumab | VEGFA, CD31, HIF1A, KDR | Anti-angiogenic therapy targets |
| Carboplatin | ERCC1, XPA, BRCA1, BRCA2 | DNA damage repair pathways |
| Paclitaxel | TUBB3, MAP2, MAPT | Microtubule stabilization |

### Biomarker-Based Genes

| Biomarker | Related Genes | Connection |
|-----------|---------------|------------|
| CA-125 (elevated) | MUC16, MSLN, HE4 | Ovarian cancer tumor markers |
| BRCA (mutation) | BRCA1, BRCA2, RAD51, PARP1 | Homologous recombination repair |

---

## Performance Metrics

### Spatial Autocorrelation Calculation

- **900 spots × 5 genes**: ~0.5 seconds
- **900 spots × 20 genes**: ~2 seconds
- **Distance threshold**: 1500 pixels (optimal for Visium, ~1000 pixel spot spacing)

### Memory Usage

- Expression data (900 spots × 31 genes): ~200 KB
- Spatial weights matrix (900 × 900): ~6.4 MB
- Peak memory: <50 MB for typical analysis

---

## Implementation Details

### Real vs. Mocked Tools

**Real Implementations (Phase 1C)**:
- ✅ `get_spatial_data_for_patient` - Clinical-spatial bridge
- ✅ `calculate_spatial_autocorrelation` - Moran's I with scipy
- ✅ `filter_quality` - QC filtering (from earlier)
- ✅ `split_by_region` - Region segmentation (from earlier)

**Still Mocked**:
- ⏳ Differential expression analysis
- ⏳ Batch correction
- ⏳ Spatial alignment
- ⏳ Cell type deconvolution

### Statistical Methods

**Moran's I Z-score Calculation**:
```python
Expected value: E[I] = -1 / (n - 1)
Variance: Var[I] = (n × S1 - S2 + 3 × W²) / (W² × (n² - 1)) - E[I]²
Z-score: z = (I - E[I]) / √Var[I]
p-value: p = 2 × (1 - Φ(|z|))  # Two-tailed test

where:
  S1 = 0.5 × Σᵢⱼ (wᵢⱼ + wⱼᵢ)²
  S2 = Σᵢ (Σⱼ wᵢⱼ + Σⱼ wⱼᵢ)²
  Φ = standard normal CDF
```

---

## Troubleshooting

### Bridge Tool Not Finding Spatial Data

**Error**: `"error": "No spatial data found for patient"`

**Fix**:
1. Check patient ID mapping in `PATIENT_SPATIAL_MAP`
2. Verify spatial data exists at expected path
3. Ensure `SPATIAL_DATA_DIR` environment variable is set correctly

### Moran's I Returns NaN

**Error**: All Moran's I values are NaN with divide-by-zero warnings

**Fix**:
1. Check distance threshold (too small → no neighbors)
2. For Visium data, use `distance_threshold=1500` pixels
3. Verify coordinate file has correct x,y values

### Gene Not Found in Expression Data

**Error**: `"message": "Gene {gene} not found in expression data"`

**Fix**:
1. Check gene naming (e.g., "Ki67" vs "MKI67")
2. Verify gene exists in expression matrix columns
3. Update gene mappings to match actual gene names in data

---

## Future Enhancements

### Phase 2 Additions
- [ ] Spatial differential expression (between regions)
- [ ] Cell type deconvolution integration
- [ ] Pathway enrichment for spatially clustered genes
- [ ] Spatial correlation networks

### Phase 3 Additions
- [ ] Multi-sample spatial comparisons
- [ ] Treatment response prediction models
- [ ] Spatial biomarker discovery
- [ ] Integration with imaging data

---

## References

### Moran's I
- Moran, P. A. P. (1950). "Notes on continuous stochastic phenomena". *Biometrika* 37(1-2): 17–23.
- Cliff, A. D. and Ord, J. K. (1981). *Spatial Processes*. London: Pion.

### Spatial Transcriptomics
- Ståhl, P. L., et al. (2016). "Visualization and analysis of gene expression in tissue sections by spatial transcriptomics". *Science* 353(6294): 78-82.

### FHIR Standards
- HL7 FHIR R4: https://hl7.org/fhir/R4/

---

## Contact

For questions or issues:
- GitHub Issues: https://github.com/lynnlangit/precision-medicine-mcp/issues
- Documentation: https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs

---

**Last Updated**: 2025-12-28
**Version**: 1.0.0 (Phase 1C Complete)
