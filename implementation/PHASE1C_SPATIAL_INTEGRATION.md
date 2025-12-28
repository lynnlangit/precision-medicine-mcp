# Phase 1C: Clinical-Spatial Integration

**Status:** In Progress
**Goal:** Connect clinical FHIR data → spatial transcriptomics analysis
**Timeline:** 1-2 days

---

## Current State Assessment

### ✅ What's Working (mcp-epic)
- Patient demographics retrieval (de-identified)
- Clinical conditions (ovarian cancer, stage IV HGSOC)
- Lab observations (CA-125: 487 U/mL, BRCA: negative)
- Medication history (Bevacizumab active, Carboplatin/Paclitaxel completed)

### 🔶 What's Partially Working (mcp-spatialtools)
**Real Implementation (40%):**
- `filter_quality` - QC filtering on spatial barcodes
- `split_by_region` - Spatial region segmentation
- `merge_tiles` - Tile combination

**Mocked (60%):**
- Alignment, differential expression, batch correction, etc.

### 📊 Available Data
- **Clinical:** Patient-001 in GCP Healthcare API
- **Spatial:** PAT001-OVC-2025 spatial data
  - `visium_gene_expression.csv` (900+ spots × 31 genes)
  - `visium_spatial_coordinates.csv` (coordinates)
  - `visium_region_annotations.csv` (tissue regions)

---

## Integration Strategy

### Phase 1C Goals (Achievable in 1-2 days)

1. **Create Clinical-Spatial Bridge Tool** ⭐
   - New MCP tool: `get_spatial_data_for_patient`
   - Links patient ID → spatial dataset
   - Enriches spatial analysis with clinical context

2. **Implement Real Spatial Autocorrelation**
   - Replace mocked version with scipy/squidpy
   - Calculate Moran's I for genes of interest
   - Focus on genes relevant to patient's condition

3. **Build Integrated Workflow**
   - Query patient clinical data (epic)
   - Retrieve relevant spatial dataset (spatialtools)
   - Filter by quality and region (spatialtools)
   - Analyze spatial patterns (autocorrelation)
   - Generate clinical-spatial report

4. **Demo with Patient 001**
   - Ovarian cancer patient
   - High CA-125 (tumor marker)
   - Platinum-resistant HGSOC
   - Analyze: Ki67 (proliferation), CD8A (immune), VIM (mesenchymal)

---

## Implementation Plan

### Task 1: Clinical-Spatial Bridge (2-3 hours)

**Create new tool in mcp-spatialtools:**

```python
@mcp.tool()
async def get_spatial_data_for_patient(
    patient_id: str,
    tissue_type: str = "tumor",
    include_clinical_context: bool = True
) -> Dict[str, Any]:
    """
    Retrieve spatial transcriptomics data for a patient with clinical context.

    Links:
    - Patient clinical data from Epic FHIR
    - Corresponding spatial dataset
    - Enriches with condition, biomarkers, treatments

    Args:
        patient_id: Patient identifier from FHIR
        tissue_type: Type of tissue ("tumor", "normal", "margin")
        include_clinical_context: Include clinical metadata

    Returns:
        spatial_data_path: Path to expression matrix
        clinical_summary: Patient condition, biomarkers, treatments
        genes_of_interest: Genes relevant to patient's condition
        suggested_analyses: Recommended spatial analyses
    """
```

**Mapping Logic:**
```python
# Map patient conditions → genes of interest
CONDITION_GENE_MAP = {
    "ovarian cancer": ["Ki67", "TP53", "BRCA1", "EPCAM"],
    "HGSOC": ["TP53", "Ki67", "FOXM1", "MYC"],
    "platinum-resistant": ["ABCB1", "ERCC1", "GSTP1"]
}

# Map treatments → biomarkers
TREATMENT_BIOMARKER_MAP = {
    "Bevacizumab": ["VEGFA", "CD31", "HIF1A"],  # Anti-angiogenic
    "Carboplatin": ["ERCC1", "XPA"],  # DNA repair
}
```

---

### Task 2: Spatial Autocorrelation (3-4 hours)

**Replace mocked implementation with real analysis:**

```python
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr

def calculate_morans_i(
    expression_values: np.ndarray,
    coordinates: np.ndarray,
    distance_threshold: float = 100.0
) -> float:
    """
    Calculate Moran's I for spatial autocorrelation.

    Moran's I ranges from -1 (dispersed) to +1 (clustered)
    0 indicates random spatial pattern
    """
    # Build spatial weights matrix (neighbors within threshold)
    distances = cdist(coordinates, coordinates)
    weights = (distances < distance_threshold).astype(float)
    np.fill_diagonal(weights, 0)

    # Normalize weights
    row_sums = weights.sum(axis=1)
    weights = weights / row_sums[:, np.newaxis]

    # Calculate Moran's I
    n = len(expression_values)
    mean_expr = expression_values.mean()
    deviations = expression_values - mean_expr

    numerator = np.sum(weights * np.outer(deviations, deviations))
    denominator = np.sum(deviations ** 2)

    morans_i = (n / weights.sum()) * (numerator / denominator)

    return morans_i
```

**Enhanced tool:**
```python
@mcp.tool()
async def calculate_spatial_autocorrelation(
    expression_file: str,
    genes: List[str],
    coordinates_file: str,
    method: str = "morans_i",
    distance_threshold: float = 100.0
) -> Dict[str, Any]:
    """
    Calculate spatial autocorrelation for genes.

    REAL IMPLEMENTATION - uses scipy for Moran's I calculation
    """
```

---

### Task 3: Integrated Workflow (2-3 hours)

**Create workflow orchestrator:**

```python
# New file: servers/mcp-spatialtools/src/mcp_spatialtools/workflows.py

async def analyze_patient_spatial_context(
    patient_id: str,
    epic_client: FHIRClient,
    spatial_data_dir: str
) -> Dict[str, Any]:
    """
    Complete clinical-spatial analysis workflow.

    Steps:
    1. Retrieve clinical data from Epic
    2. Identify genes of interest from conditions/treatments
    3. Load spatial data
    4. Filter quality (real implementation)
    5. Split by region (real implementation)
    6. Calculate spatial autocorrelation (real implementation)
    7. Generate integrated report
    """

    # Step 1: Get clinical context
    demographics = await epic_client.get("Patient/{patient_id}")
    conditions = await epic_client.search("Condition", {"patient": patient_id})
    medications = await epic_client.search("MedicationStatement", {"patient": patient_id})
    observations = await epic_client.search("Observation", {"patient": patient_id})

    # Step 2: Identify genes of interest
    genes_of_interest = identify_relevant_genes(
        conditions, medications, observations
    )

    # Step 3-6: Spatial analysis
    spatial_results = await run_spatial_analysis(
        spatial_data_dir, genes_of_interest
    )

    # Step 7: Integrate
    report = generate_clinical_spatial_report(
        clinical_data={
            "demographics": demographics,
            "conditions": conditions,
            "medications": medications,
            "observations": observations
        },
        spatial_results=spatial_results
    )

    return report
```

---

### Task 4: Demo with Patient 001 (1 hour)

**Test workflow:**

1. **Query clinical data:**
   ```
   Get patient demographics, conditions, and medications for patient-001 using epic
   ```

2. **Load spatial data:**
   ```
   Get spatial data for patient PAT001-OVC-2025 with clinical context
   ```

3. **Run spatial analysis:**
   ```
   Calculate spatial autocorrelation for Ki67, CD8A, and VIM in tumor regions
   ```

4. **Generate report:**
   ```
   Generate clinical-spatial analysis report for patient-001
   ```

**Expected Output:**
```
Clinical Context:
- Stage IV HGSOC (platinum-resistant)
- CA-125: 487 U/mL (elevated)
- Currently on Bevacizumab (anti-angiogenic)
- BRCA: negative

Spatial Analysis:
- Ki67 (proliferation): Moran's I = 0.65 (strong clustering in tumor core)
- CD8A (T-cells): Moran's I = 0.42 (moderate clustering at margins)
- VIM (mesenchymal): Moran's I = 0.78 (strong clustering, invasive front)

Clinical-Spatial Insights:
- High proliferation clusters align with platinum-resistance
- Immune infiltration (CD8A) at tumor margins suggests bevacizumab response
- Mesenchymal signature (VIM) indicates aggressive phenotype
```

---

## Success Criteria

### ✅ Phase 1C Complete When:

1. **Bridge tool works:**
   - Can retrieve spatial data for patient ID
   - Enriches with clinical context
   - Suggests relevant genes

2. **Real spatial autocorrelation:**
   - Calculates actual Moran's I (not mocked)
   - Works with real coordinates and expression
   - Returns meaningful clustering results

3. **End-to-end workflow:**
   - Epic → Spatial seamlessly
   - Can run complete analysis through Claude Desktop
   - Generates integrated clinical-spatial report

4. **Demo successful:**
   - Works with Patient 001
   - Provides actionable insights
   - Demonstrates clinical value

---

## Timeline

| Task | Duration | Dependencies |
|------|----------|--------------|
| Design bridge tool | 1 hour | - |
| Implement bridge tool | 2 hours | Design |
| Implement Moran's I | 3 hours | - |
| Test spatial autocorrelation | 1 hour | Moran's I |
| Create workflow orchestrator | 2 hours | Bridge tool |
| Test integrated workflow | 1 hour | All above |
| Demo with Patient 001 | 1 hour | Workflow |
| Documentation | 1 hour | - |

**Total: 12 hours (1.5 days)**

---

## Next Steps

1. Create clinical-spatial bridge tool
2. Implement real Moran's I calculation
3. Build workflow orchestrator
4. Test with Patient 001
5. Document the pipeline

---

**Ready to begin?** Let's start with Task 1: Clinical-Spatial Bridge Tool.
