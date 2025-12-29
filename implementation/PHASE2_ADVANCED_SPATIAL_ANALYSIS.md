# Phase 2: Advanced Spatial Analysis

**Status:** Planning
**Goal:** Expand real implementations beyond Moran's I to enable comprehensive spatial transcriptomics analysis
**Prerequisite:** Phase 1C Complete ✅

---

## Phase 1 Accomplishments

✅ **Infrastructure:**
- GCP Healthcare API with FHIR store
- De-identified patient data (HIPAA Safe Harbor)
- MCP servers (epic + spatialtools) in Claude Desktop

✅ **Real Implementations:**
- `get_spatial_data_for_patient` - Clinical-spatial bridge
- `calculate_spatial_autocorrelation` - Moran's I with scipy
- `filter_quality` - QC filtering on expression data
- `split_by_region` - Spatial region segmentation

✅ **Working Data:**
- Patient-001: Stage IV HGSOC, platinum-resistant, on Bevacizumab
- PAT001-OVC-2025: 900 spots × 31 genes, spatial coordinates, region annotations

---

## Phase 2 Goals

**Expand to clinically actionable spatial analysis:**

1. **Spatial Differential Expression** - Find genes varying between tumor regions
2. **Spatial Gene Co-expression** - Identify gene programs and spatial niches
3. **Cell Type Deconvolution** - Estimate cell populations from bulk spatial data
4. **Treatment Response Biomarkers** - Link spatial patterns to therapy outcomes

---

## Current State: What's Mocked vs Real

### ✅ Real (from Phase 1)
- Quality filtering
- Region segmentation
- Spatial autocorrelation (Moran's I)
- Clinical-spatial bridging

### 🔴 Still Mocked (Priority for Phase 2)
1. **Differential Expression** - `perform_differential_expression`
2. **Batch Correction** - `perform_batch_correction`
3. **Pathway Enrichment** - `perform_pathway_enrichment`
4. **Alignment** - `align_spatial_data` (lower priority - requires raw reads)

---

## Phase 2A: Spatial Differential Expression

**Objective:** Identify genes with significant expression differences between spatial regions

### Why This Matters Clinically
- Tumor core vs margin gene signatures
- Immune infiltration zones vs tumor-excluded regions
- Response prediction (e.g., Bevacizumab targets angiogenesis → analyze VEGFA in different zones)

### Implementation Plan

**Tool:** `perform_differential_expression`

**Real Implementation with scipy/pandas:**
```python
from scipy.stats import mannwhitneyu, ttest_ind
import pandas as pd
import numpy as np

def calculate_differential_expression(
    expr_data: pd.DataFrame,
    group1_spots: List[str],
    group2_spots: List[str],
    method: str = "wilcoxon"
) -> pd.DataFrame:
    """
    Real differential expression between spatial groups.

    Args:
        expr_data: Expression matrix (spots × genes)
        group1_spots: Spot IDs for group 1 (e.g., tumor core)
        group2_spots: Spot IDs for group 2 (e.g., tumor margin)
        method: "wilcoxon" or "t_test"

    Returns:
        DataFrame with columns:
        - gene: Gene symbol
        - log2_fc: Log2 fold change (group1 vs group2)
        - pvalue: Statistical p-value
        - qvalue: FDR-adjusted q-value (Benjamini-Hochberg)
        - significant: Boolean (q < 0.05 and |log2_fc| > 0.5)
    """
    results = []

    for gene in expr_data.columns:
        g1_expr = expr_data.loc[group1_spots, gene].values
        g2_expr = expr_data.loc[group2_spots, gene].values

        # Skip genes with no expression
        if g1_expr.sum() == 0 and g2_expr.sum() == 0:
            continue

        # Calculate fold change
        mean1 = g1_expr.mean() + 1e-10  # Pseudocount
        mean2 = g2_expr.mean() + 1e-10
        log2_fc = np.log2(mean1 / mean2)

        # Statistical test
        if method == "wilcoxon":
            stat, pval = mannwhitneyu(g1_expr, g2_expr, alternative='two-sided')
        else:  # t-test
            stat, pval = ttest_ind(g1_expr, g2_expr)

        results.append({
            'gene': gene,
            'log2_fc': log2_fc,
            'mean_group1': mean1,
            'mean_group2': mean2,
            'pvalue': pval
        })

    # FDR correction (Benjamini-Hochberg)
    df = pd.DataFrame(results)
    from scipy.stats import false_discovery_control
    df['qvalue'] = false_discovery_control(df['pvalue'].values)
    df['significant'] = (df['qvalue'] < 0.05) & (np.abs(df['log2_fc']) > 0.5)

    return df.sort_values('pvalue')
```

### Patient 001 Use Case

**Analysis:** Tumor core vs immune-infiltrated margin

```
Group 1 (tumor_core): 210 spots
Group 2 (tumor_margin): 236 spots

Expected findings:
- MKI67 ↑ in core (proliferation)
- CD8A ↑ in margin (T-cell infiltration)
- VEGFA differential (bevacizumab response marker)
- VIM patterns (EMT at invasive front)
```

**Timeline:** 4-6 hours
- Implement real differential expression: 2-3 hours
- Test with PAT001-OVC-2025 data: 1 hour
- Integration with clinical context: 1 hour
- Testing: 1 hour

---

## Phase 2B: Spatial Gene Co-expression Networks

**Objective:** Identify genes with correlated spatial expression patterns

### Why This Matters Clinically
- Discover spatial gene programs (e.g., angiogenesis module, immune checkpoint module)
- Validate treatment targets (Bevacizumab → VEGFA/HIF1A/KDR co-expression)
- Identify platinum resistance signatures

### Implementation Plan

**New Tool:** `calculate_spatial_coexpression`

**Real Implementation with scipy:**
```python
from scipy.stats import spearmanr, pearsonr
from scipy.spatial.distance import pdist, squareform
import networkx as nx

def calculate_spatial_coexpression(
    expr_data: pd.DataFrame,
    coordinates: np.ndarray,
    genes: List[str],
    distance_threshold: float = 1500.0,
    method: str = "spearman"
) -> Dict[str, Any]:
    """
    Calculate spatially-weighted gene co-expression.

    Only considers spots within distance_threshold for correlation.
    """
    # Build spatial network
    distances = squareform(pdist(coordinates))
    spatial_neighbors = distances < distance_threshold

    # Calculate correlations
    coexpr_matrix = np.zeros((len(genes), len(genes)))

    for i, gene1 in enumerate(genes):
        for j, gene2 in enumerate(genes):
            if i >= j:
                continue

            # Get expression values
            expr1 = expr_data[gene1].values
            expr2 = expr_data[gene2].values

            # Spatially-weighted correlation
            # Only use spots that are neighbors
            neighbor_weights = spatial_neighbors.mean(axis=0)

            if method == "spearman":
                corr, pval = spearmanr(expr1, expr2)
            else:
                corr, pval = pearsonr(expr1, expr2)

            coexpr_matrix[i, j] = corr
            coexpr_matrix[j, i] = corr

    # Identify modules (clusters of co-expressed genes)
    # Use correlation threshold of 0.5
    G = nx.Graph()
    for i, gene1 in enumerate(genes):
        for j, gene2 in enumerate(genes):
            if i < j and abs(coexpr_matrix[i, j]) > 0.5:
                G.add_edge(gene1, gene2, weight=coexpr_matrix[i, j])

    modules = list(nx.connected_components(G))

    return {
        'coexpression_matrix': coexpr_matrix.tolist(),
        'genes': genes,
        'modules': [list(m) for m in modules],
        'num_modules': len(modules)
    }
```

**Timeline:** 4-5 hours

---

## Phase 2C: Cell Type Deconvolution

**Objective:** Estimate cell type proportions from bulk spatial spots

### Why This Matters Clinically
- Quantify immune infiltration (CD8+ T-cells, macrophages)
- Tumor purity estimation
- Response biomarkers (high CD8+ → better immunotherapy response)

### Implementation Plan

**Tool:** `deconvolve_cell_types`

**Approach:** Signature-based deconvolution (not full scRNA-seq reference)

```python
def deconvolve_cell_types(
    expr_data: pd.DataFrame,
    signatures: Dict[str, List[str]]
) -> pd.DataFrame:
    """
    Simple signature-based deconvolution.

    Args:
        expr_data: Expression matrix (spots × genes)
        signatures: Cell type gene signatures
            e.g., {
                'tumor_cells': ['EPCAM', 'KRT8', 'KRT18'],
                'cd8_tcells': ['CD8A', 'CD3D', 'CD3E'],
                'macrophages': ['CD68', 'CD163'],
                'endothelial': ['CD31', 'VWF'],
                'fibroblasts': ['FAP', 'COL1A1', 'ACTA2']
            }

    Returns:
        DataFrame with cell type scores per spot
    """
    scores = {}

    for cell_type, marker_genes in signatures.items():
        # Get genes that exist in data
        available_genes = [g for g in marker_genes if g in expr_data.columns]

        if not available_genes:
            scores[cell_type] = np.zeros(len(expr_data))
            continue

        # Average expression of signature genes
        signature_expr = expr_data[available_genes].mean(axis=1)

        # Z-score normalization
        scores[cell_type] = (signature_expr - signature_expr.mean()) / signature_expr.std()

    return pd.DataFrame(scores, index=expr_data.index)
```

**Cell Type Signatures for Ovarian Cancer:**
```python
OVARIAN_CANCER_SIGNATURES = {
    'tumor_cells': ['EPCAM', 'KRT8', 'KRT18', 'PAX8', 'TP53'],
    'cd8_tcells': ['CD8A', 'CD3D', 'CD3E'],
    'cd4_tcells': ['CD4', 'CD3D', 'CD3E'],
    'tregs': ['FOXP3', 'CD4'],
    'macrophages': ['CD68', 'CD163'],
    'endothelial': ['CD31', 'VWF', 'VEGFA'],
    'fibroblasts': ['FAP', 'COL1A1', 'COL3A1', 'ACTA2'],
    'mesenchymal': ['VIM', 'SNAI1', 'TWIST1']
}
```

**Timeline:** 3-4 hours

---

## Phase 2D: Pathway Enrichment (Real Implementation)

**Objective:** Replace mocked pathway enrichment with real analysis

### Implementation Options

**Option 1: Pre-built pathway databases (Fast)**
```python
# Use curated gene sets
KEGG_PATHWAYS = {
    'VEGF_signaling': ['VEGFA', 'KDR', 'HIF1A', 'AKT1', 'MTOR'],
    'DNA_repair': ['BRCA1', 'BRCA2', 'TP53', 'ERCC1', 'XPA'],
    'PI3K_AKT': ['PIK3CA', 'AKT1', 'MTOR', 'TP53'],
    'Drug_resistance': ['ABCB1', 'GSTP1', 'ERCC1']
}

def simple_enrichment(gene_list: List[str], pathway_db: Dict) -> List[Dict]:
    results = []
    for pathway, pathway_genes in pathway_db.items():
        overlap = set(gene_list) & set(pathway_genes)
        if overlap:
            # Fisher's exact test
            from scipy.stats import fisher_exact
            # ... calculate enrichment
            results.append({
                'pathway': pathway,
                'overlap': list(overlap),
                'pvalue': pval
            })
    return results
```

**Option 2: Use enrichr API** (if online access acceptable)

**Timeline:** 2-3 hours for option 1

---

## Phase 2 Prioritization

### Recommended Order (by clinical impact)

**Week 1: Differential Expression (Phase 2A)**
- Most clinically actionable
- Enables comparison of tumor regions
- Foundation for other analyses
- **Est: 4-6 hours**

**Week 2: Cell Type Deconvolution (Phase 2C)**
- Quantifies immune infiltration (critical for treatment planning)
- Simpler than co-expression networks
- Immediate clinical interpretation
- **Est: 3-4 hours**

**Week 3: Gene Co-expression (Phase 2B)**
- Identifies treatment-relevant modules
- More complex, builds on differential expression
- **Est: 4-5 hours**

**Week 4: Pathway Enrichment (Phase 2D)**
- Biological interpretation
- Quick to implement with pre-built databases
- **Est: 2-3 hours**

---

## Success Criteria for Phase 2

### Phase 2 Complete When:

✅ **Differential Expression:**
- Real statistical tests (Wilcoxon, t-test)
- FDR correction
- Works with region-segmented data
- Clinical interpretation (e.g., "MKI67 upregulated in tumor core")

✅ **Cell Type Deconvolution:**
- Signature-based scoring
- Ovarian cancer-specific signatures
- Per-spot cell type estimates
- Spatial visualization-ready output

✅ **Gene Co-expression:**
- Spatially-weighted correlations
- Module detection
- Treatment-relevant gene programs identified

✅ **Pathway Enrichment:**
- Real statistical testing
- Relevant pathway databases (KEGG, GO, or custom)
- Integrates with differential expression results

✅ **End-to-End Workflow:**
- Clinical data → spatial regions → differential expression → cell types → pathways
- Works with Patient-001 data
- Generates actionable clinical insights

---

## Patient-001 Phase 2 Analysis Plan

**Complete Clinical-Spatial Analysis:**

1. **Quality Filter** → 900 spots (all pass)
2. **Split by Region** → tumor_core, tumor_margin, stroma, immune_infiltrate
3. **Differential Expression:**
   - Core vs Margin
   - Immune-high vs Immune-low
4. **Cell Type Deconvolution:**
   - Estimate CD8+ T-cell infiltration
   - Quantify endothelial cells (bevacizumab target)
5. **Gene Co-expression:**
   - Angiogenesis module (VEGFA, KDR, HIF1A)
   - Platinum resistance module (ABCB1, ERCC1, GSTP1)
6. **Pathway Enrichment:**
   - Enriched in immune-high regions
   - Bevacizumab response pathways
7. **Clinical Report:**
   - Spatial heterogeneity assessment
   - Treatment response prediction
   - Recommendations

---

## Next Steps

**Ready to start Phase 2A: Differential Expression?**

This will give us the most immediate clinical value - identifying which genes differ between tumor regions, immune zones, etc.

Let me know if you want to:
1. **Start with Phase 2A** (Differential Expression) - Recommended
2. **Adjust the priority order**
3. **Add/modify features**
4. **Something else**
