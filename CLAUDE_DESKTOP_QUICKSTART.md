# Claude Desktop Integration - Quick Start Guide

**Date:** December 29, 2025
**Status:** ✅ Fully Configured

## What's Configured

You have **10 MCP servers** configured in Claude Desktop, giving you access to spatial transcriptomics analysis, genomics tools, and FHIR data integration.

### Available MCP Servers

| Server | Status | Purpose | Tools Available |
|--------|--------|---------|-----------------|
| **spatialtools** | ✅ **ACTIVE** | Spatial transcriptomics analysis | 10 tools (Phase 2 complete) |
| **epic** | ✅ **ACTIVE** | GCP Healthcare API FHIR integration | Patient data retrieval |
| fgbio | 🟡 DRY_RUN | Genomics (fgbio toolkit) | Reference-based analysis |
| openimagedata | 🟡 DRY_RUN | Public imaging datasets | Image data access |
| seqera | 🟡 DRY_RUN | Nextflow pipelines via Seqera | Workflow execution |
| huggingface | 🟡 DRY_RUN | HuggingFace models | ML model integration |
| deepcell | 🟡 DRY_RUN | DeepCell segmentation | Cell segmentation |
| mockepic | 🟡 DRY_RUN | Mock Epic FHIR | Testing FHIR workflows |
| tcga | 🟡 DRY_RUN | TCGA cancer genomics | Public cancer data |
| multiomics | 🟡 DRY_RUN | Multi-omics integration | Data fusion |

**Legend:**
- ✅ **ACTIVE** = Real data mode, fully functional
- 🟡 **DRY_RUN** = Mock mode for testing

## Spatial Analysis Tools (spatialtools server)

### Available Tools

**Phase 2A - Differential Expression:**
- `perform_differential_expression` - Compare gene expression between groups (Mann-Whitney U + FDR)

**Phase 2B - Batch Correction:**
- `perform_batch_correction` - ComBat batch correction for multi-batch data

**Phase 2C - Cell Type Deconvolution:**
- `deconvolve_cell_types` - Signature-based cell type estimation

**Phase 2D - Pathway Enrichment:**
- `perform_pathway_enrichment` - Pathway/gene set enrichment analysis

**Phase 2E - Spatial Autocorrelation:**
- `calculate_spatial_autocorrelation` - Moran's I spatial statistics

**Data Management:**
- `filter_quality` - QC filtering by expression/coverage
- `split_by_region` - Split data by tissue regions
- `align_spatial_data` - Align expression with coordinates
- `merge_tiles` - Combine multiple Visium tiles
- `get_spatial_data_for_patient` - Bridge to FHIR patient data

## How to Use in Claude Desktop

### Step 1: Restart Claude Desktop

If you just configured the servers, restart Claude Desktop to load them:

```bash
# Quit Claude Desktop (Cmd+Q)
# Then reopen it from Applications
```

### Step 2: Verify Servers are Connected

In a new Claude conversation, ask:

```
What MCP servers are available?
```

You should see all 10 servers listed.

### Step 3: Interactive Analysis Examples

#### Example 1: Analyze Patient-001 (Complete Workflow)

```
I want to analyze Patient-001 (PAT001-OVC-2025) ovarian cancer spatial data.

Please run the complete Phase 2 workflow:
1. Load the data from /Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial/
2. Perform differential expression (tumor_core vs stroma)
3. Calculate spatial autocorrelation for all genes
4. Do cell type deconvolution
5. Summarize the findings and clinical implications
```

#### Example 2: Identify Drug Resistance Markers

```
Using the spatialtools server, analyze Patient-001 to identify drug resistance markers:

1. Run differential expression comparing tumor_core vs stroma regions
2. Filter for genes with log2FC > 2
3. Check if ABCB1, PIK3CA, or AKT1 are upregulated
4. Calculate spatial autocorrelation for resistance genes
5. Tell me if resistance markers are spatially clustered
```

#### Example 3: Find Hypoxic Zones

```
Help me identify hypoxic regions in Patient-001:

1. Load the spatial data
2. Calculate spatial autocorrelation for HIF1A and CA9
3. Show me the Moran's I values
4. Are these genes spatially clustered? What does that tell us clinically?
```

#### Example 4: Batch Correction (Multi-Site Study)

```
I have spatial data from 3 different hospitals. Run batch correction:

Files:
- /path/to/hospital1_visium.csv
- /path/to/hospital2_visium.csv
- /path/to/hospital3_visium.csv

Use ComBat to harmonize the data and tell me:
- What % of variance was due to batch effects before correction?
- How much did batch correction reduce variance?
- Is the biological signal preserved?
```

#### Example 5: Compare Cell Type Composition

```
Using deconvolve_cell_types, analyze Patient-001:

1. Deconvolve cell types using standard ovarian cancer signatures
2. Compare cell type proportions between tumor_core, stroma, and immune regions
3. Are there immune-excluded regions?
4. What are the clinical implications?
```

### Step 4: Access FHIR Data (Epic Server)

```
Using the epic MCP server, retrieve Patient-001 data from the GCP Healthcare API:

1. Get the patient demographics
2. List all conditions
3. Show CA-125 observations
4. What medications has the patient received?
```

## Advanced Usage

### Custom Analysis Workflows

You can chain multiple tools together:

```
Run a complete spatial analysis for Patient-001:

1. Use get_spatial_data_for_patient to load data
2. Filter with filter_quality (min_expression=5, min_spots=10)
3. Split into regions with split_by_region
4. Run perform_differential_expression on tumor vs normal
5. Calculate pathway enrichment for upregulated genes
6. Identify spatially variable genes
7. Generate a clinical summary with treatment recommendations
```

### Exploratory Questions

Ask Claude to explore the data:

```
What can you tell me about the spatial organization of immune cells in Patient-001?
```

```
Are there any genes that show both high differential expression and strong spatial clustering?
```

```
Based on the molecular profile, what targeted therapies would you recommend?
```

## File Paths Reference

**Patient-001 Data:**
```
Expression: /Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial/visium_gene_expression.csv
Regions: /Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial/visium_region_annotations.csv
Coordinates: /Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial/visium_spatial_coordinates.csv
```

**Reference Data:**
```
General: /Users/lynnlangit/Documents/GitHub/spatial-mcp/data/
Cache: /Users/lynnlangit/Documents/GitHub/spatial-mcp/data/cache/
```

## Tips for Best Results

### 1. Be Specific About Inputs

❌ Bad:
```
Analyze some spatial data
```

✅ Good:
```
Analyze Patient-001 spatial data located at /Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial/
Compare tumor_core vs stroma regions
```

### 2. Request Interpretation

Claude can interpret the results clinically:

```
What do these differential expression results mean for treatment?
```

```
Are there any actionable drug targets based on this pathway enrichment?
```

### 3. Iterate and Refine

Start broad, then drill down:

```
1. First: "What are the top differentially expressed genes?"
2. Then: "Tell me more about the PI3K/AKT pathway genes"
3. Finally: "Are these genes spatially clustered? What does that mean?"
```

### 4. Ask for Visualizations (Descriptions)

While Claude can't create plots directly, it can describe what to visualize:

```
Based on this spatial autocorrelation data, what plots should I create?
```

## Troubleshooting

### Issue: "No MCP servers found"

**Solution:** Restart Claude Desktop (Cmd+Q, then reopen)

### Issue: "Tool execution failed"

**Check:**
1. File paths are correct (absolute paths work best)
2. Data files exist and are readable
3. Check server logs (Claude Desktop → Settings → MCP Logs)

### Issue: "Server not responding"

**Fix:**
```bash
# Check if server can start manually
cd /Users/lynnlangit/Documents/GitHub/spatial-mcp/servers/mcp-spatialtools
uv run python -m mcp_spatialtools
```

Should output: `Server running on stdio...`

### Issue: "Permission denied" for GCP

**Fix:** Ensure service account key exists:
```bash
ls -la /Users/lynnlangit/Documents/GitHub/spatial-mcp/infrastructure/mcp-server-key.json
```

## Configuration File Location

Your MCP server configuration is stored at:

```
~/Library/Application Support/Claude/claude_desktop_config.json
```

**Current Configuration:**
- ✅ 10 servers configured
- ✅ spatialtools in REAL mode (SPATIAL_DRY_RUN=false)
- ✅ epic connected to GCP (precision-medicine-poc project)
- ✅ All file paths use absolute paths

## Example Session

Here's a complete example session you can try:

```markdown
USER: Analyze Patient-001 ovarian cancer spatial data

CLAUDE: I'll analyze the Patient-001 (PAT001-OVC-2025) spatial transcriptomics data using the spatialtools MCP server. Let me run a comprehensive analysis...

[Uses perform_differential_expression tool]
[Uses calculate_spatial_autocorrelation tool]
[Uses deconvolve_cell_types tool]

**Key Findings:**
- 17 genes differentially expressed between tumor and stroma
- TP53, ABCB1, BCL2L1 highly upregulated in tumor (drug resistance markers)
- HIF1A shows strongest spatial clustering (Moran's I = 0.141) - hypoxic zones
- Tumor regions enriched for resistance signatures
- Immune cells spatially segregated in stroma

**Clinical Implications:**
- High ABCB1 expression suggests potential chemotherapy resistance
- PI3K/AKT pathway activated (AKT1, PIK3CA upregulated)
- Consider PI3K inhibitors, hypoxia-targeting agents, BCL2 inhibitors

**Treatment Recommendations:**
1. Alpelisib or Capivasertib (PI3K/AKT inhibitors)
2. Evofosfamide (hypoxia-activated prodrug)
3. Venetoclax (BCL2 inhibitor)
4. MDR reversal agents with chemotherapy

Would you like me to drill down into any specific findings?
```

## Next Steps

**Now that Claude Desktop is configured:**

1. ✅ Try the Example Session above
2. ✅ Explore Patient-001 data interactively
3. ✅ Ask questions about the spatial organization
4. ✅ Request treatment recommendations
5. ✅ Compare different tissue regions

**Future Enhancements:**
- Add visualization generation (matplotlib/plotly outputs)
- Create PDF report generation tool
- Add more patients for cohort analysis
- Enable multi-patient comparison

## Support & Documentation

**MCP Server Docs:**
- Spatialtools: `servers/mcp-spatialtools/README.md`
- Epic FHIR: `servers/mcp-epic/README.md`

**Test Data:**
- Patient-001: `data/patient-data/PAT001-OVC-2025/`

**Phase 2 Summaries:**
- Differential Expression: `servers/mcp-spatialtools/PHASE2A_*.md`
- Batch Correction: `servers/mcp-spatialtools/PHASE2B_*.md`
- Cell Deconvolution: `servers/mcp-spatialtools/PHASE2C_*.md`
- Pathway Enrichment: `servers/mcp-spatialtools/PHASE2D_*.md`
- Spatial Autocorrelation: `servers/mcp-spatialtools/PHASE2E_*.md`
- Complete Workflow: `servers/mcp-spatialtools/COMPLETE_WORKFLOW_TEST_SUMMARY.md`

---

**🎉 You're all set! Open Claude Desktop and start analyzing spatial transcriptomics data interactively!**

**First prompt to try:**
```
List the available MCP tools from the spatialtools server and give me a brief description of what each one does.
```
