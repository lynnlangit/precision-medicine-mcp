# MCP Server Inventory

**Last Updated:** January 10, 2026

All 10 MCP servers with accurate tool counts and implementation status.

---

## Production Servers (4/9)

### 1. mcp-spatialtools
**Status:** ✅ Production (95% real)
**Deployed:** GCP Cloud Run (Jan 9, 2026)
**Revision:** mcp-spatialtools-00005-r4s
**Tools:** 14 (10 analysis + 4 visualization)

**Analysis Tools:**
1. `filter_quality` - Quality filtering of spatial barcodes
2. `split_by_region` - Segment data by tissue regions
3. `align_spatial_data` - STAR alignment (requires genome indices)
4. `merge_tiles` - Combine multiple spatial tiles
5. `calculate_spatial_autocorrelation` - Moran's I, Geary's C statistics
6. `perform_differential_expression` - Wilcoxon, t-test, DESeq2-style
7. `perform_batch_correction` - ComBat, Harmony, Scanorama
8. `perform_pathway_enrichment` - GO, KEGG, Reactome, Hallmark
9. `deconvolve_cell_types` - CIBERSORTx, EPIC, quanTIseq (synthetic)
10. `get_spatial_data_for_patient` - Bridge tool for multiomics integration

**Visualization Tools:**
11. `generate_spatial_heatmap` - Gene expression on tissue coordinates
12. `generate_gene_expression_heatmap` - Clustered heatmap (genes × regions)
13. `generate_region_composition_chart` - Bar chart of spot counts
14. `visualize_spatial_autocorrelation` - Moran's I visualization

**Real vs Mock:**
- ✅ Real: All analysis tools (1-8, 10), all visualization tools (11-14)
- ❌ Mock: Cell type deconvolution (9) - returns synthetic signatures

**URL:** https://mcp-spatialtools-ondu7mwjpa-uc.a.run.app
**README:** [mcp-spatialtools/README.md](../../servers/mcp-spatialtools/README.md)

---

### 2. mcp-fgbio
**Status:** ✅ Production (100% real)
**Deployed:** GCP Cloud Run (Dec 30, 2025)
**Tools:** 4

1. `validate_fastq` - FASTQ quality validation using FGbio
2. `extract_umis` - UMI extraction and processing
3. `filter_vcf` - VCF variant filtering
4. `get_reference_data` - Genomic reference data access

**URL:** https://mcp-fgbio-ondu7mwjpa-uc.a.run.app
**README:** [mcp-fgbio/README.md](../../servers/mcp-fgbio/README.md)

---

### 3. mcp-multiomics
**Status:** ✅ Production (100% real)
**Deployed:** GCP Cloud Run (Dec 30, 2025)
**Tools:** 9

1. `integrate_rna_protein` - HAllA-based RNA-protein integration
2. `perform_stouffer_meta_analysis` - Stouffer's Z-score meta-analysis
3. `identify_upstream_regulators` - Infer upstream regulators
4. `calculate_correlation_networks` - Build correlation networks
5. `perform_pathway_analysis` - Multi-omics pathway enrichment
6. `identify_resistance_signatures` - Drug resistance signatures
7. `compare_treatment_responses` - Compare resistant vs sensitive
8. `generate_integration_report` - Comprehensive integration report
9. `get_patient_multiomics_summary` - Bridge tool for integration

**URL:** https://mcp-multiomics-ondu7mwjpa-uc.a.run.app
**README:** [mcp-multiomics/README.md](../../servers/mcp-multiomics/README.md)

---

### 4. mcp-epic
**Status:** ✅ Production (100% real)
**Deployed:** GCP Cloud Run (Dec 30, 2025)
**Tools:** 7

**Note:** This is the REAL Epic FHIR server with OAuth 2.0 and de-identification. NOT mcp-mockepic.

1. `fetch_patient_demographics` - Epic FHIR patient resources
2. `query_lab_results` - Lab observations (CA-125, etc.)
3. `get_medication_history` - Medication statements
4. `fetch_condition_history` - Diagnosis/condition records
5. `query_procedures` - Procedure records
6. `get_allergy_intolerances` - Allergy records
7. `deidentify_fhir_resource` - HIPAA Safe Harbor de-identification

**URL:** (Private - requires Epic credentials)
**README:** [mcp-epic/README.md](../../servers/mcp-epic/README.md)

---

## Partially Real Servers (1/9)

### 5. mcp-openimagedata
**Status:** ✅ 60% Real
**Deployed:** GCP Cloud Run (Jan 9, 2026)
**Revision:** mcp-openimagedata-00004-vks
**Tools:** 5 (3 analysis + 2 visualization)

**Analysis Tools:**
1. `fetch_histology_image` - ✅ Real - Load H&E/IF images
2. `register_image_to_spatial` - ❌ Mock - Image registration (future)
3. `extract_image_features` - ❌ Mock - Feature extraction (future)

**Visualization Tools:**
4. `generate_multiplex_composite` - ✅ Real - RGB MxIF composites
5. `generate_he_annotation` - ✅ Real - Annotated H&E morphology

**Real vs Mock:**
- ✅ Real: Image loading, visualization tools (60%)
- ❌ Mock: Registration, feature extraction (40%)

**URL:** https://mcp-openimagedata-ondu7mwjpa-uc.a.run.app
**README:** [mcp-openimagedata/README.md](../../servers/mcp-openimagedata/README.md)

---

## Mock Servers (4/9)

### 6. mcp-deepcell
**Status:** ❌ Mocked (0% real)
**Deployed:** Not deployed
**Tools:** 4 (all mocked)

1. `segment_cells` - Mock: Returns synthetic segmentation masks
2. `classify_cell_states` - Mock: Random cell state classifications
3. `generate_segmentation_overlay` - Mock: Placeholder overlays
4. `generate_phenotype_visualization` - Mock: Synthetic phenotype maps

**Future:** Will use DeepCell-TF library for real MxIF segmentation
**README:** [mcp-deepcell/README.md](../../servers/mcp-deepcell/README.md)

---

### 7. mcp-tcga
**Status:** ❌ Mocked (0% real)
**Deployed:** Not deployed
**Tools:** 7 (all mocked)

1. `query_tcga_samples` - Mock: Returns synthetic TCGA metadata
2. `fetch_tcga_expression` - Mock: Random expression matrices
3. `get_mutation_data` - Mock: Synthetic mutation profiles
4. `fetch_clinical_annotations` - Mock: Random clinical data
5. `compare_to_tcga` - Mock: Placeholder comparisons
6. `get_survival_data` - Mock: Synthetic survival curves
7. `fetch_subtype_classification` - Mock: Random subtypes

**Future:** Will integrate with GDC Data Portal API
**README:** [mcp-tcga/README.md](../../servers/mcp-tcga/README.md)

---

### 8. mcp-seqera
**Status:** ❌ Mocked (0% real)
**Deployed:** Not deployed
**Tools:** 6 (all mocked)

1. `launch_nextflow_pipeline` - Mock: Simulated pipeline launches
2. `monitor_workflow_status` - Mock: Fake status updates
3. `fetch_workflow_results` - Mock: Placeholder results
4. `list_available_pipelines` - Mock: Synthetic pipeline list
5. `submit_batch_analysis` - Mock: Batch simulation
6. `get_workflow_metrics` - Mock: Random metrics

**Future:** Will integrate with Seqera Platform API
**README:** [mcp-seqera/README.md](../../servers/mcp-seqera/README.md)

---

### 9. mcp-huggingface
**Status:** ❌ Mocked (0% real)
**Deployed:** Not deployed
**Tools:** 5 (all mocked)

1. `load_genomic_model` - Mock: Fake model loading
2. `predict_cell_type` - Mock: Random cell type predictions
3. `embed_sequences` - Mock: Synthetic embeddings
4. `predict_function` - Mock: Random function predictions
5. `search_models` - Mock: Placeholder model search

**Future:** Will use transformers library + HuggingFace Hub API
**README:** [mcp-huggingface/README.md](../../servers/mcp-huggingface/README.md)

---

## Mock By Design (1/9)

### 10. mcp-mockepic
**Status:** ✅ Mock by design (intentional)
**Deployed:** GCP Cloud Run (Dec 30, 2025)
**Tools:** 6 (all synthetic FHIR)

1. `fetch_patient_demographics` - Synthetic patient data
2. `query_lab_results` - Synthetic lab values
3. `get_medication_history` - Synthetic prescriptions
4. `fetch_condition_history` - Synthetic diagnoses
5. `query_procedures` - Synthetic procedures
6. `get_allergy_intolerances` - Synthetic allergies

**Purpose:** Testing and development without real PHI
**Data:** Synthea-generated FHIR-compliant synthetic data
**README:** [mcp-mockepic/README.md](../../servers/mcp-mockepic/README.md)

---

## Summary Table

| Server | Tools | Real % | Deployed | URL Available |
|--------|-------|--------|----------|---------------|
| mcp-spatialtools | 14 | 95% | ✅ | ✅ |
| mcp-fgbio | 4 | 100% | ✅ | ✅ |
| mcp-multiomics | 9 | 100% | ✅ | ✅ |
| mcp-epic | 7 | 100% | ✅ | 🔒 Private |
| mcp-openimagedata | 5 | 60% | ✅ | ✅ |
| mcp-deepcell | 4 | 0% | ❌ | ❌ |
| mcp-tcga | 7 | 0% | ❌ | ❌ |
| mcp-seqera | 6 | 0% | ❌ | ❌ |
| mcp-huggingface | 5 | 0% | ❌ | ❌ |
| mcp-mockepic | 6 | Mock | ✅ | ✅ |

**Total Tools:** 67 tools across 10 servers
**Production Ready:** 4 servers (+ 1 partial)
**Mocked:** 4 servers
**Mock by Design:** 1 server

---

## Deployment Status

### GCP Cloud Run (Production)
- **Region:** us-central1
- **Transport:** SSE (Server-Sent Events over HTTPS)
- **Mode:** Development (public access with `--allow-unauthenticated`)
- **Resources:** 2-4 Gi memory, 1-2 CPU per server

### Deployed Servers (5 total)
1. mcp-spatialtools (Jan 9, 2026) - Revision 00005
2. mcp-openimagedata (Jan 9, 2026) - Revision 00004
3. mcp-fgbio (Dec 30, 2025)
4. mcp-multiomics (Dec 30, 2025)
5. mcp-mockepic (Dec 30, 2025)

**Epic (mcp-epic):** Deployed privately with OAuth 2.0 authentication

---

## Future Deployment Plan

### Q1 2026
- [ ] Deploy mcp-deepcell with real DeepCell-TF integration
- [ ] Deploy mcp-tcga with GDC Data Portal API

### Q2 2026
- [ ] Deploy mcp-seqera with Seqera Platform integration
- [ ] Deploy mcp-huggingface with HF transformers library

### Q3 2026
- [ ] Migrate all servers to production mode (authenticated access)
- [ ] Implement VPC networking for secure communication
- [ ] Add HIPAA compliance features

---

## See Also
- [DEPLOYMENT.md](DEPLOYMENT.md) - Deployment procedures and configuration
- [OVERVIEW.md](OVERVIEW.md) - System architecture and design
- [Main Architecture README](../README.md) - Overall precision medicine architecture
