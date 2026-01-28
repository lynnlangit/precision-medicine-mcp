# Precision Medicine MCP System - Quick Summary

**For sharing via email or quick reference**

---

## System Overview

**A production-ready AI-orchestrated platform for precision oncology research**, integrating clinical (FHIR), genomic, multi-omics, spatial transcriptomics, and imaging data through 10 specialized MCP servers (63 bioinformatics tools) deployed on GCP Cloud Run. The system uses Claude or Gemini to orchestrate complex multi-step analyses via natural language, reducing analysis time from 40 hours of manual bioinformatics work to 35 minutes of AI-orchestrated execution. Includes HIPAA compliance (de-identification, 10-year audit logs, bias detection), supports direct Google Cloud Storage file access, and achieves 95% cost savings (~$1-2 per analysis vs. $3,200 traditional). **PatientOne** demonstrates the complete workflow: Stage IV ovarian cancer analysis integrating clinical data, somatic variants, multi-omics (RNA/protein/phospho), spatial transcriptomics, and H&E imaging—with treatment response prediction via GEARS perturbation modeling.

**Technical implementation**: Multi-provider Streamlit UI supporting Claude (Anthropic's native MCP integration) and Gemini (custom SSE-based MCP client with manual tool orchestration), with 5 production servers (fgbio, multiomics, spatialtools, perturbation, epic) and 5 mock/partial servers. The system handles GCS folder URIs (loads all files automatically), supports up to 30 tool-calling iterations for complex workflows, and includes comprehensive documentation for hospitals, researchers, developers, educators, patients, and funders. Validated against actual GCP deployment with 2026 pricing: ~$0.02-0.21 per Cloud Run analysis, free tier covers ~83 hours/month of testing.

---

## MCP Servers (10 Total)

### Production Servers (5)
1. **mcp-fgbio** - Genomic reference validation
2. **mcp-multiomics** - Multi-omics data integration
3. **mcp-spatialtools** - Spatial transcriptomics analysis
4. **mcp-perturbation** - Treatment response prediction
5. **mcp-epic** - Clinical FHIR data

### Mock/Partial Servers (5)
6. **mcp-tcga** - Cancer genomics data
7. **mcp-openimagedata** - Medical image analysis
8. **mcp-seqera** - Workflow automation platform
9. **mcp-huggingface** - AI model inference
10. **mcp-deepcell** - Cell image segmentation

**Total:** 63 bioinformatics tools across all servers

---

## Key Metrics

- **Time Reduction:** 40 hours → 35 minutes (95% faster)
- **Cost Reduction:** $3,200 → $1-2 per analysis (95% savings)
- **Data Integration:** 5 modalities (clinical, genomic, multi-omics, spatial, imaging)
- **HIPAA Compliant:** De-identification, 10-year audit logs, VPC isolation
- **Deployment:** GCP Cloud Run (us-central1)
- **UI:** Streamlit with Claude & Gemini support

---

## Live Demo

- **Streamlit UI:** https://streamlit-mcp-chat-ondu7mwjpa-uc.a.run.app
- **PatientOne Scenario:** Stage IV ovarian cancer complete analysis
- **Test Data:** Available in GCS at `gs://sample-inputs-patientone/`

---

## Documentation

- **For Hospitals:** [docs/for-hospitals/](for-hospitals/)
- **For Researchers:** [docs/for-researchers/](for-researchers/)
- **For Developers:** [docs/for-developers/](for-developers/)
- **Executive Summary:** [docs/EXECUTIVE_SUMMARY.md](EXECUTIVE_SUMMARY.md)
- **Architecture:** [docs/architecture/](architecture/)

---

**Repository:** https://github.com/lynnlangit/precision-medicine-mcp

**Status:** Production-ready, validated deployment (January 2026)

**License:** Apache 2.0
