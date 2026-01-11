# Precision Medicine w/ AI-orchestrated MCP Servers

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![MCP](https://img.shields.io/badge/MCP-2025--06--18-green.svg)](https://modelcontextprotocol.io/)
[![Claude Desktop](https://img.shields.io/badge/Claude-Desktop-orange.svg)](https://claude.ai/download)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

<img src="https://github.com/lynnlangit/precision-medicine-mcp/blob/main/data/images/repo-image.png">

## What and Why
- The Problem: multi-modal precision medicine is siloed and code-heavy
- This Solution: a set of custom MCP servers for analysis and data retrieval (PatientOne use case example):
  - Makes complex analysis significantly faster and easier by presenting a customized, single interface for tools and data
  - System coordinates disparate servers and stitches results together after being given a natural-language prompt
  - "What is an MCP Server?" [(article)](https://medium.com/@elisowski/mcp-explained-the-new-standard-connecting-ai-to-everything-79c5a1c98288)
- What it is NOT: Not clinically validated yet
- Who it’s for: Researchers, Clinicians, Platform Builders, Workflow Architects
- See it / Try it: <5 minute local demo - [recording](https://www.youtube.com/watch?v=LUldOHHX5Yo) | [code](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/tests/manual_testing/PatientOne-OvarianCancer)

---

## Featured Use Case: PatientOne

<kbd><img src="https://github.com/lynnlangit/precision-medicine-mcp/blob/main/tests/manual_testing/PatientOne-OvarianCancer/architecture/patient-one-holistic.png" width=800></kbd>  

- ANALYZE complete patient profiles using **natural language**
- DEMONSTRATE end-to-end precision medicine workflows
  - Using patient example of Ovarian Cancer (Stage IV HGSOC, platinum-resistant, BRCA1 mutation)
  - Extensible for other comorbidities
- USE 10 MCP servers (9 deployed + mcp-epic local) with 55+ bioinformatics tools
  - Data Modalities: Clinical
    - (Epic FHIR) + Genomic (FGbio, TCGA) +
    - Multi-omics (RNA/Protein/Phospho) +
    - Spatial (Visium) + Imaging (H&E, multiplex IF)
  - Cost Estimates:
    - Demonstration: DRY_RUN demo in 25-35 min (~$1 tokens only) or small files in 1-3 hours ($7-29)
    - Production: Realistic hospital data in 2-4 hours ($24-92 pre-aligned) or 4-8 hours ($29-102 raw FASTQ)
    - Includes: Compute + APIs + Claude tokens (~$1-2, stays low because servers return summaries)
- LEARN More:
  - [PatientOne Documentation](architecture/patient-one/README.md)
  - [Quick Start](tests/manual_testing/PatientOne-OvarianCancer/README.md)
  - [Sample Outputs](architecture/patient-one/patient-one-outputs/)
  - [Executive Summary](docs/EXECUTIVE_SUMMARY.md) 

## Who is this For?

This repository serves multiple audiences in the precision medicine ecosystem. Find your role below:

### 🔬 Bioinformaticians
*Analyze multi-omics cancer data, build data pipelines, develop predictive models*

**What you can do:** Spatial transcriptomics • Multi-omics integration • Tumor microenvironment mapping • Drug resistance mechanisms • Reproducible data pipelines

**📖 [Full Details & Resources →](docs/README.md#bioinformaticians)**

---

### 💻 MCP Developers
*Build custom MCP servers or extend existing bioinformatics tools*

**What you'll learn:** MCP server architecture • Testing best practices • External tool integration • Real vs mocked implementation strategies

**📖 [Full Details & Resources →](docs/README.md#mcp-developers)**

---

### 🛠️ Software Engineers
*Deploy, integrate, or scale this system*

**Quick Start (5 min):** Clone repo → Install dependencies → Configure Claude Desktop → Verify servers

**Deployment options:** Local development • GCP Cloud Run (9 servers deployed ✅) • HPC clusters • Hospital production

**📖 [Full Details & Resources →](docs/README.md#software-engineers)**

---

### 🏥 Clinical Care Teams
*Understand how AI-orchestrated bioinformatics supports clinical decision-making*

**⚠️ RESEARCH USE ONLY** - Not validated for clinical decision-making or patient care

**Educational value:** Precision medicine workflows • Multi-omics integration • Pathway analysis • FHIR & de-identification

**📖 [Full Details & Resources →](docs/README.md#clinical-care-teams-oncologists-genetic-counselors)**

---

### 👥 Patients & Families
*Understand precision medicine for ovarian cancer*

**⚠️ RESEARCH DEMONSTRATION** - Always consult qualified oncologists for medical decisions

**PatientOne Story:** Named in memory of a friend who passed from HGSOC in 2025, inspiring tools to help researchers understand and combat this disease.

**📖 [Full Details & Resources →](docs/README.md#patients--families)**

---

### 🎓 Students & Educators
*Learn or teach precision medicine and bioinformatics*

**Perfect for teaching:** 100% synthetic data • Low cost (~$0.32 per analysis) • Comprehensive coverage • Well-documented

**Topics covered:** Precision oncology • Multi-omics integration • Spatial transcriptomics • AI orchestration • Statistical methods • Cloud deployment

**📖 [Full Details & Resources →](docs/README.md#students--educators)**

---

## Acknowledgments

This project is dedicated to **PatientOne** - in memory of a dear friend who passed from High-Grade Serous Ovarian Carcinoma in 2025. Her courage inspired the creation of these AI-orchestrated bioinformatics tools.

**For detailed acknowledgments, citations, and scientific references:**
- [Complete Acknowledgments](ACKNOWLEDGMENTS.md)
- [Scientific References & Publications](docs/REFERENCES.md)

---
