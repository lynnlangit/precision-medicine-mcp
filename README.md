# Precision Medicine w/ AI-orchestrated MCP Servers

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![MCP](https://img.shields.io/badge/MCP-2025--06--18-green.svg)](https://modelcontextprotocol.io/)
[![Claude Desktop](https://img.shields.io/badge/Claude-Desktop-orange.svg)](https://claude.ai/download)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

<img src="https://github.com/lynnlangit/precision-medicine-mcp/blob/main/data/images/repo-image.png">

## What and Why
- The Problem: multi-modal precision medicine is siloed and code-heavy, **too slow** for patient care treatment options
- This Solution: a set of custom MCP servers for analysis and data retrieval (PatientOne/Ovarian Cancer use case example):
  - Makes complex analysis significantly **faster and easier** by presenting a customized, single interface for tools and data accessible using **natural language** questions
  - System coordinates disparate servers and stitches results together after being given a prompt in English
  - Solution is **extensible** for other comorbidity types (i.e. other cancers, other diseases)
  - "What is an MCP Server?" [(article)](https://medium.com/@elisowski/mcp-explained-the-new-standard-connecting-ai-to-everything-79c5a1c98288)
- What it is NOT: Not clinically validated yet
- See it / Try it: <5 minute local demo (small subset of full example) - [recording](https://www.youtube.com/watch?v=LUldOHHX5Yo) | [code](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs/test-docs/patient-one-scenario)

## Featured Use Case: PatientOne

<kbd><img src="https://github.com/lynnlangit/precision-medicine-mcp/blob/main/tests/manual_testing/PatientOne-OvarianCancer/architecture/patient-one-holistic.png" width=800></kbd>  

- LEARN More:
  - [PatientOne Documentation](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/test-docs/patient-one-scenario/architecture/overview.md)
  - [Quick Start](docs/test-docs/patient-one-scenario/README.md)
  - [Sample Outputs](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs/test-docs/patient-one-scenario/architecture)
  - [Executive Summary](docs/EXECUTIVE_SUMMARY.md)

 ---

## 👥 Find Your Role

Each guide includes workflows, examples, tools, and resources tailored to your needs:

| Role | What You'll Do | Your Guide |
|------|----------------|------------|
| 🔬 **Researchers & Bioinformaticians** | Analyze multi-omics data, spatial transcriptomics, build reproducible pipelines | [Guide →](docs/guides/for-bioinformaticians.md) |
| 💻 **Developers & Engineers** | Build MCP servers, deploy to cloud, integrate bioinformatics systems | [Guide →](docs/guides/for-developers.md) |
| 🏥 **Clinical Teams & Administrators** | Understand precision medicine workflows, manage hospital deployments | [Guide →](docs/guides/for-clinicians.md) |
| 🎓 **Students & Educators** | Learn or teach precision medicine (100% synthetic data, ~$0.32/analysis) | [Guide →](docs/guides/for-researchers.md) |
| 👥 **Patients & Families** | Understand precision medicine for ovarian cancer ⚠️ Research only | [Guide →](docs/guides/for-patients.md) |

---

## Acknowledgments

This project is dedicated to **PatientOne** - in memory of a dear friend who passed from High-Grade Serous Ovarian Carcinoma in 2025. Her courage inspired the creation of these AI-orchestrated bioinformatics tools.

**For detailed acknowledgments, citations, and scientific references:**
- [Complete Acknowledgments](ACKNOWLEDGMENTS.md)
- [Scientific References & Publications](docs/architecture/references.md)

---
