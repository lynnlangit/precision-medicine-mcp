# Documentation Index

Complete documentation for the Precision Medicine MCP system.

---

## 🎯 Start Here

| Document | Purpose |
|----------|---------|
| **[Executive Summary](EXECUTIVE_SUMMARY.md)** | ROI analysis, budget, timeline |
| **[Why MCP for Healthcare?](WHY_MCP_FOR_HEALTHCARE.md)** | Understand MCP architecture advantages |
| **[Server Implementation Status](architecture/servers.md)** | All 10 servers status (167 tests ✅) |
| **[Server Details](../servers/)** | Individual server docs (10 mcp servers) |

---

## 👥 Find Your Role

Each guide includes workflows, examples, tools, and resources tailored to your needs:

| Role | What You'll Do | Your Guide |
|------|----------------|------------|
| 🔬 **Researchers & Bioinformaticians** | Analyze multi-omics data, spatial transcriptomics, build pipelines | [Guide →](for-researchers/README.md) |
| 💻 **Developers & Engineers** | Build MCP servers, deploy to cloud, integrate systems | [Guide →](for-developers/README.md) |
| 🏥 **Clinical Teams & Administrators** | Understand precision medicine workflows, manage deployments | [Guide →](for-hospitals/README.md) |
| 🎓 **Students & Educators** | Learn or teach precision medicine and bioinformatics | [Guide →](for-educators/README.md) |
| 👥 **Patients & Families** | Understand precision medicine for ovarian cancer | [Guide →](for-patients/README.md) |

---

## 📚 Documentation by Topic

### 🚀 Getting Started
- [Installation Guide](getting-started/installation.md) - Complete setup (Quick Start: 5 min)
- [Quick Test Prompts](test-docs/manual-testing/quick-test-prompts.md) - Copy-paste queries
- [Automated Patient Reports](guides/AUTOMATED_PATIENT_REPORTS.md) - Generate analysis reports

### ☁️ Deployment
- [Deployment Status](deployment/DEPLOYMENT_STATUS.md) - 9 servers on GCP Cloud Run ✅
- [GCP Testing Guide](deployment/GCP_TESTING_GUIDE.md) - Test via Claude API
- [Security](deployment/SECURITY.md) - POC security considerations

### 🏥 Hospital Production
- [Operations Manual](hospital-deployment/OPERATIONS_MANUAL.md) - System operations
- [HIPAA Compliance](hospital-deployment/compliance/hipaa.md) - De-identification, audit logs
- [Admin Guide](hospital-deployment/ADMIN_GUIDE.md) - User management, monitoring
- [User Guide](hospital-deployment/USER_GUIDE.md) - For clinicians and researchers
- [Audit Log Guide](hospital-deployment/AUDIT_LOG_GUIDE.md) - 10-year retention
- [Runbooks](hospital-deployment/RUNBOOKS/) - Incident response (server-down, SSO, Epic)

### 💰 Cost & Governance
- [Cost and Budget Management](operations/cost-and-budget.md) - Cost analysis, tracking, and optimization
- [Data Governance](operations/DATA_GOVERNANCE.md) - Privacy, retention policies
- [Risk Assessment](compliance/risk-assessment.md) - Comprehensive risk analysis and mitigation

### 🧪 Testing & QA
- [Testing Overview](test-docs/test-coverage.md) - 167 automated tests ✅
- [Quick Test Prompts](test-docs/manual-testing/quick-test-prompts.md) - Rapid testing
- [PatientOne Scenario](test-docs/patient-one-scenario/README.md) - End-to-end workflows
- [Test Code](../tests/README.md) - Python test suite

### 🔧 Technical Architecture
- [Error Handling & Retry Logic](architecture/error-handling.md) - Resilience patterns
- [Clinical-Spatial Bridge](architecture/clinical-spatial-bridge.md) - Data integration
- [Server Implementation Status](architecture/servers.md) - Production readiness

---

## 🏗️ Architecture by Modality

📋 **[See Individual Server Status →](../servers/README.md#-server-status)** - Detailed tools and documentation for all 10 servers

Detailed workflow architectures for each analysis type:

| Modality | Servers | Documentation |
|----------|---------|---------------|
| 🧬 **Clinical Data** | mcp-epic, mcp-mockepic | [Architecture →](architecture/clinical/README.md) |
| 🧪 **Genomic Cohorts** | mcp-tcga | [Architecture →](architecture/genomic/README.md) |
| 🖼️ **Imaging** | mcp-openimagedata, mcp-deepcell | [Architecture →](architecture/imaging/README.md) |
| 🔬 **Multiomics** | mcp-multiomics | [Architecture →](architecture/multiomics/README.md) |
| 📍 **Spatial Transcriptomics** | mcp-spatialtools, mcp-fgbio | [Architecture →](architecture/spatial-transcriptomics/README.md) |
| 🤖 **AI/ML Inference** | mcp-huggingface | [Architecture →](architecture/ai-ml/README.md) |
| ⚙️ **Workflow Orchestration** | mcp-seqera | [Architecture →](architecture/workflow/README.md) |

**End-to-End Example:** [PatientOne Precision Medicine Workflow](test-docs/patient-one-scenario/architecture/overview.md)

---
## Repository Structure

```
precision-medicine-mcp/
├── ACKNOWLEDGMENTS.md      # Credits & scientific references
├── FUNDING.md              # Investment opportunities & ROI
├── LICENSE                 # Apache 2.0 License
├── README.md               # This file
├── data/                   # Synthetic patient data (100% safe for demos)
├── docs/                   # Documentation organized by audience
│   ├── for-funders/        # ROI analysis, competitive landscape, grant talking points
│   ├── for-hospitals/      # Deployment checklist, security overview
│   ├── for-developers/     # Architecture, contributing guide, quick reference
│   ├── for-researchers/    # Analysis workflows, bioinformatics methods
│   ├── for-educators/      # Classroom guides, learning objectives
│   ├── for-patients/       # Patient-friendly resources
│   ├── demos/              # 90-second pitch, full PatientOne demo
│   ├── prompt-library/     # 20+ ready-to-use clinical prompts
│   ├── getting-started/    # Installation, quick start, desktop-configs
│   ├── architecture/       # System design & modality workflows
│   ├── for-hospitals/      # HIPAA compliance, operations manual
│   ├── deployment/         # GCP deployment status & guides
│   └── test-docs/          # Testing guides & PatientOne scenarios
├── infrastructure/         # Deployment, audit, environment setup
│   ├── deployment/         # GCP deployment scripts
│   ├── audit/              # Bias detection and audit tools
│   └── hospital-deployment/  # Hospital-specific infrastructure
├── servers/                # 10 MCP servers (Python)
│   ├── mcp-deepcell/       # Cell segmentation
│   ├── mcp-epic/           # Epic FHIR integration
│   ├── mcp-fgbio/          # Reference genomes, FASTQ QC
│   ├── mcp-huggingface/    # AI/ML inference
│   ├── mcp-mockepic/       # Mock Epic for testing
│   ├── mcp-multiomics/     # Multi-omics integration
│   ├── mcp-openimagedata/  # Imaging data (H&E, MxIF)
│   ├── mcp-seqera/         # Workflow orchestration
│   ├── mcp-spatialtools/   # Spatial transcriptomics
│   └── mcp-tcga/           # TCGA cohort data
├── shared/                 # Shared Python packages
│   ├── common/             # Common utilities
│   ├── models/             # Data models
│   ├── schemas/            # JSON schemas (CitL review, etc.)
│   └── utils/              # Helper functions
├── tests/                  # 167 automated tests
├── tools/                  # Automation & reporting tools
│   └── reports/            # Patient report generation, CitL submission
└── ui/                     # Streamlit chat, Jupyter notebook
```

---

