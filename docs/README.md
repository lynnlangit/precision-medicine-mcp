# Documentation Index

Complete documentation for the Precision Medicine MCP system.

---

## 🎯 Start Here

| Document | Purpose |
|----------|---------|
| **[Executive Summary](EXECUTIVE_SUMMARY.md)** | ROI analysis, budget, timeline |
| **[Server Implementation Status](architecture/servers.md)** | All 10 servers status (167 tests ✅) |
| **[Server Details](../servers/)** | Individual server docs (10 mcp servers) |

---

## 👥 Find Your Role

Each guide includes workflows, examples, tools, and resources tailored to your needs:

| Role | What You'll Do | Your Guide |
|------|----------------|------------|
| 🔬 **Researchers & Bioinformaticians** | Analyze multi-omics data, spatial transcriptomics, build pipelines | [Guide →](guides/personas/bioinformaticians.md) |
| 💻 **Developers & Engineers** | Build MCP servers, deploy to cloud, integrate systems | [Guide →](guides/personas/mcp-developers.md) • [Guide →](guides/personas/software-engineers.md) |
| 🏥 **Clinical Teams & Administrators** | Understand precision medicine workflows, manage deployments | [Guide →](guides/personas/clinical-care-teams.md) |
| 🎓 **Students & Educators** | Learn or teach precision medicine and bioinformatics | [Guide →](guides/personas/students-and-educators.md) |
| 👥 **Patients & Families** | Understand precision medicine for ovarian cancer | [Guide →](guides/personas/patients-and-families.md) |

---

## 📚 Documentation by Topic

### 🚀 Getting Started
- [Installation Guide](getting-started/installation.md) - Complete setup (Quick Start: 5 min)
- [Quick Test Prompts](../tests/manual_testing/QUICK_TEST_PROMPTS.md) - Copy-paste queries
- [Automated Patient Reports](guides/AUTOMATED_PATIENT_REPORTS.md) - Generate analysis reports

### ☁️ Deployment
- [Deployment Status](deployment/DEPLOYMENT_STATUS.md) - 9 servers on GCP Cloud Run ✅
- [GCP Testing Guide](deployment/GCP_TESTING_GUIDE.md) - Test via Claude API
- [Security](deployment/SECURITY.md) - POC security considerations

### 🏥 Hospital Production
- [Operations Manual](hospital-deployment/OPERATIONS_MANUAL.md) - System operations
- [HIPAA Compliance](hospital-deployment/HIPAA_COMPLIANCE.md) - De-identification, audit logs
- [Admin Guide](hospital-deployment/ADMIN_GUIDE.md) - User management, monitoring
- [User Guide](hospital-deployment/USER_GUIDE.md) - For clinicians and researchers
- [Audit Log Guide](hospital-deployment/AUDIT_LOG_GUIDE.md) - 10-year retention
- [Runbooks](hospital-deployment/RUNBOOKS/) - Incident response (server-down, SSO, Epic)

### 💰 Cost & Governance
- [Cost and Budget Management](operations/cost-and-budget.md) - Cost analysis, tracking, and optimization
- [Data Governance](operations/DATA_GOVERNANCE.md) - Privacy, retention policies
- [Risk Assessment](compliance/risk-assessment.md) - Comprehensive risk analysis and mitigation

### 🧪 Testing & QA
- [Testing Overview](../tests/README.md) - 167 automated tests ✅
- [Quick Test Prompts](../tests/manual_testing/QUICK_TEST_PROMPTS.md) - Rapid testing
- [PatientOne Tests](../tests/manual_testing/PatientOne-OvarianCancer/) - End-to-end workflows

### 🔧 Technical
- [Error Handling & Retry Logic](technical/ERROR_HANDLING_RETRY_LOGIC.md) - Resilience patterns
- [Clinical-Spatial Integration](guides/CLINICAL_SPATIAL_INTEGRATION.md) - Data integration

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

