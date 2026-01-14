# Documentation Index

**Welcome to the Precision Medicine MCP Documentation**

This is your central navigation hub for all documentation. Use this index to quickly find guides, references, and resources.

---

## 🚀 Getting Started

**New to the Precision Medicine MCP system? Start here:**

- **[Installation Guide](./getting-started/installation.md)** - Complete setup instructions (Quick Start: 5 minutes)
- **[README](./README.md)** - Project overview and quick links
- **[Executive Summary](./EXECUTIVE_SUMMARY.md)** - High-level overview for decision-makers

---

## 📚 User Guides

**Role-specific guides for different users:**

- **[For Funders & Grant Reviewers](./for-funders/README.md)** - ROI analysis, competitive landscape, grant talking points
- **[For Hospitals & IT Teams](./for-hospitals/README.md)** - Security, deployment, HIPAA compliance
- **[For Developers](./for-developers/README.md)** - Architecture, contributing, building new servers
- **[For Researchers](./for-researchers/README.md)** - Research workflows, analysis examples, production servers
- **[For Educators](./for-educators/README.md)** - Teaching materials, course integration, classroom activities
- **[For Patients & Families](./for-patients/README.md)** - Understanding precision medicine results (plain language)

**Specialized Guides:**
- **[Automated Patient Reports](./guides/AUTOMATED_PATIENT_REPORTS.md)** - Generate reports automatically
- **[Patient Summaries](./guides/GENERATE_PATIENT_SUMMARIES.md)** - Create patient-friendly summaries
- **[Add New Modality Server](./guides/ADD_NEW_MODALITY_SERVER.md)** - Extend the system with new data types

---

## 🏗️ Architecture

**Technical system design and implementation:**

- **[Architecture Overview](./architecture/README.md)** - System design principles
- **[Server Status](./architecture/servers.md)** - Implementation status and capabilities matrix
- **[Clinical-Spatial Bridge](./architecture/clinical-spatial-bridge.md)** - Integration between clinical and spatial data
- **[Error Handling](./architecture/error-handling.md)** - Error handling and retry logic
- **[References](./architecture/references.md)** - Technical references and citations

---

## 📋 Operations

**Day-to-day operations and management:**

- **[Cost and Budget Management](./operations/cost-and-budget.md)** - Cost estimation, tracking, and optimization
- **[Data Governance](./compliance/data-governance.md)** - Data handling policies and procedures

---

## 🔒 Compliance

**Regulatory compliance and data security:**

- **[Compliance Overview](./compliance/README.md)** - Compliance framework summary
- **[HIPAA Compliance](./compliance/hipaa.md)** - De-identification, audit logging, encryption
- **[Data Governance](./compliance/data-governance.md)** - GDPR, Common Rule, IRB requirements
- **[Risk Assessment](./compliance/risk-assessment.md)** - Risk mitigation strategies (if exists)
- **[Disclaimers](./compliance/disclaimers.md)** - Legal disclaimers and limitations

---

## 🚀 Deployment

**Deployment guides and infrastructure:**

- **[Deployment Roadmap](./deployment/roadmap.md)** - Production deployment planning
- **[Security Guide](./deployment/security.md)** - API keys, secrets management, GCP Secret Manager
- **[POC Deployment](./deployment/poc-deployment/)** - Proof-of-concept deployment guides
- **[Hospital Deployment](./deployment/hospital-deployment/)** - Enterprise deployment
  - [Operations Manual](./deployment/hospital-deployment/operations-manual.md)
  - [Admin Guide](./deployment/hospital-deployment/admin-guide.md)
  - [User Guide](./deployment/hospital-deployment/user-guide.md)
  - [Audit Log Guide](./deployment/hospital-deployment/audit-log-guide.md)
  - [Runbooks](./deployment/hospital-deployment/runbooks/) - Incident response procedures

---

## 🧪 Testing

**Test documentation, strategies, and test data:**

- **[Test Documentation Index](./test-docs/README.md)** - Overview of all test documentation
- **[Test Coverage & Guidelines](./test-docs/test-coverage.md)** - Test structure and best practices
- **[Manual Testing](./test-docs/manual-testing/)** - Quick test prompts and verification
- **[PatientOne Scenario](./test-docs/patient-one-scenario/)** - Complete end-to-end testing scenario
  - [Quick Reference](./test-docs/patient-one-scenario/quick-reference.md)
  - [CITL Quick Test](./test-docs/patient-one-scenario/citl-quick-test.md)
  - [Test Prompts](./test-docs/patient-one-scenario/test-prompts/) - Ready-to-use test prompts (6 tests)
- **[Integration Testing](./test-docs/integration-testing/)** - GCP and API testing

---

## 🏥 Clinical Workflows

**Clinical decision support and review processes:**

- **[Clinical Overview](./clinical/README.md)** - Clinical workflows overview
- **[CITL Workflow](./clinical/citl-workflow.md)** - Clinician-in-the-Loop workflow
- **[CITL Review Template](./clinical/citl-review-template.md)** - Review form template
- **[CITL Examples](./clinical/citl-examples.md)** - Example reviews

---

## ⚖️ Ethics

**Ethical considerations and bias auditing:**

- **[Ethics Overview](./ethics/README.md)** - Ethics framework
- **[Bias Framework](./ethics/bias-framework.md)** - Bias detection and mitigation
- **[Audit Checklist](./ethics/audit-checklist.md)** - Step-by-step bias audit
- **[Implementation](./ethics/implementation.md)** - Implementation guidelines

---

## 📦 Examples & Templates

**Example data, prompts, and templates:**

- **[Example Patients](./guides/examples/)** - Sample patient data and analyses
- **[Prompt Templates](./guides/prompts/)** - Reusable analysis prompts

---

## 🗄️ Archive

**Outdated documentation preserved for historical reference:**

- **[Archive Index](./archive/README.md)** - Why docs were archived
- **[2025 Q3/Q4 Archive](./archive/2025-q3-q4/)** - Pre-October 2025 docs

⚠️ **Warning:** Archived documentation is outdated and should not be used for current implementations.

---

## 📖 Additional Resources

### Quick Links

- **Installation:** [5-Minute Quick Start](./getting-started/installation.md#quick-start-5-minutes)
- **Test It:** [PatientOne Quick Test](./test-docs/patient-one-scenario/quick-reference.md)
- **Costs:** [Cost Overview](./operations/cost-and-budget.md#cost-overview)
- **HIPAA:** [HIPAA Quick Reference](./compliance/hipaa.md#executive-summary)

### By Task

**I want to:**
- **Install the system** → [Installation Guide](./getting-started/installation.md)
- **Run my first analysis** → [Quick Test Prompts](./test-docs/manual-testing/quick-test-prompts.md)
- **Understand costs** → [Cost and Budget Guide](./operations/cost-and-budget.md)
- **Deploy to production** → [Deployment Roadmap](./deployment/roadmap.md)
- **Ensure HIPAA compliance** → [HIPAA Compliance](./compliance/hipaa.md)
- **Add a new server** → [Add New Modality Server](./guides/ADD_NEW_MODALITY_SERVER.md)
- **Review patient results** → [CITL Workflow](./clinical/citl-workflow.md)
- **Test the system** → [PatientOne Scenario](./test-docs/patient-one-scenario/README.md)

### By Role

- **Funder/Grant Reviewer** → [For Funders](./for-funders/README.md) - ROI, competitive landscape, grant materials
- **Hospital Administrator** → [For Hospitals](./for-hospitals/README.md) - Security, deployment checklist, HIPAA
- **Developer** → [For Developers](./for-developers/README.md) - Architecture, contributing, building servers
- **Researcher/Bioinformatician** → [For Researchers](./for-researchers/README.md) - Workflows, analysis, production servers
- **Educator/Professor** → [For Educators](./for-educators/README.md) - Teaching materials, course integration
- **Patient/Family** → [For Patients](./for-patients/README.md) - Understanding results (plain language)
- **Clinician** → [CITL Workflow](./clinical/CITL_WORKFLOW_GUIDE.md) - Clinical workflow integration

---

## 🔍 Search Tips

**Finding documentation:**

1. **Browse this INDEX.md** - Organized by category
2. **Check the README** - Each directory has a README with local navigation
3. **Use your IDE's search** - Search across all .md files
4. **Check related docs** - Most docs have "Related Documentation" sections at the bottom

**Common searches:**
- Files ending in `-guide.md` - How-to guides
- Files starting with `for-` - Role-specific guides
- Files in `/test-docs/` - Testing documentation
- Files in `/compliance/` - Regulatory and security docs

---

## 📝 Documentation Standards

All documentation follows these standards:

- **Last Updated:** Each file has a "Last Updated" timestamp
- **Cross-Links:** Related docs are linked in "Related Documentation" sections
- **Table of Contents:** Long docs (>200 lines) have a ToC
- **Examples:** Code examples use syntax highlighting
- **File Naming:** kebab-case.md (e.g., cost-and-budget.md)

---

## 🆘 Getting Help

**If you can't find what you need:**

1. Check this INDEX.md for the right category
2. Browse the relevant directory's README.md
3. Search for keywords across all documentation
4. Create an issue on GitHub with the `documentation` label

**Documentation Feedback:**
- Found a broken link? Report it as an issue
- Documentation unclear? Suggest improvements
- Missing documentation? Request new guides

---

**Last Updated:** 2026-01-13

**Total Documentation Files:** ~40 (reduced from 47 through consolidation)

**Documentation Structure:**
```
docs/
├── INDEX.md (this file)
├── README.md
├── EXECUTIVE_SUMMARY.md
├── getting-started/
├── guides/
├── architecture/
├── operations/
├── compliance/
├── deployment/
├── clinical/
├── ethics/
├── test-docs/
└── archive/
```
