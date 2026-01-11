# Production Roadmap: Hospital Cloud Deployment

## Executive Summary

This document outlines the prioritized path to move the Precision Medicine MCP POC from development to production deployment in a hospital cloud infrastructure for testing with actual patient data.

**POC Status (2025-12-30):** All 9 servers deployed to GCP Cloud Run, tested via Claude API
**Next Goal:** HIPAA-compliant hospital deployment with real patient data
**Timeline Estimate:** 24 weeks (6 months) for Phase 1-3

**Note:** Healthcare deployments require longer timelines than typical AdTech/FinTech projects due to HIPAA compliance, clinical validation, and practitioner review requirements.

---

## POC Deployment Status (2025-12-30)

### GCP Cloud Run Deployment: ✅ COMPLETE

All 9 MCP servers successfully deployed to Google Cloud Platform and tested.

📋 **[See Server Status →](../servers/README.md#-server-status)** - Complete status, tools, and implementation details for all 10 servers

**Test Results:** 9/9 servers passed functional testing via Claude API

**POC Environment:**
- Platform: Google Cloud Run (serverless containers)
- Transport: SSE (Server-Sent Events) over HTTPS
- Authentication: Public endpoints (POC only - not for PHI)
- Cost: Minimal ($0.10-0.50/day with pay-per-use)

**POC Limitations:**
- Not HIPAA compliant (no PHI allowed)
- Public endpoints (authentication not enabled)
- No audit logging for clinical use
- Test/mock data only

**See:**
- [Deployment Status](deployment/DEPLOYMENT_STATUS.md) - Full deployment details
- [GCP Testing Guide](deployment/GCP_TESTING_GUIDE.md) - Testing with Claude API

---

## Current State Assessment (For Hospital Production Deployment)

**Note:** This assessment is for moving from POC to hospital production. All servers are deployed and functional in POC environment.

### Hospital-Ready Servers (3/9)
- ✅ **mcp-multiomics** - 95% real, 91 tests, 68% coverage, validated preprocessing workflow
- ✅ **mcp-fgbio** - 95% real, 29 tests, 77% coverage, genomic reference data access
- ✅ **mcp-spatialtools** - 95% real, 60+ tests, comprehensive spatial analysis pipeline

### Mock Servers (Require Real Implementation for Hospital Use) (5/9)
- 🔶 **mcp-tcga** - Mock (needs real TCGA GDC API integration for comparative analysis)
- 🔶 **mcp-openimagedata** - 30% real (basic features working, needs advanced imaging integration)
- 🔶 **mcp-deepcell** - Mock (needs DeepCell API or self-hosted models)
- 🔶 **mcp-huggingface** - Mock (needs real HuggingFace API integration)
- 🔶 **mcp-seqera** - Mock (needs Seqera Platform license or Nextflow integration)

### Hospital EHR Integration Required (1/9)
- 🔵 **mcp-mockepic** - Must be replaced with real EHR connector (Epic FHIR, HL7, etc.)

**Key Changes Since Initial Assessment:**
- mcp-spatialtools: 40% → 95% real (batch correction, DE, pathway enrichment implemented)
- All servers deployed to GCP Cloud Run and tested
- Ready for Phase 1: Hospital infrastructure and EHR integration

---

## Critical Prerequisites for Hospital Deployment

Before ANY server implementation work:

### 1. Infrastructure & Security (MANDATORY)
- [ ] **HIPAA-compliant cloud environment setup**
  - Azure/AWS/GCP healthcare tier
  - Encrypted storage (AES-256)
  - VPC/private networking
  - Audit logging enabled

- [ ] **Data governance implementation**
  - De-identification pipeline
  - Access control (RBAC)
  - PHI handling procedures
  - Data retention policies

- [ ] **Security hardening**
  - Remove all DRY_RUN=false paths for testing
  - Input validation on all tools
  - Rate limiting
  - Secret management (Azure Key Vault, etc.)

### 2. Regulatory Compliance (MANDATORY)
- [ ] **IRB approval** for patient data testing
- [ ] **HIPAA compliance audit**
- [ ] **Data use agreements** with hospital
- [ ] **Privacy impact assessment**

### 3. Clinical Validation (MANDATORY)
- [ ] **Clinical partner engagement**
- [ ] **Use case validation** (confirm patient-one scenario is priority)
- [ ] **Clinical workflow integration** planning

**⚠️ BLOCKER:** Do NOT proceed to server implementation until prerequisites 1-3 are complete.

---

## Phase 1: Foundation & Critical Path (Weeks 1-4)

**Goal:** Complete infrastructure and implement servers needed for basic patient-one workflow

### Priority 1A: Infrastructure (Week 1-2)
**Owner:** DevOps/IT
**Rationale:** Foundational requirement for all subsequent work

**Tasks:**
1. Set up HIPAA-compliant cloud environment
2. Configure encrypted storage and networking
3. Implement audit logging
4. Set up CI/CD with security scanning
5. Configure secret management

**Deliverables:**
- Working cloud environment
- Security hardening checklist completed
- Deployment automation

### Priority 1B: Real EHR Integration (Week 2-3)
**Current:** mcp-mockepic (mock EHR)
**Target:** Real EHR connector
**Rationale:** Clinical data is the foundation of patient-one workflow

**Tasks:**
1. Replace mcp-mockepic with real EHR connector
   - Epic FHIR API integration OR
   - HL7 interface OR
   - Custom hospital EHR connector
2. Implement FHIR resource mapping (Patient, Observation, MedicationStatement)
3. Add de-identification at ingestion
4. Write integration tests with synthetic FHIR data
5. Test with 1-2 real de-identified patient records

**Deliverables:**
- mcp-epic (or hospital-specific) server
- 50+ tests, 70%+ coverage
- Documentation for EHR team

**Estimated Effort:** 40-60 hours

### Priority 1C: Spatial Transcriptomics - ✅ MOSTLY COMPLETE (Week 3-4)
**Completed:** mcp-spatialtools (95% real, deployed to GCP)
**Remaining:** Hospital-specific integration and STAR alignment setup
**Rationale:** Core analysis for patient-one ovarian cancer case

**Completed Tasks:**
1. ✅ Implemented real file I/O (CSV/TSV parsing)
2. ✅ QC filtering (statistical thresholds, not Scanpy but functional equivalent)
3. ✅ Spatial autocorrelation (Moran's I implemented and tested)
4. ✅ Differential expression (Wilcoxon/Mann-Whitney, not DESeq2 but statistically sound)
5. ✅ Batch correction (ComBat implemented and validated, Harmony/Scanorama deferred)
6. ✅ Pathway enrichment (Fisher's exact test, 44 curated pathways)
7. ✅ Cell type deconvolution (signature scoring)
8. ✅ 60+ tests written
9. ✅ Clinical interpretation layer (bridge tool for patient→spatial mapping)
10. ✅ Deployed to GCP Cloud Run and tested

**Remaining Tasks for Hospital Deployment:**
1. STAR alignment setup (genome index download, server configuration)
2. Integration with hospital file storage (replace local file paths)
3. Performance optimization for large datasets (>10K spots)
4. Clinical validation on real de-identified spatial data

**Current Status:**
- 95% real functionality ✅
- 60+ tests, comprehensive coverage ✅
- Validated on PatientOne ovarian cancer dataset ✅
- Production-ready for POC, needs hospital integration for clinical use

**Estimated Effort for Remaining Tasks:** 20-30 hours

---

## Phase 2: Imaging & Genomics (Weeks 5-8)

**Goal:** Complete imaging and genomics servers for comprehensive patient analysis

### Priority 2A: Histology Imaging (Week 5-6)
**Current:** mcp-openimagedata (30% real)
**Target:** 90%+ real implementation
**Rationale:** Critical for spatial analysis and clinical correlation

**Tasks:**
1. Implement real image I/O (whole slide images)
2. Integrate OpenSlide for image reading
3. Add spatial registration (SIFT/ORB)
4. Implement feature extraction (HOG, texture)
5. Add tissue segmentation
6. Write 30+ tests

**Deliverables:**
- 90%+ real functionality
- 40+ tests, 70%+ coverage
- Handles TIFF, SVS, NDPI formats

**Estimated Effort:** 50-70 hours

### Priority 2B: Cell Segmentation (Week 6-7)
**Current:** mcp-deepcell (0% real - fully mocked)
**Target:** 80%+ real implementation
**Rationale:** Needed for spatial analysis and cell type identification

**Decision Point:** Build vs. Buy
- **Option A:** Integrate existing DeepCell API (faster, less control)
- **Option B:** Self-hosted DeepCell models (more control, more work)
- **Recommendation:** Option A for Phase 2, migrate to B if needed

**Tasks (Option A):**
1. Integrate DeepCell API client
2. Add image preprocessing
3. Implement result parsing
4. Add quality validation
5. Write 25+ tests
6. Add fallback for API failures

**Deliverables:**
- 80%+ real functionality
- 30+ tests, 65%+ coverage
- API key management

**Estimated Effort:** 30-40 hours (Option A), 80-100 hours (Option B)

### Priority 2C: Cancer Genomics Data (Week 7-8)
**Current:** mcp-tcga (0% real - fully mocked)
**Target:** 70%+ real implementation
**Rationale:** Needed for comparative analysis and survival data

**Tasks:**
1. Integrate real TCGA API (GDC API)
2. Implement cohort queries
3. Add expression data retrieval
4. Implement mutation data access
5. Add survival analysis
6. Write 30+ tests

**Deliverables:**
- 70%+ real functionality
- 35+ tests, 65%+ coverage
- Rate limiting and caching

**Estimated Effort:** 40-50 hours

---

## Phase 3: ML & Workflow (Weeks 9-12)

**Goal:** Complete ML model integration and workflow orchestration

### Priority 3A: Genomic ML Models (Week 9-10)
**Current:** mcp-huggingface (0% real - fully mocked)
**Target:** 60%+ real implementation
**Rationale:** Nice-to-have for advanced analysis, not critical path

**Decision Point:** Scope
- **Minimum:** Model discovery and metadata only
- **Extended:** Model inference integration
- **Recommendation:** Minimum for Phase 3, defer inference to Phase 4

**Tasks (Minimum):**
1. Integrate Hugging Face API for model search
2. Add model metadata retrieval
3. Implement model filtering (genomics/biology only)
4. Add usage documentation
5. Write 20+ tests

**Deliverables:**
- 60%+ real functionality (discovery only)
- 25+ tests, 60%+ coverage

**Estimated Effort:** 25-35 hours (Minimum), 60-80 hours (Extended)

### Priority 3B: Workflow Orchestration (Week 11-12)
**Current:** mcp-seqera (0% real - fully mocked)
**Target:** 50%+ real implementation
**Rationale:** Lower priority - can be replaced with direct pipeline execution

**Decision Point:** Build vs. Replace
- **Option A:** Integrate Seqera Platform API (requires license)
- **Option B:** Replace with Nextflow direct execution
- **Option C:** Defer to Phase 4
- **Recommendation:** Option C - defer unless hospital has Seqera license

**Tasks (if pursued):**
1. TBD based on decision

**Estimated Effort:** TBD

---

## Phase 4: Clinical Integration & Validation (Weeks 13-16)

### Priority 4A: End-to-End Testing
**Tasks:**
1. Run patient-one workflow with real de-identified data
2. Clinical validation of results
3. Performance optimization
4. Error handling hardening
5. Clinical report generation

### Priority 4B: Clinical Workflow Integration
**Tasks:**
1. Integrate with hospital LIMS
2. Add clinical decision support hooks
3. Implement results delivery to clinicians
4. Create clinical user interface (if needed)

### Priority 4C: Production Monitoring
**Tasks:**
1. Application monitoring (Prometheus/Grafana)
2. Error alerting
3. Performance dashboards
4. Cost tracking dashboards
5. Audit log review procedures

---

## Prioritization Rationale

### Why This Order?

1. **Infrastructure First** - Nothing else works without secure foundation
2. **Clinical Data (EHR) First** - Patient data is the starting point for all analysis
3. **Spatial Analysis Second** - Core of the patient-one use case (ovarian cancer)
4. **Imaging Third** - Complements spatial analysis
5. **Genomics Fourth** - Provides comparative context
6. **ML/Workflow Last** - Nice-to-have, not critical path

### What We're NOT Prioritizing (and Why)

- **mcp-seqera** - Can use direct Nextflow execution instead
- **Advanced ML inference** - Discovery is sufficient for Phase 1-3
- **Full DeepCell implementation** - API integration sufficient initially

---

## Risk Mitigation Strategy

### Technical Risks

| Risk | Mitigation | Owner |
|------|-----------|-------|
| Data privacy breach | Multi-layer security, audit logs, de-identification | Security Team |
| Integration failures | Extensive testing, fallback mechanisms | Engineering |
| Performance issues | Load testing, optimization, caching | Engineering |
| Clinical misinterpretation | Clinical validation, disclaimers, training | Clinical Team |

### Regulatory Risks

| Risk | Mitigation | Owner |
|------|-----------|-------|
| HIPAA violation | Compliance audit, legal review | Legal/Compliance |
| IRB issues | Early engagement, protocol review | Research Team |
| Data use violations | Clear agreements, access controls | Legal |

### Clinical Risks

| Risk | Mitigation | Owner |
|------|-----------|-------|
| Incorrect results | Extensive validation, clinical review | Clinical Team |
| Workflow disruption | Pilot testing, training | Operations |
| Low adoption | User engagement, feedback loops | Clinical Champions |

---

## Success Metrics

### Technical Metrics
- [ ] All production servers: 70%+ test coverage
- [ ] All production servers: 80%+ real functionality
- [ ] API response time: <5 seconds (p95)
- [ ] Zero PHI leakage in logs/errors
- [ ] 99.5%+ uptime

### Clinical Metrics
- [ ] Successful analysis of 5 patient cases
- [ ] Clinical validation: 90%+ concordance with manual analysis
- [ ] Clinician satisfaction: 4+/5
- [ ] Time to insight: <4 hours (vs. days manually)

### Regulatory Metrics
- [ ] HIPAA compliance audit: Pass
- [ ] IRB approval: Obtained
- [ ] Security audit: Pass
- [ ] Privacy impact assessment: Complete

---

## Resource Requirements

### Team Composition (Estimated)
- **2x Backend Engineers** (server implementation)
- **1x DevOps Engineer** (infrastructure, deployment)
- **1x Security Engineer** (compliance, hardening)
- **1x Clinical Informaticist** (clinical validation, workflow)
- **0.5x Product Manager** (prioritization, stakeholder management)

### Budget Estimate
- **Personnel:** 4.5 FTE × 16 weeks = 72 person-weeks
- **Cloud Infrastructure:** $2,000-5,000/month (testing tier)
- **Third-party APIs:** $500-2,000/month (DeepCell, TCGA, HuggingFace)
- **Total:** $150K-250K (Phase 1-4)

---

## Decision Points for User

Before proceeding, please confirm:

1. **Infrastructure:** Which cloud provider? (Azure/AWS/GCP)
2. **EHR System:** Epic, Cerner, other? API access available?
3. **Timeline:** Is 24 weeks (6 months) acceptable?
4. **Resources:** Team availability for work?
5. **Clinical Partner:** Hospital/clinic identified and engaged?
6. **Regulatory:** IRB submission timeline?
7. **Budget:** $150K-250K budget available?

---

## Immediate Next Steps

**If greenlit, Week 1 tasks:**

1. **Kickoff meeting** with hospital IT, security, clinical teams
2. **Cloud environment setup** (Azure/AWS/GCP healthcare tier)
3. **EHR integration scoping** (meet with Epic/Cerner team)
4. **IRB protocol draft** (research team)
5. **Security requirements gathering** (security team)

**Week 1 deliverables:**
- Cloud environment provisioned
- EHR integration plan documented
- IRB protocol submitted
- Security checklist completed
- Sprint 1 planned (spatialtools completion)

---

## Conclusion

**Recommended Path:** Execute Phase 1-3 sequentially (12 weeks), then evaluate Phase 4 based on results.

**Critical Success Factor:** Do NOT skip infrastructure/security prerequisites. Patient data security is non-negotiable.

**Quick Win:** Focus on patient-one ovarian cancer use case end-to-end, rather than trying to make all servers production-ready simultaneously.

---

**Questions or concerns? Let's discuss before starting implementation.**
