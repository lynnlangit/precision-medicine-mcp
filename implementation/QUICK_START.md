# Phase 1 Quick Start Guide

**Status:** Ready to Begin! 🚀
**Budget:** GCP free tier + minimal services (~$50-75/month)
**Timeline:** 3-4 weeks

---

## What We're Building

Moving the Precision Medicine MCP POC from local testing to GCP cloud infrastructure with real Epic EHR integration and enhanced spatial analysis capabilities.

**End Goal:** Process actual patient data (de-identified) for ovarian cancer spatial transcriptomics analysis in a HIPAA-compliant environment.

---

## Implementation Guides Created

1. ✅ **[PHASE1A_GCP_SETUP.md](PHASE1A_GCP_SETUP.md)** - GCP infrastructure (Week 1-2)
2. ✅ **[PHASE1B_EPIC_INTEGRATION.md](PHASE1B_EPIC_INTEGRATION.md)** - Epic FHIR server (Week 2-3)
3. ⏳ **PHASE1C_SPATIALTOOLS.md** - Complete spatial analysis (Week 3-4)

---

## Getting Started: Your First Steps

### Step 1: Set Up GCP Account (15 mins)

1. **Create GCP account** (if you don't have one):
   - Go to: https://cloud.google.com/free
   - Use $300 free credit for new accounts

2. **Install gcloud CLI**:
   ```bash
   # macOS
   curl https://sdk.cloud.google.com | bash
   exec -l $SHELL
   gcloud init

   # Authenticate
   gcloud auth login
   ```

3. **Get your billing account ID**:
   ```bash
   gcloud billing accounts list
   ```

   Save this for later: `BILLING_ACCOUNT_ID=YOUR_ID_HERE`

### Step 2: Review Phase 1A (30 mins)

Read through **[PHASE1A_GCP_SETUP.md](PHASE1A_GCP_SETUP.md)** to understand:
- What infrastructure will be created
- Cost estimates (~$50-75/month)
- Security features (encryption, audit logs, etc.)
- What you'll be able to do once complete

### Step 3: Execute Phase 1A (2-3 hours)

Follow the step-by-step instructions in PHASE1A_GCP_SETUP.md:

**Key Tasks:**
- [ ] Create GCP project
- [ ] Set up encrypted storage
- [ ] Configure VPC networking
- [ ] Set up Secret Manager
- [ ] Enable audit logging
- [ ] Create service accounts

**Estimated Time:** 2-3 hours (can pause and resume anytime)

**Budget:** Setup is free; ongoing costs ~$50-75/month

### Step 4: Register Epic Sandbox (30 mins)

While GCP is being set up:

1. Go to: https://fhir.epic.com/
2. Create account
3. Register new app:
   - Name: "Precision Medicine MCP"
   - FHIR Version: R4
   - Scopes: `patient/*.read`
4. Save credentials (you'll need these for Phase 1B)

### Step 5: Execute Phase 1B (1 day)

Follow **[PHASE1B_EPIC_INTEGRATION.md](PHASE1B_EPIC_INTEGRATION.md)**:

**Key Tasks:**
- [ ] Create mcp-epic server
- [ ] Implement FHIR API integration
- [ ] Add HIPAA de-identification
- [ ] Write tests (>70% coverage)
- [ ] Test in Claude Desktop

**Estimated Time:** 8-11 hours

### Step 6: Execute Phase 1C (coming soon)

Will complete mcp-spatialtools with real Scanpy integration.

---

## Cost Breakdown

### One-Time Costs
- **$0** - GCP setup (using free tier)
- **$0** - Epic sandbox access

### Monthly Costs (POC Tier)
| Service | Cost |
|---------|------|
| Cloud Storage | $1-5 |
| KMS (encryption keys) | $1 |
| Secret Manager | $1 |
| Cloud Logging | $5-10 |
| Networking (NAT) | $33 |
| Healthcare API | $5-20 |
| **Total** | **$50-75/month** |

**Budget Tip:** Use GCP $300 free credit to cover first 4-6 months!

---

## Timeline

```
Week 1: GCP Infrastructure Setup (Phase 1A)
│
├─ Day 1-2: Complete GCP setup
├─ Day 3: Register Epic sandbox
├─ Day 4-5: Test infrastructure
│
Week 2: Epic FHIR Integration (Phase 1B)
│
├─ Day 1-2: Build mcp-epic server
├─ Day 3: Write tests
├─ Day 4-5: Integration & testing
│
Week 3-4: Spatial Analysis Enhancement (Phase 1C)
│
├─ Week 3: Implement Scanpy integration
├─ Week 4: Testing & validation
```

---

## What You'll Have After Phase 1

### Infrastructure
- ✅ HIPAA-compliant GCP environment
- ✅ Encrypted storage (patient data, results)
- ✅ Private networking (VPC)
- ✅ Audit logging (10-year retention)
- ✅ Secret management
- ✅ Cost controls & budget alerts

### Servers (3/9 Production-Ready)
- ✅ **mcp-epic** - Real Epic EHR integration with de-identification
- ✅ **mcp-multiomics** - Multi-omics analysis (already done)
- ✅ **mcp-fgbio** - Genomic QC (already done)
- ✅ **mcp-spatialtools** - Complete spatial analysis with Scanpy

### Capabilities
- Retrieve real patient clinical data (de-identified)
- Process spatial transcriptomics data
- Perform multi-omics integration
- Generate results for ovarian cancer cases
- All in HIPAA-compliant environment

---

## Decision Points

Before starting, confirm:

**1. GCP Project Name**
- Default: `precision-medicine-poc`
- Custom: ________________

**2. GCP Region**
- Recommended: `us-central1` (Iowa - HIPAA eligible)
- Alternative: `us-east4` (Virginia - HIPAA eligible)
- Your choice: ________________

**3. Budget Limit**
- Set budget alert at: $______ (suggested: $100/month)

**4. Epic Credentials**
- Using sandbox (free): ✅ Recommended for Phase 1
- Using hospital Epic: ⏳ Wait until Phase 2

---

## Support & Troubleshooting

### Common Issues

**Issue: "gcloud command not found"**
```bash
# Reinstall gcloud CLI
curl https://sdk.cloud.google.com | bash
exec -l $SHELL
```

**Issue: "Permission denied" errors**
```bash
# Check you're authenticated
gcloud auth list

# Re-authenticate if needed
gcloud auth login
```

**Issue: "Billing account required"**
- Go to: https://console.cloud.google.com/billing
- Add payment method (won't be charged without your approval)

**Issue: "API not enabled"**
```bash
# Enable required APIs
gcloud services enable compute.googleapis.com storage-api.googleapis.com --project=YOUR_PROJECT_ID
```

### Getting Help

1. **Documentation Issues:** Check the detailed guides in `implementation/`
2. **GCP Issues:** GCP documentation at https://cloud.google.com/docs
3. **Epic FHIR Issues:** Epic documentation at https://fhir.epic.com/Documentation
4. **MCP Server Issues:** Create GitHub issue in this repo

---

## Next Steps

**Ready to start?**

1. ✅ Read this guide
2. ⏭️ Open `PHASE1A_GCP_SETUP.md`
3. ⏭️ Begin Task 1: GCP Project Setup

**Questions before starting?**
- Review the detailed guides
- Check cost estimates
- Confirm you have GCP access

---

## Success Metrics

You'll know Phase 1 is complete when:

- [ ] GCP infrastructure is running
- [ ] Can store encrypted patient data
- [ ] Can retrieve Epic sandbox patient data
- [ ] Data is automatically de-identified
- [ ] mcp-epic works in Claude Desktop
- [ ] Can process spatial transcriptomics data
- [ ] All tests passing (>70% coverage)
- [ ] Monthly costs < $100

---

**Let's build this! 🚀**

Start with: `implementation/PHASE1A_GCP_SETUP.md`
