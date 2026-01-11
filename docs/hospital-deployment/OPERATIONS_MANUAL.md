# Operations Manual - Precision Medicine MCP Servers
## Research Hospital Production Deployment

**Version:** 1.0
**Last Updated:** 2025-12
**Deployment:** HIPAA-Compliant Production
**Environment:** GCP Cloud Run

---

## Table of Contents

- [System Overview](#system-overview)
- [Architecture](#architecture)
- [Server Inventory](#server-inventory)
- [User Management](#user-management)
- [Monitoring & Alerting](#monitoring--alerting)
- [Backup & Disaster Recovery](#backup--disaster-recovery)
- [Incident Response](#incident-response)
- [Maintenance Procedures](#maintenance-procedures)
- [Security Operations](#security-operations)
- [Contact Information](#contact-information)

---

## System Overview

### Purpose

The Precision Medicine MCP system provides AI-powered analysis of spatial transcriptomics and multi-omics data for ovarian cancer research. The system integrates with Epic FHIR for clinical data and provides both web (Streamlit) and notebook (Jupyter) interfaces for clinicians and bioinformaticians.

### Key Components

1. **9 MCP Servers** - Bioinformatics analysis tools (Cloud Run)
2. **Streamlit UI** - Web-based chat interface
3. **Jupyter Notebooks** - Interactive data science environment
4. **OAuth2 Proxy** - Azure AD SSO authentication
5. **Epic FHIR Integration** - Clinical data access (de-identified)
6. **Cloud Logging** - 10-year audit trail (HIPAA compliant)

### Users

- **5 Testers:** 2 clinicians, 3 bioinformaticians
- **Data Scope:** 100 ovarian cancer patients during pilot
- **Access Method:** Azure AD SSO via hospital VPN

### Compliance

- HIPAA compliant
- 10-year audit log retention
- De-identification of all patient data
- Encrypted at rest and in transit

---

## Architecture

📋 **[See Complete Server Status →](../../servers/README.md#-server-status)** - All 10 MCP servers with tools, status, and documentation

### High-Level Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                      Hospital VPN                            │
│  ┌──────────────┐                   ┌──────────────┐       │
│  │   Clinician  │                   │Bioinformatician│      │
│  │   Browser    │                   │    Jupyter   │       │
│  └───────┬──────┘                   └───────┬──────┘       │
│          │                                  │              │
└──────────┼──────────────────────────────────┼──────────────┘
           │                                  │
           ▼                                  ▼
    ┌──────────────┐                  ┌──────────────┐
    │ OAuth2 Proxy │                  │ OAuth2 Proxy │
    │  (Streamlit) │                  │  (Jupyter)   │
    └──────┬───────┘                  └──────┬───────┘
           │                                  │
    ┌──────▼───────┐                  ┌──────▼───────┐
    │  Streamlit   │                  │  JupyterHub  │
    │   Chat UI    │                  │              │
    └──────┬───────┘                  └──────┬───────┘
           │                                  │
           └──────────────┬───────────────────┘
                          ▼
                  ┌───────────────┐
                  │  Claude API   │
                  │  (Anthropic)  │
                  └───────┬───────┘
                          │
        ┌─────────────────┼─────────────────┐
        │                 │                 │
    ┌───▼────┐       ┌───▼────┐       ┌───▼────┐
    │ fgbio  │       │multiomics│      │spatial │
    │        │       │          │      │ tools  │
    └───┬────┘       └───┬─────┘      └───┬────┘
        │                │                 │
        └────────────────┼─────────────────┘
                         ▼
                  ┌─────────────┐
                  │ Epic FHIR   │
                  │  (Research) │
                  └─────────────┘
```

### Network Architecture

- **VPC:** Hospital VPC network (existing)
- **VPC Connector:** `mcp-connector` for Cloud Run private networking
- **Ingress:** Internal and Cloud Load Balancing (no public access)
- **Egress:** All traffic routes through VPC

### Data Flow

1. User logs in via Azure AD SSO (OAuth2 Proxy)
2. Query submitted to Streamlit/Jupyter
3. Forwarded to Claude API with MCP server context
4. MCP servers called as needed:
   - Epic FHIR for clinical data (de-identified)
   - Spatial analysis tools for genomics
   - Multi-omics integration
5. Results returned to user
6. All actions logged to Cloud Logging (10-year retention)

---

## Server Inventory

### Production Servers (Real Analysis)

| Server | Status | URL | Service Account | Memory | CPU | Min Instances |
|--------|--------|-----|-----------------|--------|-----|---------------|
| **mcp-fgbio** | Production | `mcp-fgbio-{hash}.run.app` | `mcp-fgbio-sa` | 2Gi | 2 | 0 |
| **mcp-multiomics** | Production | `mcp-multiomics-{hash}.run.app` | `mcp-multiomics-sa` | 4Gi | 2 | 0 |
| **mcp-spatialtools** | Production | `mcp-spatialtools-{hash}.run.app` | `mcp-spatialtools-sa` | 4Gi | 2 | 0 |
| **mcp-epic** | Production | `mcp-epic-{hash}.run.app` | `mcp-epic-sa` | 2Gi | 2 | 0 |

### Mock/Demo Servers (Limited Functionality)

| Server | Status | URL | Purpose |
|--------|--------|-----|---------|
| **mcp-tcga** | Mock | `mcp-tcga-{hash}.run.app` | TCGA data demo |
| **mcp-openimagedata** | Partial | `mcp-openimagedata-{hash}.run.app` | Imaging data (30% real) |
| **mcp-seqera** | Mock | `mcp-seqera-{hash}.run.app` | Workflow demo |
| **mcp-huggingface** | Mock | `mcp-huggingface-{hash}.run.app` | Model integration demo |
| **mcp-deepcell** | Mock | `mcp-deepcell-{hash}.run.app` | Segmentation demo |

### User Interfaces

| Interface | URL | Authentication | Min Instances |
|-----------|-----|----------------|---------------|
| **Streamlit UI** | `streamlit-mcp-chat-{hash}.run.app` | OAuth2 Proxy | 1 |
| **Jupyter Notebook** | `jupyter-mcp-notebook-{hash}.run.app` | OAuth2 Proxy | 1 |

### Authentication Services

| Service | URL | Purpose |
|---------|-----|---------|
| **OAuth2 Proxy (Streamlit)** | `oauth2-proxy-streamlit-{hash}.run.app` | Azure AD SSO for Streamlit |
| **OAuth2 Proxy (Jupyter)** | `oauth2-proxy-jupyter-{hash}.run.app` | Azure AD SSO for Jupyter |

### Health Check Commands

```bash
# Check all servers are running
for server in mcp-fgbio mcp-multiomics mcp-spatialtools mcp-epic; do
  echo "Checking $server..."
  gcloud run services describe $server --region=us-central1 --format="value(status.url)"
done

# Check OAuth2 Proxy health
curl https://oauth2-proxy-streamlit-{hash}.run.app/ping
curl https://oauth2-proxy-jupyter-{hash}.run.app/ping
```

---

## User Management

### Current Users

- **Clinician 1:** Dr. Sarah Johnson (sarah.johnson@hospital.org)
- **Clinician 2:** Dr. Michael Chen (michael.chen@hospital.org)
- **Bioinformatician 1:** Dr. Emily Rodriguez (emily.rodriguez@hospital.org)
- **Bioinformatician 2:** Alex Kim (alex.kim@hospital.org)
- **Bioinformatician 3:** Jordan Taylor (jordan.taylor@hospital.org)

### Adding a New User

1. **Add to Azure AD Group:**
   ```
   Azure Portal -> Azure Active Directory -> Groups
   -> "precision-medicine-users" -> Add Member
   ```

2. **Verify Access:**
   - User logs in to Streamlit UI or Jupyter
   - Check audit logs for successful login event

3. **Training:**
   - Schedule training session (see USER_GUIDE.md)
   - Provide access to documentation
   - Assign test patient IDs for practice

### Removing a User

1. **Remove from Azure AD Group:**
   ```
   Azure Portal -> Azure Active Directory -> Groups
   -> "precision-medicine-users" -> Remove Member
   ```

2. **Verify Removal:**
   - User should be unable to log in
   - Check audit logs to confirm no new access

3. **Offboarding:**
   - Document any in-progress analyses
   - Transfer work to another user if needed

### Access Levels

**Current Configuration (Phase 1 - Simplified):**
- All users have identical access (no RBAC)
- Can use all 9 MCP servers
- Can analyze all 100 test patients

**Future Enhancements (Phase 2):**
- Role-based access control (RBAC)
- Clinician role: Read-only analysis results
- Bioinformatician role: Full analysis capabilities
- Admin role: User management, monitoring

---

## Monitoring & Alerting

### Cloud Monitoring Dashboard

**Access:** [Cloud Console Monitoring](https://console.cloud.google.com/monitoring/dashboards)

**Key Metrics:**
- Server request rate
- Server error rate (5xx responses)
- User query volume by server
- De-identification success rate
- Cost per analysis
- Token usage trends

### Alert Policies

| Alert | Condition | Notification |
|-------|-----------|--------------|
| **Server Down** | 10+ 5xx errors in 5 minutes | Hospital IT + Dev Team |
| **High Error Rate** | >5% error rate | Hospital IT + Dev Team |
| **Budget Alert** | 50%, 75%, 90%, 100% of $1,000/month | PI + Finance |
| **Epic FHIR Failures** | 5+ failures in 10 minutes | Dev Team + Hospital IT |
| **De-ID Failures** | Any de-identification failure | Dev Team + Privacy Officer |

### Log-Based Metrics

```bash
# View de-identification success rate
gcloud logging read 'jsonPayload.event="deidentification"' \
  --limit=100 --format=json

# View Epic FHIR failures
gcloud logging read 'jsonPayload.event="epic_fhir_call" AND jsonPayload.status="error"' \
  --limit=50

# View user access events
gcloud logging read 'jsonPayload.event="user_login" OR jsonPayload.event="mcp_query"' \
  --limit=100
```

### Daily Monitoring Checklist

- [ ] Check Cloud Monitoring dashboard for anomalies
- [ ] Review error logs (5xx responses)
- [ ] Check budget spend vs. forecast
- [ ] Verify all services are healthy
- [ ] Review Epic FHIR connection status
- [ ] Check user access patterns for anomalies

---

## Backup & Disaster Recovery

### What is Backed Up

1. **User Data:** None (system is stateless, all data in hospital GCS)
2. **Configuration:** Secret Manager, Cloud Run service configs
3. **Audit Logs:** 10-year retention in Cloud Logging
4. **Container Images:** Stored in Container Registry

### Backup Procedures

**Configuration Backup (Weekly):**
```bash
# Export all service configurations
for server in mcp-fgbio mcp-multiomics mcp-spatialtools mcp-epic; do
  gcloud run services describe $server --region=us-central1 \
    --format=yaml > backups/$server-config-$(date +%Y%m%d).yaml
done

# Export secrets list (not values)
gcloud secrets list --format=yaml > backups/secrets-list-$(date +%Y%m%d).yaml
```

**Audit Log Export (Monthly):**
```bash
# Export audit logs for compliance
gcloud logging read 'resource.type="cloud_run_revision"' \
  --format=json --limit=10000 > audit-logs-$(date +%Y%m).json
```

### Disaster Recovery

**RTO (Recovery Time Objective):** 4 hours
**RPO (Recovery Point Objective):** 24 hours

**Recovery Procedure:**

1. **Infrastructure Failure:**
   - Cloud Run auto-heals failed containers
   - If entire region fails, redeploy to secondary region:
     ```bash
     ./infrastructure/hospital-deployment/setup-project.sh
     REGION=us-east1 ./infrastructure/hospital-deployment/setup-vpc.sh
     # ... redeploy all services
     ```

2. **Epic FHIR Connection Failure:**
   - Switch to mcp-mockepic server temporarily
   - Contact hospital IT to restore Epic connection
   - See RUNBOOKS/epic-connection-failure.md

3. **OAuth2 Proxy Failure:**
   - Users cannot log in
   - Redeploy OAuth2 Proxy:
     ```bash
     ./infrastructure/hospital-deployment/deploy-oauth2-proxy.sh
     ```

4. **Complete System Failure:**
   - Estimated recovery: 2-4 hours
   - Follow deployment plan from scratch
   - Restore configuration from backups

### Testing DR Procedures

**Quarterly DR Test:**
- Simulate server failure
- Test OAuth2 Proxy failover
- Verify Epic FHIR fallback to mock
- Document lessons learned

---

## Incident Response

### Severity Levels

| Severity | Description | Response Time | Escalation |
|----------|-------------|---------------|------------|
| **P0 - Critical** | System down, data breach | 15 minutes | Immediate to CTO + Privacy Officer |
| **P1 - High** | Major functionality broken | 1 hour | Hospital IT Lead |
| **P2 - Medium** | Partial functionality issue | 4 hours | Dev Team Lead |
| **P3 - Low** | Minor issue, workaround available | 24 hours | Dev Team |

### Incident Response Procedure

1. **Detection:**
   - Alert notification
   - User report
   - Monitoring dashboard anomaly

2. **Triage:**
   - Assess severity (P0-P3)
   - Identify affected users
   - Determine root cause

3. **Communication:**
   - P0/P1: Immediately notify all users
   - P2/P3: Email update within 4 hours
   - Use template: `docs/hospital-deployment/templates/incident-notification.md`

4. **Resolution:**
   - Follow appropriate runbook
   - Document all actions taken
   - Test fix before declaring resolved

5. **Post-Mortem:**
   - Write incident report within 48 hours
   - Identify root cause
   - Document preventive measures
   - Update runbooks if needed

### Security Incident Response

**If PHI Exposure Suspected:**

1. **IMMEDIATE:** Shut down affected service
2. **NOTIFY:** Privacy Officer + Hospital IT Security + Dev Team
3. **INVESTIGATE:** Review audit logs for extent of exposure
4. **CONTAIN:** Revoke access, rotate credentials
5. **DOCUMENT:** Complete incident report for compliance
6. **REPORT:** Follow hospital's breach notification procedures

### Recent Incidents

| Date | Severity | Description | Resolution | Time to Resolve |
|------|----------|-------------|------------|-----------------|
| 2025-XX-XX | P2 | Epic FHIR timeout errors | Increased Cloud Run timeout to 300s | 2 hours |
| 2025-XX-XX | P3 | Slow query response times | Set min-instances=1 for core servers | 1 hour |

---

## Maintenance Procedures

### Scheduled Maintenance

**Weekly (Sundays 2-4 AM):**
- Review and archive old logs
- Check for security updates
- Review cost reports
- Backup configurations

**Monthly:**
- Update dependencies
- Security patch review
- Performance optimization review
- User feedback review

**Quarterly:**
- Disaster recovery test
- Security audit
- HIPAA compliance review
- User satisfaction survey

### Deploying Updates

**Server Updates:**
```bash
# Build new container version
cd servers/mcp-fgbio
docker build -t gcr.io/{PROJECT_ID}/mcp-fgbio:v1.1 .
docker push gcr.io/{PROJECT_ID}/mcp-fgbio:v1.1

# Deploy with zero downtime
gcloud run deploy mcp-fgbio \
  --image=gcr.io/{PROJECT_ID}/mcp-fgbio:v1.1 \
  --region=us-central1

# Monitor for errors
gcloud run services logs read mcp-fgbio --limit=50
```

**Rollback Procedure:**
```bash
# List previous revisions
gcloud run revisions list --service=mcp-fgbio --region=us-central1

# Rollback to previous revision
gcloud run services update-traffic mcp-fgbio \
  --to-revisions=mcp-fgbio-00001-abc=100 \
  --region=us-central1
```

---

## Security Operations

### Access Control

- **Azure AD Groups:** `precision-medicine-users` (all 5 testers)
- **GCP IAM:** Service accounts with least-privilege access
- **Network:** VPC isolation, no public access
- **VPN:** Required for all access

### Secret Management

**Secrets in Secret Manager:**
- `anthropic-api-key`
- `epic-fhir-endpoint`
- `epic-client-id`
- `epic-client-secret`
- `azure-ad-client-id`
- `azure-ad-client-secret`
- `azure-ad-tenant-id`

**Secret Rotation:**
- Epic FHIR credentials: Annually (hospital IT manages)
- Azure AD credentials: Annually (hospital IT manages)
- Anthropic API key: As needed (dev team manages)

**Access Audit:**
```bash
# List who can access secrets
for secret in anthropic-api-key epic-fhir-endpoint; do
  echo "=== $secret ==="
  gcloud secrets get-iam-policy $secret
done
```

### De-identification Validation

**Monthly Validation:**
```bash
# Check de-identification logs
gcloud logging read 'jsonPayload.event="deidentification"' --limit=100

# Verify no PHI in logs
gcloud logging read 'jsonPayload.prompt~"patient"' --limit=10
# Should show truncated prompts only
```

### Security Scanning

**Container Security:**
```bash
# Scan for vulnerabilities
gcloud container images scan gcr.io/{PROJECT_ID}/mcp-fgbio:latest
gcloud container images describe gcr.io/{PROJECT_ID}/mcp-fgbio:latest \
  --show-package-vulnerability
```

---

## Contact Information

### Internal Team

| Role | Name | Email | Phone | Escalation |
|------|------|-------|-------|------------|
| **PI (Ovarian Cancer)** | Dr. Jennifer Martinez | j.martinez@hospital.org | (555) 123-4567 | Budget, research questions |
| **Hospital IT Lead** | Robert Kim | r.kim@hospital.org | (555) 234-5678 | Infrastructure, networking |
| **Privacy Officer** | Lisa Thompson | l.thompson@hospital.org | (555) 345-6789 | HIPAA compliance, PHI issues |
| **IT Security Lead** | David Brown | d.brown@hospital.org | (555) 456-7890 | Security incidents |

### Development Team

| Role | Contact | Hours | SLA |
|------|---------|-------|-----|
| **Tier 2 Support** | mcp-support@devteam.com | Mon-Fri 9-5 PT | 4 hours |
| **On-Call Engineer** | oncall-mcp@devteam.com | 24/7 | 1 hour (P0/P1) |
| **Project Lead** | lead@devteam.com | Mon-Fri 9-5 PT | 24 hours |

### Vendor Support

| Vendor | Service | Contact | Hours |
|--------|---------|---------|-------|
| **Anthropic** | Claude API | support@anthropic.com | 24/7 |
| **Google Cloud** | GCP Platform | Enterprise support portal | 24/7 |
| **Epic** | FHIR API | Hospital IT manages | Business hours |

### Emergency Contacts

**P0 Incident (System Down / Data Breach):**
1. On-Call Engineer (immediate)
2. Hospital IT Security (immediate)
3. Privacy Officer (within 1 hour if PHI involved)
4. Development Team Lead (within 1 hour)

**P1 Incident (Major Functionality Broken):**
1. Hospital IT Lead
2. Development Team Support
3. PI (if research impacted)

---

**Document History:**
- v1.0 (2025-12-30): Initial operations manual for production deployment
- Next Review: 2026-03-30 (quarterly)

**Related Documents:**
- [User Guide](USER_GUIDE.md) - For end users
- [Admin Guide](ADMIN_GUIDE.md) - For administrators
- [HIPAA Compliance](HIPAA_COMPLIANCE.md) - Compliance documentation
- [Runbooks](RUNBOOKS/) - Incident response procedures
