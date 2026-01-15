# HIPAA Compliance Documentation
## Precision Medicine MCP System

**Version:** 1.0
**Last Updated:** 2025-12
**Compliance Officer:** [Hospital Privacy Officer Name]
**System:** Precision Medicine MCP Servers

---

## Table of Contents

- [Executive Summary](#executive-summary)
- [HIPAA Requirements](#hipaa-requirements)
- [De-identification Validation](#de-identification-validation)
- [Access Control](#access-control)
- [Audit Controls](#audit-controls)
- [Encryption](#encryption)
- [PHI Handling Procedures](#phi-handling-procedures)
- [Compliance Validation](#compliance-validation)
- [Incident Response](#incident-response)

---

## Executive Summary

### Compliance Status

✅ **COMPLIANT** - The Precision Medicine MCP system meets all applicable HIPAA requirements for a research system handling de-identified patient data.

### Key Compliance Features

1. **De-identification:** All patient data de-identified using HIPAA Safe Harbor method
2. **Access Control:** Azure AD SSO with MFA, VPN required
3. **Audit Logging:** 10-year retention in encrypted Cloud Logging
4. **Encryption:** All data encrypted at rest and in transit
5. **BAA:** Business Associate Agreement with Google Cloud in place

### Scope

**What is covered:**
- Epic FHIR clinical data (de-identified)
- Spatial transcriptomics data
- Multi-omics data
- User queries and system logs

**What is NOT covered:**
- Patient identifiable information (removed via de-identification)
- Non-research patient data (not in system)

---

## HIPAA Requirements

### Administrative Safeguards

| Requirement | Implementation | Status |
|-------------|----------------|--------|
| **Security Management** | Incident response procedures, risk assessment | ✅ Implemented |
| **Workforce Security** | Azure AD user management, training required | ✅ Implemented |
| **Information Access** | Role-based access (planned), current: all users same access | ⚠️ Simplified (Phase 1) |
| **Security Awareness** | User training on PHI handling, cost management | ✅ Implemented |
| **Contingency Plan** | Disaster recovery procedures, backups | ✅ Implemented |

### Physical Safeguards

| Requirement | Implementation | Status |
|-------------|----------------|--------|
| **Facility Access** | GCP data centers (hospital has BAA) | ✅ GCP responsibility |
| **Workstation Security** | Hospital VPN required, device management | ✅ Hospital IT manages |
| **Device Controls** | No PHI on local devices (web-based system) | ✅ Implemented |

### Technical Safeguards

| Requirement | Implementation | Status |
|-------------|----------------|--------|
| **Access Control** | Azure AD SSO, MFA, unique user identification | ✅ Implemented |
| **Audit Controls** | 10-year audit logging, monthly reviews | ✅ Implemented |
| **Integrity Controls** | De-identification validation, error logging | ✅ Implemented |
| **Transmission Security** | TLS 1.3, VPC isolation, no public access | ✅ Implemented |

---

## De-identification Validation

### Safe Harbor Method

The system implements **HIPAA Safe Harbor** de-identification (§164.514(b)(2)), removing all 18 identifiers:

| Identifier | Removal Method | Validation |
|------------|----------------|------------|
| 1. Names | Removed | ✅ Verified in Epic FHIR client |
| 2. Geographic subdivisions | Zip codes → 3-digit only | ✅ Verified |
| 3. Dates | Dates → Year only | ✅ Verified |
| 4. Telephone numbers | Removed | ✅ Verified |
| 5. Fax numbers | Removed | ✅ Verified |
| 6. Email addresses | Removed | ✅ Verified |
| 7. Social Security numbers | Removed | ✅ Verified |
| 8. Medical record numbers | Research ID substituted | ✅ Verified |
| 9. Health plan numbers | Removed | ✅ Verified |
| 10. Account numbers | Removed | ✅ Verified |
| 11. Certificate/license numbers | Removed | ✅ Verified |
| 12. Vehicle identifiers | Removed | ✅ N/A (not in data) |
| 13. Device identifiers | Removed | ✅ Verified |
| 14. Web URLs | Removed | ✅ N/A (not in data) |
| 15. IP addresses | Removed | ✅ Verified |
| 16. Biometric identifiers | Hashed (SHA-256) | ✅ Verified |
| 17. Full-face photos | Removed | ✅ N/A (not in data) |
| 18. Other unique identifiers | Research ID only | ✅ Verified |

### Implementation Location

**File:** `/servers/mcp-epic/src/mcp_epic/deidentify.py`

**Key Functions:**
```python
def deidentify_patient(patient: dict) -> dict:
    """Remove all 18 HIPAA identifiers from patient resource."""
    # Implementation removes names, addresses, contact info, etc.
    # Keeps only: research ID, year of birth, gender, de-identified diagnosis codes

def deidentify_observation(observation: dict) -> dict:
    """Remove identifiers from observation (lab results, vitals)."""
    # Keeps: de-identified patient reference, observation date (year only), values

def hash_identifier(value: str) -> str:
    """Create deterministic hash for tracking (e.g., for linking records)."""
    # Uses SHA-256 for one-way hashing
```

### Validation Tests

**Automated tests:** `/servers/mcp-epic/tests/test_deidentify.py`

```python
def test_deidentification_removes_names():
    """Verify names are removed from patient data."""

def test_deidentification_removes_addresses():
    """Verify addresses are removed."""

def test_deidentification_removes_dates():
    """Verify dates reduced to year only."""

def test_deidentification_removes_contact_info():
    """Verify phone, email, etc. removed."""
```

**Manual validation:**
```bash
# Review sample de-identified data
gcloud logging read \
  'jsonPayload.event="deidentification"' \
  --limit=10 \
  --format=json

# Verify no PHI in logs
gcloud logging read \
  'jsonPayload.prompt~"name" OR jsonPayload.prompt~"phone"' \
  --limit=10
# Should return truncated prompts only, no actual PHI
```

### De-identification Success Rate

**Target:** 100% of Epic FHIR data de-identified
**Current:** Monitored via Cloud Logging metric

```bash
# Check de-identification failures (should be zero)
gcloud logging read \
  'jsonPayload.event="deidentification" AND jsonPayload.success=false' \
  --limit=50
```

---

## Access Control

### User Authentication

**Method:** Azure Active Directory (Azure AD) SSO with OAuth 2.0

**Components:**
1. **Azure AD:** Hospital's identity provider
2. **OAuth2 Proxy:** Authentication layer (deployed on Cloud Run)
3. **Multi-Factor Authentication:** Required by hospital Azure AD policy

**Login Flow:**
```
User → OAuth2 Proxy → Azure AD Login → MFA → OAuth2 Proxy → Application
```

### User Authorization

**Current (Phase 1 - Simplified):**
- All authenticated users have identical access
- Access controlled via Azure AD group: `precision-medicine-users`
- 5 authorized users (2 clinicians, 3 bioinformaticians)

**Future (Phase 2 - RBAC):**
- Clinician role: Read-only analysis results
- Bioinformatician role: Full analysis capabilities
- Admin role: User management, system configuration

### User Identification

Every action is logged with user identity:

```python
# From audit_logger.py
{
    "event": "mcp_query",
    "user_email_hash": "abc123def456",  # SHA-256 hash of email
    "user_id": "abc123def456",
    "timestamp": "2025-12-30T10:30:00Z",
    "servers": ["epic", "spatialtools"],
    "prompt_length": 245
}
```

### Network Access Control

**Requirements:**
- Hospital VPN (for remote access)
- Cloud Run ingress: Internal and Cloud Load Balancing (no public access)
- VPC isolation

**IP Whitelisting:**
- Configured at hospital VPN level
- Cloud Run accessible only via VPC connector

---

## Audit Controls

### Audit Logging Requirements

**HIPAA §164.312(b):** "Implement hardware, software, and/or procedural mechanisms that record and examine activity in information systems that contain or use electronic protected health information."

### What is Logged

1. **User Access Events:**
   - Login (user_login)
   - Logout (user_logout)
   - Session start/end

2. **Data Access Events:**
   - MCP query (mcp_query) - with truncated prompt
   - MCP response (mcp_response) - with metadata only
   - Epic FHIR call (epic_fhir_call)
   - De-identification operation (deidentification)

3. **System Events:**
   - Server errors
   - Configuration changes
   - Epic FHIR connection failures

4. **Security Events:**
   - Failed login attempts (Azure AD logs)
   - Unauthorized access attempts
   - De-identification failures

### Log Retention

**Requirement:** HIPAA minimum 6 years

**Implementation:** 10 years (3,650 days)

```bash
# Configured in setup-audit-logging.sh
gcloud logging buckets create hipaa-audit-logs \
  --location=us-central1 \
  --retention-days=3650 \
  --description="HIPAA audit logs - 10 year retention"
```

### Log Review Procedures

**Daily:**
- Automated alerts for critical issues (server errors, de-ID failures)

**Weekly:**
- Review error logs
- Check for unusual access patterns

**Monthly:**
- Comprehensive audit log review
- User access audit
- De-identification validation
- Export logs for compliance records

### Audit Log Access

**Who can access:**
- Hospital Privacy Officer
- IT Security team
- Authorized auditors
- System administrators (read-only)

**How to access:**
```bash
# Via gcloud CLI
gcloud logging read 'resource.type="cloud_run_revision"' \
  --limit=100 \
  --project=<PROJECT_ID>

# Via Cloud Console
https://console.cloud.google.com/logs/query?project=<PROJECT_ID>
```

### Log Integrity

**Protection mechanisms:**
- Cloud Logging automatically signs logs (tamper-evident)
- 10-year retention prevents deletion
- Access logged (audit of audits)
- Exported monthly for offline storage

---

## Encryption

### Encryption at Rest

**All data encrypted using AES-256:**

| Data Type | Storage | Encryption |
|-----------|---------|------------|
| **Secrets** | Secret Manager | Google-managed keys (AES-256) |
| **Audit Logs** | Cloud Logging | Google-managed keys (AES-256) |
| **Container Images** | Container Registry | Google-managed keys (AES-256) |
| **Patient Data** | Hospital GCS buckets | Google-managed keys (customer-managed optional) |

**Validation:**
```bash
# Verify encryption enabled (default for all GCP services)
gcloud storage buckets describe gs://<HOSPITAL_BUCKET> \
  --format='value(encryption)'
```

### Encryption in Transit

**All network traffic encrypted with TLS 1.3:**

| Connection | Encryption | Certificate |
|------------|------------|-------------|
| **User → OAuth2 Proxy** | TLS 1.3 | Google-managed |
| **OAuth2 Proxy → Azure AD** | TLS 1.3 | Microsoft-managed |
| **OAuth2 Proxy → Streamlit** | TLS 1.3 | Google-managed |
| **Streamlit → Claude API** | TLS 1.3 | Anthropic-managed |
| **Claude API → MCP Servers** | TLS 1.3 | Google-managed |
| **MCP Servers → Epic FHIR** | TLS 1.3 | Hospital-managed |

**Validation:**
```bash
# Test TLS configuration
curl -vI https://oauth2-proxy-streamlit-{hash}.run.app 2>&1 | grep "TLS"
# Should show TLS 1.3
```

### Key Management

**Secret Manager:**
- Keys automatically rotated by Google
- Access controlled via IAM
- All access logged

**Manual key rotation:** Annually for Epic FHIR and Azure AD credentials

---

## PHI Handling Procedures

### What Qualifies as PHI

**Protected Health Information (PHI)** includes:
- Patient names, addresses, contact information
- Medical record numbers
- Diagnosis and treatment information **when linked to patient identity**

### What is NOT PHI (De-identified Data)

After Safe Harbor de-identification, data is NOT PHI:
- Research patient ID (e.g., "RESEARCH-PAT001")
- Year of birth (not exact date)
- De-identified diagnosis codes
- Lab values without identifiers
- Spatial transcriptomics data

### User Responsibilities

**DO:**
- ✅ Use research patient IDs only
- ✅ Keep queries focused on scientific analysis
- ✅ Report any suspected PHI exposure immediately
- ✅ Complete HIPAA training annually

**DO NOT:**
- ❌ Include patient names in queries
- ❌ Try to re-identify patients from de-identified data
- ❌ Share de-identified data outside authorized users
- ❌ Export data to unsecured locations

### PHI Breach Procedure

**If PHI is suspected to be exposed:**

1. **IMMEDIATE:** Notify Privacy Officer and IT Security
   - Privacy Officer: l.thompson@hospital.org / (555) 345-6789
   - IT Security: d.brown@hospital.org / (555) 456-7890

2. **Within 1 hour:** Document incident
   - What data was exposed
   - How exposure occurred
   - Who had access
   - When exposure was discovered

3. **Within 4 hours:** Contain breach
   - Shut down affected system if needed
   - Revoke access
   - Secure exposed data

4. **Within 24 hours:** Investigate
   - Review audit logs
   - Determine extent of exposure
   - Identify affected patients (if any)

5. **Within 60 days:** Report (if required)
   - HIPAA breach notification rules apply if PHI confirmed
   - Hospital Privacy Officer coordinates reporting

### De-identification Failure Response

**If de-identification fails:**

```bash
# Alert is automatically triggered
# Check failure reason
gcloud logging read \
  'jsonPayload.event="deidentification" AND jsonPayload.success=false' \
  --limit=1 \
  --format=json
```

**Response:**
1. System automatically blocks request (fails closed)
2. Alert sent to development team
3. Epic FHIR connection paused until fixed
4. Incident logged for review

---

## Compliance Validation

### Monthly Compliance Checklist

**Performed by:** System Administrator + Privacy Officer
**Date:** Last business day of each month

```markdown
## Monthly HIPAA Compliance Review - [Month Year]

### Access Control
- [ ] Review Azure AD group members (should be 5 users)
- [ ] Verify all users have completed HIPAA training (current year)
- [ ] Check for unauthorized access attempts in audit logs
- [ ] Verify VPN access is working correctly

### De-identification
- [ ] Run de-identification test suite (100% pass required)
- [ ] Review de-identification failure logs (should be zero)
- [ ] Spot-check sample de-identified data for PHI
- [ ] Verify Epic FHIR connection using correct endpoint

### Audit Logging
- [ ] Verify audit logs are being created (check last 24 hours)
- [ ] Check log retention is set to 10 years (3650 days)
- [ ] Export logs for offline archive
- [ ] Review unusual activity or errors

### Encryption
- [ ] Verify TLS 1.3 is enforced on all services
- [ ] Check Secret Manager access (no unauthorized access)
- [ ] Review GCS bucket encryption settings

### Incidents
- [ ] No PHI breaches this month (check incident log)
- [ ] No de-identification failures (check failure log)
- [ ] No security incidents (check security log)
- [ ] All incidents properly documented and resolved

### Training & Awareness
- [ ] All users completed annual HIPAA training
- [ ] User guide updated with any changes
- [ ] No user policy violations

### Signature
Reviewed by: _________________ Date: _________
Privacy Officer: _________________ Date: _________
```

### Quarterly Audits

**Additional items (quarterly):**
- Penetration testing (by hospital IT security)
- Disaster recovery test
- User satisfaction survey
- Cost and usage analysis

### Annual Certification

**Requirements:**
- Privacy Officer certifies HIPAA compliance
- IT Security certifies technical safeguards
- PI certifies research use is appropriate
- External audit (if required by hospital)

---

## Incident Response

### Incident Classification

| Severity | Description | Examples |
|----------|-------------|----------|
| **P0 - Critical** | PHI breach or suspected breach | Identifiable data exposed, unauthorized access to Epic |
| **P1 - High** | De-identification failure | Safe Harbor method failed, system exposed PHI |
| **P2 - Medium** | Security control failure | MFA bypassed, audit logging failed |
| **P3 - Low** | Minor compliance issue | User missed training, incorrect access level |

### Response Procedures

**P0 - Critical (PHI Breach):**
1. **0-15 min:** Contain breach, notify Privacy Officer + IT Security + CTO
2. **15-60 min:** Document incident, assess scope
3. **1-4 hours:** Investigate, review audit logs
4. **4-24 hours:** Remediate, implement fixes
5. **24-60 days:** Report per HIPAA breach notification rules

**P1 - High (De-ID Failure):**
1. **0-1 hour:** System auto-blocks, notify development team
2. **1-4 hours:** Investigate failure, test fix
3. **4-24 hours:** Deploy fix, validate de-identification
4. **24-48 hours:** Document root cause, update procedures

### Breach Notification

**Thresholds:**
- Breach affecting <500 individuals: Annual reporting to HHS
- Breach affecting ≥500 individuals: Immediate notification (within 60 days)
- Breach of unsecured PHI: Individual notification required

**Coordination:** Hospital Privacy Officer manages all breach notifications per hospital policy

---

## Compliance Documentation

### Required Documentation (On File)

- [x] Business Associate Agreement (BAA) with Google Cloud
- [x] HIPAA Security Risk Assessment
- [x] Policies and Procedures Manual (this document)
- [x] User Training Materials
- [x] Incident Response Plan
- [x] Disaster Recovery Plan
- [x] Audit Logging Procedures
- [x] De-identification Validation Report

### External References

- **HIPAA Security Rule:** 45 CFR §164.308-318
- **HIPAA Privacy Rule:** 45 CFR §164.500-534
- **De-identification Standard:** 45 CFR §164.514(b)
- **HHS Guidance:** https://www.hhs.gov/hipaa/

---

**Document History:**
- v1.0 (2025-12-30): Initial HIPAA compliance documentation
- Next Review: 2026-01-30 (monthly), 2026-03-30 (quarterly)

**Approval:**
- Privacy Officer: _________________ Date: _________
- IT Security Lead: _________________ Date: _________
- Principal Investigator: _________________ Date: _________

**Related Documents:**
- [Operations Manual](../hospital-deployment/OPERATIONS_MANUAL.md)
- [Admin Guide](../hospital-deployment/ADMIN_GUIDE.md)
- [Audit Log Guide](../hospital-deployment/AUDIT_LOG_GUIDE.md)
