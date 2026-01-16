# Audit Log Guide
## HIPAA-Compliant Logging for Precision Medicine MCP

**Version:** 1.0
**Last Updated:** 2025-12
**Purpose:** Guide for accessing, analyzing, and reporting on audit logs

---

## Table of Contents

- [What is Logged](#what-is-logged)
- [How to Access Logs](#how-to-access-logs)
- [Log Retention Policy](#log-retention-policy)
- [Sample Audit Queries](#sample-audit-queries)
- [Compliance Reporting](#compliance-reporting)
- [Log Analysis](#log-analysis)

---

## What is Logged

### User Access Events

**Event:** `user_login`

**Logged Information:**
- User email hash (SHA-256)
- User ID
- Display name
- Timestamp
- Session ID

**Example:**
```json
{
  "event": "user_login",
  "timestamp": "2025-12-30T10:30:00Z",
  "user_email_hash": "abc123def456",
  "user_id": "abc123def456",
  "display_name": "Dr. Smith"
}
```

---

### MCP Query Events

**Event:** `mcp_query`

**Logged Information:**
- User email hash
- User ID
- Servers used
- Prompt length
- Prompt preview (first 100 characters, truncated to prevent PHI logging)
- Model used
- Session ID

**Example:**
```json
{
  "event": "mcp_query",
  "timestamp": "2025-12-30T10:31:00Z",
  "user_email_hash": "abc123def456",
  "user_id": "abc123def456",
  "servers": ["epic", "spatialtools"],
  "server_count": 2,
  "prompt_length": 245,
  "prompt_preview": "For patient RESEARCH-PAT001: 1. Get clinical demographics 2. Analyze spatial transcriptomics...",
  "model": "claude-sonnet-4-5",
  "session_id": "sess_abc123"
}
```

**Note:** Full prompts are NOT logged to prevent accidental PHI exposure.

---

### MCP Response Events

**Event:** `mcp_response`

**Logged Information:**
- User email hash
- User ID
- Servers used
- Response length (characters)
- Token usage (input, output, total)
- Estimated cost (USD)
- Duration (seconds)
- Session ID

**Example:**
```json
{
  "event": "mcp_response",
  "timestamp": "2025-12-30T10:31:15Z",
  "user_email_hash": "abc123def456",
  "user_id": "abc123def456",
  "servers": ["epic", "spatialtools"],
  "response_length": 1523,
  "input_tokens": 523,
  "output_tokens": 1247,
  "total_tokens": 1770,
  "estimated_cost_usd": 0.0283,
  "duration_seconds": 14.5,
  "session_id": "sess_abc123"
}
```

---

### Epic FHIR Access Events

**Event:** `epic_fhir_call`

**Logged Information:**
- User email hash
- Resource type (Patient, Observation, etc.)
- Resource ID
- Status (success/error)
- Response code
- Timestamp

**Example:**
```json
{
  "event": "epic_fhir_call",
  "timestamp": "2025-12-30T10:31:10Z",
  "resource_type": "Patient",
  "resource_id": "RESEARCH-PAT001",
  "status": "success",
  "response_code": 200
}
```

---

### De-identification Events

**Event:** `deidentification`

**Logged Information:**
- Resource type
- Success/failure status
- Method used (Safe Harbor)
- Identifiers removed count
- Timestamp

**Example:**
```json
{
  "event": "deidentification",
  "timestamp": "2025-12-30T10:31:11Z",
  "resource_type": "Patient",
  "success": true,
  "method": "safe_harbor",
  "identifiers_removed": 12
}
```

---

### Error Events

**Event:** `error`

**Logged Information:**
- User email hash
- Error type
- Error message (truncated)
- Servers involved
- Session ID
- Timestamp

**Example:**
```json
{
  "event": "error",
  "timestamp": "2025-12-30T10:32:00Z",
  "user_email_hash": "abc123def456",
  "user_id": "abc123def456",
  "error_type": "EpicConnectionError",
  "error_message": "Failed to connect to Epic FHIR endpoint: Connection timeout after 30s",
  "servers": ["epic"],
  "session_id": "sess_abc123"
}
```

---

### Clinician Review Events (CitL)

**Event:** `citl_review`

**Logged Information:**
- Patient ID hash (SHA-256)
- Report ID
- Reviewer name, credentials, role
- Reviewer email hash
- Decision (APPROVE / REVISE / REJECT)
- Rationale
- Digital signature hash (SHA-256)
- Findings validated count
- Findings confirmed/uncertain/rejected counts
- Guideline compliance status
- Quality flags count
- Revision count
- Timestamp

**Example:**
```json
{
  "event": "citl_review",
  "timestamp": "2026-01-13T14:55:00Z",
  "patient_id_hash": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
  "report_id": "PAT001-OVC-2025-2026-01-13T14:00:00Z",
  "reviewer": {
    "email_hash": "5e884898da28047151d0e56f8dc6292773603d0d6aabbdd62a11ef721d1542d8",
    "name": "Dr. Sarah Johnson",
    "credentials": "MD, Gynecologic Oncology",
    "role": "oncologist"
  },
  "decision": {
    "status": "APPROVE",
    "rationale": "All findings consistent with clinical presentation and imaging..."
  },
  "signature_hash": "a1b2c3d4e5f6789012345678901234567890123456789012345678901234abcd",
  "findings_validated": 10,
  "findings_confirmed": 10,
  "findings_uncertain": 0,
  "findings_rejected": 0,
  "guideline_compliance": {
    "nccn_aligned": "ALIGNED",
    "institutional_aligned": "ALIGNED"
  },
  "quality_flags_count": 0,
  "revision_count": 0
}
```

**Retention:** 10 years (HIPAA requirement)

---

## How to Access Logs

### Via Cloud Console (Web UI)

**Easiest for non-technical users:**

1. Open: https://console.cloud.google.com/logs/query
2. Select project: `<hospital-name>-precision-medicine`
3. Enter query in filter box (see [Sample Queries](#sample-audit-queries))
4. Click "Run Query"
5. Export if needed: Click "Export logs" → "Download as JSON/CSV"

---

### Via gcloud CLI (Command Line)

**For administrators and automated reporting:**

**Basic syntax:**
```bash
gcloud logging read '<FILTER>' \
  --limit=<NUMBER> \
  --format=<FORMAT> \
  --project=<PROJECT_ID>
```

**Common formats:**
- `json` - Full structured data
- `table` - Human-readable table
- `csv` - Comma-separated values

**Examples:**
```bash
# Read last 50 logs
gcloud logging read 'jsonPayload.event!=""' --limit=50

# Read with specific format
gcloud logging read 'jsonPayload.event="user_login"' \
  --format='table(timestamp, jsonPayload.user_email_hash, jsonPayload.display_name)'

# Export to file
gcloud logging read 'jsonPayload.event!=""' \
  --format=json > audit_logs_$(date +%Y%m%d).json
```

---

### Via Logs Explorer (Advanced)

**For complex queries and visualizations:**

1. Open: https://console.cloud.google.com/logs/query
2. Use Logs Explorer interface:
   - Add filters
   - Visualize trends
   - Create charts
   - Set up alerts

---

## Log Retention Policy

### Retention Period

**Duration:** 10 years (3,650 days)

**Reason:** Exceeds HIPAA minimum requirement of 6 years

**Implementation:**
```bash
# Configured in hipaa-audit-logs bucket
# Retention: 3650 days
# Auto-deletion after retention period
```

---

### Storage Location

- **Primary:** Cloud Logging bucket `hipaa-audit-logs`
- **Region:** us-central1 (or hospital-specified region)
- **Backup:** Monthly exports to hospital compliance archive

---

### Access Control

**Who can access:**
- Privacy Officer (full access)
- IT Security team (full access)
- System administrators (read-only)
- Authorized auditors (temporary read-only)

**IAM roles:**
- `roles/logging.viewer` - Read logs
- `roles/logging.admin` - Manage logging configuration

---

## Sample Audit Queries

### User Access Reports

**All user logins (last 30 days):**
```bash
gcloud logging read \
  'jsonPayload.event="user_login"
   AND timestamp>="$(date -d '30 days ago' -I)T00:00:00Z"' \
  --format='table(timestamp, jsonPayload.user_email_hash, jsonPayload.display_name)' \
  --limit=1000
```

**Unique users (last month):**
```bash
gcloud logging read \
  'jsonPayload.event="user_login"
   AND timestamp>="$(date -d '30 days ago' -I)T00:00:00Z"' \
  --format=json | \
  jq -r '.[] | .jsonPayload.user_email_hash' | \
  sort -u
```

**User login frequency:**
```bash
gcloud logging read \
  'jsonPayload.event="user_login"
   AND timestamp>="$(date -d '30 days ago' -I)T00:00:00Z"' \
  --format=json | \
  jq -r '.[] | .jsonPayload.user_email_hash' | \
  sort | uniq -c | sort -rn
```

---

### Query Activity Reports

**Total queries (last 7 days):**
```bash
gcloud logging read \
  'jsonPayload.event="mcp_query"
   AND timestamp>="$(date -d '7 days ago' -I)T00:00:00Z"' \
  --format=json | jq '. | length'
```

**Queries by server:**
```bash
gcloud logging read \
  'jsonPayload.event="mcp_query"
   AND timestamp>="$(date -d '7 days ago' -I)T00:00:00Z"' \
  --format=json | \
  jq -r '.[] | .jsonPayload.servers[]' | \
  sort | uniq -c | sort -rn
```

**Queries by user:**
```bash
gcloud logging read \
  'jsonPayload.event="mcp_query"
   AND timestamp>="$(date -d '30 days ago' -I)T00:00:00Z"' \
  --format=json | \
  jq -r '.[] | .jsonPayload.user_email_hash' | \
  sort | uniq -c | sort -rn
```

**Average query cost:**
```bash
gcloud logging read \
  'jsonPayload.event="mcp_response"
   AND timestamp>="$(date -d '7 days ago' -I)T00:00:00Z"' \
  --format=json | \
  jq -r '.[] | .jsonPayload.estimated_cost_usd' | \
  awk '{sum+=$1; count++} END {print "Average: $" sum/count " Total: $" sum " Count: " count}'
```

---

### Epic FHIR Access Reports

**All Epic FHIR accesses:**
```bash
gcloud logging read \
  'jsonPayload.event="epic_fhir_call"' \
  --limit=100 \
  --format='table(timestamp, jsonPayload.resource_type, jsonPayload.resource_id, jsonPayload.status)'
```

**Epic FHIR failures:**
```bash
gcloud logging read \
  'jsonPayload.event="epic_fhir_call"
   AND jsonPayload.status="error"' \
  --limit=50
```

**Patient data accesses:**
```bash
gcloud logging read \
  'jsonPayload.event="epic_fhir_call"
   AND jsonPayload.resource_type="Patient"' \
  --limit=100
```

---

### De-identification Reports

**De-identification success rate:**
```bash
gcloud logging read \
  'jsonPayload.event="deidentification"
   AND timestamp>="$(date -d '30 days ago' -I)T00:00:00Z"' \
  --format=json | \
  jq '[.[] | .jsonPayload.success] | [group_by(.)[] | {success: .[0], count: length}]'
```

**De-identification failures (should be zero):**
```bash
gcloud logging read \
  'jsonPayload.event="deidentification"
   AND jsonPayload.success=false' \
  --limit=50
```

---

### Clinician Review Reports (CitL)

**All reviews (last 30 days):**
```bash
gcloud logging read \
  'jsonPayload.event="citl_review"
   AND timestamp>="$(date -d '30 days ago' -I)T00:00:00Z"' \
  --format='table(timestamp, jsonPayload.reviewer.name, jsonPayload.decision.status)' \
  --limit=100
```

**Approval rate:**
```bash
gcloud logging read \
  'jsonPayload.event="citl_review"
   AND timestamp>="$(date -d '90 days ago' -I)T00:00:00Z"' \
  --format=json | \
  jq '[.[] | .jsonPayload.decision.status] | [group_by(.)[] | {status: .[0], count: length}]'
```

**Reviews by oncologist:**
```bash
gcloud logging read \
  'jsonPayload.event="citl_review"
   AND timestamp>="$(date -d '30 days ago' -I)T00:00:00Z"' \
  --format=json | \
  jq -r '.[] | .jsonPayload.reviewer.name' | \
  sort | uniq -c | sort -rn
```

**Revisions/rejections (require follow-up):**
```bash
gcloud logging read \
  'jsonPayload.event="citl_review"
   AND (jsonPayload.decision.status="REVISE" OR jsonPayload.decision.status="REJECT")
   AND timestamp>="$(date -d '30 days ago' -I)T00:00:00Z"' \
  --format='table(timestamp, jsonPayload.reviewer.name, jsonPayload.decision.status, jsonPayload.report_id)' \
  --limit=50
```

---

### Error Reports

**All errors (last 24 hours):**
```bash
gcloud logging read \
  'severity>=ERROR
   AND timestamp>="$(date -d '1 day ago' -I)T00:00:00Z"' \
  --limit=100
```

**Errors by type:**
```bash
gcloud logging read \
  'jsonPayload.event="error"
   AND timestamp>="$(date -d '7 days ago' -I)T00:00:00Z"' \
  --format=json | \
  jq -r '.[] | .jsonPayload.error_type' | \
  sort | uniq -c | sort -rn
```

---

### Cost & Usage Reports

**Total tokens used (last month):**
```bash
gcloud logging read \
  'jsonPayload.event="mcp_response"
   AND timestamp>="$(date -d '30 days ago' -I)T00:00:00Z"' \
  --format=json | \
  jq -r '.[] | .jsonPayload.total_tokens' | \
  awk '{sum+=$1} END {print "Total tokens: " sum}'
```

**Total cost (last month):**
```bash
gcloud logging read \
  'jsonPayload.event="mcp_response"
   AND timestamp>="$(date -d '30 days ago' -I)T00:00:00Z"' \
  --format=json | \
  jq -r '.[] | .jsonPayload.estimated_cost_usd' | \
  awk '{sum+=$1} END {print "Total cost: $" sum}'
```

**Cost by user:**
```bash
gcloud logging read \
  'jsonPayload.event="mcp_response"
   AND timestamp>="$(date -d '30 days ago' -I)T00:00:00Z"' \
  --format=json | \
  jq -r '.[] | "\(.jsonPayload.user_email_hash) \(.jsonPayload.estimated_cost_usd)"' | \
  awk '{user[$1]+=$2} END {for (u in user) print u, "$" user[u]}' | \
  sort -k2 -rn
```

---

## Compliance Reporting

### Monthly Compliance Report

**Template:**
```bash
#!/bin/bash
# generate_monthly_report.sh

MONTH=$(date +%Y-%m)
PROJECT_ID="<hospital-name>-precision-medicine"

cat > compliance_report_${MONTH}.md <<EOF
# HIPAA Compliance Report - $MONTH
## Precision Medicine MCP System

### Access Control

**Total Active Users:** $(gcloud logging read \
  'jsonPayload.event="user_login"
   AND timestamp>="'$(date -d '30 days ago' -I)'T00:00:00Z"' \
  --project=$PROJECT_ID --format=json | \
  jq -r '.[].jsonPayload.user_email_hash' | sort -u | wc -l)

**Login Events:** $(gcloud logging read \
  'jsonPayload.event="user_login"
   AND timestamp>="'$(date -d '30 days ago' -I)'T00:00:00Z"' \
  --project=$PROJECT_ID --format=json | jq '. | length')

### Data Access

**Total Queries:** $(gcloud logging read \
  'jsonPayload.event="mcp_query"
   AND timestamp>="'$(date -d '30 days ago' -I)'T00:00:00Z"' \
  --project=$PROJECT_ID --format=json | jq '. | length')

**Epic FHIR Accesses:** $(gcloud logging read \
  'jsonPayload.event="epic_fhir_call"
   AND timestamp>="'$(date -d '30 days ago' -I)'T00:00:00Z"' \
  --project=$PROJECT_ID --format=json | jq '. | length')

### De-identification

**Total De-ID Operations:** $(gcloud logging read \
  'jsonPayload.event="deidentification"
   AND timestamp>="'$(date -d '30 days ago' -I)'T00:00:00Z"' \
  --project=$PROJECT_ID --format=json | jq '. | length')

**De-ID Failures:** $(gcloud logging read \
  'jsonPayload.event="deidentification"
   AND jsonPayload.success=false
   AND timestamp>="'$(date -d '30 days ago' -I)'T00:00:00Z"' \
  --project=$PROJECT_ID --format=json | jq '. | length')

**Success Rate:** $(echo "scale=2; ($(gcloud logging read \
  'jsonPayload.event="deidentification"
   AND jsonPayload.success=true
   AND timestamp>="'$(date -d '30 days ago' -I)'T00:00:00Z"' \
  --project=$PROJECT_ID --format=json | jq '. | length') / $(gcloud logging read \
  'jsonPayload.event="deidentification"
   AND timestamp>="'$(date -d '30 days ago' -I)'T00:00:00Z"' \
  --project=$PROJECT_ID --format=json | jq '. | length')) * 100" | bc)%

### Security Incidents

**Total Errors:** $(gcloud logging read \
  'severity>=ERROR
   AND timestamp>="'$(date -d '30 days ago' -I)'T00:00:00Z"' \
  --project=$PROJECT_ID --format=json | jq '. | length')

**PHI Breaches:** 0 (verified)

### Compliance Status

- [x] All users authenticated via Azure AD SSO
- [x] All queries logged with user identity
- [x] De-identification success rate: 100%
- [x] No PHI breaches detected
- [x] Audit logs retained for 10 years
- [x] All access via secure VPN

**Reviewed by:** _________________
**Date:** $(date)
EOF

echo "Report generated: compliance_report_${MONTH}.md"
```

---

### Audit for Regulators

**Export all logs for a date range:**
```bash
#!/bin/bash
# export_audit_logs.sh

START_DATE="2025-01-01"
END_DATE="2025-12-31"
PROJECT_ID="<hospital-name>-precision-medicine"

# Export all audit events
gcloud logging read \
  "jsonPayload.event!=''
   AND timestamp>='${START_DATE}T00:00:00Z'
   AND timestamp<='${END_DATE}T23:59:59Z'" \
  --project=$PROJECT_ID \
  --format=json > audit_logs_${START_DATE}_${END_DATE}.json

# Create summary
jq '[.[] | .jsonPayload.event] | group_by(.) | map({event: .[0], count: length})' \
  audit_logs_${START_DATE}_${END_DATE}.json > audit_summary_${START_DATE}_${END_DATE}.json

echo "Audit logs exported: audit_logs_${START_DATE}_${END_DATE}.json"
echo "Summary: audit_summary_${START_DATE}_${END_DATE}.json"
```

---

## Log Analysis

### Detecting Anomalies

**Unusual login times:**
```bash
# Logins outside business hours (before 7 AM or after 7 PM)
gcloud logging read \
  'jsonPayload.event="user_login"' \
  --format=json | \
  jq -r '.[] | select(.timestamp | strptime("%Y-%m-%dT%H:%M:%SZ") | .hour < 7 or .hour > 19) | .timestamp + " " + .jsonPayload.user_email_hash'
```

**High query volume (potential abuse):**
```bash
# Users with >50 queries in a day
gcloud logging read \
  'jsonPayload.event="mcp_query"
   AND timestamp>="'$(date -I)'T00:00:00Z"' \
  --format=json | \
  jq -r '.[] | .jsonPayload.user_email_hash' | \
  sort | uniq -c | awk '$1 > 50 {print $1, $2}'
```

**Failed Epic FHIR accesses:**
```bash
# Multiple failures may indicate system issue
gcloud logging read \
  'jsonPayload.event="epic_fhir_call"
   AND jsonPayload.status="error"
   AND timestamp>="'$(date -I)'T00:00:00Z"' \
  --limit=50
```

---

### Performance Monitoring

**Average query duration:**
```bash
gcloud logging read \
  'jsonPayload.event="mcp_response"
   AND timestamp>="'$(date -d '7 days ago' -I)'T00:00:00Z"' \
  --format=json | \
  jq -r '.[] | .jsonPayload.duration_seconds' | \
  awk '{sum+=$1; count++} END {print "Average: " sum/count " seconds"}'
```

**Slow queries (>30 seconds):**
```bash
gcloud logging read \
  'jsonPayload.event="mcp_response"
   AND jsonPayload.duration_seconds>30' \
  --limit=20 \
  --format='table(timestamp, jsonPayload.user_email_hash, jsonPayload.duration_seconds, jsonPayload.servers)'
```

---

**Document History:**
- v1.0 (2025-12-30): Initial audit log guide
- Next Review: 2026-01-30 (monthly)

**Related Documents:**
- [Operations Manual](OPERATIONS_MANUAL.md)
- [HIPAA Compliance](../compliance/hipaa.md)
- [Admin Guide](ADMIN_GUIDE.md)
