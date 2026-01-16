# Runbook: MCP Server Down

**Severity:** P1 (High)
**Response Time:** 1 hour
**Owner:** Development Team + Hospital IT

---

## Symptoms

- Users report "Server timeout" or "Connection refused" errors
- Cloud Monitoring alert: "MCP Server Down"
- High 5xx error rate in logs
- Server not responding to health checks

---

## Diagnosis

### Step 1: Identify Affected Server

```bash
# Check all server statuses
for server in mcp-fgbio mcp-multiomics mcp-spatialtools mcp-epic \
              mcp-tcga mcp-openimagedata mcp-seqera mcp-huggingface \
              mcp-deepcell mcp-mockepic; do
  echo "=== $server ==="
  gcloud run services describe $server \
    --region=us-central1 \
    --project=<PROJECT_ID> \
    --format='value(status.conditions[0].status,status.conditions[0].message)'
done
```

**Expected Output:**
- Healthy server: `True` with no error message
- Down server: `False` or `Unknown` with error message

### Step 2: Check Recent Logs

```bash
# View last 50 logs for affected server
gcloud run services logs read mcp-<SERVER_NAME> \
  --region=us-central1 \
  --limit=50 \
  | grep -i error

# Check for crash loops
gcloud logging read \
  'resource.type="cloud_run_revision"
   AND resource.labels.service_name="mcp-<SERVER_NAME>"
   AND severity>=ERROR
   AND timestamp>="'$(date -d '1 hour ago' -I)'T00:00:00Z"' \
  --limit=20
```

**Common Error Patterns:**
- `ImportError:` Missing Python dependencies
- `PermissionError:` Service account lacks access to Secret Manager/GCS
- `ConnectionError:` Cannot reach Epic FHIR / upstream service
- `MemoryError:` Container out of memory (OOM)

### Step 3: Check Service Configuration

```bash
# View current configuration
gcloud run services describe mcp-<SERVER_NAME> \
  --region=us-central1 \
  --format=yaml > /tmp/server_config.yaml

# Check for recent changes
gcloud run revisions list \
  --service=mcp-<SERVER_NAME> \
  --region=us-central1 \
  --format='table(name,creationTimestamp,status)'
```

---

## Resolution

### Option 1: Restart Container (Simple Issues)

```bash
# Force new revision by updating environment variable
gcloud run services update mcp-<SERVER_NAME> \
  --update-env-vars=RESTART_TIMESTAMP=$(date +%s) \
  --region=us-central1

# Monitor restart
gcloud run services logs read mcp-<SERVER_NAME> \
  --region=us-central1 \
  --follow
```

**Success criteria:** Server starts responding within 2 minutes

---

### Option 2: Rollback to Previous Revision

**Use when:** Recent deployment caused the issue

```bash
# 1. List revisions
gcloud run revisions list \
  --service=mcp-<SERVER_NAME> \
  --region=us-central1 \
  --format='table(name,creationTimestamp,status)'

# 2. Identify last working revision (usually 2nd in list)
LAST_WORKING_REVISION="mcp-<SERVER_NAME>-00042-xyz"

# 3. Rollback
gcloud run services update-traffic mcp-<SERVER_NAME> \
  --to-revisions=$LAST_WORKING_REVISION=100 \
  --region=us-central1

# 4. Verify
gcloud run services describe mcp-<SERVER_NAME> \
  --region=us-central1 \
  --format='value(status.traffic[0].revisionName)'

# 5. Test
# Send test query via Streamlit or curl
```

**Success criteria:** Server responds with expected behavior

---

### Option 3: Fix Permission Issues

**Use when:** Logs show `PermissionError`

```bash
# Check service account
SERVICE_ACCOUNT=$(gcloud run services describe mcp-<SERVER_NAME> \
  --region=us-central1 \
  --format='value(spec.template.spec.serviceAccountName)')

echo "Service Account: $SERVICE_ACCOUNT"

# Check if service account can access secrets
for secret in anthropic-api-key epic-client-secret; do
  echo "Checking $secret..."
  gcloud secrets get-iam-policy $secret \
    --project=<PROJECT_ID> \
    | grep "$SERVICE_ACCOUNT"
done

# Grant access if missing
gcloud secrets add-iam-policy-binding epic-client-secret \
  --member=serviceAccount:$SERVICE_ACCOUNT \
  --role=roles/secretmanager.secretAccessor
```

---

### Option 4: Fix Memory Issues

**Use when:** Logs show `MemoryError` or OOM kills

```bash
# Increase memory allocation
gcloud run services update mcp-<SERVER_NAME> \
  --memory=8Gi \
  --region=us-central1

# May also need to increase CPU
gcloud run services update mcp-<SERVER_NAME> \
  --cpu=4 \
  --region=us-central1
```

---

### Option 5: Fix Epic Connection Issues

**Use when:** Epic FHIR server is unreachable

```bash
# Test Epic endpoint manually
EPIC_ENDPOINT=$(gcloud secrets versions access latest \
  --secret=epic-fhir-endpoint \
  --project=<PROJECT_ID>)

curl -I "$EPIC_ENDPOINT/metadata"

# If Epic is down, temporarily switch to mock
# See: epic-connection-failure.md runbook
```

---

## Verification

### 1. Health Check

```bash
# Get service URL
SERVICE_URL=$(gcloud run services describe mcp-<SERVER_NAME> \
  --region=us-central1 \
  --format='value(status.url)')

# Test health endpoint
curl -I "$SERVICE_URL/health" || curl -I "$SERVICE_URL/"
```

**Expected:** HTTP 200 response

### 2. Test Query

**Via Streamlit:**
1. Log in to Streamlit UI
2. Select affected server
3. Send simple test query: "List available tools"
4. Verify response is received

**Via curl:**
```bash
# Get auth token (requires gcloud auth)
TOKEN=$(gcloud auth print-access-token)

# Test MCP server
curl -X POST "$SERVICE_URL/sse" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"method": "tools/list", "params": {}}'
```

### 3. Monitor for Errors

```bash
# Watch logs for 5 minutes
gcloud run services logs read mcp-<SERVER_NAME> \
  --region=us-central1 \
  --follow &

sleep 300
kill %1

# No errors should appear
```

---

## Communication

### During Incident

**Send to users (if P1):**
```
Subject: [Incident] MCP Server Temporarily Unavailable

We're experiencing an issue with the <SERVER_NAME> server.

Impact: Queries using this server may fail or timeout
Status: Investigating
ETA: Resolution within 1 hour

We'll send an update when the issue is resolved.

For urgent analysis needs, try using alternative servers if applicable.

- Hospital IT Team
```

### After Resolution

```
Subject: [Resolved] MCP Server Back Online

The <SERVER_NAME> server has been restored to full functionality.

Root cause: [brief description]
Resolution: [what was done]
Downtime: [duration]

No further action needed. System is operating normally.

Thank you for your patience.

- Hospital IT Team
```

---

## Post-Incident

### 1. Root Cause Analysis

**Template:**
```markdown
## Incident Report: MCP Server Down

**Date:** YYYY-MM-DD
**Duration:** HH:MM
**Severity:** P1
**Affected Server:** mcp-<SERVER_NAME>
**Users Impacted:** X

### Timeline
- HH:MM - Incident detected via alert
- HH:MM - Investigation started
- HH:MM - Root cause identified
- HH:MM - Fix applied
- HH:MM - Service restored
- HH:MM - Verification complete

### Root Cause
[Detailed explanation]

### Contributing Factors
- [Factor 1]
- [Factor 2]

### Resolution
[What was done to fix]

### Preventive Measures
- [ ] Action item 1
- [ ] Action item 2
- [ ] Update monitoring/alerting
- [ ] Update runbook if needed
```

### 2. Update Monitoring

**If this issue should have been caught earlier:**

```bash
# Create alert for this specific error
gcloud alpha monitoring policies create \
  --notification-channels=<CHANNEL_ID> \
  --display-name="<SERVER_NAME> - <ERROR_TYPE>" \
  --condition-filter='...' \
  --project=<PROJECT_ID>
```

### 3. Document Lessons Learned

- What went well?
- What could be improved?
- Any recurring patterns?

---

## Related Runbooks

- [Epic Connection Failure](epic-connection-failure.md)
- [SSO Issues](sso-issues.md)

---

**Document History:**
- v1.0 (2025-12-30): Initial runbook
- Last Incident: N/A
- Last Updated: 2025-12-30
