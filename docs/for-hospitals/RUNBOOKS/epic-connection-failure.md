# Runbook: Epic FHIR Connection Failure

**Severity:** P1 (High)
**Response Time:** 1 hour
**Owner:** Development Team + Hospital IT

---

## Symptoms

- Users report "Cannot retrieve patient data from Epic"
- Queries involving `epic` server fail
- Logs show Epic FHIR connection errors
- Alert: "High Epic FHIR Failures"

---

## Diagnosis

### Step 1: Check Epic FHIR Server Status

```bash
# Get Epic endpoint
EPIC_ENDPOINT=$(gcloud secrets versions access latest \
  --secret=epic-fhir-endpoint)

echo "Epic FHIR Endpoint: $EPIC_ENDPOINT"

# Test connectivity
curl -I "$EPIC_ENDPOINT/metadata"
```

**Expected:** HTTP 200 response
**If fails:** Epic server may be down (contact hospital IT)

### Step 2: Check Epic FHIR Logs

```bash
# View recent Epic FHIR calls
gcloud logging read \
  'jsonPayload.event="epic_fhir_call"
   AND timestamp>="'$(date -d '1 hour ago' -I)'T00:00:00Z"' \
  --limit=50 \
  --format='table(timestamp, jsonPayload.status, jsonPayload.response_code, jsonPayload.resource_type)'

# Count failures
gcloud logging read \
  'jsonPayload.event="epic_fhir_call"
   AND jsonPayload.status="error"
   AND timestamp>="'$(date -d '1 hour ago' -I)'T00:00:00Z"' \
  --format=json | jq '. | length'
```

### Step 3: Test Epic Credentials

```bash
# Get credentials
EPIC_CLIENT_ID=$(gcloud secrets versions access latest --secret=epic-client-id)
EPIC_CLIENT_SECRET=$(gcloud secrets versions access latest --secret=epic-client-secret)

# Test OAuth token
curl -X POST "$EPIC_ENDPOINT/oauth2/token" \
  -d "grant_type=client_credentials" \
  -d "client_id=$EPIC_CLIENT_ID" \
  -d "client_secret=$EPIC_CLIENT_SECRET" \
  -d "scope=patient/*.read"
```

**Expected:** JSON with `access_token`
**If fails:** Credentials may be invalid or expired

---

## Resolution

### Option 1: Epic Server is Down (Temporary)

**Use mcp-mockepic as fallback:**

**In Streamlit UI (instruct users):**
1. Deselect `epic` server
2. Select `mockepic` server instead
3. Continue analysis with synthetic data

**Note:** Mock data will not match real patient IDs

---

### Option 2: Epic Credentials Expired

```bash
# Contact hospital IT to get new credentials
# Once obtained, update Secret Manager:

echo -n "new-client-id" | gcloud secrets versions add epic-client-id --data-file=-
echo -n "new-client-secret" | gcloud secrets versions add epic-client-secret --data-file=-

# Restart mcp-epic server to pick up new credentials
gcloud run services update mcp-epic \
  --update-env-vars=RESTART_TIMESTAMP=$(date +%s) \
  --region=us-central1

# Verify
sleep 30
gcloud run services logs read mcp-epic --limit=20 | grep -i "epic"
```

---

### Option 3: Epic Endpoint Changed

**If hospital IT notifies Epic endpoint URL has changed:**

```bash
# Update Epic endpoint
echo -n "https://new-epic-endpoint.hospital.org/api/FHIR/R4/" | \
  gcloud secrets versions add epic-fhir-endpoint --data-file=-

# Restart mcp-epic
gcloud run services update mcp-epic \
  --update-env-vars=RESTART_TIMESTAMP=$(date +%s) \
  --region=us-central1

# Test
EPIC_ENDPOINT=$(gcloud secrets versions access latest --secret=epic-fhir-endpoint)
curl -I "$EPIC_ENDPOINT/metadata"
```

---

### Option 4: Network/Firewall Issue

**If Epic is up but unreachable:**

```bash
# Check if VPC connector can reach Epic
# May need hospital IT to whitelist Cloud Run IP ranges

# Get Cloud Run NAT IPs
gcloud compute networks vpc-access connectors describe mcp-connector \
  --region=us-central1 \
  --format='value(ipCidrRange)'

# Provide to hospital IT for whitelisting in Epic firewall
```

---

### Option 5: Rate Limiting

**If Epic returns 429 Too Many Requests:**

```bash
# Check request frequency
gcloud logging read \
  'jsonPayload.event="epic_fhir_call"
   AND timestamp>="'$(date -d '1 hour ago' -I)'T00:00:00Z"' \
  --format=json | jq '. | length'

# If >100/hour, may need to:
# 1. Implement caching in mcp-epic server
# 2. Request higher rate limit from hospital IT
# 3. Add retry logic with exponential backoff
```

---

## Verification

### 1. Test Epic Connection

```bash
# Get Epic test patient
TEST_PATIENT_ID="RESEARCH-PAT001"

# Test via mcp-epic server
# Use Streamlit or curl to query:
# "Get demographics for patient RESEARCH-PAT001"

# Check logs
gcloud logging read \
  'jsonPayload.event="epic_fhir_call"
   AND jsonPayload.resource_id="'$TEST_PATIENT_ID'"' \
  --limit=5
```

**Expected:** `status: "success"`, `response_code: 200`

### 2. Verify De-identification

```bash
# Ensure de-identification still working
gcloud logging read \
  'jsonPayload.event="deidentification"
   AND timestamp>="'$(date -d '5 minutes ago' -I)'T00:00:00Z"' \
  --limit=10 \
  --format='table(timestamp, jsonPayload.success)'
```

**Expected:** All `success: true`

---

## Communication

### To Users

```
Subject: Epic FHIR Connection Issue

We're experiencing connectivity issues with the Epic FHIR server.

Impact: Queries requiring Epic clinical data may fail
Workaround: Use 'mockepic' server for demonstration/testing
Status: [Working with hospital IT / Restored]

Updates will be sent as we learn more.

- Hospital IT Team
```

### To Hospital IT

```
Subject: Epic FHIR API Connection Issue - Precision Medicine MCP

We're unable to connect to the Epic FHIR API from our Precision Medicine MCP system.

Epic Endpoint: [endpoint URL]
Client ID: [client ID]
Error: [error message from logs]

Possible causes:
- Epic server maintenance/downtime
- Credentials expired
- Firewall blocking Cloud Run IPs
- Rate limiting

Cloud Run NAT IP range: [IP range]

Can you please check:
1. Epic FHIR API status
2. Our client credentials validity
3. Firewall rules for our IP range

Thank you,
Precision Medicine MCP Team
```

---

## Post-Incident

### Update Documentation

**If Epic endpoint or credentials changed:**

Update `/docs/for-hospitals/OPERATIONS_MANUAL.md`:
- Epic FHIR endpoint URL
- Contact for Epic credentials
- Any new requirements

### Monitor Epic Health

**Create proactive monitoring:**

```bash
# Cron job to test Epic connectivity every 15 minutes
# Alert if 3 consecutive failures

# Add to Cloud Scheduler:
gcloud scheduler jobs create http epic-health-check \
  --location=us-central1 \
  --schedule="*/15 * * * *" \
  --uri="https://monitoring.hospital.internal/epic-health" \
  --http-method=GET
```

---

## Related Runbooks

- [Server Down](server-down.md)
- [SSO Issues](sso-issues.md)

---

**Document History:**
- v1.0 (2025-12-30): Initial runbook
- Last Incident: N/A
