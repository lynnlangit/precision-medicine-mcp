# Runbook: Azure AD SSO Login Issues

**Severity:** P2 (Medium) - P1 if all users affected
**Response Time:** 4 hours (1 hour if P1)
**Owner:** Hospital IT + Development Team

---

## Symptoms

- Users cannot log in to Streamlit or Jupyter
- "Access Denied" error after Azure AD login
- OAuth2 Proxy returns errors
- Users report being redirected in a loop

---

## Diagnosis

### Step 1: Identify Scope of Issue

**Questions to ask:**
- Is it affecting all users or just one user?
- Can any user log in successfully?
- Did it start after a recent change?
- Are users on VPN?

### Step 2: Check OAuth2 Proxy Status

```bash
# Check OAuth2 Proxy services
for service in oauth2-proxy-streamlit oauth2-proxy-jupyter; do
  echo "=== $service ==="
  gcloud run services describe $service \
    --region=us-central1 \
    --format='value(status.conditions[0].status,status.url)'
done

# Check OAuth2 Proxy logs
gcloud run services logs read oauth2-proxy-streamlit \
  --region=us-central1 \
  --limit=50 \
  | grep -i "error\|denied\|failed"
```

**Common errors:**
- `invalid_client` - Azure AD client ID/secret wrong
- `redirect_uri_mismatch` - Redirect URI not configured in Azure AD
- `unauthorized_client` - User not in authorized group
- `invalid_token` - Token expired or invalid

### Step 3: Test OAuth2 Flow

```bash
# Get OAuth2 Proxy URL
OAUTH_URL=$(gcloud run services describe oauth2-proxy-streamlit \
  --region=us-central1 \
  --format='value(status.url)')

# Test health endpoint
curl -I "$OAUTH_URL/ping"

# Expected: HTTP 200

# Test redirect
curl -I "$OAUTH_URL/" | grep Location
# Should redirect to Azure AD login
```

---

## Resolution

### Issue 1: User Not in Azure AD Group

**Symptom:** "Access Denied" after successful Azure AD login

**Resolution:**

```powershell
# In PowerShell (Azure AD admin required)
Connect-AzureAD

# Get group ID
$group = Get-AzureADGroup -SearchString "precision-medicine-users"

# Check if user is member
Get-AzureADGroupMember -ObjectId $group.ObjectId |
  Where-Object {$_.UserPrincipalName -eq "user@hospital.org"}

# Add user if missing
$user = Get-AzureADUser -ObjectId "user@hospital.org"
Add-AzureADGroupMember -ObjectId $group.ObjectId -RefObjectId $user.ObjectId

# Verify
Get-AzureADGroupMember -ObjectId $group.ObjectId |
  Where-Object {$_.UserPrincipalName -eq "user@hospital.org"}
```

**Test:** User should be able to log in within 5 minutes

---

### Issue 2: Redirect URI Not Configured

**Symptom:** `redirect_uri_mismatch` error in OAuth2 Proxy logs

**Resolution:**

1. **Get current OAuth2 Proxy URL:**
```bash
OAUTH_URL=$(gcloud run services describe oauth2-proxy-streamlit \
  --region=us-central1 \
  --format='value(status.url)')

REDIRECT_URI="$OAUTH_URL/oauth2/callback"
echo "Redirect URI should be: $REDIRECT_URI"
```

2. **Update Azure AD App Registration:**
   - Go to Azure Portal → App Registrations
   - Select "Precision Medicine MCP"
   - Click "Authentication"
   - Under "Redirect URIs", click "Add URI"
   - Add: `https://oauth2-proxy-streamlit-{hash}.run.app/oauth2/callback`
   - Add: `https://oauth2-proxy-jupyter-{hash}.run.app/oauth2/callback`
   - Click "Save"

**Test:** Try logging in again

---

### Issue 3: Expired Azure AD Client Secret

**Symptom:** `invalid_client` error

**Resolution:**

1. **Generate new client secret in Azure Portal:**
   - Azure Portal → App Registrations → "Precision Medicine MCP"
   - Certificates & secrets → New client secret
   - Description: "OAuth2 Proxy - 2025"
   - Expires: 12 months (or per hospital policy)
   - Click "Add"
   - **Copy secret immediately** (can't view again)

2. **Update Secret Manager:**
```bash
# Update secret
echo -n "new-client-secret-value" | \
  gcloud secrets versions add azure-ad-client-secret --data-file=-

# Restart OAuth2 Proxy services
gcloud run services update oauth2-proxy-streamlit \
  --update-env-vars=RESTART_TIMESTAMP=$(date +%s) \
  --region=us-central1

gcloud run services update oauth2-proxy-jupyter \
  --update-env-vars=RESTART_TIMESTAMP=$(date +%s) \
  --region=us-central1
```

**Test:** Should be able to log in within 2 minutes

---

### Issue 4: VPN Not Connected

**Symptom:** User cannot reach OAuth2 Proxy URL at all

**Resolution:**

**Verify URL whitelisting:**
```bash
# Get OAuth2 Proxy URLs
gcloud run services list \
  --filter="metadata.name~'oauth2-proxy'" \
  --format='value(status.url)'
```

**Instruct user:**
1. Connect to hospital VPN
2. Verify connection status
3. Try accessing OAuth2 Proxy URL
4. If still fails, check VPN whitelist includes Cloud Run domains

**Contact hospital IT if:**
- VPN whitelist needs Cloud Run domains added
- Firewall blocking *.run.app domains

---

### Issue 5: Cookie Issues / Browser Cache

**Symptom:** Login loop, inconsistent behavior

**Resolution:**

**Instruct user:**
1. Clear browser cookies for *.run.app
2. Clear browser cache
3. Try incognito/private browsing mode
4. Try different browser

**If persistent:**
```bash
# Check cookie configuration in OAuth2 Proxy
gcloud run services describe oauth2-proxy-streamlit \
  --region=us-central1 \
  --format='value(spec.template.spec.containers[0].env)'

# Verify:
# OAUTH2_PROXY_COOKIE_SECURE=true
# OAUTH2_PROXY_COOKIE_HTTPONLY=true
# OAUTH2_PROXY_COOKIE_SAMESITE=lax
```

---

### Issue 6: Azure AD Outage

**Symptom:** All users affected, Azure AD login page doesn't load

**Resolution:**

1. **Check Azure AD status:**
   - https://status.azure.com/en-us/status
   - Or ask hospital IT if Azure AD is down

2. **If Azure AD is down:**
   - No immediate fix available
   - Wait for Microsoft to restore service
   - Notify users of expected downtime

3. **Workaround (emergency only):**
   - Temporarily disable authentication for testing
   - **DO NOT DO THIS IN PRODUCTION** without Privacy Officer approval
   - Only for development/emergency troubleshooting

---

## Verification

### 1. Test User Login Flow

**Complete login test:**
1. Open OAuth2 Proxy URL in browser
2. Should redirect to Azure AD login
3. Enter hospital credentials
4. Complete MFA if prompted
5. Should redirect back to application
6. Should see application interface (Streamlit/Jupyter)

### 2. Check Audit Logs

```bash
# Verify login event was logged
gcloud logging read \
  'jsonPayload.event="user_login"
   AND timestamp>="'$(date -d '5 minutes ago' -I)'T00:00:00Z"' \
  --limit=5 \
  --format='table(timestamp, jsonPayload.user_email_hash, jsonPayload.display_name)'
```

**Expected:** Recent login event for test user

### 3. Test Multiple Users

- Ask 2-3 users to test login
- Verify all can access
- Check for any error patterns

---

## Communication

### To Affected User

```
Subject: Login Issue - Precision Medicine MCP

We've resolved the login issue you reported.

Issue: [brief description]
Resolution: [what was done]

Please try logging in again:
- Streamlit: https://oauth2-proxy-streamlit-{hash}.run.app
- Jupyter: https://oauth2-proxy-jupyter-{hash}.run.app

If you still have issues:
1. Clear browser cookies and cache
2. Make sure you're on hospital VPN
3. Try a different browser

Contact help@hospital.org if problems persist.

- Hospital IT Team
```

### To All Users (if widespread issue)

```
Subject: [Resolved] Login System Restored

The login system for Precision Medicine MCP has been restored.

Affected time: [start] - [end]
Cause: [brief description]
Resolution: [what was done]

All users should now be able to log in normally.

No action needed unless you experience issues.

- Hospital IT Team
```

---

## Post-Incident

### Update Azure AD Configuration

**Document current settings:**
```markdown
## Azure AD App Registration: Precision Medicine MCP

**Application (client) ID:** [copy from Azure Portal]
**Directory (tenant) ID:** [copy from Azure Portal]

**Redirect URIs:**
- https://oauth2-proxy-streamlit-{hash}.run.app/oauth2/callback
- https://oauth2-proxy-jupyter-{hash}.run.app/oauth2/callback

**API Permissions:**
- Microsoft Graph: User.Read, User.ReadBasic.All
- Azure AD Graph: Directory.Read.All

**Client Secret:**
- Current: Created [date], Expires [date]
- Stored in: Secret Manager `azure-ad-client-secret`

**Group Assignment:**
- Required: Yes
- Group: precision-medicine-users ([group-id])
```

### Set Calendar Reminder

**For client secret expiration:**
```bash
# Get expiration date from Azure Portal
# Set reminder 30 days before expiration
# Update Secret Manager before expiration
```

### Review Access Logs

```bash
# Check for failed login attempts
gcloud logging read \
  'resource.type="cloud_run_revision"
   AND resource.labels.service_name~"oauth2-proxy"
   AND severity>=WARNING
   AND timestamp>="'$(date -d '7 days ago' -I)'T00:00:00Z"' \
  --limit=100
```

---

## Preventive Measures

### 1. Monitoring

**Create alert for SSO failures:**
```bash
gcloud alpha monitoring policies create \
  --notification-channels=<CHANNEL_ID> \
  --display-name="OAuth2 Proxy High Failure Rate" \
  --condition-filter='resource.type="cloud_run_revision"
    AND resource.labels.service_name~"oauth2-proxy"
    AND severity>=ERROR' \
  --condition-threshold-value=10 \
  --condition-threshold-duration=300s
```

### 2. Documentation

**Update user guide with troubleshooting:**
- Common login issues
- How to clear cookies
- VPN requirements
- Who to contact

### 3. Secret Expiration Tracking

**Create calendar reminders:**
- 90 days before expiration: Plan renewal
- 30 days before: Execute renewal
- 7 days before: Verify new secret working

---

## Related Runbooks

- [Server Down](server-down.md)
- [Epic Connection Failure](epic-connection-failure.md)

---

**Document History:**
- v1.0 (2025-12-30): Initial runbook
- Last Incident: N/A
