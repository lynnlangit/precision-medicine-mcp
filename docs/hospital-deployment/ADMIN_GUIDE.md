# Administrator Guide - Precision Medicine MCP System
## For Hospital IT and System Administrators

**Version:** 1.0
**Last Updated:** December 2025
**Target Audience:** Hospital IT staff, GCP administrators

---

## Table of Contents

- [Administrator Responsibilities](#administrator-responsibilities)
- [User Management](#user-management)
- [Server Management](#server-management)
- [Monitoring & Alerts](#monitoring--alerts)
- [Cost Management](#cost-management)
- [Security Administration](#security-administration)
- [Configuration Management](#configuration-management)
- [Troubleshooting Guide](#troubleshooting-guide)

---

## Administrator Responsibilities

### Primary Responsibilities

1. **User Access Management**
   - Add/remove users from Azure AD groups
   - Verify user access and permissions
   - Troubleshoot login issues

2. **System Monitoring**
   - Check Cloud Monitoring dashboard daily
   - Respond to alerts
   - Review error logs weekly

3. **Cost Management**
   - Track spending vs. budget ($1,000/month)
   - Optimize resource usage
   - Report overages to PI

4. **Security & Compliance**
   - Ensure HIPAA compliance
   - Review audit logs monthly
   - Rotate credentials annually
   - Respond to security incidents

5. **Maintenance**
   - Deploy updates
   - Backup configurations
   - Test disaster recovery procedures

### Access Requirements

**You should have:**
- GCP Project Owner or Editor role
- Azure AD User Administrator role (or delegated to Hospital IT)
- Access to hospital VPN admin console
- Access to Secret Manager

---

## User Management

### Adding a New User

**Prerequisites:**
- User has hospital email address (@hospital.org)
- User has completed IRB training
- PI approval obtained

**Step 1: Add to Azure AD Group**

```bash
# Via Azure Portal (recommended)
1. Login to Azure Portal (portal.azure.com)
2. Navigate to Azure Active Directory
3. Go to Groups
4. Find "precision-medicine-users"
5. Click "Members" → "Add members"
6. Search for user by name or email
7. Select user and click "Select"
8. Click "Add"

# Via PowerShell (alternative)
Connect-AzureAD
Add-AzureADGroupMember `
    -ObjectId <group-object-id> `
    -RefObjectId <user-object-id>
```

**Step 2: Verify Access**

```bash
# Test login
1. User should log in to Streamlit UI
2. Check audit logs for successful login event:

gcloud logging read \
  'jsonPayload.event="user_login" AND jsonPayload.user_email_hash~"<user-hash>"' \
  --limit=10 \
  --project=<PROJECT_ID>
```

**Step 3: Notify User**

Email template:
```
Subject: Access Granted - Precision Medicine MCP System

Hi [Name],

Your access to the Precision Medicine MCP system has been granted.

Login URLs:
- Streamlit UI: https://oauth2-proxy-streamlit-{hash}.run.app
- Jupyter Notebook: https://oauth2-proxy-jupyter-{hash}.run.app

Credentials: Use your hospital email and Azure AD password

Resources:
- User Guide: [link to USER_GUIDE.md]
- Training Session: [schedule if applicable]

Support: help@hospital.org

Best regards,
Hospital IT
```

### Removing a User

**When to remove:**
- User leaves research team
- User violates policies
- Access no longer needed

**Step 1: Remove from Azure AD Group**

```bash
# Via Azure Portal
1. Azure Portal → Azure Active Directory → Groups
2. "precision-medicine-users" → Members
3. Find user → Select → "Remove member"

# Via PowerShell
Remove-AzureADGroupMember `
    -ObjectId <group-object-id> `
    -MemberId <user-object-id>
```

**Step 2: Verify Removal**

```bash
# User should not be able to log in
# Check they get "Access Denied" error
```

**Step 3: Document Removal**

```bash
# Log the removal for compliance
echo "[$(date)] Removed user: <email> - Reason: <reason>" >> user_access.log
```

### Listing Current Users

```bash
# Via Azure Portal
Azure Portal → Azure Active Directory → Groups
→ "precision-medicine-users" → Members

# Via PowerShell
Get-AzureADGroupMember -ObjectId <group-object-id> | Select DisplayName, UserPrincipalName

# Check active users from audit logs (last 30 days)
gcloud logging read \
  'jsonPayload.event="user_login"' \
  --limit=100 \
  --format='table(jsonPayload.user_email_hash, timestamp)' \
  --project=<PROJECT_ID>
```

### User Access Audit

**Monthly procedure:**

```bash
# 1. Export Azure AD group members
Get-AzureADGroupMember -ObjectId <group-object-id> |
  Export-Csv -Path "ad_users_$(date +%Y%m).csv"

# 2. Export active users from logs
gcloud logging read \
  'jsonPayload.event="user_login" AND timestamp>"2025-12-01T00:00:00Z"' \
  --format=json \
  --project=<PROJECT_ID> > active_users_$(date +%Y%m).json

# 3. Compare lists
# Users in AD but not in logs = inactive users (consider removing)
# Users in logs but not in AD = error (investigate)

# 4. Document findings
echo "Access audit completed: $(date)" >> audit_reports/access_audit_$(date +%Y%m).log
```

---

## Server Management

### Viewing Server Status

```bash
# List all Cloud Run services
gcloud run services list \
  --region=us-central1 \
  --project=<PROJECT_ID> \
  --format='table(name,status.url,status.latestCreatedRevision)'

# Get detailed status for specific server
gcloud run services describe mcp-fgbio \
  --region=us-central1 \
  --project=<PROJECT_ID> \
  --format=yaml

# Check server health
for server in mcp-fgbio mcp-multiomics mcp-spatialtools mcp-epic; do
  echo "=== $server ==="
  gcloud run services describe $server \
    --region=us-central1 \
    --project=<PROJECT_ID> \
    --format='value(status.conditions[0].status)'
done
```

### Viewing Server Logs

```bash
# Real-time logs for a server
gcloud run services logs read mcp-fgbio \
  --region=us-central1 \
  --project=<PROJECT_ID> \
  --follow

# Recent errors only
gcloud logging read \
  'resource.type="cloud_run_revision"
   AND resource.labels.service_name="mcp-fgbio"
   AND severity>=ERROR' \
  --limit=50 \
  --project=<PROJECT_ID>

# Query-specific logs
gcloud logging read \
  'resource.type="cloud_run_revision"
   AND resource.labels.service_name="mcp-fgbio"
   AND timestamp>="2025-12-30T10:00:00Z"' \
  --limit=100 \
  --project=<PROJECT_ID>
```

### Deploying Server Updates

**Before deploying:**
1. Review changelog
2. Test in development environment (if available)
3. Schedule maintenance window
4. Notify users of potential downtime

**Deployment procedure:**

```bash
# 1. Build new container version
cd servers/mcp-fgbio
docker build -t gcr.io/<PROJECT_ID>/mcp-fgbio:v1.1 .

# 2. Push to Container Registry
docker push gcr.io/<PROJECT_ID>/mcp-fgbio:v1.1

# 3. Deploy to Cloud Run (zero downtime)
gcloud run deploy mcp-fgbio \
  --image=gcr.io/<PROJECT_ID>/mcp-fgbio:v1.1 \
  --region=us-central1 \
  --project=<PROJECT_ID>

# 4. Monitor for errors
gcloud run services logs read mcp-fgbio \
  --region=us-central1 \
  --limit=50

# 5. Test functionality
# Send test query via Streamlit or Jupyter

# 6. If successful, update version tag
gcloud run services update mcp-fgbio \
  --update-labels=version=v1.1 \
  --region=us-central1
```

### Rolling Back a Deployment

```bash
# 1. List previous revisions
gcloud run revisions list \
  --service=mcp-fgbio \
  --region=us-central1 \
  --format='table(name,creationTimestamp,status)'

# 2. Rollback to previous revision
gcloud run services update-traffic mcp-fgbio \
  --to-revisions=mcp-fgbio-00042-xyz=100 \
  --region=us-central1

# 3. Verify rollback
gcloud run services describe mcp-fgbio \
  --region=us-central1 \
  --format='value(status.traffic[0].revisionName)'

# 4. Delete failed revision (optional)
gcloud run revisions delete mcp-fgbio-00043-bad \
  --region=us-central1
```

### Scaling Configuration

```bash
# Update min/max instances
gcloud run services update mcp-spatialtools \
  --min-instances=1 \
  --max-instances=10 \
  --region=us-central1

# Update memory/CPU
gcloud run services update mcp-multiomics \
  --memory=8Gi \
  --cpu=4 \
  --region=us-central1

# Update timeout
gcloud run services update mcp-fgbio \
  --timeout=600 \
  --region=us-central1
```

---

## Monitoring & Alerts

### Cloud Monitoring Dashboard

**Access:** https://console.cloud.google.com/monitoring/dashboards

**Key metrics to watch:**

1. **Request Rate**
   - Normal: 10-50 requests/hour during business hours
   - Alert if: Sudden spike (>100/hour) or drop to zero

2. **Error Rate**
   - Normal: <1% errors
   - Alert if: >5% errors

3. **Latency**
   - Normal: P50 < 5s, P95 < 30s
   - Alert if: P95 > 60s consistently

4. **Cost per Day**
   - Normal: $30-40/day
   - Alert if: >$50/day (approaching budget limit)

### Viewing Alerts

```bash
# List alert policies
gcloud alpha monitoring policies list \
  --project=<PROJECT_ID>

# View alert incidents
gcloud alpha monitoring alerting incidents list \
  --policy=<POLICY_NAME> \
  --project=<PROJECT_ID>
```

### Creating Custom Alerts

**Example: Alert on high de-identification failures**

```bash
gcloud alpha monitoring policies create \
  --notification-channels=<CHANNEL_ID> \
  --display-name="High De-identification Failures" \
  --condition-filter='metric.type="logging.googleapis.com/user/deidentification_failures"
                     resource.type="global"' \
  --condition-threshold-value=5 \
  --condition-threshold-duration=300s \
  --project=<PROJECT_ID>
```

### Notification Channels

```bash
# List notification channels
gcloud alpha monitoring channels list \
  --project=<PROJECT_ID>

# Create email notification channel
gcloud alpha monitoring channels create \
  --display-name="Hospital IT Team" \
  --type=email \
  --channel-labels=email_address=it-team@hospital.org \
  --project=<PROJECT_ID>

# Create SMS notification channel (for P0 incidents)
gcloud alpha monitoring channels create \
  --display-name="On-Call Pager" \
  --type=sms \
  --channel-labels=number=+15551234567 \
  --project=<PROJECT_ID>
```

---

## Cost Management

### Viewing Current Spend

```bash
# Via gcloud (requires billing permissions)
gcloud billing accounts list

# Get project billing info
gcloud beta billing projects describe <PROJECT_ID>

# Export billing data
gcloud beta billing export describe

# Via Cloud Console (easier)
https://console.cloud.google.com/billing/<BILLING_ACCOUNT_ID>/reports
```

### Cost Breakdown

**Monthly Budget:** $1,000

**Typical Breakdown:**
- Cloud Run (servers): $400-450
- VPC/Networking: $50
- Storage (GCS): $100-150
- Anthropic API: $400-500
- Monitoring/Logging: $50

### Optimizing Costs

**1. Reduce Min Instances for Mock Servers:**

```bash
# Mock servers don't need to be always-on
for server in mcp-tcga mcp-seqera mcp-huggingface mcp-deepcell mcp-mockepic; do
  gcloud run services update $server \
    --min-instances=0 \
    --region=us-central1
done

# Save ~$100/month
```

**2. Use Haiku Model for Simple Queries:**

Encourage users to use Haiku instead of Sonnet for simple queries:
- Sonnet: $3/M input, $15/M output
- Haiku: $0.25/M input, $1.25/M output
- 10x cheaper!

**3. Implement Query Caching:**

For identical queries, cache results to avoid duplicate API calls:
```python
# Future enhancement - not yet implemented
# Would save ~20% on API costs
```

**4. Set Up Budget Alerts:**

```bash
# Alert at 50%, 75%, 90%, 100% of budget
# Already configured in setup-project.sh
```

**5. Monitor Token Usage:**

```bash
# Find users with highest token usage
gcloud logging read \
  'jsonPayload.event="mcp_response"' \
  --format=json \
  --limit=1000 | \
  jq -r '.[] | [.jsonPayload.user_email_hash, .jsonPayload.total_tokens] | @csv' | \
  sort | uniq -c | sort -rn | head -10

# Educate high-usage users on cost-saving practices
```

### Monthly Cost Report

**Generate report for PI:**

```bash
# Export billing data
gcloud billing projects describe <PROJECT_ID> --format=json

# Create summary report
cat > monthly_cost_report_$(date +%Y%m).md <<EOF
# Precision Medicine MCP - Monthly Cost Report
## $(date +"%B %Y")

**Total Spend:** \$XXX.XX / \$1,000.00 budget

**Breakdown:**
- Cloud Run: \$XXX.XX (XX%)
- Anthropic API: \$XXX.XX (XX%)
- Networking: \$XXX.XX (XX%)
- Storage: \$XXX.XX (XX%)
- Other: \$XXX.XX (XX%)

**Usage Statistics:**
- Total Queries: XXX
- Active Users: X/5
- Average Cost per Query: \$X.XX

**Optimization Opportunities:**
- [List any cost-saving recommendations]

**Next Month Forecast:** \$XXX.XX
EOF
```

---

## Security Administration

### Secret Management

**View secrets:**

```bash
# List all secrets
gcloud secrets list --project=<PROJECT_ID>

# View secret metadata (not value)
gcloud secrets describe anthropic-api-key --project=<PROJECT_ID>

# View secret access permissions
gcloud secrets get-iam-policy anthropic-api-key --project=<PROJECT_ID>
```

**Rotate secret:**

```bash
# 1. Create new secret version
echo -n "new-secret-value" | gcloud secrets versions add anthropic-api-key \
  --data-file=- \
  --project=<PROJECT_ID>

# 2. Update Cloud Run service to use new version
gcloud run services update streamlit-mcp-chat \
  --update-secrets=ANTHROPIC_API_KEY=anthropic-api-key:latest \
  --region=us-central1

# 3. Test that service still works

# 4. Destroy old secret version
gcloud secrets versions destroy <VERSION_NUMBER> \
  --secret=anthropic-api-key \
  --project=<PROJECT_ID>
```

**Grant service account access:**

```bash
# Allow service account to read specific secret
gcloud secrets add-iam-policy-binding epic-client-secret \
  --member=serviceAccount:mcp-epic-sa@<PROJECT_ID>.iam.gserviceaccount.com \
  --role=roles/secretmanager.secretAccessor \
  --project=<PROJECT_ID>
```

### Audit Log Review

**Monthly audit checklist:**

```bash
# 1. Review all user logins
gcloud logging read \
  'jsonPayload.event="user_login" AND timestamp>"$(date -d '30 days ago' -I)T00:00:00Z"' \
  --project=<PROJECT_ID> \
  --format='table(timestamp, jsonPayload.user_email_hash, jsonPayload.display_name)'

# 2. Check for failed de-identification attempts
gcloud logging read \
  'jsonPayload.event="deidentification" AND jsonPayload.success=false' \
  --project=<PROJECT_ID> \
  --limit=50

# 3. Review Epic FHIR access
gcloud logging read \
  'jsonPayload.event="epic_fhir_call"' \
  --project=<PROJECT_ID> \
  --limit=100

# 4. Check for errors
gcloud logging read \
  'severity>=ERROR AND timestamp>"$(date -d '30 days ago' -I)T00:00:00Z"' \
  --project=<PROJECT_ID> \
  --limit=100

# 5. Export for compliance
gcloud logging read \
  'timestamp>"$(date -d '30 days ago' -I)T00:00:00Z"' \
  --format=json \
  --project=<PROJECT_ID> > audit_logs_$(date +%Y%m).json
```

### Security Scanning

**Scan containers for vulnerabilities:**

```bash
# Enable Container Scanning API
gcloud services enable containerscanning.googleapis.com

# Scan container image
gcloud container images scan gcr.io/<PROJECT_ID>/mcp-fgbio:latest

# View scan results
gcloud container images describe gcr.io/<PROJECT_ID>/mcp-fgbio:latest \
  --show-package-vulnerability
```

### VPC Security

**Review firewall rules:**

```bash
# List firewall rules
gcloud compute firewall-rules list --project=<PROJECT_ID>

# Check VPC connector configuration
gcloud compute networks vpc-access connectors describe mcp-connector \
  --region=us-central1 \
  --project=<PROJECT_ID>
```

---

## Configuration Management

### Backing Up Configuration

**Weekly backup procedure:**

```bash
#!/bin/bash
# backup_configuration.sh

BACKUP_DIR="backups/$(date +%Y%m%d)"
mkdir -p $BACKUP_DIR

# 1. Export Cloud Run service configs
for service in mcp-fgbio mcp-multiomics mcp-spatialtools mcp-epic \
               streamlit-mcp-chat jupyter-mcp-notebook \
               oauth2-proxy-streamlit oauth2-proxy-jupyter; do
  gcloud run services describe $service \
    --region=us-central1 \
    --format=yaml > $BACKUP_DIR/$service.yaml
done

# 2. Export IAM policies
gcloud projects get-iam-policy <PROJECT_ID> > $BACKUP_DIR/iam_policy.yaml

# 3. Export secret list (not values)
gcloud secrets list --format=yaml > $BACKUP_DIR/secrets_list.yaml

# 4. Export VPC configuration
gcloud compute networks describe <VPC_NETWORK> \
  --format=yaml > $BACKUP_DIR/vpc_network.yaml

# 5. Export firewall rules
gcloud compute firewall-rules list --format=yaml > $BACKUP_DIR/firewall_rules.yaml

# 6. Compress backup
tar -czf $BACKUP_DIR.tar.gz $BACKUP_DIR
rm -rf $BACKUP_DIR

echo "Backup completed: $BACKUP_DIR.tar.gz"
```

### Updating Environment Variables

```bash
# Update single environment variable
gcloud run services update mcp-spatialtools \
  --update-env-vars=SPATIAL_DATA_DIR=gs://new-bucket/spatial \
  --region=us-central1

# Update multiple environment variables
gcloud run services update mcp-multiomics \
  --update-env-vars=MULTIOMICS_DRY_RUN=false,MULTIOMICS_LOG_LEVEL=DEBUG \
  --region=us-central1

# Remove environment variable
gcloud run services update mcp-fgbio \
  --remove-env-vars=OLD_VAR \
  --region=us-central1
```

### Managing Service Accounts

```bash
# Create new service account
gcloud iam service-accounts create new-server-sa \
  --display-name="New Server Service Account" \
  --project=<PROJECT_ID>

# Grant permissions
gcloud projects add-iam-policy-binding <PROJECT_ID> \
  --member=serviceAccount:new-server-sa@<PROJECT_ID>.iam.gserviceaccount.com \
  --role=roles/storage.objectViewer

# List service accounts
gcloud iam service-accounts list --project=<PROJECT_ID>

# View service account permissions
gcloud projects get-iam-policy <PROJECT_ID> \
  --flatten="bindings[].members" \
  --filter="bindings.members:serviceAccount:mcp-fgbio-sa@<PROJECT_ID>.iam.gserviceaccount.com"
```

---

## Troubleshooting Guide

### Common Issues

**Issue: OAuth2 Proxy returns "Access Denied"**

**Diagnosis:**
```bash
# Check OAuth2 Proxy logs
gcloud run services logs read oauth2-proxy-streamlit \
  --region=us-central1 \
  --limit=50

# Check user is in Azure AD group
# (Use Azure Portal or PowerShell)
```

**Resolution:**
1. Verify user is in `precision-medicine-users` AD group
2. Check redirect URI is configured in Azure AD app
3. Verify Azure AD secrets in Secret Manager are correct

---

**Issue: Server returns 5xx errors**

**Diagnosis:**
```bash
# View server logs
gcloud run services logs read mcp-fgbio \
  --region=us-central1 \
  --limit=100 \
  | grep -i error

# Check server status
gcloud run services describe mcp-fgbio \
  --region=us-central1 \
  --format='value(status.conditions[0].message)'
```

**Resolution:**
1. Check if secret is accessible (Epic credentials, API keys)
2. Verify service account has required permissions
3. Check VPC connectivity
4. Review recent deployments (may need rollback)

---

**Issue: Epic FHIR connection fails**

**Diagnosis:**
```bash
# Check Epic FHIR logs
gcloud logging read \
  'jsonPayload.event="epic_fhir_call" AND jsonPayload.status="error"' \
  --limit=20
```

**Resolution:**
1. Test Epic FHIR endpoint manually: `curl https://epic-endpoint/api/FHIR/R4/metadata`
2. Verify Epic credentials in Secret Manager
3. Check Epic OAuth token is valid
4. Contact hospital IT if Epic is down
5. Temporary: Switch to mcp-mockepic server

---

**Issue: Budget overrun**

**Diagnosis:**
```bash
# Check spending
# View in Cloud Console: Billing → Reports

# Find high token usage users
gcloud logging read 'jsonPayload.event="mcp_response"' \
  --format=json --limit=1000 | \
  jq -r '.[] | [.jsonPayload.user_email_hash, .jsonPayload.total_tokens] | @csv'
```

**Resolution:**
1. Identify users with highest usage
2. Educate on cost-saving practices (Haiku model, concise prompts)
3. Consider setting per-user quotas
4. Request budget increase from PI if needed

---

**Issue: Slow query responses**

**Diagnosis:**
```bash
# Check server metrics
gcloud monitoring time-series list \
  --filter='metric.type="run.googleapis.com/request_latencies"' \
  --project=<PROJECT_ID>
```

**Resolution:**
1. Increase Cloud Run CPU/memory allocation
2. Set min-instances=1 to avoid cold starts
3. Optimize MCP server code if needed
4. Consider using Haiku model for faster responses

---

For additional troubleshooting, see:
- [Runbooks](RUNBOOKS/) - Detailed procedures for specific issues
- [Operations Manual](OPERATIONS_MANUAL.md) - System architecture and incident response
- Development Team Support: mcp-support@devteam.com

---

**Document History:**
- v1.0 (2025-12-30): Initial admin guide for production deployment
- Next Review: 2026-01-30 (monthly)

**Related Documents:**
- [Operations Manual](OPERATIONS_MANUAL.md)
- [User Guide](USER_GUIDE.md)
- [HIPAA Compliance](HIPAA_COMPLIANCE.md)
- [Runbooks](RUNBOOKS/)
