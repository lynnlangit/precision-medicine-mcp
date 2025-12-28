# Phase 1A: GCP Infrastructure Setup

**Timeline:** Week 1-2
**Budget Focus:** Use GCP free tier + minimal services
**Goal:** HIPAA-compliant foundation for patient data testing

---

## Overview

Set up Google Cloud Platform infrastructure to support HIPAA-compliant precision medicine analysis with Epic FHIR integration.

---

## Prerequisites

- [ ] GCP account with billing enabled
- [ ] GCP Organization (or use personal account for POC)
- [ ] gcloud CLI installed locally
- [ ] Terraform installed (optional but recommended)

---

## Task 1: GCP Project Setup (30 mins)

### 1.1 Create HIPAA-Compliant Project

```bash
# Set variables
export PROJECT_ID="precision-medicine-poc"
export PROJECT_NAME="Precision Medicine POC"
export BILLING_ACCOUNT_ID="YOUR_BILLING_ACCOUNT_ID"
export REGION="us-central1"  # HIPAA-eligible region

# Create project
gcloud projects create $PROJECT_ID \
  --name="$PROJECT_NAME" \
  --set-as-default

# Link billing
gcloud billing projects link $PROJECT_ID \
  --billing-account=$BILLING_ACCOUNT_ID

# Enable required APIs
gcloud services enable \
  compute.googleapis.com \
  storage-api.googleapis.com \
  cloudkms.googleapis.com \
  logging.googleapis.com \
  monitoring.googleapis.com \
  secretmanager.googleapis.com \
  healthcare.googleapis.com \
  --project=$PROJECT_ID
```

### 1.2 Enable HIPAA Compliance Features

**Sign BAA (Business Associate Agreement) with Google:**
- Go to: https://cloud.google.com/security/compliance/hipaa
- Request BAA through GCP support or sales
- **Note:** Required for production, can defer for POC testing

**Enable Compliance Features:**
```bash
# Enable audit logging (MANDATORY for HIPAA)
gcloud logging sinks create hipaa-audit-logs \
  storage.googleapis.com/precision-medicine-audit-logs \
  --log-filter='protoPayload.serviceName="healthcare.googleapis.com"' \
  --project=$PROJECT_ID

# Enable Security Command Center (if available)
gcloud services enable securitycenter.googleapis.com --project=$PROJECT_ID
```

---

## Task 2: Encrypted Storage Setup (45 mins)

### 2.1 Create Customer-Managed Encryption Keys (CMEK)

```bash
# Create keyring
gcloud kms keyrings create precision-medicine-keyring \
  --location=$REGION \
  --project=$PROJECT_ID

# Create encryption key
gcloud kms keys create data-encryption-key \
  --location=$REGION \
  --keyring=precision-medicine-keyring \
  --purpose=encryption \
  --rotation-period=90d \
  --next-rotation-time=$(date -u -d '+90 days' +%Y-%m-%dT%H:%M:%SZ) \
  --project=$PROJECT_ID

# Grant service account access to key
export SERVICE_ACCOUNT="cloudservices@${PROJECT_ID}.iam.gserviceaccount.com"
gcloud kms keys add-iam-policy-binding data-encryption-key \
  --location=$REGION \
  --keyring=precision-medicine-keyring \
  --member="serviceAccount:${SERVICE_ACCOUNT}" \
  --role=roles/cloudkms.cryptoKeyEncrypterDecrypter \
  --project=$PROJECT_ID
```

### 2.2 Create Encrypted Storage Buckets

```bash
# Get encryption key resource name
export KEY_NAME=$(gcloud kms keys describe data-encryption-key \
  --location=$REGION \
  --keyring=precision-medicine-keyring \
  --project=$PROJECT_ID \
  --format="value(name)")

# Create bucket for patient data (CMEK encrypted)
gsutil mb \
  -p $PROJECT_ID \
  -c STANDARD \
  -l $REGION \
  -b on \
  gs://precision-medicine-patient-data/

# Enable encryption
gsutil kms encryption \
  -k $KEY_NAME \
  gs://precision-medicine-patient-data/

# Create bucket for results
gsutil mb \
  -p $PROJECT_ID \
  -c STANDARD \
  -l $REGION \
  -b on \
  gs://precision-medicine-results/

gsutil kms encryption \
  -k $KEY_NAME \
  gs://precision-medicine-results/

# Set lifecycle rules (delete after 7 years per HIPAA)
cat > lifecycle.json <<EOF
{
  "lifecycle": {
    "rule": [
      {
        "action": {"type": "Delete"},
        "condition": {"age": 2555}
      }
    ]
  }
}
EOF

gsutil lifecycle set lifecycle.json gs://precision-medicine-patient-data/
gsutil lifecycle set lifecycle.json gs://precision-medicine-results/
```

---

## Task 3: Networking & VPC Setup (30 mins)

### 3.1 Create Private VPC

```bash
# Create VPC
gcloud compute networks create precision-medicine-vpc \
  --subnet-mode=custom \
  --bgp-routing-mode=regional \
  --project=$PROJECT_ID

# Create subnet
gcloud compute networks subnets create precision-medicine-subnet \
  --network=precision-medicine-vpc \
  --region=$REGION \
  --range=10.0.0.0/24 \
  --enable-private-ip-google-access \
  --project=$PROJECT_ID

# Create firewall rule - deny all ingress by default
gcloud compute firewall-rules create deny-all-ingress \
  --network=precision-medicine-vpc \
  --action=deny \
  --direction=ingress \
  --rules=all \
  --priority=65534 \
  --project=$PROJECT_ID

# Allow internal communication
gcloud compute firewall-rules create allow-internal \
  --network=precision-medicine-vpc \
  --action=allow \
  --direction=ingress \
  --rules=all \
  --source-ranges=10.0.0.0/24 \
  --priority=1000 \
  --project=$PROJECT_ID
```

### 3.2 Set Up Cloud NAT (for outbound internet access)

```bash
# Create Cloud Router
gcloud compute routers create precision-medicine-router \
  --network=precision-medicine-vpc \
  --region=$REGION \
  --project=$PROJECT_ID

# Create Cloud NAT
gcloud compute routers nats create precision-medicine-nat \
  --router=precision-medicine-router \
  --region=$REGION \
  --nat-all-subnet-ip-ranges \
  --auto-allocate-nat-external-ips \
  --project=$PROJECT_ID
```

---

## Task 4: Secret Management (20 mins)

### 4.1 Set Up Secret Manager

```bash
# Create secrets for Epic FHIR credentials (placeholders for now)
echo -n "YOUR_EPIC_CLIENT_ID" | \
  gcloud secrets create epic-client-id \
  --data-file=- \
  --replication-policy=automatic \
  --project=$PROJECT_ID

echo -n "YOUR_EPIC_CLIENT_SECRET" | \
  gcloud secrets create epic-client-secret \
  --data-file=- \
  --replication-policy=automatic \
  --project=$PROJECT_ID

echo -n "https://fhir.epic.com/interconnect-fhir-oauth/api/FHIR/R4/" | \
  gcloud secrets create epic-fhir-endpoint \
  --data-file=- \
  --replication-policy=automatic \
  --project=$PROJECT_ID

# Create service account for accessing secrets
gcloud iam service-accounts create mcp-server-sa \
  --display-name="MCP Server Service Account" \
  --project=$PROJECT_ID

# Grant secret access
gcloud secrets add-iam-policy-binding epic-client-id \
  --member="serviceAccount:mcp-server-sa@${PROJECT_ID}.iam.gserviceaccount.com" \
  --role="roles/secretmanager.secretAccessor" \
  --project=$PROJECT_ID

gcloud secrets add-iam-policy-binding epic-client-secret \
  --member="serviceAccount:mcp-server-sa@${PROJECT_ID}.iam.gserviceaccount.com" \
  --role="roles/secretmanager.secretAccessor" \
  --project=$PROJECT_ID

gcloud secrets add-iam-policy-binding epic-fhir-endpoint \
  --member="serviceAccount:mcp-server-sa@${PROJECT_ID}.iam.gserviceaccount.com" \
  --role="roles/secretmanager.secretAccessor" \
  --project=$PROJECT_ID
```

---

## Task 5: Audit Logging & Monitoring (30 mins)

### 5.1 Configure Comprehensive Audit Logging

```bash
# Enable Data Access audit logs for all services
cat > audit-config.yaml <<EOF
auditConfigs:
- auditLogConfigs:
  - logType: ADMIN_READ
  - logType: DATA_READ
  - logType: DATA_WRITE
  service: allServices
EOF

gcloud projects set-iam-policy $PROJECT_ID audit-config.yaml

# Create log bucket for long-term retention (10 years for HIPAA)
gcloud logging buckets create hipaa-logs \
  --location=$REGION \
  --retention-days=3650 \
  --project=$PROJECT_ID

# Create log sink for all audit logs
gcloud logging sinks create all-audit-logs \
  logging.googleapis.com/projects/$PROJECT_ID/locations/$REGION/buckets/hipaa-logs \
  --log-filter='logName:"cloudaudit.googleapis.com"' \
  --project=$PROJECT_ID
```

### 5.2 Set Up Monitoring Alerts

```bash
# Create notification channel (email)
gcloud alpha monitoring channels create \
  --display-name="Security Alerts" \
  --type=email \
  --channel-labels=email_address=YOUR_EMAIL@example.com \
  --project=$PROJECT_ID

# Get channel ID
export CHANNEL_ID=$(gcloud alpha monitoring channels list \
  --filter="displayName='Security Alerts'" \
  --format="value(name)" \
  --project=$PROJECT_ID)

# Create alert for unusual data access
cat > alert-policy.yaml <<EOF
displayName: "Unusual Data Access Pattern"
conditions:
  - displayName: "High volume of GCS reads"
    conditionThreshold:
      filter: 'resource.type="gcs_bucket" AND metric.type="storage.googleapis.com/api/request_count"'
      comparison: COMPARISON_GT
      thresholdValue: 1000
      duration: 60s
      aggregations:
        - alignmentPeriod: 60s
          perSeriesAligner: ALIGN_RATE
notificationChannels:
  - ${CHANNEL_ID}
alertStrategy:
  autoClose: 604800s
EOF

gcloud alpha monitoring policies create --policy-from-file=alert-policy.yaml --project=$PROJECT_ID
```

---

## Task 6: Service Account & IAM Setup (20 mins)

### 6.1 Create Service Accounts with Least Privilege

```bash
# Service account for MCP servers
gcloud iam service-accounts create mcp-server-sa \
  --display-name="MCP Server Service Account" \
  --project=$PROJECT_ID

# Grant minimal permissions
gcloud projects add-iam-policy-binding $PROJECT_ID \
  --member="serviceAccount:mcp-server-sa@${PROJECT_ID}.iam.gserviceaccount.com" \
  --role="roles/storage.objectViewer" \
  --condition=None

gcloud projects add-iam-policy-binding $PROJECT_ID \
  --member="serviceAccount:mcp-server-sa@${PROJECT_ID}.iam.gserviceaccount.com" \
  --role="roles/storage.objectCreator" \
  --condition=None

# Create key for local development (TEMPORARY - use Workload Identity in production)
gcloud iam service-accounts keys create ./mcp-server-key.json \
  --iam-account=mcp-server-sa@${PROJECT_ID}.iam.gserviceaccount.com \
  --project=$PROJECT_ID

# IMPORTANT: Store this key securely and add to .gitignore
echo "mcp-server-key.json" >> .gitignore
```

---

## Task 7: De-identification Service Setup (45 mins)

### 7.1 Set Up Cloud Healthcare API for De-identification

```bash
# Create Healthcare dataset
gcloud healthcare datasets create precision-medicine-dataset \
  --location=$REGION \
  --project=$PROJECT_ID

# Create FHIR store for de-identified data
gcloud healthcare fhir-stores create deidentified-fhir-store \
  --dataset=precision-medicine-dataset \
  --location=$REGION \
  --version=R4 \
  --enable-update-create \
  --project=$PROJECT_ID

# Create FHIR store for identified data (temporary, for de-identification)
gcloud healthcare fhir-stores create identified-fhir-store \
  --dataset=precision-medicine-dataset \
  --location=$REGION \
  --version=R4 \
  --enable-update-create \
  --project=$PROJECT_ID
```

### 7.2 Configure De-identification Template

Create a de-identification config file:

```bash
cat > deidentify-config.json <<'EOF'
{
  "deidentifyConfig": {
    "fhir": {
      "fieldMetadataList": [
        {
          "paths": ["Patient.name", "Patient.telecom", "Patient.address"],
          "action": "REMOVE"
        },
        {
          "paths": ["Patient.birthDate"],
          "action": "DATE_SHIFT"
        },
        {
          "paths": ["Patient.identifier"],
          "action": "HASH"
        }
      ]
    }
  }
}
EOF
```

---

## Task 8: Cost Controls & Budgets (15 mins)

### 8.1 Set Up Budget Alerts

```bash
# Create budget (adjust amount as needed)
gcloud billing budgets create \
  --billing-account=$BILLING_ACCOUNT_ID \
  --display-name="Precision Medicine POC Budget" \
  --budget-amount=500 \
  --threshold-rule=percent=50 \
  --threshold-rule=percent=75 \
  --threshold-rule=percent=90 \
  --threshold-rule=percent=100
```

### 8.2 Enable Cost Optimization

```bash
# Create custom cost report
# Go to: https://console.cloud.google.com/billing/reports
# Filter by project: precision-medicine-poc
# Set up weekly email reports

# Enable idle resource recommendations
gcloud recommender recommendations list \
  --project=$PROJECT_ID \
  --location=$REGION \
  --recommender=google.compute.instance.IdleResourceRecommender
```

---

## Task 9: Security Hardening Checklist

### 9.1 Security Configuration Review

- [ ] **VPC Service Controls** (if available in your GCP tier)
  ```bash
  # Create security perimeter (Enterprise tier only)
  # gcloud access-context-manager perimeters create ...
  ```

- [ ] **Enable OS Login** (for VM access if using Compute Engine)
  ```bash
  gcloud compute project-info add-metadata \
    --metadata enable-oslogin=TRUE \
    --project=$PROJECT_ID
  ```

- [ ] **Disable default service account** (security best practice)
  ```bash
  gcloud projects remove-iam-policy-binding $PROJECT_ID \
    --member="serviceAccount:${PROJECT_NUMBER}-compute@developer.gserviceaccount.com" \
    --role="roles/editor"
  ```

- [ ] **Enable Security Scanner** (for web applications if deploying UI)

- [ ] **Review IAM permissions** - ensure least privilege
  ```bash
  gcloud projects get-iam-policy $PROJECT_ID
  ```

---

## Task 10: Documentation & Handoff (30 mins)

### 10.1 Create Infrastructure Documentation

Create `infrastructure/GCP_SETUP.md`:

```markdown
# GCP Infrastructure Documentation

## Project Details
- Project ID: precision-medicine-poc
- Region: us-central1
- Environment: POC/Development

## Resources Created
- VPC: precision-medicine-vpc
- Subnets: precision-medicine-subnet (10.0.0.0/24)
- Storage Buckets:
  - gs://precision-medicine-patient-data (encrypted)
  - gs://precision-medicine-results (encrypted)
- KMS Keys: data-encryption-key
- Service Accounts: mcp-server-sa
- Healthcare Dataset: precision-medicine-dataset

## Access Credentials
- Service Account Key: mcp-server-key.json (LOCAL ONLY - DO NOT COMMIT)
- Epic FHIR Secrets: stored in Secret Manager

## Audit & Compliance
- Audit logs: 10-year retention in hipaa-logs bucket
- Data retention: 7-year automatic deletion
- Encryption: Customer-managed keys (CMEK)
```

### 10.2 Update Environment Configuration

Create `.env.gcp` file:

```bash
cat > .env.gcp <<EOF
# GCP Configuration
export GCP_PROJECT_ID="precision-medicine-poc"
export GCP_REGION="us-central1"
export GOOGLE_APPLICATION_CREDENTIALS="./mcp-server-key.json"

# Storage
export PATIENT_DATA_BUCKET="gs://precision-medicine-patient-data"
export RESULTS_BUCKET="gs://precision-medicine-results"

# Healthcare API
export HEALTHCARE_DATASET="precision-medicine-dataset"
export FHIR_STORE="deidentified-fhir-store"
export HEALTHCARE_LOCATION="us-central1"

# Epic FHIR (secrets in Secret Manager)
export EPIC_CLIENT_ID_SECRET="projects/precision-medicine-poc/secrets/epic-client-id"
export EPIC_CLIENT_SECRET_SECRET="projects/precision-medicine-poc/secrets/epic-client-secret"
export EPIC_FHIR_ENDPOINT_SECRET="projects/precision-medicine-poc/secrets/epic-fhir-endpoint"

# Security
export DRY_RUN="false"  # Set to false for real data processing
EOF

# Add to .gitignore
echo ".env.gcp" >> .gitignore
echo "mcp-server-key.json" >> .gitignore
```

---

## Validation Checklist

Before proceeding to Phase 1B, verify:

- [ ] GCP project created and billing enabled
- [ ] All required APIs enabled
- [ ] Encrypted storage buckets created and accessible
- [ ] VPC and networking configured
- [ ] Secret Manager configured with placeholder secrets
- [ ] Service accounts created with proper IAM roles
- [ ] Audit logging enabled and tested
- [ ] Budget alerts configured
- [ ] Security hardening checklist completed
- [ ] Documentation created
- [ ] .env.gcp file created and tested

---

## Testing the Infrastructure

```bash
# Test storage access
echo "test" | gsutil cp - gs://precision-medicine-patient-data/test.txt
gsutil cat gs://precision-medicine-patient-data/test.txt
gsutil rm gs://precision-medicine-patient-data/test.txt

# Test secret access
gcloud secrets versions access latest --secret=epic-fhir-endpoint --project=$PROJECT_ID

# Test service account
gcloud auth activate-service-account \
  --key-file=./mcp-server-key.json

gsutil ls gs://precision-medicine-patient-data/

# Revert to your user account
gcloud config set account YOUR_EMAIL@example.com
```

---

## Cost Estimate (Monthly)

**Infrastructure Costs (POC tier):**
- Cloud Storage: $0.02/GB = ~$1-5/month (assuming <250GB)
- KMS: $0.06/key/month = $1/month
- Secret Manager: $0.06/secret/month = $1/month
- Cloud Logging: $0.50/GB = $5-10/month
- Networking (NAT): $0.045/hour = ~$33/month
- Healthcare API: Pay per use (~$5-20/month for POC)

**Total Estimated: $50-75/month**

---

## Troubleshooting

### Common Issues

**Issue: Permission denied when accessing buckets**
```bash
# Check service account permissions
gcloud projects get-iam-policy $PROJECT_ID \
  --flatten="bindings[].members" \
  --filter="bindings.members:mcp-server-sa@${PROJECT_ID}.iam.gserviceaccount.com"
```

**Issue: Cannot access secrets**
```bash
# Verify secret exists
gcloud secrets describe epic-client-id --project=$PROJECT_ID

# Check IAM bindings
gcloud secrets get-iam-policy epic-client-id --project=$PROJECT_ID
```

**Issue: VPC connection issues**
```bash
# Check firewall rules
gcloud compute firewall-rules list --project=$PROJECT_ID

# Test connectivity
gcloud compute ssh test-vm --zone=us-central1-a --project=$PROJECT_ID
```

---

## Next Steps

Once Phase 1A is complete:
1. Proceed to [Phase 1B: Epic FHIR Integration](PHASE1B_EPIC_INTEGRATION.md)
2. Test de-identification pipeline with sample FHIR data
3. Set up CI/CD for MCP server deployment

---

**Questions or Issues?** Document in: `infrastructure/ISSUES.md`
