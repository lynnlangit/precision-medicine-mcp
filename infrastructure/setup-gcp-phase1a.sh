#!/bin/bash

# Phase 1A: GCP Infrastructure Setup Script
# Run this script step-by-step, reviewing each section before executing

set -e  # Exit on error

echo "================================================================================"
echo "Phase 1A: GCP Infrastructure Setup"
echo "================================================================================"
echo ""
echo "⚠️  This script will create GCP resources that incur costs (~$50-75/month)"
echo "⚠️  Review each section before running"
echo ""

# ============================================================================
# CONFIGURATION - UPDATE THESE VALUES
# ============================================================================

read -p "Enter your GCP Project ID (e.g., precision-medicine-poc): " PROJECT_ID
read -p "Enter your GCP Billing Account ID: " BILLING_ACCOUNT_ID
read -p "Enter your preferred region (default: us-central1): " REGION
REGION=${REGION:-us-central1}

export PROJECT_ID
export BILLING_ACCOUNT_ID
export REGION
export PROJECT_NAME="Precision Medicine POC"

echo ""
echo "Configuration:"
echo "  Project ID: $PROJECT_ID"
echo "  Billing Account: $BILLING_ACCOUNT_ID"
echo "  Region: $REGION"
echo ""
read -p "Is this correct? (y/n): " CONFIRM
if [ "$CONFIRM" != "y" ]; then
    echo "Setup cancelled. Please run again with correct values."
    exit 1
fi

# ============================================================================
# TASK 1: Create Project and Enable APIs
# ============================================================================

echo ""
echo "============================================================================"
echo "TASK 1: Creating GCP Project and Enabling APIs"
echo "============================================================================"
read -p "Press Enter to continue (or Ctrl+C to cancel)..."

# Create project
echo "Creating project: $PROJECT_ID..."
gcloud projects create $PROJECT_ID \
  --name="$PROJECT_NAME" \
  --set-as-default || echo "Project may already exist, continuing..."

# Link billing
echo "Linking billing account..."
gcloud billing projects link $PROJECT_ID \
  --billing-account=$BILLING_ACCOUNT_ID

# Enable required APIs
echo "Enabling required APIs (this may take 2-3 minutes)..."
gcloud services enable \
  compute.googleapis.com \
  storage-api.googleapis.com \
  cloudkms.googleapis.com \
  logging.googleapis.com \
  monitoring.googleapis.com \
  secretmanager.googleapis.com \
  healthcare.googleapis.com \
  --project=$PROJECT_ID

echo "✅ Task 1 Complete: Project created and APIs enabled"

# ============================================================================
# TASK 2: Create Encrypted Storage
# ============================================================================

echo ""
echo "============================================================================"
echo "TASK 2: Setting Up Encrypted Storage"
echo "============================================================================"
read -p "Press Enter to continue..."

# Create keyring
echo "Creating KMS keyring..."
gcloud kms keyrings create precision-medicine-keyring \
  --location=$REGION \
  --project=$PROJECT_ID || echo "Keyring may already exist, continuing..."

# Create encryption key
echo "Creating encryption key with 90-day rotation..."
gcloud kms keys create data-encryption-key \
  --location=$REGION \
  --keyring=precision-medicine-keyring \
  --purpose=encryption \
  --rotation-period=90d \
  --next-rotation-time=$(date -u -v+90d +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || date -u -d '+90 days' +%Y-%m-%dT%H:%M:%SZ) \
  --project=$PROJECT_ID || echo "Key may already exist, continuing..."

# Get key name
KEY_NAME=$(gcloud kms keys describe data-encryption-key \
  --location=$REGION \
  --keyring=precision-medicine-keyring \
  --project=$PROJECT_ID \
  --format="value(name)")

echo "Encryption key created: $KEY_NAME"

# Create storage buckets
echo "Creating encrypted storage buckets..."

# Patient data bucket
gsutil mb \
  -p $PROJECT_ID \
  -c STANDARD \
  -l $REGION \
  -b on \
  gs://${PROJECT_ID}-patient-data/ 2>/dev/null || echo "Bucket may already exist, continuing..."

gsutil kms encryption \
  -k $KEY_NAME \
  gs://${PROJECT_ID}-patient-data/

# Results bucket
gsutil mb \
  -p $PROJECT_ID \
  -c STANDARD \
  -l $REGION \
  -b on \
  gs://${PROJECT_ID}-results/ 2>/dev/null || echo "Bucket may already exist, continuing..."

gsutil kms encryption \
  -k $KEY_NAME \
  gs://${PROJECT_ID}-results/

# Set lifecycle rules (7-year retention)
cat > /tmp/lifecycle.json <<EOF
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

gsutil lifecycle set /tmp/lifecycle.json gs://${PROJECT_ID}-patient-data/
gsutil lifecycle set /tmp/lifecycle.json gs://${PROJECT_ID}-results/

echo "✅ Task 2 Complete: Encrypted storage created"

# ============================================================================
# TASK 3: VPC Networking
# ============================================================================

echo ""
echo "============================================================================"
echo "TASK 3: Configuring VPC Networking"
echo "============================================================================"
read -p "Press Enter to continue..."

# Create VPC
echo "Creating private VPC..."
gcloud compute networks create precision-medicine-vpc \
  --subnet-mode=custom \
  --bgp-routing-mode=regional \
  --project=$PROJECT_ID || echo "VPC may already exist, continuing..."

# Create subnet
echo "Creating subnet..."
gcloud compute networks subnets create precision-medicine-subnet \
  --network=precision-medicine-vpc \
  --region=$REGION \
  --range=10.0.0.0/24 \
  --enable-private-ip-google-access \
  --project=$PROJECT_ID || echo "Subnet may already exist, continuing..."

# Create firewall rules
echo "Creating firewall rules..."
gcloud compute firewall-rules create deny-all-ingress \
  --network=precision-medicine-vpc \
  --action=deny \
  --direction=ingress \
  --rules=all \
  --priority=65534 \
  --project=$PROJECT_ID 2>/dev/null || echo "Firewall rule may already exist, continuing..."

gcloud compute firewall-rules create allow-internal \
  --network=precision-medicine-vpc \
  --action=allow \
  --direction=ingress \
  --rules=all \
  --source-ranges=10.0.0.0/24 \
  --priority=1000 \
  --project=$PROJECT_ID 2>/dev/null || echo "Firewall rule may already exist, continuing..."

# Create Cloud Router and NAT
echo "Creating Cloud Router and NAT..."
gcloud compute routers create precision-medicine-router \
  --network=precision-medicine-vpc \
  --region=$REGION \
  --project=$PROJECT_ID || echo "Router may already exist, continuing..."

gcloud compute routers nats create precision-medicine-nat \
  --router=precision-medicine-router \
  --region=$REGION \
  --nat-all-subnet-ip-ranges \
  --auto-allocate-nat-external-ips \
  --project=$PROJECT_ID 2>/dev/null || echo "NAT may already exist, continuing..."

echo "✅ Task 3 Complete: VPC networking configured"

# ============================================================================
# TASK 4: Secret Manager
# ============================================================================

echo ""
echo "============================================================================"
echo "TASK 4: Setting Up Secret Manager"
echo "============================================================================"
read -p "Press Enter to continue..."

# Create secrets with placeholder values
echo "Creating secrets for Epic FHIR credentials..."

echo -n "PLACEHOLDER_CLIENT_ID" | \
  gcloud secrets create epic-client-id \
  --data-file=- \
  --replication-policy=automatic \
  --project=$PROJECT_ID 2>/dev/null || \
  echo "Secret epic-client-id already exists, skipping..."

echo -n "PLACEHOLDER_CLIENT_SECRET" | \
  gcloud secrets create epic-client-secret \
  --data-file=- \
  --replication-policy=automatic \
  --project=$PROJECT_ID 2>/dev/null || \
  echo "Secret epic-client-secret already exists, skipping..."

echo -n "https://fhir.epic.com/interconnect-fhir-oauth/api/FHIR/R4/" | \
  gcloud secrets create epic-fhir-endpoint \
  --data-file=- \
  --replication-policy=automatic \
  --project=$PROJECT_ID 2>/dev/null || \
  echo "Secret epic-fhir-endpoint already exists, skipping..."

echo "✅ Task 4 Complete: Secret Manager configured"

# ============================================================================
# TASK 5: Audit Logging
# ============================================================================

echo ""
echo "============================================================================"
echo "TASK 5: Configuring Audit Logging"
echo "============================================================================"
read -p "Press Enter to continue..."

# Create log bucket
echo "Creating audit log bucket with 10-year retention..."
gcloud logging buckets create hipaa-logs \
  --location=$REGION \
  --retention-days=3650 \
  --project=$PROJECT_ID 2>/dev/null || echo "Log bucket may already exist, continuing..."

# Create log sink
echo "Creating log sink for audit logs..."
gcloud logging sinks create all-audit-logs \
  logging.googleapis.com/projects/$PROJECT_ID/locations/$REGION/buckets/hipaa-logs \
  --log-filter='logName:"cloudaudit.googleapis.com"' \
  --project=$PROJECT_ID 2>/dev/null || echo "Log sink may already exist, continuing..."

echo "✅ Task 5 Complete: Audit logging configured"

# ============================================================================
# TASK 6: Service Account
# ============================================================================

echo ""
echo "============================================================================"
echo "TASK 6: Creating Service Account"
echo "============================================================================"
read -p "Press Enter to continue..."

# Create service account
echo "Creating MCP server service account..."
gcloud iam service-accounts create mcp-server-sa \
  --display-name="MCP Server Service Account" \
  --project=$PROJECT_ID 2>/dev/null || echo "Service account may already exist, continuing..."

# Grant permissions
echo "Granting permissions..."
gcloud projects add-iam-policy-binding $PROJECT_ID \
  --member="serviceAccount:mcp-server-sa@${PROJECT_ID}.iam.gserviceaccount.com" \
  --role="roles/storage.objectViewer" \
  --condition=None

gcloud projects add-iam-policy-binding $PROJECT_ID \
  --member="serviceAccount:mcp-server-sa@${PROJECT_ID}.iam.gserviceaccount.com" \
  --role="roles/storage.objectCreator" \
  --condition=None

# Grant secret access
for SECRET in epic-client-id epic-client-secret epic-fhir-endpoint; do
  gcloud secrets add-iam-policy-binding $SECRET \
    --member="serviceAccount:mcp-server-sa@${PROJECT_ID}.iam.gserviceaccount.com" \
    --role="roles/secretmanager.secretAccessor" \
    --project=$PROJECT_ID 2>/dev/null || true
done

# Create service account key
echo "Creating service account key..."
KEY_FILE="./mcp-server-key.json"
gcloud iam service-accounts keys create $KEY_FILE \
  --iam-account=mcp-server-sa@${PROJECT_ID}.iam.gserviceaccount.com \
  --project=$PROJECT_ID

echo "⚠️  Service account key created: $KEY_FILE"
echo "⚠️  Keep this file secure and DO NOT commit to git!"

echo "✅ Task 6 Complete: Service account created"

# ============================================================================
# TASK 7: Healthcare API Dataset
# ============================================================================

echo ""
echo "============================================================================"
echo "TASK 7: Setting Up Healthcare API for De-identification"
echo "============================================================================"
read -p "Press Enter to continue..."

# Create Healthcare dataset
echo "Creating Healthcare API dataset..."
gcloud healthcare datasets create precision-medicine-dataset \
  --location=$REGION \
  --project=$PROJECT_ID 2>/dev/null || echo "Dataset may already exist, continuing..."

# Create FHIR stores
echo "Creating FHIR stores..."
gcloud healthcare fhir-stores create deidentified-fhir-store \
  --dataset=precision-medicine-dataset \
  --location=$REGION \
  --version=R4 \
  --enable-update-create \
  --project=$PROJECT_ID 2>/dev/null || echo "FHIR store may already exist, continuing..."

gcloud healthcare fhir-stores create identified-fhir-store \
  --dataset=precision-medicine-dataset \
  --location=$REGION \
  --version=R4 \
  --enable-update-create \
  --project=$PROJECT_ID 2>/dev/null || echo "FHIR store may already exist, continuing..."

echo "✅ Task 7 Complete: Healthcare API configured"

# ============================================================================
# TASK 8: Cost Controls
# ============================================================================

echo ""
echo "============================================================================"
echo "TASK 8: Setting Up Cost Controls"
echo "============================================================================"
read -p "Press Enter to continue..."

echo "Creating budget with $100 threshold..."
gcloud billing budgets create \
  --billing-account=$BILLING_ACCOUNT_ID \
  --display-name="Precision Medicine POC Budget" \
  --budget-amount=100 \
  --threshold-rule=percent=50 \
  --threshold-rule=percent=75 \
  --threshold-rule=percent=90 \
  --threshold-rule=percent=100 2>/dev/null || echo "Budget may already exist, continuing..."

echo "✅ Task 8 Complete: Cost controls configured"

# ============================================================================
# FINAL SUMMARY
# ============================================================================

echo ""
echo "================================================================================"
echo "✅ PHASE 1A COMPLETE!"
echo "================================================================================"
echo ""
echo "Infrastructure Created:"
echo "  ✅ GCP Project: $PROJECT_ID"
echo "  ✅ Region: $REGION"
echo "  ✅ Encrypted Storage:"
echo "     - gs://${PROJECT_ID}-patient-data/"
echo "     - gs://${PROJECT_ID}-results/"
echo "  ✅ VPC: precision-medicine-vpc"
echo "  ✅ Secret Manager: epic-* secrets created"
echo "  ✅ Service Account: mcp-server-sa"
echo "  ✅ Healthcare Dataset: precision-medicine-dataset"
echo "  ✅ Audit Logging: 10-year retention"
echo "  ✅ Budget Alerts: $100 threshold"
echo ""
echo "Service Account Key: $KEY_FILE"
echo "⚠️  IMPORTANT: Keep this key secure!"
echo ""
echo "Next Steps:"
echo "  1. Update Epic credentials in Secret Manager (when you have them)"
echo "  2. Review the infrastructure in GCP Console"
echo "  3. Proceed to Phase 1B: Epic FHIR Integration"
echo ""
echo "Estimated Monthly Cost: $50-75"
echo "================================================================================"

# Create environment file
cat > .env.gcp <<EOF
# GCP Configuration
export GCP_PROJECT_ID="$PROJECT_ID"
export GCP_REGION="$REGION"
export GOOGLE_APPLICATION_CREDENTIALS="$(pwd)/mcp-server-key.json"

# Storage
export PATIENT_DATA_BUCKET="gs://${PROJECT_ID}-patient-data"
export RESULTS_BUCKET="gs://${PROJECT_ID}-results"

# Healthcare API
export HEALTHCARE_DATASET="precision-medicine-dataset"
export FHIR_STORE="deidentified-fhir-store"
export HEALTHCARE_LOCATION="$REGION"

# Epic FHIR (secrets in Secret Manager)
export EPIC_CLIENT_ID_SECRET="projects/$PROJECT_ID/secrets/epic-client-id"
export EPIC_CLIENT_SECRET_SECRET="projects/$PROJECT_ID/secrets/epic-client-secret"
export EPIC_FHIR_ENDPOINT_SECRET="projects/$PROJECT_ID/secrets/epic-fhir-endpoint"
EOF

echo "Environment file created: .env.gcp"
echo "Source it with: source .env.gcp"
echo ""
