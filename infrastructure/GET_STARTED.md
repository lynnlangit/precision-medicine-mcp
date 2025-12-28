# Phase 1A: Let's Get Started! 🚀

**Time Required:** 30 minutes - 1 hour
**Cost:** Free (using GCP $300 credit)

---

## Prerequisites Check

Before we start, make sure you have:

### 1. GCP Account Setup ✅

Check if you have gcloud installed:
```bash
gcloud --version
```

**If not installed:**
```bash
# macOS
curl https://sdk.cloud.google.com | bash
exec -l $SHELL

# Then authenticate
gcloud auth login
```

### 2. Get Your Billing Account ID

```bash
gcloud billing accounts list
```

You should see something like:
```
ACCOUNT_ID            NAME                OPEN  MASTER_ACCOUNT_ID
01234-567890-ABCDEF  My Billing Account  True
```

**Copy your ACCOUNT_ID** - you'll need it!

---

## Option 1: Automated Setup (Recommended)

I've created an automated script that will guide you through each step.

### Run the Setup Script

```bash
cd /Users/lynnlangit/Documents/GitHub/spatial-mcp/infrastructure

# Run the setup script
./setup-gcp-phase1a.sh
```

**The script will ask you for:**
1. **Project ID** - Choose a unique ID (e.g., `precision-medicine-poc-yourname`)
2. **Billing Account ID** - The ID you got from the command above
3. **Region** - Press Enter for default (`us-central1`) or type your preferred region

Then it will:
- Create your GCP project
- Set up encrypted storage
- Configure networking
- Set up secrets
- Enable audit logging
- Create service accounts
- Configure Healthcare API
- Set budget alerts

**Time:** ~10-15 minutes

---

## Option 2: Manual Step-by-Step

If you prefer to understand each step, follow the detailed guide:

Open: `implementation/PHASE1A_GCP_SETUP.md`

This guide explains every command and what it does.

---

## During Setup

The script will:
- ✅ Pause at each major step
- ✅ Ask you to confirm before proceeding
- ✅ Show you what it's doing
- ✅ Handle errors gracefully (you can re-run if something fails)

**You can stop at any time** (Ctrl+C) and resume later!

---

## After Setup

When complete, you'll have:

1. **Infrastructure Created:**
   - GCP project
   - Encrypted storage buckets
   - Private VPC network
   - Secret Manager (for Epic credentials)
   - Service account with key
   - Healthcare API dataset
   - Audit logging (10-year retention)
   - Budget alerts

2. **Files Created:**
   - `mcp-server-key.json` - Service account key (keep secure!)
   - `.env.gcp` - Environment variables

3. **Cost Estimate:**
   - ~$50-75/month (covered by $300 free credit for 4-6 months)

---

## Troubleshooting

### "Permission denied" Error

Make sure you're authenticated:
```bash
gcloud auth login
```

### "Project ID already exists" Error

Choose a different project ID - they must be globally unique:
```bash
# Try adding your initials or a number
precision-medicine-poc-jd
precision-medicine-poc-2024
```

### "Billing account required" Error

Add a payment method:
1. Go to: https://console.cloud.google.com/billing
2. Add payment method
3. Note: You won't be charged without your approval

### "API not enabled" Error

The script should enable all APIs automatically. If you see this error:
```bash
gcloud services enable healthcare.googleapis.com --project=YOUR_PROJECT_ID
```

### Script Fails Midway

You can safely re-run the script - it will skip resources that already exist.

---

## Validation

After setup completes, verify everything works:

```bash
# Source the environment file
source .env.gcp

# Test storage access
echo "test" | gsutil cp - ${PATIENT_DATA_BUCKET}/test.txt
gsutil cat ${PATIENT_DATA_BUCKET}/test.txt
gsutil rm ${PATIENT_DATA_BUCKET}/test.txt

# Test secret access
gcloud secrets versions access latest --secret=epic-fhir-endpoint --project=$GCP_PROJECT_ID

# Check your infrastructure
gcloud projects describe $GCP_PROJECT_ID
```

If all commands succeed, you're ready for Phase 1B! ✅

---

## Next Steps

Once Phase 1A is complete:

1. ✅ Register for Epic Sandbox (while waiting)
   - Go to: https://fhir.epic.com/
   - Create account and register app
   - Save credentials

2. ✅ Review what was created
   - Go to: https://console.cloud.google.com
   - Explore your new infrastructure

3. ✅ Proceed to Phase 1B
   - Open: `implementation/PHASE1B_EPIC_INTEGRATION.md`
   - Build the mcp-epic server

---

## Questions?

- **What did we just create?** Read `PHASE1A_GCP_SETUP.md` for details
- **Why do we need this?** It's HIPAA-compliant infrastructure for patient data
- **Can I delete this later?** Yes, delete the entire project: `gcloud projects delete PROJECT_ID`
- **How much will this cost?** ~$50-75/month, covered by $300 credit

---

## Ready?

Let's do this! Run:

```bash
cd /Users/lynnlangit/Documents/GitHub/spatial-mcp/infrastructure
./setup-gcp-phase1a.sh
```

**Estimated time:** 10-15 minutes ⏱️

Good luck! 🚀
