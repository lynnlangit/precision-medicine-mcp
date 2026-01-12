# Deployment Guide

For comprehensive deployment information, see:

## 📋 [Deployment Status](../../docs/deployment/DEPLOYMENT_STATUS.md)
Complete deployment history with revision numbers, new features, and configuration details.

## 🧪 [GCP Testing Guide](../../docs/deployment/GCP_TESTING_GUIDE.md)
How to test deployed servers using Claude API (HTTP/SSE remote access).

## 🔒 [Security Guide](../../docs/deployment/SECURITY.md)
API key and secrets management, GCP Secret Manager integration.

---

## Quick Reference

**Deployment Script:** `/scripts/deployment/deploy_to_gcp.sh`

```bash
# Deploy single server (development)
./scripts/deployment/deploy_to_gcp.sh --development --server mcp-spatialtools

# Deploy all servers
./scripts/deployment/deploy_to_gcp.sh --development
```

**Current Status (as of Jan 10, 2026):**
- **Deployed:** 6/10 servers (spatialtools, openimagedata, fgbio, multiomics, epic, mockepic)
- **Mode:** Development (public) + Production (epic - private)
- **Region:** us-central1

---

For detailed deployment procedures, rollback instructions, monitoring setup, and troubleshooting, see the comprehensive [Deployment Status](../../docs/deployment/DEPLOYMENT_STATUS.md) documentation.
