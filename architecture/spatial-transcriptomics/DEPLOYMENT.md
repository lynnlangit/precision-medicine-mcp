# Deployment Guide

GCP Cloud Run deployment procedures and status.

---

## Current Deployment Status

**Last Deployment:** January 10, 2026
**Servers Deployed:** 6/10 (spatialtools, openimagedata, fgbio, multiomics, epic, mockepic)
**Mode:** Development (public access) + Production (epic - private)
**Region:** us-central1

---

## Deployed Servers

### mcp-spatialtools
- **Revision:** mcp-spatialtools-00005-r4s
- **Date:** Jan 9, 2026
- **Tools:** 14 (10 analysis + 4 visualization)
- **Resources:** 4Gi memory, 2 CPU
- **URL:** https://mcp-spatialtools-ondu7mwjpa-uc.a.run.app

### mcp-openimagedata
- **Revision:** mcp-openimagedata-00004-vks
- **Date:** Jan 9, 2026
- **Tools:** 5 (3 analysis + 2 visualization)
- **Resources:** 2Gi memory, 2 CPU
- **URL:** https://mcp-openimagedata-ondu7mwjpa-uc.a.run.app

---

## Deployment Script

**Location:** `/scripts/deployment/deploy_to_gcp.sh`

**Usage:**
```bash
# Deploy single server (development mode)
./deploy_to_gcp.sh --development --server mcp-spatialtools

# Deploy all servers
./deploy_to_gcp.sh --development

# Production mode (authenticated)
./deploy_to_gcp.sh --production --server mcp-spatialtools
```

**Pre-deployment checklist:**
1. ✅ GCP project configured: precision-medicine-poc
2. ✅ Cloud Run API enabled
3. ✅ Authenticated: `gcloud auth login`
4. ✅ Docker images build successfully
5. ✅ Shared utilities staged correctly

---

## Configuration

### Environment Variables
```bash
# Spatialtools
SPATIAL_DRY_RUN=false
SPATIAL_DATA_DIR=/app/data/spatial
SPATIAL_CACHE_DIR=/app/data/cache/spatial
MCP_TRANSPORT=sse
MCP_PORT=3002

# OpenImageData
IMAGE_DRY_RUN=false
IMAGE_DATA_DIR=/app/data/images
IMAGE_OUTPUT_DIR=/app/data/output
MCP_TRANSPORT=sse
MCP_PORT=3004
```

### Resource Allocation
- **spatialtools:** 4Gi memory, 2 CPU (analysis + visualization)
- **openimagedata:** 2Gi memory, 2 CPU (image processing)
- **fgbio:** 2Gi memory, 2 CPU (QC tools)
- **multiomics:** 4Gi memory, 2 CPU (HAllA integration)

---

## Rollback Procedure

```bash
# Get previous revision
PREV_REVISION=$(gcloud run services describe mcp-spatialtools \
  --region=us-central1 \
  --format="value(status.latestReadyRevisionName)")

# Rollback to previous
gcloud run services update-traffic mcp-spatialtools \
  --region=us-central1 \
  --to-revisions=${PREV_REVISION}=100
```

**Rollback revisions saved:**
- spatialtools: mcp-spatialtools-00004-xr9
- openimagedata: mcp-openimagedata-00003-f2m

---

## Monitoring

**Cloud Run Logs:**
```bash
gcloud run services logs read mcp-spatialtools \
  --region=us-central1 \
  --limit=100
```

**Health Check:**
```bash
curl https://mcp-spatialtools-ondu7mwjpa-uc.a.run.app/health
```

---

See [DEPLOYMENT_STATUS.md](../../docs/deployment/DEPLOYMENT_STATUS.md) for detailed deployment history.
