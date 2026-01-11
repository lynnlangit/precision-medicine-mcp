# GCP Cloud Run Deployment Status

---

## Latest Update: January 9, 2026

**Date:** 2026-01-09
**Scope:** Redeployment of 2 servers with new visualization capabilities
**Result:** ✅ **SUCCESS**

### Updated Servers

| Server | Previous Revision | New Revision | Status |
|--------|-------------------|--------------|--------|
| **mcp-spatialtools** | mcp-spatialtools-00004-xr9 | mcp-spatialtools-00005-r4s | ✅ Deployed |
| **mcp-openimagedata** | mcp-openimagedata-00003-f2m | mcp-openimagedata-00004-vks | ✅ Deployed |

### New Features Added

**mcp-spatialtools (4 new visualization tools):**
- `generate_spatial_heatmap` - Creates spatial coordinate heatmaps
- `generate_gene_expression_heatmap` - Generates clustered gene expression heatmaps
- `generate_region_composition_chart` - Produces region composition bar charts
- `visualize_spatial_autocorrelation` - Visualizes Moran's I spatial autocorrelation

**mcp-openimagedata (2 new visualization tools):**
- `generate_multiplex_composite` - Creates RGB composite images from multiplex channels
- `generate_he_annotation` - Generates annotated H&E images with region overlays

### Dependencies Added

**mcp-spatialtools:**
- `matplotlib>=3.8.0` - Core plotting library
- `seaborn>=0.13.0` - Statistical visualization

**mcp-openimagedata:**
- `matplotlib>=3.8.0` - Image rendering and annotation
- `pandas>=2.2.0` - Timestamp generation for output files

### Deployment Details

**Configuration:**
- Mode: Development (public access with `--allow-unauthenticated`)
- Region: us-central1
- Transport: SSE
- Build context: Repository root with staged shared utilities

**Resource Allocation:**
- mcp-spatialtools: 4Gi memory, 2 CPU (unchanged)
- mcp-openimagedata: 2Gi memory, 2 CPU (unchanged)

**Deployment URLs:**
- mcp-spatialtools: https://mcp-spatialtools-ondu7mwjpa-uc.a.run.app
- mcp-openimagedata: https://mcp-openimagedata-ondu7mwjpa-uc.a.run.app

### Issues Resolved

**Issue: Docker Build Context for Shared Utilities**

**Problem:**
- Dockerfiles expected `_shared_temp/utils/` directory for shared utilities
- Deployment script didn't stage these files before building
- Build failed: `COPY _shared_temp/utils/ /app/shared/utils/` - directory not found

**Solution:**
Updated deployment script to stage shared utilities before build:
```bash
# Stage shared utilities for Docker build
print_info "Staging shared utilities..."
mkdir -p "${server_path}/_shared_temp"
cp -r "${REPO_ROOT}/shared/utils" "${server_path}/_shared_temp/"
```

Added cleanup after deployment:
```bash
# Cleanup staged files
rm -rf "${server_path}/_shared_temp"
```

### Files Modified

**Deployment Script:**
- `scripts/deployment/deploy_to_gcp.sh`
  - Added staging logic before `gcloud run deploy`
  - Added cleanup logic after deployment (success or failure)

**No server code changes required** - existing implementations already complete

### Verification

Both servers verified running with new revisions:
```bash
$ gcloud run services describe mcp-spatialtools --region us-central1
REVISION: mcp-spatialtools-00005-r4s
TRAFFIC: 100%

$ gcloud run services describe mcp-openimagedata --region us-central1
REVISION: mcp-openimagedata-00004-vks
TRAFFIC: 100%
```

### Tool Counts

- **mcp-spatialtools**: 10 → 14 tools (added 4 visualization tools)
- **mcp-openimagedata**: 3 → 5 tools (added 2 visualization tools)

---

## Initial Deployment: December 30, 2025

**Date:** 2025-12-30
**Result:** ✅ **COMPLETE SUCCESS - All 9 Servers Deployed**

---

## Final Deployment Results

📋 **[See Server Status →](../../servers/README.md#-server-status)** - Complete server status, tools, and implementation details

### ✅ All 9 Servers Successfully Deployed and Running

| Server | URL | Status |
|--------|-----|--------|
| **mcp-deepcell** | https://mcp-deepcell-ondu7mwjpa-uc.a.run.app | ✅ Running |
| **mcp-fgbio** | https://mcp-fgbio-ondu7mwjpa-uc.a.run.app | ✅ Running |
| **mcp-huggingface** | https://mcp-huggingface-ondu7mwjpa-uc.a.run.app | ✅ Running |
| **mcp-mockepic** | https://mcp-mockepic-ondu7mwjpa-uc.a.run.app | ✅ Running |
| **mcp-multiomics** | https://mcp-multiomics-ondu7mwjpa-uc.a.run.app | ✅ Running |
| **mcp-openimagedata** | https://mcp-openimagedata-ondu7mwjpa-uc.a.run.app | ✅ Running |
| **mcp-seqera** | https://mcp-seqera-ondu7mwjpa-uc.a.run.app | ✅ Running |
| **mcp-spatialtools** | https://mcp-spatialtools-ondu7mwjpa-uc.a.run.app | ✅ Running |
| **mcp-tcga** | https://mcp-tcga-ondu7mwjpa-uc.a.run.app | ✅ Running |

---

## Issues Identified and Resolved

### Issue 1: Shared Utilities Import Error (5 servers)
**Affected:** mcp-fgbio, mcp-multiomics, mcp-tcga, mcp-seqera, mcp-huggingface

**Problem:**
- Servers used hardcoded path `Path(__file__).resolve().parents[4] / "shared" / "utils"`
- This path doesn't exist in container directory structure
- Caused `IndexError: 4` preventing server startup

**Solution:**
```python
# Try direct import first (works if PYTHONPATH is set)
try:
    from api_retry import retry_with_backoff, optional_api_call
except ImportError:
    # Fallback to parents[4] for local development
    _shared_utils_path = Path(__file__).resolve().parents[4] / "shared" / "utils"
    if str(_shared_utils_path) not in sys.path:
        sys.path.insert(0, str(_shared_utils_path))
    from api_retry import retry_with_backoff, optional_api_call
```

**Dockerfile update:**
```dockerfile
ENV PYTHONPATH=/app/shared/utils:${PYTHONPATH}
```

---

### Issue 2: Environment Variable Caching (all servers)
**Affected:** All 9 servers

**Problem:**
- Cloud Run was caching old `MCP_TRANSPORT=http` from previous deployments
- Dockerfile `ENV` statements alone don't override cached Cloud Run configurations
- Servers started with HTTP transport on 127.0.0.1:8000 instead of SSE on 0.0.0.0:PORT

**Evidence:**
- Dockerfile: `ENV MCP_TRANSPORT=sse`
- Deployed container: `MCP_TRANSPORT=http` (cached old value)

**Solution:**
```bash
# Added to deployment script
gcloud run deploy "${server_name}" \
  --update-env-vars MCP_TRANSPORT=sse \
  # ... other flags
```

This explicitly sets the environment variable, overriding any cached values.

---

### Issue 3: RPY2 Build Dependency (mcp-multiomics only)
**Affected:** mcp-multiomics

**Problem:**
- Dependency `rpy2` requires R (statistical programming language) to build
- R not available in `python:3.11-slim` base image
- Build failed with: `Error: rpy2 in API mode cannot be built without R in the PATH`

**Solution:**
```dockerfile
# Set RPY2_CFFI_MODE=ABI to build rpy2 without R
ENV RPY2_CFFI_MODE=ABI
RUN pip install --no-cache-dir -e .
```

**Note:** The code has a Python fallback for HAllA analysis, so rpy2 can be installed but doesn't need to be functional.

---

## Files Modified

### Server Code (5 files)
- `servers/mcp-fgbio/src/mcp_fgbio/server.py`
- `servers/mcp-multiomics/src/mcp_multiomics/server.py`
- `servers/mcp-tcga/src/mcp_tcga/server.py`
- `servers/mcp-seqera/src/mcp_seqera/server.py`
- `servers/mcp-huggingface/src/mcp_huggingface/server.py`

**Change:** Wrapped shared utilities imports in try/except pattern

### Dockerfiles (6 files)
- All 5 servers above + `servers/mcp-fgbio/Dockerfile`
- Added `ENV PYTHONPATH=/app/shared/utils:${PYTHONPATH}`
- mcp-multiomics: Added `ENV RPY2_CFFI_MODE=ABI`

### Deployment Script
- `scripts/deployment/deploy_to_gcp.sh`
- Added `--update-env-vars MCP_TRANSPORT=sse` to gcloud deploy command

---

## Deployment Process

**Methodical one-at-a-time deployment:**
1. ✅ mcp-fgbio (tested fix)
2. ✅ mcp-multiomics (added RPY2 fix)
3. ✅ mcp-tcga, mcp-seqera, mcp-huggingface (sequential)
4. ✅ mcp-spatialtools, mcp-openimagedata (final batch)

**Already working:**
- ✅ mcp-deepcell
- ✅ mcp-mockepic

**Time:** ~2 hours total deployment time

---

## Verification

All servers verified running with correct configuration:
- ✅ Transport: SSE (not http)
- ✅ Port: Cloud Run assigned port (read from PORT env var)
- ✅ Host: 0.0.0.0 (not 127.0.0.1)
- ✅ Health checks: Passing

**Example log (mcp-fgbio):**
```
INFO: Starting MCP server 'fgbio' with transport 'sse' on http://0.0.0.0:3000/sse
INFO: Default STARTUP TCP probe succeeded after 1 attempt for container "mcp-fgbio-1" on port 3000
```

---

## Testing the Servers

Use the automated test script:
```bash
export ANTHROPIC_API_KEY=your_key_here
python tests/integration/test_all_gcp_servers.py
```

Or test individual servers:
```python
import anthropic

client = anthropic.Anthropic()

response = client.beta.messages.create(
    model="claude-sonnet-4-5",
    max_tokens=512,
    messages=[{
        "role": "user",
        "content": "List the available tools from the fgbio MCP server."
    }],
    mcp_servers=[{
        "type": "url",
        "url": "https://mcp-fgbio-ondu7mwjpa-uc.a.run.app/sse",
        "name": "fgbio",
    }],
    tools=[{"type": "mcp_toolset", "mcp_server_name": "fgbio"}],
    betas=["mcp-client-2025-11-20"]
)

print(response.content[0].text)
```

---

## Commits Made

1. **Fix shared utilities import for container deployment** (9f35320)
   - Try/except import pattern for all 5 servers
   - Added PYTHONPATH to Dockerfiles

2. **Force MCP_TRANSPORT=sse in deployment script** (f6606da)
   - Added --update-env-vars to deployment command

3. **Fix mcp-multiomics build: add RPY2_CFFI_MODE=ABI** (b01c7fe)
   - Allows rpy2 to build without R installed

---

## Summary

**Final Status:** 🎉 **100% Success - All 9 Servers Deployed**

**Key Learnings:**
1. Docker ENV vars can be overridden by Cloud Run cached configs - always use `--update-env-vars` for critical settings
2. Container path structure differs from local development - use PYTHONPATH or try/except imports
3. Optional dependencies (like rpy2) can be built in ABI mode without full installation
4. Methodical one-at-a-time deployment is more efficient for debugging than batch deployments

**Next Steps:**
- ✅ Run automated test suite (`test_all_gcp_servers.py`)
- ✅ Configure Claude Desktop to use deployed servers
- ✅ Document usage examples for each server
