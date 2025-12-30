# Testing GCP-Deployed MCP Servers

## ⚠️ **CRITICAL: Claude Desktop vs Claude API**

### The Important Distinction:

| Client | Transport | Use Case | GCP Servers? |
|--------|-----------|----------|--------------|
| **Claude Desktop** | STDIO only | Local development | ❌ **NO** - Cannot connect to remote servers |
| **Claude API** | HTTP/SSE | Production, remote servers | ✅ **YES** - Designed for remote MCP servers |

---

## **Answer to Your Question:**

**Q: Can I still use Claude Desktop to test the 9 MCP servers deployed on GCP?**

**A: NO** - Claude Desktop only supports STDIO transport (local process spawning). It cannot connect to remote HTTP/SSE endpoints.

**To test GCP-deployed servers, you must use:**
1. **Claude API** (via Python/TypeScript SDK)
2. **Direct HTTP requests** (curl, Postman)
3. **Custom test scripts** (provided below)

---

## Why Claude Desktop Can't Connect to GCP

### Claude Desktop Configuration (STDIO):
```json
{
  "command": "python",
  "args": ["-m", "mcp_multiomics"]
}
```
→ Spawns **local process** with stdin/stdout pipes

### GCP Cloud Run (HTTP/SSE):
```
https://mcp-multiomics-xxxxx.run.app/sse
```
→ **Remote HTTPS endpoint**, not a local process

**Claude Desktop has no mechanism to connect to HTTP/SSE endpoints.**

---

## Testing Options for GCP-Deployed Servers

### **Option 1: Claude API (Python SDK)** ✅ RECOMMENDED

**Install SDK:**
```bash
pip install anthropic
```

**Python Test Script:**
```python
import anthropic
import os

# Initialize Claude API client
client = anthropic.Anthropic(
    api_key=os.environ.get("ANTHROPIC_API_KEY")
)

# Configure MCP servers (after deployment)
response = client.beta.messages.create(
    model="claude-opus-4-5",
    max_tokens=2000,
    messages=[{
        "role": "user",
        "content": "Use the multiomics server to show available HAllA analysis capabilities."
    }],
    mcp_servers=[
        {
            "type": "url",
            "url": "https://mcp-fgbio-xxxxx.run.app/sse",
            "name": "fgbio",
        },
        {
            "type": "url",
            "url": "https://mcp-multiomics-xxxxx.run.app/sse",
            "name": "multiomics",
        },
        {
            "type": "url",
            "url": "https://mcp-spatialtools-xxxxx.run.app/sse",
            "name": "spatialtools",
        },
        # ... add all 9 servers
    ],
    tools=[{
        "type": "mcp_toolset",
        "mcp_server_name": "multiomics"
    }],
    betas=["mcp-client-2025-11-20"]
)

print(response.content)
```

---

### **Option 2: Direct HTTP Health Checks**

Test each server's health endpoint:

```bash
# After deployment, use the URLs from deployment_urls.txt

# Test mcp-fgbio
curl https://mcp-fgbio-xxxxx.run.app/health

# Test mcp-multiomics
curl https://mcp-multiomics-xxxxx.run.app/health

# Test all servers
while IFS='=' read -r server url; do
  echo "Testing $server..."
  curl -s "$url/health" && echo " ✓" || echo " ✗"
done < deployment_urls.txt
```

---

### **Option 3: Automated Test Script**

I'll create a comprehensive test script that verifies all deployed servers:

**File:** `test_gcp_servers.py` (see below)

---

## Authentication for GCP-Deployed Servers

### Current Deployment (Public):
```bash
gcloud run deploy mcp-multiomics \
  --allow-unauthenticated  # ← Anyone can access
```

### Production (Authenticated):
```bash
gcloud run deploy mcp-multiomics \
  --no-allow-unauthenticated  # ← Requires authentication

# Generate auth token
gcloud auth print-identity-token
```

**Claude API with authentication:**
```python
mcp_servers=[{
    "type": "url",
    "url": "https://mcp-multiomics-xxxxx.run.app/sse",
    "name": "multiomics",
    "authorization_token": "Bearer YOUR_GCP_TOKEN_HERE"
}]
```

---

## Migration Path: Local → GCP

### Phase 1: Local Development (Current)
```
Claude Desktop (STDIO)
  ↓ spawns local processes
  ├── mcp-fgbio (local)
  ├── mcp-multiomics (local)
  └── ... (7 more servers)
```

**Use Case:** Development, testing, debugging

---

### Phase 2: GCP Deployment (Future)
```
Claude API (HTTP/SSE)
  ↓ HTTP requests
  ├── https://mcp-fgbio.run.app
  ├── https://mcp-multiomics.run.app
  └── ... (7 more servers)
```

**Use Case:** Production, multi-user, scaling

---

## Hybrid Approach (Best of Both Worlds)

### Keep Both Environments:

**Local (Claude Desktop):**
- Fast iteration
- No API costs
- Easy debugging
- File system access

**GCP (Claude API):**
- Production analysis
- Team collaboration
- Scalable compute
- Centralized logging

**Workflow:**
1. Develop/test locally with Claude Desktop (STDIO)
2. Deploy to GCP when ready for production
3. Use Claude API to access GCP-deployed servers
4. Keep local setup for rapid development

---

## Cost Comparison

| Environment | Cost | Speed | Use Case |
|-------------|------|-------|----------|
| **Claude Desktop (Local)** | Free (your hardware) | Fastest | Development |
| **Claude API + GCP** | $0.05-0.20/request + $0.10/hr compute | Fast | Production |

**Recommendation:** Use both - local for dev, GCP for production.

---

## Verification Workflow

### After Deployment:

**Step 1: Verify Health**
```bash
./tests/integration/test_gcp_servers.sh
```

**Step 2: Test with Claude API**
```python
python tests/integration/test_claude_api_integration.py
```

**Step 3: Run PatientOne Workflow**
```python
python run_patientone_gcp.py
```

---

## FAQ

### Q: Can I use my existing Claude Desktop config with GCP servers?
**A:** No. Claude Desktop only supports STDIO (local). You must use Claude API for GCP servers.

### Q: Do I need to rewrite my servers for GCP?
**A:** No! FastMCP automatically supports both transports. Just set `MCP_TRANSPORT=http` env var.

### Q: Will my local tests still work after deploying to GCP?
**A:** Yes! Your local Claude Desktop setup is unchanged. You can use both simultaneously.

### Q: How do I switch between local and GCP testing?
**A:**
- **Local:** Use Claude Desktop (STDIO)
- **GCP:** Use Claude API (HTTP/SSE)
- No code changes needed

### Q: What's the best way to test GCP servers during development?
**A:**
1. Test locally first with Claude Desktop (fast iteration)
2. Deploy to GCP when stable
3. Use Claude API for final validation
4. Keep local setup for continued development

---

## Next Steps

1. ✅ **Deploy servers:** `./scripts/deployment/deploy_to_gcp.sh`
2. ✅ **Verify health:** `./tests/integration/test_gcp_servers.sh`
3. ✅ **Test with Claude API:** Use Python SDK (examples above)
4. ✅ **Set up authentication:** For production use

---

**Key Takeaway:**
- **Claude Desktop = Local STDIO only**
- **Claude API = Remote HTTP/SSE required for GCP**
- **Keep both environments for optimal workflow**

**Last Updated:** December 29, 2025
