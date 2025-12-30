# Quick Server Verification

## Single Command Verification ⚡

Verify all 9 MCP servers are working in **2-3 seconds**:

```bash
python3 tests/verify_all_servers.py
```

### Expected Output:

```
============================================================
MCP Server Verification - Quick Import Test
============================================================
Repo: /Users/lynnlangit/Documents/GitHub/spatial-mcp

Testing mcp-fgbio... ✓
Testing mcp-multiomics... ✓
Testing mcp-spatialtools... ✓
Testing mcp-tcga... ✓
Testing mcp-openimagedata... ✓
Testing mcp-seqera... ✓
Testing mcp-huggingface... ✓
Testing mcp-deepcell... ✓
Testing mcp-mockepic... ✓

============================================================
Summary
============================================================

  mcp-fgbio                 ✓ PASS
  mcp-multiomics            ✓ PASS
  mcp-spatialtools          ✓ PASS
  mcp-tcga                  ✓ PASS
  mcp-openimagedata         ✓ PASS
  mcp-seqera                ✓ PASS
  mcp-huggingface           ✓ PASS
  mcp-deepcell              ✓ PASS
  mcp-mockepic              ✓ PASS

Results: 9/9 servers passed
Time: 2.5s

✅ All servers verified successfully!
Ready for Claude Desktop testing.
```

---

## What This Tests

For each server, it verifies:
- ✅ Virtual environment exists
- ✅ Server module can be imported
- ✅ DRY_RUN mode is configured
- ✅ No import errors or missing dependencies

---

## If Any Server Fails

**Example failure:**
```
Testing mcp-fgbio... ✗ venv not found
```

**Fix:**
```bash
cd servers/mcp-fgbio
python3 -m venv venv
venv/bin/pip install -e .
```

**Then re-run:**
```bash
python3 tests/verify_all_servers.py
```

---

## Next Steps

### After Verification Passes:

**1. Update Claude Desktop Config:**
```bash
cp desktop-configs/claude_desktop_config.json \
   ~/Library/Application\ Support/Claude/claude_desktop_config.json
```

**2. Restart Claude Desktop**

**3. Test in Claude Desktop:**

Open Claude Desktop and paste this prompt:

```
Please verify all MCP servers are connected by listing available tools from each:

1. fgbio - show available reference genomes
2. multiomics - show HAllA capabilities
3. spatialtools - show spatial analysis tools
4. tcga - show TCGA projects
5. openimagedata - show image datasets
6. seqera - show Nextflow workflows
7. huggingface - show genomic models
8. deepcell - show cell segmentation models
9. mockepic - get patient-001 demographics

Just confirm each server is connected and responding.
```

**Expected:** All 9 servers respond with ✓

---

## Troubleshooting

### Script says "venv not found"
```bash
# Install all server dependencies
cd tests/manual_testing/Solution-Testing
./install_dependencies.sh
```

### Script says "Import failed"
```bash
# Check server logs
cd servers/mcp-<server-name>
venv/bin/python -m pytest tests/ -v
```

### Claude Desktop doesn't see servers
1. Check config location: `~/Library/Application Support/Claude/claude_desktop_config.json`
2. Restart Claude Desktop
3. Check logs: `~/Library/Logs/Claude/mcp*.log`

---

**Time:** 2-3 seconds for verification
**Last Updated:** December 29, 2025
