# Claude Desktop Configuration Files

This folder contains configuration files for connecting Claude Desktop to the Precision Medicine MCP servers.

---

## 📁 Files

### `claude_desktop_config.json` ✅ USE THIS FILE

**Status:** ✅ **Ready to use** - Pre-configured for this system

**Purpose:** Production-ready configuration with full absolute paths to all 9 MCP servers.

**Contains:**
- All 9 MCP servers (fgbio, spatialtools, openimagedata, seqera, huggingface, deepcell, epic, tcga, multiomics)
- Full absolute paths to Python 3.11 virtual environments
- All required environment variables
- DRY_RUN mode enabled by default for quick & safe testing

**How to use:**
```bash
# Copy to Claude Desktop config location
cp claude_desktop_config.json ~/Library/Application\ Support/Claude/claude_desktop_config.json

# Restart Claude Desktop
# All 9 servers should now be available
```

---

### `claude_desktop_config.template.json` 📝 TEMPLATE

**Purpose:** Template for other users or different installations.

**Contains:**
- Same structure as the main config
- Placeholder paths: `/ABSOLUTE/PATH/TO/precision-medicine-mcp/...`

**How to use:**
1. Copy template to a new file
2. Replace all `/ABSOLUTE/PATH/TO/` placeholders with your actual installation path
3. Verify all venv Python executables exist
4. Copy to Claude Desktop config location

**Example replacements:**
```
Before: "/ABSOLUTE/PATH/TO/precision-medicine-mcp/servers/mcp-fgbio/venv/bin/python"
After:  "/Users/yourname/projects/precision-medicine-mcp/servers/mcp-fgbio/venv/bin/python"
```

---


## 🔧 Configuration Details

### Why Full Venv Paths?

Claude Desktop runs in a sandboxed environment and cannot resolve `python` or `python3` commands. You **must** provide the full absolute path to the Python executable in each server's virtual environment.

**❌ This will NOT work:**
```json
"command": "python"
```

**✅ This WILL work:**
```json
"command": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-fgbio/venv/bin/python"
```

### Environment Variables

📋 **[See Complete Server Status →](../servers/README.md#-server-status)** - All 10 servers with tools, status, and documentation

Each server requires specific environment variables:

| Server | Required Variables | Purpose |
|--------|-------------------|---------|
| **fgbio** | `FGBIO_REFERENCE_DATA_DIR`<br>`FGBIO_CACHE_DIR`<br>`FGBIO_DRY_RUN` | Reference genome paths<br>Cache location<br>Mock execution mode |
| **spatialtools** | `SPATIAL_DATA_DIR`<br>`SPATIAL_DRY_RUN` | Data directory<br>Mock execution mode |
| **openimagedata** | `IMAGE_DATA_DIR`<br>`IMAGE_DRY_RUN` | Image storage<br>Mock execution mode |
| **seqera** | `SEQERA_DRY_RUN` | Mock execution mode |
| **huggingface** | `HF_DRY_RUN` | Mock execution mode |
| **deepcell** | `DEEPCELL_DRY_RUN` | Mock execution mode |
| **epic** | `EPIC_DRY_RUN` | Mock execution mode |
| **tcga** | `TCGA_DRY_RUN` | Mock execution mode |
| **multiomics** | `MULTIOMICS_DATA_DIR`<br>`MULTIOMICS_CACHE_DIR`<br>`MULTIOMICS_DRY_RUN` | Multi-omics data directory<br>Cache location<br>Mock execution mode |

### DRY_RUN Mode

All servers are configured with `DRY_RUN=true`. This means:
- ✅ Tools simulate execution without external dependencies
- ✅ No actual downloads or intensive computations
- ✅ Perfect for testing workflow orchestration
- ✅ Responses include expected results format

To disable DRY_RUN mode (use real tools):
```json
"FGBIO_DRY_RUN": "false"
```

**Warning:** Disabling DRY_RUN requires:
- Actual bioinformatics tools installed (STAR, samtools, bedtools, etc.)
- Reference genomes downloaded (~30GB for hg38)
- Significant compute resources (16GB+ RAM)

---

## 📍 Installation Path

The current configuration is set up for:
```
/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/
```

If you cloned the repository to a different location, you'll need to:
1. Copy `claude_desktop_config.template.json` to a new file
2. Replace all placeholder paths with your installation path
3. Use that customized config file

---

## ✅ Verification

### Verify Config is Installed
```bash
cat ~/Library/Application\ Support/Claude/claude_desktop_config.json | grep "mcpServers"
```

Expected output: Should show the config JSON

### Verify All Server Paths Exist
```bash
for server in fgbio spatialtools openimagedata seqera huggingface deepcell epic tcga multiomics; do
  echo "Checking mcp-$server..."
  ls /Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-$server/venv/bin/python
done
```

Expected: All 9 Python executables should be found

### Verify in Claude Desktop

After installing the config and restarting Claude Desktop:

```
What MCP servers are available?
```

Expected response: List of all 9 servers with their tools

---

## 🚀 Quick Start

1. **Install server dependencies** (if not done):
   ```bash
   cd ../manual_testing
   ./install_dependencies.sh
   ```

2. **Install config**:
   ```bash
   cp claude_desktop_config.json ~/Library/Application\ Support/Claude/claude_desktop_config.json
   ```

3. **Restart Claude Desktop**

4. **Verify servers loaded**:
   ```
   What MCP servers are available?
   ```

5. **Run first test**:
   - Open `../manual_testing/PatientOne-OvarianCancer/implementation/TEST_1_CLINICAL_GENOMIC.txt`
   - Copy and paste test into Claude Desktop

---

## 🔧 Troubleshooting

### Issue: "Could not connect to MCP server"

**Cause:** Config file not in correct location or invalid paths

**Solutions:**
1. Verify config location:
   ```bash
   ls -l ~/Library/Application\ Support/Claude/claude_desktop_config.json
   ```

2. Check JSON syntax is valid:
   ```bash
   python3 -m json.tool ~/Library/Application\ Support/Claude/claude_desktop_config.json
   ```

3. Verify all venv paths exist (see verification section above)

4. Restart Claude Desktop completely (Quit and relaunch)

### Issue: Some servers load, others don't

**Cause:** Individual server venv not properly installed

**Solution:** Reinstall specific server:
```bash
cd ../servers/mcp-fgbio  # Change to problematic server
rm -rf venv
python3.11 -m venv venv
source venv/bin/activate
pip install -e ".[dev]"
```

### Issue: Config file gets overwritten

**Cause:** Claude Desktop may reset config on updates

**Solution:** Keep a backup and re-apply:
```bash
# Create backup
cp claude_desktop_config.json ~/Desktop/claude_desktop_config_backup.json

# Restore if needed
cp ~/Desktop/claude_desktop_config_backup.json ~/Library/Application\ Support/Claude/claude_desktop_config.json
```

---



