# Installation Guide

Complete guide for setting up the Precision Medicine MCP system with Claude Desktop.

## Table of Contents

- [Quick Start (5 Minutes)](#quick-start-5-minutes)
- [Complete Setup](#complete-setup)
- [Verification](#verification)
- [Usage Examples](#usage-examples)
- [Troubleshooting](#troubleshooting)

---

## Quick Start (5 Minutes)

If you already have Claude Desktop installed and want to get started quickly:

### 1. Verify Prerequisites

```bash
python --version  # Should be 3.11 or higher
git --version
```

### 2. Clone and Setup

```bash
git clone https://github.com/lynnlangit/spatial-mcp.git
cd spatial-mcp
chmod +x scripts/setup_environment.sh
./scripts/setup_environment.sh
```

### 3. Configure Claude Desktop

```bash
# macOS
cp configs/claude_desktop_config.json \
   ~/Library/Application\ Support/Claude/claude_desktop_config.json

# Update paths in the config file to match your installation
```

### 4. Restart Claude Desktop

Quit Claude Desktop completely (Cmd+Q on macOS), then relaunch it.

### 5. Test It

Open Claude Desktop and try:

```
What MCP servers are available?
```

You should see all 10 servers listed.

**🎉 Done! Jump to [Usage Examples](#usage-examples) to start analyzing spatial transcriptomics data.**

---

## Complete Setup

For first-time installation or if you need detailed instructions:

### System Requirements

- **Operating System:** Ubuntu 24.04 LTS, macOS 13+, or Windows 11 with WSL2
- **Python:** 3.11 or higher
- **Memory:** 16GB RAM minimum (32GB recommended for production)
- **Disk Space:** 50GB free space
- **Network:** Internet connection for downloading reference data

### 1. Install Claude Desktop

Download and install from [claude.ai/download](https://claude.ai/download).

### 2. Clone Repository

```bash
git clone https://github.com/lynnlangit/spatial-mcp.git
cd spatial-mcp
```

### 3. Environment Setup

The setup script creates the directory structure and installs dependencies:

```bash
chmod +x scripts/setup_environment.sh
./scripts/setup_environment.sh
```

**What this script does:**
- Creates data directories (reference, patient-data, cache)
- Sets up Python virtual environments for each server
- Installs dependencies
- Configures environment variables

**Directory Structure Created:**

```
spatial-mcp/
├── data/
│   ├── reference/          # Reference genomes
│   ├── patient-data/       # Patient data (PatientOne synthetic data)
│   ├── cache/              # Temporary files and caches
│   └── multiomics/         # Multi-omics data
├── servers/
│   ├── mcp-fgbio/          # Genomic reference data server
│   ├── mcp-spatialtools/   # Spatial processing server
│   ├── mcp-multiomics/     # Multi-omics integration server
│   ├── mcp-tcga/           # TCGA data server
│   └── ... (6 more servers)
└── configs/
    └── claude_desktop_config.json
```

### 4. Install MCP Servers

Each server has its own virtual environment:

```bash
# Example: Install mcp-spatialtools
cd servers/mcp-spatialtools
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -e ".[dev]"
pytest  # Run tests to verify
deactivate
```

Repeat for other servers you need.

### 5. Configure Environment Variables

Each server can use a `.env` file for configuration:

**Example: servers/mcp-spatialtools/.env**

```bash
# Data directories
SPATIAL_DATA_DIR=/absolute/path/to/spatial-mcp/data
SPATIAL_CACHE_DIR=/absolute/path/to/spatial-mcp/data/cache

# Mode
SPATIAL_DRY_RUN=false  # Set to true for testing

# Performance
SPATIAL_TIMEOUT_SECONDS=300
```

### 6. Configure Claude Desktop

**Configuration File Location:**

- **macOS:** `~/Library/Application Support/Claude/claude_desktop_config.json`
- **Linux:** `~/.config/Claude/claude_desktop_config.json`
- **Windows:** `%APPDATA%\Claude\claude_desktop_config.json`

**Configuration Template:**

```json
{
  "mcpServers": {
    "spatialtools": {
      "command": "python",
      "args": ["-m", "mcp_spatialtools"],
      "cwd": "/absolute/path/to/spatial-mcp/servers/mcp-spatialtools",
      "env": {
        "PYTHONPATH": "/absolute/path/to/spatial-mcp/servers/mcp-spatialtools/src",
        "SPATIAL_DATA_DIR": "/absolute/path/to/spatial-mcp/data",
        "SPATIAL_CACHE_DIR": "/absolute/path/to/spatial-mcp/data/cache",
        "SPATIAL_DRY_RUN": "false"
      }
    },
    "epic": {
      "command": "python",
      "args": ["-m", "mcp_epic"],
      "cwd": "/absolute/path/to/spatial-mcp/servers/mcp-epic",
      "env": {
        "PYTHONPATH": "/absolute/path/to/spatial-mcp/servers/mcp-epic/src",
        "GCP_PROJECT_ID": "precision-medicine-poc",
        "GCP_REGION": "us-central1",
        "GOOGLE_APPLICATION_CREDENTIALS": "/absolute/path/to/spatial-mcp/infrastructure/mcp-server-key.json"
      }
    }
  }
}
```

**⚠️ Important:**
- Replace all `/absolute/path/to/spatial-mcp` with your actual installation path
- Use absolute paths, not relative paths
- Verify JSON syntax before saving

### 7. Restart Claude Desktop

After editing configuration:
1. Quit Claude Desktop completely (not just close window)
2. Relaunch Claude Desktop
3. Servers will load automatically

---

## Verification

### 1. Check Server Installation

```bash
cd servers/mcp-spatialtools
source venv/bin/activate
python -c "import mcp_spatialtools; print('✓ mcp-spatialtools installed')"
```

### 2. Run Tests

```bash
cd servers/mcp-spatialtools
pytest -v

# Expected output:
# test_server.py::test_differential_expression PASSED
# test_server.py::test_spatial_autocorrelation PASSED
# ... (more tests)
```

### 3. Test with Claude Desktop

Open Claude Desktop and try:

**Check Available Servers:**
```
What MCP servers are available?
```

Expected: Should list all configured servers (spatialtools, epic, fgbio, etc.)

**List Server Tools:**
```
What tools does the spatialtools server provide?
```

Expected: Should list 14 tools including perform_differential_expression, calculate_spatial_autocorrelation, etc.

**Test a Tool:**
```
Using the spatialtools server, list the available test data files.
```

Expected: Should return information about available patient data.

---

## Usage Examples

### Available MCP Servers

After installation, you have access to **10 MCP servers**:

📋 **[Complete Server Status](/docs/architecture/servers.md)** - Detailed capabilities matrix

**Quick Reference:**

- ✅ **spatialtools** - Spatial transcriptomics analysis (14 tools, 95% real implementation)
- ✅ **epic** - GCP Healthcare API FHIR integration (4 tools, 100% real)
- 🟡 **fgbio, multiomics, tcga, openimagedata, seqera, huggingface, deepcell, mockepic** - Available in DRY_RUN mode

**Legend:**
- ✅ **ACTIVE** = Real data mode, fully functional
- 🟡 **DRY_RUN** = Mock mode for testing

### Example 1: Analyze Patient-001 (Complete Workflow)

```
I want to analyze Patient-001 (PAT001-OVC-2025) ovarian cancer spatial data.

Please run the complete analysis:
1. Load the data from /path/to/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial/
2. Perform differential expression (tumor_core vs stroma)
3. Calculate spatial autocorrelation for all genes
4. Do cell type deconvolution
5. Summarize findings and clinical implications
```

### Example 2: Identify Drug Resistance Markers

```
Using the spatialtools server, analyze Patient-001 to identify drug resistance markers:

1. Run differential expression comparing tumor_core vs stroma regions
2. Filter for genes with log2FC > 2
3. Check if ABCB1, PIK3CA, or AKT1 are upregulated
4. Calculate spatial autocorrelation for resistance genes
5. Tell me if resistance markers are spatially clustered
```

### Example 3: Access FHIR Clinical Data

```
Using the epic MCP server, retrieve Patient-001 data from the GCP Healthcare API:

1. Get patient demographics
2. List all conditions
3. Show CA-125 observations
4. What medications has the patient received?
```

### Example 4: Batch Correction (Multi-Site Study)

```
I have spatial data from 3 different hospitals. Run batch correction:

Files:
- /path/to/hospital1_visium.csv
- /path/to/hospital2_visium.csv
- /path/to/hospital3_visium.csv

Use ComBat to harmonize the data and tell me:
- What % of variance was due to batch effects?
- How much did batch correction reduce variance?
- Is the biological signal preserved?
```

### Tips for Best Results

**1. Be Specific About Inputs**

❌ Bad:
```
Analyze some spatial data
```

✅ Good:
```
Analyze Patient-001 spatial data located at /Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial/
Compare tumor_core vs stroma regions
```

**2. Request Interpretation**

```
What do these differential expression results mean for treatment?
```

```
Are there any actionable drug targets based on this pathway enrichment?
```

**3. Iterate and Refine**

```
1. First: "What are the top differentially expressed genes?"
2. Then: "Tell me more about the PI3K/AKT pathway genes"
3. Finally: "Are these genes spatially clustered? What does that mean?"
```

---

## Troubleshooting

### Issue: "No MCP servers found"

**Solution:** Restart Claude Desktop

```bash
# macOS: Cmd+Q, then relaunch
# Verify config file exists
cat ~/Library/Application\ Support/Claude/claude_desktop_config.json
```

### Issue: "Tool execution failed"

**Check:**
1. File paths are correct (use absolute paths)
2. Data files exist and are readable
3. Server logs (Claude Desktop → Settings → MCP Logs)

**Enable Debug Mode:**

```json
"env": {
  "SPATIAL_LOG_LEVEL": "DEBUG"
}
```

### Issue: "Server not responding"

**Test server manually:**

```bash
cd /path/to/spatial-mcp/servers/mcp-spatialtools
source venv/bin/activate
python -m mcp_spatialtools

# Expected output: "Server running on stdio..."
```

### Issue: "Permission denied" for GCP

**Fix:** Ensure service account key exists:

```bash
ls -la /path/to/spatial-mcp/infrastructure/mcp-server-key.json
```

### Issue: "ModuleNotFoundError"

**Solutions:**

1. **Verify PYTHONPATH in config:**
   ```json
   "env": {
     "PYTHONPATH": "/absolute/path/to/spatial-mcp/servers/mcp-spatialtools/src"
   }
   ```

2. **Reinstall in development mode:**
   ```bash
   cd servers/mcp-spatialtools
   source venv/bin/activate
   pip install -e .
   ```

### Issue: Performance/Timeouts

**Solutions:**

1. **Increase timeout:**
   ```json
   "env": {
     "SPATIAL_TIMEOUT_SECONDS": "600"
   }
   ```

2. **Check disk space:**
   ```bash
   df -h /path/to/spatial-mcp/data/
   ```

3. **Use local cache on SSD**

### Issue: Validate JSON Configuration

**Check syntax:**

```bash
python -m json.tool < ~/Library/Application\ Support/Claude/claude_desktop_config.json
```

---

## Next Steps

**Now that you're set up:**

1. ✅ Try the [Example Prompts](#usage-examples) above
2. ✅ Explore [PatientOne Testing Scenario](/docs/test-docs/patient-one-scenario/README.md)
3. ✅ Read [Test Prompts](/docs/test-docs/patient-one-scenario/test-prompts/) for comprehensive workflows
4. ✅ Review [Architecture Documentation](/docs/architecture/README.md)

**For Development:**
- See [Developer Guide](../for-developers/README.md) - Architecture, contributing, building servers
- Add new tools following existing patterns
- Run tests before submitting changes

---

## Related Documentation

- **[Architecture Overview](/docs/architecture/README.md)** - System design
- **[Server Status](/docs/architecture/servers.md)** - Implementation status and capabilities
- **[User Guides](/docs/guides/README.md)** - Role-specific guides
- **[Test Documentation](/docs/test-docs/README.md)** - Testing strategies and prompts

---

**Last Updated:** 2026-01-13

**Support:**
- Documentation: `/docs` directory
- Issues: GitHub Issues
- MCP Specification: https://modelcontextprotocol.io/
