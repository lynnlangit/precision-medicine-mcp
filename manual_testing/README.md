# Manual Testing Resources

This folder contains scripts and documentation for manually testing the Spatial MCP POC.

## 📁 Contents

### Scripts (Executable)

| File | Purpose | Usage |
|------|---------|-------|
| `install_dependencies.sh` | Install all dependencies for 9 MCP servers | `./install_dependencies.sh` |
| `verify_servers.sh` | Verify all servers can be imported | `./verify_servers.sh` |
| `setup_and_test_servers.sh` | Combined setup and verification | `./setup_and_test_servers.sh --install` |
| `test_all_servers.py` | Python-based server verification | `python3 test_all_servers.py` |

### Documentation

| File | Purpose |
|------|---------|
| `MANUAL_TESTING_GUIDE.md` | **Comprehensive step-by-step testing guide** |
| `TESTING_SUMMARY.md` | Quick reference summary |
| `TESTING_STATUS.md` | Complete testing status and verification checklist |
| `CLAUDE_DESKTOP_TEST_PROMPTS.md` | **8 ready-to-paste test prompts with absolute paths for Claude Desktop** |

### Testing Directories

| Directory | Purpose |
|-----------|---------|
| `Solution-Testing/` | Test results and verification documentation |
| `PatientOne-OvarianCancer/` | **Comprehensive end-to-end patient scenario with synthetic data** |

**See:** `PatientOne-OvarianCancer/implementation/` for complete test suite

---

## 🚀 Quick Start

### 1. Install All Server Dependencies

```bash
cd manual_testing
./install_dependencies.sh
```

This will:
- Create Python 3.11 virtual environments for each server
- Install FastMCP and all dependencies
- Set up 9 servers in development mode

**Time:** ~5-10 minutes

### 2. Verify All Servers

```bash
./verify_servers.sh
```

Expected output:
```
Servers working: 9/9
Total tools: 36
🎉 All MCP servers are operational!
```

**Time:** ~10-30 seconds

---

## 📖 Full Testing Guide

For complete testing instructions, see:
- **`MANUAL_TESTING_GUIDE.md`** - Detailed step-by-step guide
- **`TESTING_SUMMARY.md`** - Quick reference

---

## ⚠️ Important Notes

### Claude Code vs Claude Desktop

**These scripts run in Claude Code** (VSCode extension):
- ✅ Can install dependencies
- ✅ Can verify server code
- ❌ Cannot orchestrate MCP protocol

**To test MCP workflows:**
- ✅ Use Claude Desktop (standalone app)
- ✅ Configure with `../configs/claude_desktop_config.json`
- ✅ Test with prompts from `CLAUDE_DESKTOP_TEST_PROMPTS.md`

### Python Version Requirement

All servers require **Python 3.11+**. The install script automatically uses `python3.11`.

Check your version:
```bash
python3.11 --version
# Should show: Python 3.11.x
```

---

## 🔧 Troubleshooting

### Issue: Scripts won't execute

**Solution:** Make them executable
```bash
chmod +x *.sh
```

### Issue: Dependencies fail to install

**Solution:** Check Python version
```bash
which python3.11
python3.11 --version
```

If Python 3.11 not found, install it:
```bash
brew install python@3.11  # macOS with Homebrew
```

### Issue: Servers won't import

**Solution:** Reinstall in specific server
```bash
cd ../servers/mcp-fgbio
rm -rf venv
python3.11 -m venv venv
source venv/bin/activate
pip install -e ".[dev]"
```

---

## 📊 What Gets Tested

| Server | Tools | Status |
|--------|-------|--------|
| mcp-fgbio | 4 | ✅ |
| mcp-spatialtools | 8 | ✅ |
| mcp-openimagedata | 3 | ✅ |
| mcp-seqera | 3 | ✅ |
| mcp-huggingface | 3 | ✅ |
| mcp-deepcell | 2 | ✅ |
| mcp-mockepic | 3 | ✅ |
| mcp-tcga | 5 | ✅ |
| mcp-multiomics | 5 | ✅ |
| **TOTAL** | **36** | **✅** |

---

## 🧪 End-to-End Patient Testing

### PatientOne-OvarianCancer Scenario

A complete synthetic patient dataset for testing all 9 MCP servers:

**Patient:** PAT001-OVC-2025 (Sarah Anderson, 58yo, Stage IV HGSOC)
**Data:** 17 synthetic files (clinical, genomic, multi-omics, spatial, imaging)
**Servers Tested:** All 9 (mockepic, fgbio, tcga, multiomics, spatialtools, openimagedata, deepcell, seqera, huggingface)

**Location:** `PatientOne-OvarianCancer/implementation/`

**Test Approach:**
- 5 focused tests (avoid context limits)
- Each test runs independently
- Total time: 25-45 minutes

**Quick Start:**
```bash
cd PatientOne-OvarianCancer/implementation/
cat QUICK_TEST_REFERENCE.md  # Fast overview
cat TEST_1_CLINICAL_GENOMIC.txt  # Copy and paste into Claude Desktop
```

**See:** `TESTING_STRATEGY.md` for complete testing guide

---

## 🎯 Next Steps After Installation

1. ✅ Run installation scripts (this folder)
2. ✅ Verify all servers working
3. 📋 Configure Claude Desktop (`../configs/`)
4. 📋 Test with example prompts (`../docs/spatial/MCP_POC_Example_Prompts.md`)
5. 📋 **Test end-to-end patient scenario** (`PatientOne-OvarianCancer/`)

---

## 📚 Related Documentation

- `../docs/spatial/MCP_POC_Example_Prompts.md` - 18 test prompts
- `../synthetic_data/README.md` - Synthetic test data
- `../architecture/spatial/README.md` - Spatial MCP architecture
- `../architecture/patient-one/README.md` - Patient-one scenario architecture
- `../README.md` - Main project README

---

**Last Updated:** November 12, 2025
**Status:** ✅ Ready for Testing
**New:** PatientOne-OvarianCancer end-to-end testing scenario with 17 synthetic data files
