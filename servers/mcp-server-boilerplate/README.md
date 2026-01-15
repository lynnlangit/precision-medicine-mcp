# MCP Server Boilerplate Template

**Purpose:** Reusable template for creating new MCP servers for precision medicine modalities

**Usage:** Copy this directory and customize for your modality

---

## Quick Start

```bash
# 1. Copy template
cp -r docs/mcp-server-boilerplate servers/mcp-{your-modality}

# 2. Find and replace placeholders
cd servers/mcp-{your-modality}
find . -type f -exec sed -i '' 's/{{SERVER_NAME}}/{your-modality}/g' {} +
find . -type f -exec sed -i '' 's/{{MODALITY_DESCRIPTION}}/Your description/g' {} +

# 3. Rename directories
mv src/mcp_template src/mcp_{your-modality}

# 4. Install dependencies
python3.11 -m venv venv
source venv/bin/activate
pip install -e ".[dev]"

# 5. Customize server.py with your tools
# Edit src/mcp_{your-modality}/server.py

# 6. Test
export {{SERVER_NAME}}_DRY_RUN=true
pytest tests/ -v
```

---

## Placeholders to Replace

| Placeholder | Replace With | Example |
|-------------|--------------|---------|
| `{{SERVER_NAME}}` | Your server name (lowercase) | `metabolomics` |
| `{{SERVER_NAME_UPPER}}` | Uppercase version | `METABOLOMICS` |
| `{{MODALITY_DESCRIPTION}}` | Brief description | `Metabolomics data analysis` |
| `{{TOOL_EXAMPLE}}` | Example tool name | `identify_differential_metabolites` |

---

## Directory Structure

```
mcp-server-boilerplate/
├── README.md                          # This file
├── pyproject.toml                     # Python dependencies
├── src/
│   └── mcp_template/
│       ├── __init__.py
│       ├── __main__.py                # Entry point
│       └── server.py                  # Main server implementation
└── tests/
    ├── test_server.py                 # Smoke tests
    └── fixtures/                      # Test data (create as needed)
        └── .gitkeep
```

---

## Customization Guide

### 1. Update pyproject.toml
- Change `name = "mcp-template"` to your server name
- Add modality-specific dependencies
- Update description

### 2. Update server.py
- Replace example tools with your tools
- Update DRY_RUN logic for your data types
- Add modality-specific imports

### 3. Add Tests
- Copy patterns from `tests/unit/mcp-multiomics/`
- Create fixtures for your data types
- Aim for >35% coverage

### 4. Write Documentation
- Create `servers/mcp-{name}/README.md` with tool descriptions
- Create `architecture/{modality}/README.md` with workflow diagram
- Add examples for Claude Desktop usage

---

## Integration Checklist

After customizing, complete the [Integration Checklist](../../docs/for-developers/automation-guides/ADD_NEW_MODALITY_SERVER.md#integration-checklist) in the main guide.

---

## Next Steps

1. **Read the full guide:** [docs/for-developers/automation-guides/ADD_NEW_MODALITY_SERVER.md](../../docs/for-developers/automation-guides/ADD_NEW_MODALITY_SERVER.md)
2. **Study existing servers:** [servers/](../../servers/)
3. **Test locally:** Use Claude Desktop to verify tools work
4. **Deploy to GCP:** Follow [deployment guide](../../docs/deployment/DEPLOYMENT_STATUS.md)

---

**Questions?** See [Getting Help](../../docs/for-developers/automation-guides/ADD_NEW_MODALITY_SERVER.md#getting-help) in the main guide.
