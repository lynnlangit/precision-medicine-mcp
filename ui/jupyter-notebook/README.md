# Jupyter Notebook MCP Client - Precision Medicine

Interactive Jupyter notebook for testing MCP servers via Claude API. Perfect for data scientists and researchers who want to analyze spatial transcriptomics, multi-omics, and genomics data.

🌐 **Live JupyterLab:** https://jupyter-mcp-notebook-305650208648.us-central1.run.app

<img src="../../data/images/jupyter-preview.png" width="800" alt="Jupyter Notebook Preview" style="display:none">

## Features

- 📓 **Interactive Notebook** - Run analysis step-by-step with Python
- 🔧 **9 MCP Servers** - Access all deployed bioinformatics tools
- 💬 **MCPClient Helper Class** - Simple API for calling MCP servers
- 📊 **Token Usage Tracking** - Monitor API costs per query
- 📈 **Built-in Visualizations** - matplotlib, seaborn, plotly examples
- 💾 **Reproducible Analysis** - Save and share workflows
- 🎯 **Example Workflows** - Pre-built analyses for common tasks

## Quick Start (2 minutes)

### Prerequisites

- Python 3.11+
- Anthropic API key ([get one here](https://console.anthropic.com/))

### Local Installation

```bash
# 1. Navigate to notebook directory
cd ui/jupyter-notebook

# 2. Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. Set your API key
export ANTHROPIC_API_KEY=your_key_here  # On Windows: set ANTHROPIC_API_KEY=your_key_here

# Or create a .env file:
cp .env.example .env
# Edit .env and add your API key

# 5. Start Jupyter
jupyter notebook mcp_client.ipynb
```

The notebook will open in your browser at http://localhost:8888

### Cloud Access (No Installation)

Simply open the live JupyterLab instance:

👉 **https://jupyter-mcp-notebook-305650208648.us-central1.run.app**

**Note:** The cloud instance is public and resets on each deployment. For persistent work, use local installation or save your notebooks elsewhere.

## Notebook Contents

The `mcp_client.ipynb` notebook includes:

### 1. Setup and Configuration
- Import dependencies
- Configure API key
- Define MCP server endpoints

### 2. MCPClient Helper Class
- Simple API for calling MCP servers
- Conversation history tracking
- Token usage and cost estimation

### 3. Quick Test
List all available tools from any MCP server:

```python
client = MCPClient()
result = client.call_servers(
    prompt="List all available tools from the spatialtools server.",
    servers=["spatialtools"]
)
print(result["response"])
```

### 4. Example Workflows

**Spatial Transcriptomics Analysis:**
```python
result = client.call_servers(
    prompt="Analyze spatial transcriptomics data for Patient-001. Perform cell type deconvolution.",
    servers=["spatialtools"]
)
```

**Pathway Enrichment:**
```python
result = client.call_servers(
    prompt="For genes [TP53, BRCA1, MYC], perform pathway enrichment using GO_BP.",
    servers=["spatialtools"]
)
```

**Multi-omics Integration:**
```python
result = client.call_servers(
    prompt="Integrate RNA, protein, and phospho data. Run HAllA analysis.",
    servers=["multiomics"]
)
```

**Complete PatientOne Workflow:**
```python
result = client.call_servers(
    prompt="""For Patient-001 (ovarian cancer):
    1. Get clinical data
    2. Retrieve spatial transcriptomics
    3. Perform cell type deconvolution
    4. Run differential expression
    5. Generate treatment recommendations""",
    servers=["spatialtools", "multiomics"]
)
```

### 5. Token Usage Tracking

Every query returns usage statistics:

```python
{
    "response": "...",
    "usage": {
        "input_tokens": 5124,
        "output_tokens": 483,
        "total_tokens": 5607,
        "estimated_cost_usd": 0.0226
    },
    "model": "claude-sonnet-4-5",
    "servers_used": ["spatialtools"]
}
```

### 6. Visualization Examples

Plot token usage over time:

```python
import matplotlib.pyplot as plt

# Track usage across multiple queries
usage_data = []
for query in queries:
    result = client.call_servers(query["prompt"], query["servers"])
    usage_data.append(result["usage"])

# Plot
plt.plot([u["total_tokens"] for u in usage_data])
plt.title("Token Usage Over Time")
plt.show()
```

## Available MCP Servers

All 9 servers are pre-configured:

### Production Servers (Real Analysis)
- **fgbio** - Genomic reference data, FASTQ validation
- **multiomics** - Multi-omics integration (RNA/Protein/Phospho)
- **spatialtools** - Spatial transcriptomics analysis

### Mock Servers (Demo Only)
- **tcga** - TCGA genomic data
- **openimagedata** - Imaging data access
- **seqera** - Nextflow workflow execution
- **huggingface** - Model integration
- **deepcell** - Image segmentation
- **mockepic** - FHIR clinical data

## MCPClient API Reference

### Class: `MCPClient`

**Constructor:**
```python
client = MCPClient(api_key: str = None)
```
- `api_key`: Anthropic API key (reads from ANTHROPIC_API_KEY env var if not provided)

**Methods:**

**`call_servers()`**
```python
result = client.call_servers(
    prompt: str,                    # Your query
    servers: List[str],             # Server names to use
    model: str = "claude-sonnet-4-5",  # Claude model
    max_tokens: int = 4096,         # Max response length
    clear_history: bool = False     # Reset conversation
)
```

Returns:
```python
{
    "response": str,      # Claude's response text
    "usage": dict,        # Token counts and cost
    "model": str,         # Model used
    "servers_used": list  # Servers called
}
```

## Configuration

### API Key Security

**The ANTHROPIC_API_KEY is stored differently depending on environment:**

| Environment | Storage Method | Security |
|-------------|---------------|----------|
| **Local Development** | `.env` file (gitignored) | Not committed to git, local machine only |
| **Cloud JupyterLab** | Environment variable (encrypted) | Encrypted at rest, managed by Google Cloud |
| **Browser/Client** | Never exposed | Key stays server-side, never sent to browser |

**Important Security Notes:**
- ✅ `.env` file is in `.gitignore` - never committed to git
- ✅ Cloud Run environment variables are encrypted at rest
- ✅ API key is only used server-side in container
- ✅ Use separate API keys for development vs production
- ❌ Never hardcode API keys in notebooks
- ❌ Never commit `.env` files to git

### Environment Variables

**For Local Development:**

Create a `.env` file (from `.env.example`):

```bash
# Required
ANTHROPIC_API_KEY=your_key_here
```

**For Cloud Deployment:**

API key is passed as environment variable during deployment:

```bash
export ANTHROPIC_API_KEY=your_key_here
./deploy.sh
```

The deployment script automatically sets the key as a Cloud Run environment variable (encrypted).

### MCP Server URLs

All server URLs are pre-configured in the notebook. To add/modify:

```python
MCP_SERVERS = {
    "your_server": {
        "url": "https://your-server.run.app/sse",
        "description": "Server description"
    }
}
```

## Cost Estimates

**Per Query:**
- Input: ~500-2000 tokens ($0.003-0.012 with Sonnet)
- Output: ~1000-4000 tokens ($0.015-0.060 with Sonnet)
- **Total: ~$0.02-0.08 per query**

**Typical Analysis Session (10 queries):**
- ~$0.20-0.80 total

**Cloud JupyterLab (Cloud Run):**
- ~$0.50-2.00 per hour of active use
- Scales to zero when not in use
- No cost when idle

**See:** [Cost Analysis](../../docs/for-hospitals/operations/COST_ANALYSIS.md) for detailed breakdowns

## Troubleshooting

### "API Key Missing" Error

```bash
# Set the environment variable
export ANTHROPIC_API_KEY=your_key_here

# Or create .env file
echo "ANTHROPIC_API_KEY=your_key_here" > .env
```

### MCP Servers Not Responding

1. Check server status: [GCP Cloud Run Console](https://console.cloud.google.com/run)
2. Verify URLs in notebook match deployed URLs
3. Test individual server: `curl https://mcp-spatialtools-ondu7mwjpa-uc.a.run.app/sse`

### Slow Responses

- Use **claude-haiku-4** for faster responses
- Reduce **max_tokens** parameter
- Select fewer MCP servers per query

### Kernel Crashes

```bash
# Restart Jupyter kernel
# In notebook: Kernel → Restart

# Or restart entire Jupyter server
jupyter notebook stop
jupyter notebook mcp_client.ipynb
```

### Module Import Errors

```bash
# Reinstall dependencies
pip install -r requirements.txt

# Or upgrade specific package
pip install --upgrade anthropic
```

## Development

### Project Structure

```
ui/jupyter-notebook/
├── mcp_client.ipynb       # Main Jupyter notebook
├── requirements.txt       # Python dependencies
├── .env.example          # Environment variable template
├── .gitignore            # Git ignore rules
├── Dockerfile            # Container for Cloud Run
├── .dockerignore         # Docker ignore rules
├── deploy.sh             # Cloud Run deployment script
├── deploy_to_vertex_ai.sh # Vertex AI Workbench deployment
├── test_notebook.py      # Local testing script
└── README.md             # This file
```

### Adding New Analysis Workflows

Create a new cell in the notebook:

```python
# Example: Custom pathway enrichment
result = client.call_servers(
    prompt="""
    For the gene list [GENE1, GENE2, GENE3]:
    1. Perform pathway enrichment
    2. Filter for p-value < 0.05
    3. Return top 10 pathways
    """,
    servers=["spatialtools"],
    max_tokens=2048
)

print(result["response"])
```

### Adding Visualizations

```python
import pandas as pd
import matplotlib.pyplot as plt

# Parse results (example)
data = pd.DataFrame({
    "pathway": ["Path1", "Path2", "Path3"],
    "p_value": [0.001, 0.01, 0.05]
})

# Plot
plt.bar(data["pathway"], -np.log10(data["p_value"]))
plt.ylabel("-log10(p-value)")
plt.title("Pathway Enrichment")
plt.show()
```

### Saving Analysis Results

```python
# Save conversation history
import json

with open("analysis_results.json", "w") as f:
    json.dump(client.conversation_history, f, indent=2)

# Save specific results
result = client.call_servers(...)
with open("pathway_results.txt", "w") as f:
    f.write(result["response"])
```

## Deployment

### Local Development

```bash
jupyter notebook mcp_client.ipynb
```

### GCP Cloud Run (Current Deployment)

```bash
# Set API key
export ANTHROPIC_API_KEY=your_key_here

# Deploy
./deploy.sh
```

**Live URL:** https://jupyter-mcp-notebook-305650208648.us-central1.run.app

**Cost:** ~$0.50-2.00 per hour of active use (pay-per-use)

### GCP Vertex AI Workbench (Alternative)

For a dedicated managed Jupyter instance:

```bash
# Set API key
export ANTHROPIC_API_KEY=your_key_here

# Deploy
./deploy_to_vertex_ai.sh
```

**Cost:** ~$150-200/month (always-on VM)

**When to use:**
- Need persistent environment
- Working with large datasets
- Long-running analyses
- Need GPUs for deep learning

## Roadmap

**Planned Features:**
- [ ] Interactive widgets for parameter tuning
- [ ] Automated report generation (PDF/HTML)
- [ ] Workflow templates library
- [ ] Integration with Google Cloud Storage
- [ ] Multi-user authentication
- [ ] Notebook scheduling (cron jobs)
- [ ] GPU support for image analysis
- [ ] R kernel support

**Want to contribute?** Open an issue or pull request on GitHub!

## Support

- **Issues:** [GitHub Issues](https://github.com/lynnlangit/precision-medicine-mcp/issues)
- **Documentation:** [Main README](../../README.md)
- **MCP Spec:** [Model Context Protocol](https://modelcontextprotocol.io/)
- **Anthropic Docs:** [Claude API](https://docs.anthropic.com/)

## License

See the main repository [LICENSE](../../LICENSE) file.

---

**Built with:**
- [Jupyter](https://jupyter.org/) - Interactive computing
- [Anthropic Claude API](https://www.anthropic.com/) - AI model
- [Model Context Protocol](https://modelcontextprotocol.io/) - MCP standard
- [GCP Cloud Run](https://cloud.google.com/run) - Container hosting

**Part of the Precision Medicine MCP suite** - Enabling AI-driven bioinformatics for cancer research.
