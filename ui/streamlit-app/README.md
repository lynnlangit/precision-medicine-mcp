# MCP Chat - Streamlit UI for Precision Medicine MCP Servers

A visual chat interface for testing deployed MCP servers on GCP Cloud Run. Provides a Claude Desktop-like experience for bioinformatics workflows.

🌐 **Live Demo:** https://streamlit-mcp-chat-ondu7mwjpa-uc.a.run.app

<img src="../../data/images/streamlit-ui-preview.png" width="800" alt="Streamlit Chat UI Preview">

## Features

- 💬 **Chat Interface** - Natural language interaction with MCP servers
- 🔧 **Server Selection** - Choose which of the 9 MCP servers to use
- 🎯 **Example Prompts** - Quick-start templates for common workflows
- 📊 **Token Usage** - Track API usage per message
- 🎨 **Clean UI** - Simple, Claude Desktop-like interface
- ⚡ **Real-time** - Instant responses from deployed servers

## Quick Start (2 minutes)

### Prerequisites

- Python 3.11+
- Anthropic API key ([get one here](https://console.anthropic.com/))

### Installation

```bash
# 1. Navigate to the UI directory
cd ui/streamlit-app

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

# 5. Run the app
streamlit run app.py
```

The app will open in your browser at http://localhost:8501

## Usage

### 1. Select MCP Servers

Use the sidebar to select which servers to enable:

**Production Servers (Real Analysis):**
- **fgbio** - Genomic reference data and FASTQ validation
- **multiomics** - Multi-omics integration (RNA/Protein/Phospho)
- **spatialtools** - Spatial transcriptomics analysis

**Mock Servers (Demo Only):**
- tcga, openimagedata, seqera, huggingface, deepcell, mockepic

### 2. Choose a Model

Select from:
- **claude-sonnet-4-5** (recommended) - Best balance of speed/quality
- **claude-opus-4-5** - Highest quality, slower
- **claude-haiku-4** - Fastest, most cost-effective

### 3. Start Chatting

Type your question or:
- Click "Example Prompts" in sidebar to load pre-built queries
- Try spatial analysis, pathway enrichment, multi-omics integration

### 4. View Responses

- Chat history shows full conversation
- Token usage displayed per message
- Server status cards show active servers

## Example Workflows

### Quick Tests

**Spatial Analysis:**
```
Analyze the spatial transcriptomics data for Patient-001.
Perform cell type deconvolution and identify key cell populations.
```

**Pathway Enrichment:**
```
For the upregulated genes [TP53, BRCA1, MYC, KRAS],
perform pathway enrichment analysis using GO_BP database.
```

**Multi-omics Integration:**
```
Integrate RNA, protein, and phosphorylation data.
Run HAllA association analysis and identify significant correlations.
```

### Complete PatientOne Workflow

```
For Patient-001 (ovarian cancer):
1. Get clinical data from FHIR
2. Retrieve spatial transcriptomics data
3. Perform cell type deconvolution
4. Run differential expression between tumor core and margin
5. Generate treatment recommendations
```

## Configuration

### API Key Security

**The ANTHROPIC_API_KEY is stored differently depending on environment:**

| Environment | Storage Method | Security |
|-------------|---------------|----------|
| **Local Development** | `.env` file (gitignored) | Not committed to git, local machine only |
| **GCP Cloud Run** | Environment variable (encrypted) | Encrypted at rest, managed by Google Cloud |
| **Browser/Client** | Never exposed | Key stays on server, never sent to browser |

**Important Security Notes:**
- ✅ `.env` file is in `.gitignore` - never committed to git
- ✅ Cloud Run environment variables are encrypted at rest
- ✅ API key is only used server-side, never exposed to browser
- ✅ Use separate API keys for development vs production
- ❌ Never hardcode API keys in source code
- ❌ Never commit `.env` files to git

### Environment Variables

**For Local Development:**

Create a `.env` file (from `.env.example`):

```bash
# Required
ANTHROPIC_API_KEY=your_key_here

# Optional
DEFAULT_MODEL=claude-sonnet-4-5
DEFAULT_MAX_TOKENS=4096
```

**For Cloud Run Deployment:**

API key is passed as environment variable during deployment:

```bash
export ANTHROPIC_API_KEY=your_key_here
./deploy.sh
```

The deployment script automatically sets the key as a Cloud Run environment variable (encrypted).

### MCP Server Configuration

Server URLs are configured in `utils/mcp_config.py`. All 9 servers are pre-configured with GCP Cloud Run URLs.

To add/modify servers:
```python
# Edit utils/mcp_config.py
MCP_SERVERS = {
    "your_server": {
        "name": "your_server",
        "url": "https://your-server.run.app/sse",
        "description": "Server description",
        "status": "production",  # or "mock"
        "tools_count": 5
    }
}
```

## Architecture

```
Streamlit UI (Browser)
    ↓
Python Chat Handler
    ↓
Anthropic Claude API (with MCP support)
    ↓
GCP Cloud Run MCP Servers (9 servers)
    ↓
Bioinformatics Tools (STAR, ComBat, HAllA, etc.)
```

## Cost Estimates

**Per Message:**
- Input: ~500-2000 tokens ($0.003-0.012 with Sonnet)
- Output: ~1000-4000 tokens ($0.015-0.060 with Sonnet)
- **Total: ~$0.02-0.08 per exchange**

**Typical Session (10 messages):**
- ~$0.20-0.80 total

**See:** [Cost Analysis](../../docs/hospital-deployment/operations/COST_ANALYSIS.md) for detailed breakdowns

## Troubleshooting

### "API Key Missing" Error

```bash
# Set the environment variable
export ANTHROPIC_API_KEY=your_key_here

# Or create .env file
echo "ANTHROPIC_API_KEY=your_key_here" > .env
```

### Servers Not Responding

1. Check server status: [GCP Cloud Run Console](https://console.cloud.google.com/run)
2. Verify URLs in `utils/mcp_config.py` match deployed URLs
3. Test individual server: `curl https://mcp-spatialtools-ondu7mwjpa-uc.a.run.app/sse`

### Slow Responses

- Use **claude-haiku-4** for faster responses
- Reduce **max_tokens** slider
- Select fewer MCP servers
- Check GCP Cloud Run logs for server performance

### Connection Errors

```python
# Error: "Connection refused"
# Fix: Verify server URLs are correct (check deployment docs)

# Error: "Rate limit exceeded"
# Fix: Wait 60 seconds, then retry (Anthropic API rate limits)

# Error: "Invalid API key"
# Fix: Check ANTHROPIC_API_KEY is set correctly
```

## Development

### Project Structure

```
ui/streamlit-app/
├── app.py                 # Main Streamlit application
├── requirements.txt       # Python dependencies
├── .env.example          # Environment variable template
├── .gitignore            # Git ignore rules
├── README.md             # This file
└── utils/
    ├── __init__.py       # Package init
    ├── mcp_config.py     # MCP server configurations
    └── chat_handler.py   # Claude API integration
```

### Adding New Features

**Add a new example prompt:**
```python
# Edit utils/mcp_config.py
EXAMPLE_PROMPTS["Your New Prompt"] = "Your prompt text here..."
```

**Add custom styling:**
```python
# Edit app.py, add to st.markdown() CSS block
st.markdown("""
<style>
    .your-custom-class {
        /* your styles */
    }
</style>
""", unsafe_allow_html=True)
```

**Add response visualization:**
```python
# In app.py, after displaying response:
if "spatial_data" in response_text:
    import pandas as pd
    # Create visualization
    st.plotly_chart(your_plot)
```

## Deployment

### Local Development

```bash
streamlit run app.py
```

### Streamlit Cloud (Free)

1. Push code to GitHub
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Connect your repo
4. Add secrets (ANTHROPIC_API_KEY)
5. Deploy!

**Cost:** Free for public apps

### GCP Cloud Run

```bash
# Create Dockerfile
cat > Dockerfile << 'EOF'
FROM python:3.11-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY . .
EXPOSE 8501
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
EOF

# Deploy
gcloud run deploy streamlit-mcp-chat \
  --source . \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated \
  --set-env-vars ANTHROPIC_API_KEY=your_key_here
```

**Cost:** ~$5-20/month (Cloud Run pay-per-use)

## Roadmap

**Planned Features:**
- [ ] Data visualization (spatial plots, pathway networks)
- [ ] File upload (FASTQ, VCF, spatial data)
- [ ] Export conversation to PDF/Markdown
- [ ] Workflow templates (save/load common workflows)
- [ ] Multi-user support (authentication)
- [ ] Response streaming (real-time token display)
- [ ] Server health monitoring
- [ ] Cost tracking dashboard

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
- [Streamlit](https://streamlit.io/) - Web UI framework
- [Anthropic Claude API](https://www.anthropic.com/) - AI model
- [Model Context Protocol](https://modelcontextprotocol.io/) - MCP standard
- [GCP Cloud Run](https://cloud.google.com/run) - Server hosting

**Part of the Precision Medicine MCP suite** - Enabling AI-driven bioinformatics for cancer research.
