# MCP Cost & Performance Observatory Dashboard

**Real-time observability for MCP server token usage, latency, and costs**

Monitor and optimize your precision medicine MCP server architecture for production deployment. Track "cost per insight" and identify optimization opportunities across all 9 deployed servers.

---

## Quick Start

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Run dashboard
streamlit run streamlit_app.py

# 3. Open browser
# Dashboard automatically opens at http://localhost:8501
```

---

## What This Dashboard Does

### For Platform Builders
- **Cost Visibility**: See exactly how much each tool call costs in real-time
- **Performance Monitoring**: Track latency across all 9 MCP servers
- **Resource Planning**: Project monthly/annual costs based on usage patterns
- **Optimization Guidance**: Get specific recommendations to reduce costs

### For Workflow Architects
- **Cost Per Insight**: Understand the ROI of each analytical step
- **Server Comparison**: Identify bottlenecks and expensive operations
- **Model Selection**: Compare costs across Claude Sonnet/Opus/Haiku
- **Scaling Analysis**: Estimate costs as you scale from demo to production

---

## Dashboard Views

### 📈 Overview
- High-level metrics (total cost, tool calls, duration)
- Cost breakdown by server
- Token distribution (input vs output)
- Workflow execution summary

![Overview](https://via.placeholder.com/800x400?text=Overview+Dashboard)

### 💰 Cost Analysis
- Detailed cost breakdown by server and tool
- Top 10 most expensive tool calls
- Model cost comparison (Sonnet vs Opus vs Haiku)
- Monthly/annual cost projections

### ⚡ Performance
- Latency analysis by server
- Execution timeline (cost/tokens over time)
- Token processing efficiency
- Bottleneck identification

### 🔧 Optimization
- Automated cost-reduction recommendations
- Potential savings calculator
- Strategy comparison (caching, batching, model switching)
- Export reports for planning

---

## Sample Data

The dashboard includes sample metrics from **PatientOne Complete Workflow**:
- **Patient**: PAT001-OVC-2025 (Stage IV HGSOC)
- **Analysis**: End-to-end precision medicine workflow
- **Servers**: All 9 MCP servers (mcp-epic, mcp-spatialtools, mcp-multiomics, etc.)
- **Tool Calls**: 27 tool invocations
- **Total Cost**: $1.45 (Claude Sonnet 4.5)
- **Duration**: 127 seconds

See [sample_data/patientone_workflow.yaml](sample_data/patientone_workflow.yaml) for full metrics.

---

## Architecture

```
ui/dashboard/
├── streamlit_app.py              # Main Streamlit dashboard
├── cost_calculator.py            # Anthropic API pricing calculations
├── metrics_aggregator.py         # Load and process metrics
├── requirements.txt              # Python dependencies
├── sample_data/
│   └── patientone_workflow.yaml  # Sample PatientOne metrics
└── README.md                     # This file
```

### Data Flow

```
Workflow Metrics (YAML)
        ↓
metrics_aggregator.py
        ↓
cost_calculator.py
        ↓
streamlit_app.py
        ↓
Interactive Dashboard
```

---

## Using Your Own Data

### Option 1: YAML File (Recommended for Quick Win)

Create a YAML file with your workflow metrics:

```yaml
workflow_metadata:
  name: "My Analysis Workflow"
  patient_id: "patient-123"
  analysis_date: "2026-01-12T10:00:00Z"
  model: "claude-sonnet-4-5-20250929"

tool_calls:
  - server: "mcp-spatialtools"
    tool: "perform_differential_expression"
    timestamp: "2026-01-12T10:05:00Z"
    input_tokens: 8500
    output_tokens: 2100
    latency_ms: 3500
    status: "success"

  # ... more tool calls ...
```

Load in dashboard:
```python
from metrics_aggregator import MetricsAggregator

aggregator = MetricsAggregator("path/to/your/metrics.yaml")
```

### Option 2: Parse GCP Cloud Logging (Advanced)

Extract metrics from GCP Cloud Run logs:

```bash
# Export logs for a specific timeframe
gcloud logging read \
  "resource.type=cloud_run_revision AND resource.labels.service_name=mcp-spatialtools" \
  --limit=1000 \
  --format=json > mcp_logs.json

# Parse logs and convert to YAML format
python parse_gcp_logs.py mcp_logs.json > my_metrics.yaml
```

*(Log parser script not included - create custom parser based on your log format)*

### Option 3: Live Instrumentation (Future Work)

Add metrics endpoints to MCP servers:

```python
# In each MCP server
@app.get("/metrics")
def get_metrics():
    return {
        "tool_calls": metrics_collector.get_tool_calls(),
        "token_usage": metrics_collector.get_token_usage(),
        "latency": metrics_collector.get_latency_stats()
    }
```

Poll metrics from dashboard for real-time updates.

---

## Cost Calculator Usage

### Standalone Cost Calculations

```python
from cost_calculator import calculate_cost, estimate_monthly_cost, compare_model_costs

# Calculate cost for single tool call
cost = calculate_cost(input_tokens=8000, output_tokens=2000)
print(f"Total: ${cost['total_cost']:.4f}")
# Output: Total: $0.0540

# Compare models
comparison = compare_model_costs(10000, 2000)
for model, data in comparison.items():
    print(f"{model}: ${data['total_cost']:.4f}")
# Output:
# Claude Sonnet 4.5: $0.0600
# Claude Opus 4.5: $0.3000
# Claude Haiku 4: $0.0050

# Project monthly costs
monthly = estimate_monthly_cost(
    daily_tool_calls=50,
    avg_input_tokens=8000,
    avg_output_tokens=2000
)
print(f"Monthly: ${monthly['monthly_cost']:.2f}")
# Output: Monthly: $81.00
```

### Anthropic API Pricing (as of Jan 2025)

| Model | Input (per 1M tokens) | Output (per 1M tokens) |
|-------|----------------------|------------------------|
| Claude Sonnet 4.5 | $3.00 | $15.00 |
| Claude Opus 4.5 | $15.00 | $75.00 |
| Claude Haiku 4 | $0.25 | $1.25 |

**Key Insight**: Haiku is ~92% cheaper than Sonnet. Use for simple lookups!

---

## Optimization Strategies

### Strategy 1: Use Haiku for Simple Queries

**Savings**: 90%+ for low-complexity tools

**Good candidates**:
- `mcp-fgbio.get_reference_genome_info` (reference lookups)
- `mcp-tcga.query_mutation_frequency` (database queries)
- `mcp-epic.get_patient_demographics` (FHIR data retrieval)

**Implementation**:
```python
# In MCP server, add model parameter
@mcp.tool()
async def get_reference_genome_info(genome: str, model: str = "haiku"):
    # Use Haiku for simple reference lookups
    pass
```

### Strategy 2: Cache Expensive Computations

**Savings**: 30-50% for repeated analyses

**Target**: `mcp-multiomics.perform_halla_analysis` (18K input, 6.8K output tokens)

**Implementation**:
```python
import hashlib
import json

def cache_key(data):
    return hashlib.md5(json.dumps(data).encode()).hexdigest()

@mcp.tool()
async def perform_halla_analysis(data):
    key = cache_key(data)
    if key in cache:
        return cache[key]

    result = expensive_halla_computation(data)
    cache[key] = result
    return result
```

### Strategy 3: Batch Visualization Generation

**Savings**: 15-20% by reducing overhead

**Target**: `mcp-spatialtools` (6 separate visualization calls)

**Implementation**:
```python
@mcp.tool()
async def generate_all_visualizations(data):
    # Generate all 6 visualizations in one call
    # Reduces per-call overhead from 6× to 1×
    return {
        "spatial_heatmap": generate_spatial_heatmap(data),
        "gene_expression": generate_gene_expression(data),
        # ...
    }
```

### Strategy 4: Response Streaming

**Savings**: 10-15% by enabling early cancellation

**Implementation**: Use Claude API's streaming mode to stop generation when sufficient information is received.

---

## Cost Projections

### Example: Scaling PatientOne Workflow

**Current Demo Cost**: $1.45 per patient (Sonnet 4.5)

| Workload | Daily Cost | Monthly Cost | Annual Cost |
|----------|-----------|--------------|-------------|
| 10 patients/day | $14.50 | $319 (22 days) | $3,828 |
| 50 patients/day | $72.50 | $1,595 | $19,140 |
| 100 patients/day | $145.00 | $3,190 | $38,280 |

**With Optimizations** (Haiku for simple queries + caching):

| Workload | Daily Cost | Monthly Cost | Annual Cost | Savings |
|----------|-----------|--------------|-------------|---------|
| 10 patients/day | $8.70 | $191 | $2,293 | 40% |
| 50 patients/day | $43.50 | $957 | $11,466 | 40% |
| 100 patients/day | $87.00 | $1,914 | $22,932 | 40% |

---

## Integration with CI/CD

### Automated Cost Tracking

Add dashboard to your CI/CD pipeline:

```yaml
# .github/workflows/cost_check.yml
name: Cost Check

on: [pull_request]

jobs:
  cost-analysis:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3

      - name: Run test workflow
        run: python run_test_workflow.py --metrics-output test_metrics.yaml

      - name: Analyze costs
        run: |
          pip install -r ui/dashboard/requirements.txt
          python -c "
          from metrics_aggregator import MetricsAggregator
          agg = MetricsAggregator('test_metrics.yaml')
          cost = agg.get_cost_breakdown()['total_cost']
          print(f'Workflow cost: ${cost:.4f}')
          if cost > 2.0:
              print('⚠️ Cost exceeds threshold!')
              exit(1)
          "

      - name: Comment on PR
        uses: actions/github-script@v6
        with:
          script: |
            github.rest.issues.createComment({
              issue_number: context.issue.number,
              body: '✅ Cost analysis passed: $1.45'
            })
```

---

## Troubleshooting

### Dashboard Won't Start

**Issue**: `ModuleNotFoundError: No module named 'streamlit'`

**Solution**:
```bash
pip install -r requirements.txt
```

### No Metrics Data

**Issue**: `FileNotFoundError: Sample metrics file not found`

**Solution**: Ensure you're running from the `ui/dashboard/` directory:
```bash
cd ui/dashboard
streamlit run streamlit_app.py
```

### Costs Don't Match Actual API Bills

**Issue**: Dashboard costs are estimates only

**Solution**: Pricing in `cost_calculator.py` reflects Anthropic's published rates (Jan 2025). Verify with:
```bash
python cost_calculator.py
```

Update pricing constants if rates change:
```python
# In cost_calculator.py
PRICING = {
    "claude-sonnet-4-5-20250929": ModelPricing(
        input_price=3.00,  # Update these values
        output_price=15.00,
        name="Claude Sonnet 4.5"
    ),
    # ...
}
```

---

## Deployment Options

### Option 1: Local Development

```bash
streamlit run streamlit_app.py
```

**Use for**: Testing, ad-hoc analysis, local demos

### Option 2: GCP Cloud Run

```bash
# Create Dockerfile
cat > Dockerfile <<EOF
FROM python:3.11-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY . .
EXPOSE 8080
CMD streamlit run streamlit_app.py --server.port=8080 --server.address=0.0.0.0
EOF

# Build and deploy
gcloud builds submit --tag gcr.io/precision-medicine-poc/mcp-dashboard
gcloud run deploy mcp-dashboard \
  --image gcr.io/precision-medicine-poc/mcp-dashboard \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated
```

**Use for**: Team access, production monitoring, stakeholder demos

### Option 3: Embedded in Jupyter Notebook

```python
# In Jupyter notebook
from metrics_aggregator import load_sample_metrics
import plotly.express as px

aggregator = load_sample_metrics("patientone_workflow")
server_summary = aggregator.get_server_summary()

fig = px.bar(server_summary, x='server', y='total_cost')
fig.show()
```

**Use for**: Research notebooks, interactive analysis, reports

---

## Roadmap

### v1.0 (Current - Quick Win)
- ✅ Sample data from PatientOne workflow
- ✅ Cost calculator with Anthropic API pricing
- ✅ Streamlit dashboard with 4 views
- ✅ Optimization recommendations
- ✅ Export reports (TXT, CSV)

### v1.1 (Planned)
- [ ] GCP Cloud Logging integration
- [ ] Real-time metrics polling
- [ ] Multi-workflow comparison
- [ ] Historical trend analysis
- [ ] Alert thresholds

### v2.0 (Future)
- [ ] Live instrumentation of MCP servers
- [ ] WebSocket real-time updates
- [ ] BigQuery data storage
- [ ] Advanced ML-based optimization suggestions
- [ ] Cost anomaly detection

---

## API Reference

### MetricsAggregator

```python
class MetricsAggregator:
    """Aggregates and analyzes MCP server metrics."""

    def __init__(self, metrics_file: str):
        """Load metrics from YAML file."""

    def get_summary_stats(self) -> Dict[str, Any]:
        """Get high-level summary statistics."""

    def get_server_summary(self) -> pd.DataFrame:
        """Get per-server metrics aggregation."""

    def get_cost_breakdown(self) -> Dict[str, float]:
        """Get detailed cost breakdown."""

    def get_top_cost_tools(self, n: int = 10) -> pd.DataFrame:
        """Get top N most expensive tool calls."""

    def get_optimization_recommendations(self) -> List[str]:
        """Get cost optimization recommendations."""
```

### CostCalculator

```python
def calculate_cost(input_tokens: int, output_tokens: int, model: str) -> Dict:
    """Calculate cost in USD for given token usage."""

def compare_model_costs(input_tokens: int, output_tokens: int) -> Dict:
    """Compare costs across all Claude models."""

def estimate_monthly_cost(daily_tool_calls: int, avg_input_tokens: int,
                         avg_output_tokens: int, model: str) -> Dict:
    """Estimate monthly costs based on usage patterns."""

def get_optimization_recommendations(server_metrics: Dict) -> List[str]:
    """Generate cost optimization recommendations."""
```

---

## Related Documentation

- **[Cost Calculator](cost_calculator.py)** - Standalone cost calculations
- **[Metrics Aggregator](metrics_aggregator.py)** - Data processing and analysis
- **[Sample Metrics](sample_data/patientone_workflow.yaml)** - PatientOne workflow data
- **[MCP Servers](../servers/)** - Individual server READMEs
- **[Architecture](../architecture/)** - System design documentation

---

## Support

**Questions?**
- Check [troubleshooting](#troubleshooting) section
- Review [sample_data/patientone_workflow.yaml](sample_data/patientone_workflow.yaml) for data format
- Open issue on GitHub

**Contributing:**
- Add new visualization views to `streamlit_app.py`
- Update pricing in `cost_calculator.py` if Anthropic rates change
- Create parsers for other log formats (GCP, AWS, local logs)

---

**Last Updated:** 2026-01-12
**Version:** 1.0 (Quick Win)
**Maintainer:** precision-medicine-mcp project

**🎉 Monitor your MCP costs. Optimize your architecture. Scale with confidence.**
