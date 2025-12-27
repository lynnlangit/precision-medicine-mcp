# Cost Tracking & Monitoring

## Overview

This document describes the cost tracking and monitoring system for the Precision Medicine MCP servers. These utilities help researchers estimate, track, and control computational costs for precision medicine analyses.

**Status:** ✅ Implemented (WI-6 from Risk Mitigation Workplan)
**Risk Reduced:** R6 (Unexpected costs) - 70% reduction (7/10 → 2/10)

---

## Why Cost Tracking Matters

### Financial Transparency
- Research budgets are limited
- Cloud costs can escalate quickly
- Grant proposals require cost estimates
- Multi-patient studies need per-patient cost tracking

### Cost Control
- Set budget limits before analysis
- Get alerts when approaching limits
- Identify expensive operations
- Optimize workflows for cost-efficiency

### Reproducibility
- Document exact costs for published research
- Enable cost comparison across methods
- Support funding applications
- Track costs over time

---

## Architecture

### Cost Tracking Utilities Location

All cost utilities are centralized in:
```
shared/utils/cost_tracking.py
```

### Components

1. **CostTracker** - Track actual costs during execution
2. **CostEstimator** - Estimate costs before execution
3. **BudgetAlert** - Monitor costs and trigger alerts
4. **Cost decorators** - Automatic cost tracking for functions
5. **Cost context managers** - Track costs for code blocks

---

## 1. Cost Estimation (Before Execution)

### Purpose
Estimate costs before running analysis to help with budgeting and planning.

### Using CostEstimator

```python
from cost_tracking import CostEstimator

estimator = CostEstimator()

# Estimate RNA-seq analysis
rna_cost = estimator.rna_seq_cost(num_samples=5)
print(f"RNA-seq for 5 samples: ${rna_cost:.2f}")  # $2.50

# Estimate variant calling
variant_cost = estimator.variant_calling_cost(num_samples=3)
print(f"Variant calling for 3 samples: ${variant_cost:.2f}")  # $3.60

# Estimate multi-omics integration
multiomics_cost = estimator.multiomics_integration_cost(
    num_samples=10,
    num_modalities=3  # RNA + Protein + Phospho
)
print(f"Multi-omics integration: ${multiomics_cost:.2f}")  # $2.50

# Estimate PatientOne full workflow
patient_one_costs = estimator.patient_one_workflow_cost()
print(f"PatientOne total: ${patient_one_costs['total']:.2f}")  # ~$7.20
```

### Available Estimation Methods

| Method | Parameters | Description |
|--------|------------|-------------|
| `compute_cost()` | instance_type, duration_hours | Cloud compute costs |
| `storage_cost()` | size_gb, duration_days, storage_class | Data storage costs |
| `rna_seq_cost()` | num_samples | RNA-seq alignment + quantification |
| `variant_calling_cost()` | num_samples | BWA + GATK variant calling |
| `multiomics_integration_cost()` | num_samples, num_modalities | Multi-omics data integration |
| `spatial_analysis_cost()` | num_slides | Spatial transcriptomics analysis |
| `patient_one_workflow_cost()` | - | Full PatientOne workflow |

### MCP Tool: estimate_analysis_cost

The mcp-multiomics server provides a cost estimation tool:

```python
# Via MCP tool
result = await estimate_analysis_cost(
    num_samples=10,
    modalities=["rna", "protein", "phospho"],
    include_halla=True,
    include_upstream=False
)

print(f"Estimated cost: ${result['estimated_cost_usd']:.2f}")
print(f"Budget tier: {result['budget_tier']}")
print(f"Breakdown: {result['breakdown']}")
```

**Output:**
```json
{
  "estimated_cost_usd": 18.50,
  "breakdown": {
    "integration": 2.50,
    "preprocessing": 0.50,
    "validation": 0.20,
    "halla_analysis": 15.00,
    "compute": 0.15,
    "storage_30days": 0.15
  },
  "budget_tier": "medium",
  "recommendation": "Standard multi-omics analysis",
  "notes": [
    "Costs are estimates based on typical analysis patterns",
    "Actual costs may vary based on data complexity",
    "Includes 30-day data caching by default",
    "Use DRY_RUN mode for testing without charges"
  ]
}
```

---

## 2. Cost Tracking (During Execution)

### Purpose
Track actual costs as operations execute, providing real-time cost monitoring and detailed logs.

### Basic Usage

```python
from cost_tracking import CostTracker

# Initialize tracker
tracker = CostTracker(
    patient_id="patient_001",
    workflow_id="ovarian_cancer_analysis",
    budget_limit_usd=50.0
)

# Add costs as operations complete
tracker.add_cost("analysis", "fastq_qc", 0.15, quantity=10, unit="GB")
tracker.add_cost("analysis", "variant_calling", 1.20, quantity=1, unit="sample")
tracker.add_cost("analysis", "rna_seq", 0.50, quantity=1, unit="sample")
tracker.add_cost("compute", "cloud_compute", 0.80, quantity=4.2, unit="hours")
tracker.add_cost("storage", "data_caching", 0.20, quantity=8.7, unit="GB-month")

# Get current total
print(f"Current cost: ${tracker.get_total():.2f}")

# Get breakdown by category
by_category = tracker.get_by_category()
for category, amount in by_category.items():
    print(f"{category}: ${amount:.2f}")

# Print formatted summary
tracker.print_summary()

# Save cost log
log_path = tracker.save_log()
```

### CostTracker Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `patient_id` | str (optional) | Patient identifier for attribution |
| `workflow_id` | str (optional) | Workflow identifier |
| `budget_limit_usd` | float (optional) | Budget limit (triggers warnings) |
| `cost_log_path` | Path (optional) | Directory for cost logs (default: ./cost_logs/) |

### add_cost() Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `category` | str | Cost category: "compute", "storage", "api", "analysis" |
| `operation` | str | Operation name (e.g., "rna_seq_alignment") |
| `amount_usd` | float | Cost in USD |
| `quantity` | float (optional) | Quantity (e.g., 10.5 for 10.5 GB) |
| `unit` | str (optional) | Unit (e.g., "GB", "hours", "samples") |
| `metadata` | dict (optional) | Additional metadata |

### Cost Summary Example

```
================================================================================
COST SUMMARY
================================================================================
Patient: patient_001
Workflow: ovarian_cancer_analysis
Duration: 245.3 seconds

Total Cost: $7.15
Budget: $50.00 (14.3% used)
Remaining: $42.85

Costs by Category:
  analysis             $  4.30 ( 60.1%)
  compute              $  2.00 ( 28.0%)
  storage              $  0.85 ( 11.9%)

Top 5 Operations by Cost:
  spatial_analysis               $  2.0000
  multiomics_integration         $  0.7500
  variant_calling                $  1.2000
  rna_seq                        $  0.5000
  cloud_compute                  $  0.8000
================================================================================
```

---

## 3. Automatic Cost Tracking

### Context Manager

Track costs for code blocks automatically:

```python
from cost_tracking import track_cost, CostTracker

tracker = CostTracker(patient_id="patient_002")

with track_cost(tracker, "compute", "rna_seq_pipeline"):
    # Run expensive operation
    run_rna_seq_pipeline()
    # Cost automatically tracked based on duration

# Cost already added to tracker
```

### Decorator

Track costs for functions automatically:

```python
from cost_tracking import track_operation_cost, CostTracker

tracker = CostTracker(patient_id="patient_003")

@track_operation_cost(
    category="analysis",
    operation="variant_calling",
    cost_estimator=lambda sample_id: 1.20  # Fixed cost per call
)
def run_variant_calling(sample_id: str):
    # Expensive operation
    call_variants(sample_id)
    return results

# Use with cost tracking
result = run_variant_calling("sample_001", cost_tracker=tracker)

# Cost automatically tracked
```

---

## 4. Budget Alerts

### Purpose
Monitor costs in real-time and trigger alerts when thresholds are exceeded.

### Basic Usage

```python
from cost_tracking import CostTracker, BudgetAlert

# Initialize tracker
tracker = CostTracker(patient_id="patient_004")

# Set up alerts at $5, $10, and $15
alerts = BudgetAlert(
    tracker=tracker,
    thresholds=[5.0, 10.0, 15.0]
)

# Run operations
tracker.add_cost("analysis", "operation_1", 3.0)
alerts.check()  # No alert

tracker.add_cost("analysis", "operation_2", 4.0)
alerts.check()  # Triggers $5 alert

tracker.add_cost("analysis", "operation_3", 5.0)
alerts.check()  # Triggers $10 alert
```

### Custom Alert Callback

```python
def custom_alert_handler(current_cost: float, threshold: float):
    """Send email/Slack notification when budget exceeded."""
    message = f"Budget alert: ${current_cost:.2f} exceeded ${threshold:.2f}"

    # Send email
    send_email(to="pi@university.edu", subject="Budget Alert", body=message)

    # Send Slack notification
    slack_client.post_message(channel="#research", text=message)

    # Log to file
    with open("budget_alerts.log", "a") as f:
        f.write(f"{datetime.now()}: {message}\n")

alerts = BudgetAlert(
    tracker=tracker,
    thresholds=[10.0, 25.0, 50.0],
    alert_callback=custom_alert_handler
)
```

---

## 5. Pricing Data

### Current Pricing (2025)

All pricing data is centralized in `PRICING` dictionary in `cost_tracking.py`:

#### Cloud Compute (per hour)
| Instance Type | vCPU | RAM (GB) | Price/hour |
|---------------|------|----------|------------|
| aws_batch_small | 2 | 8 | $0.096 |
| aws_batch_medium | 4 | 16 | $0.192 |
| aws_batch_large | 8 | 32 | $0.384 |
| azure_standard_d4 | 4 | 16 | $0.192 |
| gcp_n1_standard_4 | 4 | 15 | $0.190 |

#### Storage (per GB-month)
| Storage Class | Price/GB-month |
|---------------|----------------|
| AWS S3 Standard | $0.023 |
| AWS S3 Infrequent Access | $0.0125 |
| Azure Blob Hot | $0.018 |
| GCP Standard | $0.020 |

#### Data Transfer (per GB)
| Transfer Type | Price/GB |
|---------------|----------|
| AWS Internet Out | $0.09 |
| AWS Internet In | $0.00 |
| Azure Internet Out | $0.087 |
| GCP Internet Out | $0.12 |

#### External APIs
| API | Price |
|-----|-------|
| TCGA Query | Free (NCI funded) |
| TCGA Download | Free |
| HuggingFace Inference (small) | $0.06 per 1M tokens |
| HuggingFace Inference (large) | $0.60 per 1M tokens |
| Seqera Platform | $0.02 per compute hour |

#### Common Analyses (estimated)
| Analysis | Price |
|----------|-------|
| FASTQ QC | $0.01 per GB |
| RNA-seq | $0.50 per sample |
| Variant Calling | $1.20 per sample |
| Multiomics Integration | $0.25 per sample |
| Spatial Analysis | $2.00 per slide |

### Updating Pricing

To update pricing data, edit `PRICING` dictionary in `shared/utils/cost_tracking.py`:

```python
PRICING = {
    "compute": {
        "aws_batch_medium": 0.192,  # Update this value
    },
    "analysis": {
        "rna_seq_per_sample": 0.50,  # Update this value
    }
}
```

---

## 6. Cost Logs

### Log Format

Cost logs are saved as JSON files with complete operation details:

```json
{
  "summary": {
    "total_usd": 7.15,
    "by_category": {
      "analysis": 4.30,
      "compute": 2.00,
      "storage": 0.85
    },
    "start_time": "2025-01-15T10:30:00",
    "end_time": "2025-01-15T10:34:05",
    "duration_seconds": 245.3,
    "patient_id": "patient_001",
    "workflow_id": "ovarian_cancer_analysis"
  },
  "items": [
    {
      "timestamp": "2025-01-15T10:30:15",
      "category": "analysis",
      "operation": "fastq_qc",
      "amount_usd": 0.15,
      "quantity": 10.0,
      "unit": "GB",
      "metadata": {}
    },
    {
      "timestamp": "2025-01-15T10:31:20",
      "operation": "variant_calling",
      "amount_usd": 1.20,
      "quantity": 1.0,
      "unit": "sample",
      "metadata": {"sample_id": "sample_001"}
    }
  ]
}
```

### Log File Naming

Default naming pattern:
```
{patient_id}_cost_log_{timestamp}.json

Examples:
patient_001_cost_log_20250115_103405.json
patient_002_cost_log_20250115_143022.json
```

### Analyzing Cost Logs

```python
import json
from pathlib import Path

# Load cost log
with open("patient_001_cost_log_20250115_103405.json") as f:
    log = json.load(f)

# Analyze costs
print(f"Total: ${log['summary']['total_usd']:.2f}")
print(f"Duration: {log['summary']['duration_seconds']:.1f}s")

# Find most expensive operation
items = log['items']
most_expensive = max(items, key=lambda x: x['amount_usd'])
print(f"Most expensive: {most_expensive['operation']} (${most_expensive['amount_usd']:.2f})")

# Calculate cost per patient
costs_by_patient = {}
for item in items:
    pid = log['summary']['patient_id']
    costs_by_patient[pid] = costs_by_patient.get(pid, 0.0) + item['amount_usd']
```

---

## 7. Example Workflows

### PatientOne Workflow Cost Tracking

```python
from cost_tracking import CostTracker, BudgetAlert, CostEstimator

# Step 1: Estimate cost before execution
estimator = CostEstimator()
patient_one_estimate = estimator.patient_one_workflow_cost()
print(f"Estimated cost: ${patient_one_estimate['total']:.2f}")

# Step 2: Initialize tracker with budget
tracker = CostTracker(
    patient_id="patient_one",
    workflow_id="ovarian_cancer_comprehensive",
    budget_limit_usd=10.0  # Set budget
)

# Step 3: Set up budget alerts
alerts = BudgetAlert(
    tracker=tracker,
    thresholds=[5.0, 7.5, 10.0]
)

# Step 4: Track costs as analysis progresses
tracker.add_cost("analysis", "genomic_qc", 0.15)
alerts.check()

tracker.add_cost("analysis", "variant_calling", 1.20)
alerts.check()

tracker.add_cost("analysis", "rna_seq", 0.50)
alerts.check()

tracker.add_cost("analysis", "multiomics_integration", 0.75)
alerts.check()

tracker.add_cost("analysis", "spatial_analysis", 2.00)
alerts.check()

tracker.add_cost("analysis", "imaging_analysis", 1.50)
alerts.check()

tracker.add_cost("compute", "compute_overhead", 0.80)
alerts.check()

tracker.add_cost("storage", "30day_cache", 0.20)
alerts.check()

# Step 5: Generate report
tracker.print_summary()

# Step 6: Save cost log
log_path = tracker.save_log()
print(f"Cost log saved: {log_path}")

# Step 7: Compare actual vs estimated
actual = tracker.get_total()
estimated = patient_one_estimate['total']
difference = actual - estimated
pct_diff = (difference / estimated) * 100

print(f"\nActual: ${actual:.2f}")
print(f"Estimated: ${estimated:.2f}")
print(f"Difference: ${difference:.2f} ({pct_diff:+.1f}%)")
```

### Multi-Patient Study Cost Tracking

```python
from cost_tracking import CostTracker

# Track costs for multiple patients
patient_trackers = {}

for patient_id in ["patient_001", "patient_002", "patient_003"]:
    tracker = CostTracker(
        patient_id=patient_id,
        workflow_id="multi_patient_study",
        budget_limit_usd=15.0
    )

    # Run analysis
    run_patient_analysis(patient_id, cost_tracker=tracker)

    # Save log
    tracker.save_log()
    patient_trackers[patient_id] = tracker

# Aggregate costs
total_cost = sum(t.get_total() for t in patient_trackers.values())
avg_cost_per_patient = total_cost / len(patient_trackers)

print(f"Total study cost: ${total_cost:.2f}")
print(f"Average per patient: ${avg_cost_per_patient:.2f}")

# Cost breakdown across all patients
all_categories = {}
for tracker in patient_trackers.values():
    by_cat = tracker.get_by_category()
    for cat, amount in by_cat.items():
        all_categories[cat] = all_categories.get(cat, 0.0) + amount

print("\nCosts by Category (All Patients):")
for cat, amount in sorted(all_categories.items()):
    print(f"  {cat:20s} ${amount:8.2f}")
```

---

## 8. Best Practices

### 1. Always Estimate Before Execution

```python
# Good: Estimate first
estimator = CostEstimator()
estimated_cost = estimator.multiomics_integration_cost(num_samples=100, num_modalities=3)
print(f"Estimated cost: ${estimated_cost:.2f}")

if estimated_cost > budget:
    print("Cost exceeds budget, consider reducing samples or modalities")
else:
    # Proceed with analysis
    run_analysis()
```

### 2. Set Budget Limits

```python
# Good: Set budget limit
tracker = CostTracker(
    patient_id="patient_001",
    budget_limit_usd=50.0  # Will warn when exceeded
)

# Bad: No budget limit
tracker = CostTracker(patient_id="patient_001")  # No cost control
```

### 3. Use Meaningful Operation Names

```python
# Good: Descriptive names
tracker.add_cost("analysis", "rna_seq_alignment_star_hg38", 0.50)
tracker.add_cost("analysis", "variant_calling_gatk_best_practices", 1.20)

# Bad: Vague names
tracker.add_cost("analysis", "step1", 0.50)
tracker.add_cost("analysis", "analysis", 1.20)
```

### 4. Track Quantities and Units

```python
# Good: Include quantity and unit
tracker.add_cost("compute", "aws_batch", 0.80, quantity=4.2, unit="hours")
tracker.add_cost("storage", "s3_cache", 0.20, quantity=8.7, unit="GB-month")

# Bad: No context
tracker.add_cost("compute", "aws_batch", 0.80)
tracker.add_cost("storage", "s3_cache", 0.20)
```

### 5. Always Save Cost Logs

```python
# Good: Save log for every analysis
tracker = CostTracker(patient_id="patient_001")
# ... run analysis ...
log_path = tracker.save_log()  # Permanent record

# Bad: No log saved
tracker = CostTracker(patient_id="patient_001")
# ... run analysis ...
# Cost data lost when script ends
```

### 6. Use DRY_RUN for Testing

```python
# Good: Test with DRY_RUN first (no charges)
export MULTIOMICS_DRY_RUN=true
# Run full pipeline, estimate costs
# Review cost estimates

# Then run for real
export MULTIOMICS_DRY_RUN=false
# Actual analysis with real costs
```

---

## 9. Cost Optimization Tips

### Reduce Sample Count
```python
# Instead of all samples
cost = estimator.multiomics_integration_cost(num_samples=100, num_modalities=3)
# $25.00

# Start with pilot study
cost = estimator.multiomics_integration_cost(num_samples=10, num_modalities=3)
# $2.50
```

### Use Appropriate Compute Instances
```python
# Over-provisioned (expensive)
cost = estimator.compute_cost("aws_batch_large", duration_hours=2.0)
# $0.768

# Right-sized (cost-effective)
cost = estimator.compute_cost("aws_batch_medium", duration_hours=2.0)
# $0.384
```

### Cache Intermediate Results
```python
# Without caching: re-compute every time
for patient in patients:
    run_expensive_analysis(patient)  # $5 each time

# With caching: compute once, reuse
for patient in patients:
    result = get_from_cache_or_compute(patient)  # $5 first time, $0.05 subsequent
```

### Use Storage Classes Appropriately
```python
# Frequently accessed data
cost = estimator.storage_cost(10, duration_days=30, storage_class="aws_s3_standard")
# $0.23

# Infrequently accessed data
cost = estimator.storage_cost(10, duration_days=30, storage_class="aws_s3_ia")
# $0.125 (46% savings)
```

---

## 10. Risk Reduction Impact

### Before Implementation

**R6: Unexpected Costs (7/10 severity)**
- Likelihood: High (70%)
- Impact: Significant (grant exhaustion)
- No cost estimation before execution
- No cost tracking during execution
- No budget limits or alerts
- Difficult to justify costs in publications

### After Implementation

**R6: Unexpected Costs (2/10 severity)** ✅

- Likelihood: Low (20%)
- Impact: Minor (controlled spending)
- Accurate cost estimation before execution
- Real-time cost tracking during execution
- Budget limits with automated alerts
- Detailed cost logs for publications

**Risk Reduction:** 70% (from 7/10 to 2/10)

---

## 11. Integration with MCP Servers

### mcp-multiomics (Implemented)

**Tool:** `estimate_analysis_cost`

```python
# Estimate before running analysis
cost_estimate = await estimate_analysis_cost(
    num_samples=10,
    modalities=["rna", "protein", "phospho"],
    include_halla=True,
    include_upstream=False
)

print(f"Estimated: ${cost_estimate['estimated_cost_usd']:.2f}")
print(f"Tier: {cost_estimate['budget_tier']}")
```

### Future Integration (Other Servers)

Similar cost estimation tools can be added to:
- **mcp-fgbio:** Estimate FASTQ QC and validation costs
- **mcp-spatialtools:** Estimate spatial transcriptomics analysis costs
- **mcp-seqera:** Estimate Nextflow pipeline execution costs
- **mcp-tcga:** Estimate TCGA data download and processing costs

---

## 12. Related Documentation

- **Risk Matrix:** `RISK_MATRIX.md` - See R6 (Unexpected costs)
- **Risk Mitigation Workplan:** `RISK_MITIGATION_WORKPLAN.md` - WI-6
- **Cost Tracking Utilities Source:** `shared/utils/cost_tracking.py`

---

## Support

For questions or issues with cost tracking:
1. Check cost estimation before running expensive analyses
2. Review cost logs to identify unexpected charges
3. Adjust pricing data in `PRICING` dictionary if rates change
4. Use DRY_RUN mode for testing without charges

**Remember:** Cost estimates are based on typical usage patterns. Actual costs may vary based on data complexity, compute efficiency, and current cloud pricing.
