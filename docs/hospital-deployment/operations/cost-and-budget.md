# Cost and Budget Management

Comprehensive guide to understanding, estimating, tracking, and controlling costs for the Precision Medicine MCP system.

## Table of Contents

1. [Cost Overview](#cost-overview)
2. [Cost Analysis by Mode](#cost-analysis-by-mode)
3. [Cost Estimation Tools](#cost-estimation-tools)
4. [Cost Tracking During Execution](#cost-tracking-during-execution)
5. [Budget Alerts and Monitoring](#budget-alerts-and-monitoring)
6. [Infrastructure Requirements](#infrastructure-requirements)
7. [Cost Optimization Strategies](#cost-optimization-strategies)
8. [Use Case Examples](#use-case-examples)
9. [Best Practices](#best-practices)

---

## Cost Overview

### Why Cost Management Matters

**Financial Transparency:**
- Research budgets are limited
- Cloud costs can escalate quickly
- Grant proposals require accurate cost estimates
- Multi-patient studies need per-patient cost tracking

**Cost Control:**
- Set budget limits before analysis
- Get alerts when approaching limits
- Identify expensive operations
- Optimize workflows for cost-efficiency

**Reproducibility:**
- Document exact costs for published research
- Enable cost comparison across methods
- Support funding applications
- Track costs over time

### Cost Summary by Mode

| Mode | Total Time | Total Cost | Data Size | Best For |
|------|-----------|------------|-----------|----------|
| **DRY_RUN** | 25-35 min | ~$1 | ~5 MB | Demo, learning, CI/CD testing |
| **Automated Report** | ~12 seconds | ~$1 | Pre-aligned data | Quick analysis with pre-aligned data |
| **Real Data (Small Files)** | 1-2 hours | $7-19 | ~5 MB | Workflow testing with small synthetic data |
| **Real Data (Small Files + STAR)** | 1.5-3 hours | $12-29 | ~5 MB | Testing from raw FASTQ (small files) |
| **Production (Pre-aligned)** | 2-4 hours | **$25-75** | 3-8 GB | Production analysis with Space Ranger output |
| **Production (Raw FASTQ)** | 4-8 hours | **$50-120** | 12-35 GB | Production analysis from raw sequencing data |

**Key differences across modes:**
- **Spatial data:** 315 KB demo → 100-500 MB production (300-1500× larger)
- **Multi-omics data:** 505 KB demo → 2.7 GB raw or 15-20 MB processed (400-500× larger)
- **Processing time:** 25 min (DRY_RUN) → 4-8 hours (production)
- **Compute cost:** ~$1 (DRY_RUN) → $50-120 (production)

---

## Cost Analysis by Mode

### 1. DRY_RUN Mode (Demonstration - $1)

**Overview:**
- Uses synthetic responses from all MCP servers
- No external API calls or computational processing
- Ideal for demonstration and workflow validation

**Per-Test Breakdown:**

| Test | Servers Used | Est. Tokens (In/Out) | Time | Cost |
|------|--------------|---------------------|------|------|
| **TEST_1: Clinical + Genomic** | MockEpic, FGbio, TCGA | 2,000 / 3,500 | 3-5 min | ~$1 |
| **TEST_2: Multi-Omics** | MultiOmics | 2,500 / 4,000 | 5-8 min | ~$1 |
| **TEST_3: Spatial** | SpatialTools, DeepCell | 2,000 / 3,500 | 4-6 min | ~$1 |
| **TEST_4: Imaging** | OpenImageData, DeepCell | 1,800 / 3,000 | 3-5 min | ~$1 |
| **TEST_5: Integration** | All 9 servers | 3,000 / 5,000 | 5-7 min | ~$1 |
| **TOTAL** | - | **11,300 / 19,000** | **25-35 min** | **~$1** |

**Cost Calculation (Claude Sonnet 4 pricing):**
- Input: 11,300 tokens × $3/M = ~$0.03
- Output: 19,000 tokens × $15/M = ~$0.29
- **Total: ~$0.32** → rounded to **~$1** for simplicity

### 2. Automated Patient Report Generator ($1, 12 seconds)

**Overview:**
- Standalone Python script for rapid patient analysis
- Requires pre-aligned spatial transcriptomics data
- Integrates FHIR clinical data from GCP Healthcare API

**Script:** `tools/reports/generate_patient_report.py`

**Performance:**
- Runtime: ~12 seconds per patient
- Compute cost: ~$0.50 (local CPU processing)
- API cost: ~$0.50 (FHIR data retrieval)
- **Total: ~$1 per patient**

**Output Files (10 total, ~3.4 MB):**
- 5 data files (CSV/JSON/TXT)
- 5 visualizations (PNG, 300 DPI)

**Best for:**
- Quick patient analysis with pre-aligned data
- Batch processing of multiple patients
- Automated report generation pipelines
- Research studies requiring standardized analysis

### 3. Real Patient Data Mode - Small Synthetic Files ($7-29)

**Overview:**
- Processes actual patient genomic, transcriptomic, and imaging data
- Executes computational workflows (alignment, segmentation, etc.)
- **Data size:** Small synthetic files (~4.9 MB total)

**Per-Test Breakdown:**

| Test | Processing Time (Pre-aligned) | Processing Time (with STAR) | Total Cost (Pre-aligned) | Total Cost (with STAR) |
|------|-------------------------------|------------------------------|--------------------------|------------------------|
| **TEST_1: Clinical + Genomic** | 10-15 min | 10-15 min | ~$1 | ~$1 |
| **TEST_2: Multi-Omics** | 15-25 min | 15-25 min | $2-4 | $2-4 |
| **TEST_3: Spatial** | 10-20 min | 40-80 min | $1-3 | $6-13 |
| **TEST_4: Imaging** | 20-40 min | 20-40 min | $3-7 | $3-7 |
| **TEST_5: Integration** | 5-10 min | 5-10 min | ~$1 | ~$1 |
| **TOTAL** | **1-2 hours** | **1.5-3 hours** | **$7-19** | **$12-29** |

### 4. Production Patient Data Mode ($25-120)

**Overview:**
- Processes realistic hospital data volumes
- **Data size per patient:** ~3-8 GB total (pre-aligned) or 12-35 GB (raw FASTQ)
- Higher memory requirements (64-128 GB RAM)
- Larger Cloud Run instances needed

**Per-Test Breakdown:**

| Test | Processing Time (Pre-aligned) | Processing Time (Raw FASTQ) | Per-Test Cost (Pre-aligned) | Per-Test Cost (Raw FASTQ) |
|------|-------------------------------|------------------------------|----------------------------|--------------------------|
| **TEST_1: Clinical + Genomic** | 15-30 min | 15-30 min | ~$1-2 | ~$1-2 |
| **TEST_2: Multi-Omics** | 30-60 min | 30-60 min | $6-20 | $6-20 |
| **TEST_3: Spatial** | 45-120 min | 90-240 min | $5-30 | $10-40 |
| **TEST_4: Imaging** | 40-90 min | 40-90 min | $10-36 | $10-36 |
| **TEST_5: Integration** | 10-20 min | 10-20 min | ~$1-2 | ~$1-2 |
| **TOTAL** | **2-4 hours** | **4-8 hours** | **$24-92** | **$29-102** |

**Key differences from small files:**
- **Data volumes:** 300-1500× larger spatial, 400-500× larger multi-omics
- **Processing time:** 2-6× longer (optimized algorithms prevent 300× slowdown)
- **Compute costs:** 3-5× higher due to longer processing, larger memory
- **Token costs:** Only ~2× higher despite massive data increase (servers return summaries!)
- **Total cost:** $7-29 → **$24-102** (3-4× more expensive for realistic hospital data)

---

## Cost Estimation Tools

### Before You Run: Estimate Costs

Always estimate costs before running expensive analyses to help with budgeting and planning.

### Using CostEstimator

```python
from cost_tracking import CostEstimator

estimator = CostEstimator()

# Estimate RNA-seq analysis
rna_cost = estimator.rna_seq_cost(num_samples=5)
print(f"RNA-seq for 5 samples: ${rna_cost:.2f}")

# Estimate variant calling
variant_cost = estimator.variant_calling_cost(num_samples=3)
print(f"Variant calling for 3 samples: ${variant_cost:.2f}")

# Estimate multi-omics integration
multiomics_cost = estimator.multiomics_integration_cost(
    num_samples=10,
    num_modalities=3  # RNA + Protein + Phospho
)
print(f"Multi-omics integration: ${multiomics_cost:.2f}")

# Estimate PatientOne full workflow
patient_one_costs = estimator.patient_one_workflow_cost()
print(f"PatientOne total: ${patient_one_costs['total']:.2f}")
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
```

**Output:**
```json
{
  "estimated_cost_usd": 18.50,
  "breakdown": {
    "integration": 2.50,
    "preprocessing": 0.50,
    "halla_analysis": 15.00,
    "compute": 0.15,
    "storage_30days": 0.15
  },
  "budget_tier": "medium",
  "recommendation": "Standard multi-omics analysis"
}
```

---

## Cost Tracking During Execution

### Real-Time Cost Monitoring

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

### Automatic Cost Tracking

**Context Manager:**

```python
from cost_tracking import track_cost, CostTracker

tracker = CostTracker(patient_id="patient_002")

with track_cost(tracker, "compute", "rna_seq_pipeline"):
    run_rna_seq_pipeline()
    # Cost automatically tracked based on duration
```

**Decorator:**

```python
from cost_tracking import track_operation_cost

@track_operation_cost(
    category="analysis",
    operation="variant_calling",
    cost_estimator=lambda sample_id: 1.20
)
def run_variant_calling(sample_id: str):
    call_variants(sample_id)
    return results

# Use with cost tracking
result = run_variant_calling("sample_001", cost_tracker=tracker)
```

---

## Budget Alerts and Monitoring

### Real-Time Budget Alerts

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

alerts = BudgetAlert(
    tracker=tracker,
    thresholds=[10.0, 25.0, 50.0],
    alert_callback=custom_alert_handler
)
```

---

## Infrastructure Requirements

### DRY_RUN Mode (Demonstration)
- **CPU:** Any modern laptop (2+ cores)
- **RAM:** 4GB minimum, 8GB recommended
- **Storage:** 5GB for MCP servers + Python dependencies
- **Data size:** ~5 MB total

### Real Patient Data - Small Synthetic Files

**Without STAR alignment (pre-aligned data):**
- **CPU:** 4-8 cores
- **RAM:** 16GB minimum, 32GB recommended
- **Storage:** 20GB minimum
- **Data size:** ~5 MB per patient

**With STAR alignment (from raw FASTQ, small files):**
- **CPU:** 8-16 cores recommended
- **RAM:** 32GB minimum, 64GB recommended
- **Storage:** 100GB minimum
- **Data size:** ~5 MB per patient (synthetic small FASTQ)

### Real Patient Data - Production Volumes

**Without STAR alignment (pre-aligned from Space Ranger):**
- **CPU:** 8-16 cores recommended
- **RAM:** 32GB minimum, **64GB recommended** for production
- **Storage:** **300-800 GB for 100 patients**
- **Data size:** 3-8 GB per patient

**With STAR alignment (from raw FASTQ, production):**
- **CPU:** 16-32 cores recommended
- **RAM:** 64GB minimum, **128GB recommended** for production
- **Storage:** **500 GB - 1.5 TB for 100 patients**
- **Data size:** 12-35 GB per patient

---

## Cost Optimization Strategies

### 1. Batch Processing
- Run multiple patients in sequence to amortize setup costs
- Cache reference genomes and model weights
- **Estimated savings:** 20-30% for 10+ patients

### 2. Pre-computed Data Reuse
- Store intermediate results (aligned BAMs, cell segmentation masks)
- Reuse for hypothesis testing without re-processing
- **Estimated savings:** 50-70% for follow-up analyses

### 3. Selective Analysis
- Run only relevant tests (skip imaging if not available)
- Use DRY_RUN mode for workflow validation
- **Estimated savings:** 40-60% when only 2-3 tests needed

### 4. Use Appropriate Compute Instances

```python
# Over-provisioned (expensive)
cost = estimator.compute_cost("aws_batch_large", duration_hours=2.0)
# $0.768

# Right-sized (cost-effective)
cost = estimator.compute_cost("aws_batch_medium", duration_hours=2.0)
# $0.384 (50% savings)
```

### 5. Cache Intermediate Results

```python
# Without caching: re-compute every time
for patient in patients:
    run_expensive_analysis(patient)  # $5 each time

# With caching: compute once, reuse
for patient in patients:
    result = get_from_cache_or_compute(patient)  # $5 first time, $0.05 subsequent
```

### 6. Use Storage Classes Appropriately

```python
# Frequently accessed data
cost = estimator.storage_cost(10, duration_days=30, storage_class="aws_s3_standard")
# $0.23

# Infrequently accessed data
cost = estimator.storage_cost(10, duration_days=30, storage_class="aws_s3_ia")
# $0.125 (46% savings)
```

---

## Use Case Examples

### Academic Research Lab (Demonstration Data)

**Scenario:** 50 patient analyses per year with small synthetic files

| Mode | Cost per Patient | Annual Cost | Use Case |
|------|-----------------|-------------|----------|
| DRY_RUN | ~$1 | $50 | Workflow development, testing |
| Automated Report | ~$1 | $50 | Quick analysis (pre-aligned data) |
| Real Data (Small Files) | $13 (avg) | $650 | Full analysis with MCP orchestration |

**ROI:** Replaces ~40 hours of manual bioinformatics work per patient ($80/hr × 40 = $3,200)
**Annual Savings:** $159,350 using automated reports

### Academic Research Lab (Production Data)

**Scenario:** 50 patient analyses per year with realistic hospital data volumes

| Mode | Cost per Patient | Annual Cost |
|------|-----------------|-------------|
| Real Data (Pre-aligned) | $50 (avg) | $2,500 |
| Real Data (Raw FASTQ) | $85 (avg) | $4,250 |

**ROI:** Replaces ~40 hours of manual work per patient ($3,200)
**Annual Savings:** $157,500 (pre-aligned) or $155,750 (raw FASTQ)

### Clinical Genomics Center (Production Data)

**Scenario:** 500 patient analyses per year

| Mode | Cost per Patient | Annual Cost |
|------|-----------------|-------------|
| Automated Report | ~$1 | $500 |
| Real Data (Pre-aligned) | $50 (avg) | $25,000 |
| Real Data (Raw FASTQ) | $85 (avg) | $42,500 |

**ROI:** Reduces time-to-result from 2-3 weeks to 4-6 hours
**Annual Savings:** ~$1.6M using automated reports (vs manual $3,200/patient)

### Pharmaceutical R&D (Production Data)

**Scenario:** 200 PDX model analyses per year

| Mode | Cost per Analysis | Annual Cost |
|------|------------------|-------------|
| Automated Report | ~$1 | $200 |
| Real Data (Pre-aligned) | $50 (avg) | $10,000 |
| Real Data (Raw FASTQ) | $85 (avg) | $17,000 |

**ROI:** Accelerates target identification from weeks → hours
**Cost avoidance:** Reduces need for manual bioinformatics staff (~$500K/year)

---

## Best Practices

### 1. Always Estimate Before Execution

```python
# Good: Estimate first
estimator = CostEstimator()
estimated_cost = estimator.multiomics_integration_cost(num_samples=100, num_modalities=3)

if estimated_cost > budget:
    print("Cost exceeds budget, consider reducing samples")
else:
    run_analysis()
```

### 2. Set Budget Limits

```python
# Good: Set budget limit
tracker = CostTracker(
    patient_id="patient_001",
    budget_limit_usd=50.0  # Will warn when exceeded
)
```

### 3. Use Meaningful Operation Names

```python
# Good: Descriptive names
tracker.add_cost("analysis", "rna_seq_alignment_star_hg38", 0.50)

# Bad: Vague names
tracker.add_cost("analysis", "step1", 0.50)
```

### 4. Track Quantities and Units

```python
# Good: Include quantity and unit
tracker.add_cost("compute", "aws_batch", 0.80, quantity=4.2, unit="hours")

# Bad: No context
tracker.add_cost("compute", "aws_batch", 0.80)
```

### 5. Always Save Cost Logs

```python
# Good: Save log for every analysis
tracker = CostTracker(patient_id="patient_001")
# ... run analysis ...
log_path = tracker.save_log()  # Permanent record
```

### 6. Use DRY_RUN for Testing

```bash
# Test with DRY_RUN first (no charges)
export MULTIOMICS_DRY_RUN=true
# Run full pipeline, estimate costs

# Then run for real
export MULTIOMICS_DRY_RUN=false
```

---

## Frequently Asked Questions

### Q: Why is Real Data mode 25-120× more expensive than DRY_RUN?

**A:** Real data requires actual computational processing rather than synthetic responses:
- **Scientific computation** (not just LLM orchestration)
- **Data volume**: Production data is 300-1500× larger
- **Processing time**: 2-8 hours vs 25-35 minutes
- **Memory requirements**: 64-128 GB RAM vs 4-8 GB

### Q: Can I reduce costs by using smaller data files?

**A:** Yes! Costs scale with:
- **Number of spatial spots:** 900 demo → 3,000-5,000 production
- **Number of genes:** 31 demo → 18,000-30,000 production
- **Multi-omics volume:** 505 KB demo → 2.7 GB production
- **Image resolution:** 4 MB demo → 500 MB - 2 GB production

**Demo costs:** ~$1 (DRY_RUN) to $7-29 (small files)
**Production costs:** $25-75 (pre-aligned) to $50-120 (raw FASTQ)

### Q: What happens if a test fails midway?

**A:**
- DRY_RUN: No cost beyond Claude tokens consumed
- Real Data: You pay for compute time used, but intermediate results are cached for retry

### Q: Are there free tiers for external APIs?

**A:**
- TCGA: Always free (NIH public data)
- HuggingFace: Free tier with rate limits (30 requests/hour)
- Seqera: 10 compute-hours/month free

---

## Related Documentation

- [Executive Summary](../EXECUTIVE_SUMMARY.md) - High-level cost overview
- [Risk Assessment](../compliance/risk-assessment.md) - R6 (Unexpected costs)
- **Cost Utilities Source:** `shared/utils/cost_tracking.py`

---

**Last Updated:** 2026-01-13

**Pricing basis:** Claude Sonnet 4 ($3/M input, $15/M output), AWS EC2, Google Cloud Run

**Status:** ✅ Implemented (WI-6 from Risk Mitigation Workplan)
**Risk Reduced:** R6 (Unexpected costs) - 70% reduction (7/10 → 2/10)
