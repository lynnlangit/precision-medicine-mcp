# Risk Mitigation Work Plan - Next 1-2 Weeks

## Executive Summary

**Timeline:** 1-2 weeks (solo developer)
**Priority:** Full clinical integration (all 9 servers)
**Data Type:** Research data (PDX models, not clinical patients)
**Cost Target:** $15-45/patient (acceptable)

**Focus:** Code quality and documentation improvements to reduce critical risks before research deployment

---

## Work Items Prioritized by Risk Reduction

### 🔴 CRITICAL PRIORITY (Week 1 - Must Complete)

These items address the highest-severity risks and should be completed before any pilot deployment.

---

#### **WI-1: Server Implementation Status Documentation**
**Addresses:** R3 (Mocked Servers in Production) - Severity 9/10
**Effort:** 4-6 hours
**Impact:** 🔴 CRITICAL - Prevents catastrophic deployment of mocked servers

**Problem:**
- 5 servers are 0-30% real implementation but not clearly documented
- Users/developers cannot easily tell what's mocked vs real
- High risk of using synthetic results for research decisions

**Solution - Create Implementation Status Matrix:**

**File to create:** `docs/SERVER_IMPLEMENTATION_STATUS.md`

```markdown
# MCP Server Implementation Status Matrix

**Last Updated:** [DATE]
**Purpose:** Clearly document what's real vs mocked in each server to prevent accidental use of synthetic data

## Quick Reference

| Server | Real % | Status | Safe for Production? | Notes |
|--------|--------|--------|---------------------|-------|
| mcp-multiomics | 85% | ✅ Production Ready | YES | Preprocessing, Stouffer's, HAllA, upstream regulators all real |
| mcp-fgbio | 65% | ✅ Production Ready | YES | FASTQ validation, UMI extraction real |
| mcp-spatialtools | 40% | ⚠️ Partial | CONDITIONAL | Quality filtering real, STAR alignment mocked |
| mcp-openimagedata | 30% | ❌ Not Ready | NO | Basic PIL operations only, registration mocked |
| mcp-tcga | 0% | ❌ Fully Mocked | NO | **ALL RESULTS SYNTHETIC** |
| mcp-deepcell | 0% | ❌ Fully Mocked | NO | **ALL RESULTS SYNTHETIC** |
| mcp-huggingface | 0% | ❌ Fully Mocked | NO | **ALL RESULTS SYNTHETIC** |
| mcp-seqera | 0% | ❌ Fully Mocked | NO | **ALL RESULTS SYNTHETIC** |
| mcp-mockepic | 0% | ✅ Intentional Mock | N/A | Mock EHR by design |

## Detailed Implementation Status

### mcp-multiomics (85% Real) ✅

**Real Capabilities:**
- ✅ `validate_multiomics_data` - Full pandas/numpy validation
- ✅ `preprocess_multiomics_data` - Real batch correction (limma), imputation (KNN)
- ✅ `visualize_data_quality` - Real matplotlib PCA plots
- ✅ `calculate_stouffer_meta` - 100% real statistical implementation
- ✅ `predict_upstream_regulators` - Real kinase/TF prediction with databases
- ✅ `run_halla_analysis` - Real HAllA with chunking
- ✅ `integrate_omics_data` - Real pandas merging
- ✅ `create_multiomics_heatmap` - Real seaborn heatmaps
- ✅ `run_multiomics_pca` - Real PCA computation

**Mocked Capabilities:**
- 🔶 R integration for advanced HAllA features (uses Python fallback)

**DRY_RUN Behavior:**
- When `MULTIOMICS_DRY_RUN=true`: Returns synthetic results immediately
- When `MULTIOMICS_DRY_RUN=false`: Executes real computations (tested with 91 automated tests)

**Production Readiness:** ✅ YES - Extensively tested, 68% code coverage

---

### mcp-tcga (0% Real) ❌

**Status:** **FULLY MOCKED - DO NOT USE FOR RESEARCH DECISIONS**

**All Tools Return Synthetic Data:**
- ❌ `query_tcga_cohorts` - Returns hardcoded cohort list
- ❌ `fetch_expression_data` - Returns random numbers
- ❌ `compare_to_cohort` - Returns synthetic z-scores
- ❌ `get_survival_data` - Returns fake survival curves
- ❌ `get_mutation_data` - Returns mock mutation frequencies

**Why Mocked:**
- TCGA GDC API client not yet implemented
- Requires authentication setup
- Estimated 1 week to implement real API integration

**Production Readiness:** ❌ NO - Not safe for any research use

**Recommended Action:** Either implement real TCGA API or remove from production config
```

**Acceptance Criteria:**
- [ ] Document created with all 9 servers
- [ ] Each server has % real implementation clearly stated
- [ ] "Safe for Production?" column with YES/NO/CONDITIONAL
- [ ] Link from main README to this document
- [ ] Add prominent warning at top of README: "Check SERVER_IMPLEMENTATION_STATUS.md before production use"

---

#### **WI-2: Runtime Warnings for DRY_RUN Mode**
**Addresses:** R3 (Mocked Servers in Production) - Severity 9/10
**Effort:** 3-4 hours
**Impact:** 🔴 CRITICAL - Prevents accidental use of synthetic data

**Problem:**
- DRY_RUN mode returns synthetic data but doesn't clearly warn users
- Developers might forget to set `DRY_RUN=false` in production
- No visual indication in Claude Desktop that results are mocked

**Solution - Add Warnings to All Server Responses:**

**Code changes needed in each server's `__main__.py`:**

```python
# Example: servers/mcp-multiomics/src/mcp_multiomics/__main__.py

import os
import logging

# Add at module level
DRY_RUN = os.getenv("MULTIOMICS_DRY_RUN", "true").lower() == "true"

# Add startup warning
if DRY_RUN:
    logging.warning("=" * 80)
    logging.warning("⚠️  DRY_RUN MODE ENABLED - RETURNING SYNTHETIC DATA")
    logging.warning("⚠️  Set MULTIOMICS_DRY_RUN=false for real analysis")
    logging.warning("=" * 80)

# Modify each tool response to include warning
@app.tool()
async def integrate_omics_data(rna_file, protein_file, phospho_file):
    """Integrate multi-omics data."""

    if DRY_RUN:
        warning = """
⚠️  SYNTHETIC DATA WARNING ⚠️
This result was generated in DRY_RUN mode and does NOT represent real analysis.
Set MULTIOMICS_DRY_RUN=false for production use.
"""
        return warning + "\n" + generate_synthetic_result()
    else:
        # Real implementation
        return perform_real_integration(rna_file, protein_file, phospho_file)
```

**Files to modify (9 total):**
1. `servers/mcp-multiomics/src/mcp_multiomics/__main__.py`
2. `servers/mcp-fgbio/src/mcp_fgbio/__main__.py`
3. `servers/mcp-spatialtools/src/mcp_spatialtools/__main__.py`
4. `servers/mcp-openimagedata/src/mcp_openimagedata/__main__.py`
5. `servers/mcp-tcga/src/mcp_tcga/__main__.py`
6. `servers/mcp-deepcell/src/mcp_deepcell/__main__.py`
7. `servers/mcp-huggingface/src/mcp_huggingface/__main__.py`
8. `servers/mcp-seqera/src/mcp_seqera/__main__.py`
9. `servers/mcp-mockepic/src/mcp_mockepic/__main__.py`

**Acceptance Criteria:**
- [ ] All 9 servers log warning on startup if DRY_RUN=true
- [ ] All tool responses include "SYNTHETIC DATA WARNING" when DRY_RUN=true
- [ ] Warning is visible in Claude Desktop output
- [ ] Test with DRY_RUN=true and DRY_RUN=false to verify

**Estimated Time:** 3-4 hours (30 min per server)

---

#### **WI-3: Input Validation & Error Messages**
**Addresses:** R5 (Poor Data Quality) - Severity 5/10
**Effort:** 6-8 hours
**Impact:** 🔴 HIGH - Prevents pipeline failures from malformed data

**Problem:**
- Servers assume well-formed input files (like synthetic test data)
- Real research data may have format variations, missing fields, encoding issues
- Error messages like "FileNotFoundError" don't help users fix issues

**Solution - Add Input Validation Layer:**

**Create validation utilities:**

**File to create:** `servers/mcp-multiomics/src/mcp_multiomics/validation.py`

```python
"""Input validation utilities for multi-omics data."""

import pandas as pd
from pathlib import Path
from typing import Tuple, List

class ValidationError(Exception):
    """Raised when input data fails validation."""
    pass

def validate_multiomics_file(
    file_path: str,
    required_columns: List[str],
    file_type: str = "multi-omics"
) -> Tuple[bool, List[str]]:
    """
    Validate multi-omics input file format.

    Returns:
        (is_valid, error_messages)
    """
    errors = []

    # Check file exists
    if not Path(file_path).exists():
        errors.append(f"❌ File not found: {file_path}")
        errors.append(f"💡 Make sure the file path is absolute, not relative")
        return False, errors

    # Check file not empty
    if Path(file_path).stat().st_size == 0:
        errors.append(f"❌ File is empty: {file_path}")
        return False, errors

    # Try to read as TSV/CSV
    try:
        df = pd.read_csv(file_path, sep='\t', nrows=5)
    except Exception as e:
        errors.append(f"❌ Cannot parse file as tab-separated: {file_path}")
        errors.append(f"💡 Error: {str(e)}")
        errors.append(f"💡 Try: Excel → Save As → Tab delimited text (.txt)")
        return False, errors

    # Check required columns exist
    missing_cols = set(required_columns) - set(df.columns)
    if missing_cols:
        errors.append(f"❌ Missing required columns: {missing_cols}")
        errors.append(f"💡 Found columns: {list(df.columns)[:10]}")
        errors.append(f"💡 Expected columns must include: {required_columns}")
        return False, errors

    # Check for common issues
    if df.isnull().all().any():
        empty_cols = df.columns[df.isnull().all()].tolist()
        errors.append(f"⚠️  Warning: Columns with all missing values: {empty_cols}")

    # Success!
    return True, []

def validate_fastq_file(file_path: str) -> Tuple[bool, List[str]]:
    """Validate FASTQ file format."""
    errors = []

    if not Path(file_path).exists():
        errors.append(f"❌ FASTQ file not found: {file_path}")
        return False, errors

    # Check gzip compression
    is_gzipped = file_path.endswith('.gz')

    # Read first few lines
    try:
        if is_gzipped:
            import gzip
            with gzip.open(file_path, 'rt') as f:
                lines = [f.readline() for _ in range(8)]
        else:
            with open(file_path, 'r') as f:
                lines = [f.readline() for _ in range(8)]
    except Exception as e:
        errors.append(f"❌ Cannot read FASTQ file: {str(e)}")
        return False, errors

    # Validate FASTQ format (4 lines per read)
    if not lines[0].startswith('@'):
        errors.append(f"❌ Invalid FASTQ format: First line should start with '@'")
        errors.append(f"💡 Found: {lines[0][:50]}")
        return False, errors

    if not lines[2].startswith('+'):
        errors.append(f"❌ Invalid FASTQ format: Third line should start with '+'")
        return False, errors

    # Check quality score encoding
    qual = lines[3].strip()
    if qual:
        qual_min = ord(min(qual))
        qual_max = ord(max(qual))

        # Phred+33 (standard): 33-126
        # Phred+64 (old Illumina): 64-126
        if qual_min < 33:
            errors.append(f"⚠️  Warning: Unusual quality score encoding detected")
            errors.append(f"💡 Min quality character: ASCII {qual_min}")

        if qual_max > 126:
            errors.append(f"❌ Invalid quality scores (ASCII > 126)")
            return False, errors

    return True, []

# Add similar validators for VCF, spatial data, images...
```

**Then use in tools:**

```python
@app.tool()
async def integrate_omics_data(rna_file: str, protein_file: str, phospho_file: str):
    """Integrate RNA, protein, phospho data."""

    # Validate inputs BEFORE processing
    from mcp_multiomics.validation import validate_multiomics_file, ValidationError

    required_cols = ["gene_id", "log2FC", "pvalue"]

    # Validate RNA file
    is_valid, errors = validate_multiomics_file(rna_file, required_cols, "RNA")
    if not is_valid:
        error_msg = "RNA File Validation Failed:\n" + "\n".join(errors)
        raise ValidationError(error_msg)

    # Validate protein file
    is_valid, errors = validate_multiomics_file(protein_file, required_cols, "Protein")
    if not is_valid:
        error_msg = "Protein File Validation Failed:\n" + "\n".join(errors)
        raise ValidationError(error_msg)

    # ... rest of implementation
```

**Files to create/modify:**
1. Create `servers/mcp-multiomics/src/mcp_multiomics/validation.py`
2. Create `servers/mcp-fgbio/src/mcp_fgbio/validation.py`
3. Create `servers/mcp-spatialtools/src/mcp_spatialtools/validation.py`
4. Modify tool implementations to use validation (3-4 tools per server)

**Acceptance Criteria:**
- [ ] Validation utilities created for mcp-multiomics, mcp-fgbio, mcp-spatialtools
- [ ] All tools validate inputs before processing
- [ ] Error messages include:
  - ❌ Clear statement of what's wrong
  - 💡 Actionable suggestion to fix it
  - Example of expected format
- [ ] Test with intentionally malformed files (missing columns, wrong encoding, empty files)

**Estimated Time:** 6-8 hours

---

#### **WI-4: Research Use Disclaimers**
**Addresses:** R7 (Incorrect Clinical Recommendations) - Severity 6/10
**Effort:** 2-3 hours
**Impact:** 🔴 HIGH - Critical for patient safety and liability

**Problem:**
- Outputs could be misinterpreted as clinical recommendations
- No disclaimers about research-only use
- No guidance on human review requirements

**Solution - Add Disclaimers to All Outputs:**

**Create disclaimer templates:**

**File to create:** `docs/DISCLAIMERS.md`

```markdown
# Output Disclaimers & Safety Warnings

## Research Use Only Disclaimer

**Add to ALL analysis outputs (especially mcp-multiomics, PatientOne integration):**

```
⚠️  RESEARCH USE ONLY - NOT FOR CLINICAL DECISION-MAKING ⚠️

This analysis is provided for RESEARCH PURPOSES ONLY and has not been
clinically validated. Results should NOT be used for:
- Patient diagnosis
- Treatment selection
- Clinical decision-making

All findings must be:
1. Reviewed by qualified bioinformatics personnel
2. Validated with orthogonal methods (e.g., qPCR, IHC)
3. Interpreted by board-certified oncologist before any clinical use

The developers assume NO liability for clinical decisions based on this output.
```

## Uncertainty Quantification

**For statistical results (p-values, fold changes):**

```
📊 Statistical Confidence Indicators:

- ⭐⭐⭐ HIGH CONFIDENCE: p < 0.001, concordant across modalities
- ⭐⭐ MODERATE CONFIDENCE: p < 0.05, some modality disagreement
- ⭐ LOW CONFIDENCE: p < 0.1, exploratory finding only

Always report:
- p-values (raw and FDR-corrected)
- Effect sizes (log2FC, not just significance)
- Sample size and statistical power
- Missing data percentage
```

## Data Quality Warnings

**When data quality is suboptimal:**

```
⚠️  DATA QUALITY WARNING:
- Missing values: 15% (imputed using KNN)
- Batch effects detected: Corrected using limma
- Low sequencing depth: <10M reads in 3/15 samples

Interpretation should account for these limitations.
Consider validating top findings with higher-quality samples.
```
```

**Then add to tool outputs:**

```python
@app.tool()
async def predict_upstream_regulators(diff_expr_file: str):
    """Predict kinases, TFs, drug targets."""

    # ... perform analysis ...

    # Add disclaimer to output
    disclaimer = """
⚠️  RESEARCH USE ONLY - NOT FOR CLINICAL DECISION-MAKING ⚠️
This upstream regulator prediction is for research hypothesis generation.
Not validated for clinical use. Requires experimental validation.
"""

    result = perform_prediction(diff_expr_file)

    return disclaimer + "\n\n" + result
```

**Files to modify:**
1. `servers/mcp-multiomics/src/mcp_multiomics/__main__.py` - Add disclaimers to all 9 tools
2. `tests/manual_testing/PatientOne-OvarianCancer/README.md` - Add prominent disclaimer
3. `architecture/patient-one/README.md` - Add disclaimer section

**Acceptance Criteria:**
- [ ] `docs/DISCLAIMERS.md` created with templates
- [ ] All mcp-multiomics tools output research use disclaimer
- [ ] PatientOne documentation includes disclaimer prominently
- [ ] Outputs include uncertainty indicators (high/medium/low confidence)

**Estimated Time:** 2-3 hours

---

### 🟡 HIGH PRIORITY (Week 1-2 - Should Complete)

These items significantly improve robustness and user experience.

---

#### **WI-5: Error Handling & Retry Logic**
**Addresses:** R4 (External API Failures) - Severity 6/10
**Effort:** 4-6 hours
**Impact:** 🟡 HIGH - Improves reliability

**Problem:**
- API calls to TCGA, HuggingFace, Seqera have no retry logic
- Transient failures cause entire analysis to fail
- No graceful degradation

**Solution - Add Retry Utilities:**

**File to create:** `shared/utils/api_retry.py`

```python
"""Retry utilities for external API calls."""

import time
import logging
from typing import Callable, Any, Optional
from functools import wraps

logger = logging.getLogger(__name__)

def retry_with_backoff(
    max_retries: int = 3,
    base_delay: float = 1.0,
    max_delay: float = 30.0,
    exceptions: tuple = (Exception,)
):
    """
    Decorator for retrying functions with exponential backoff.

    Args:
        max_retries: Maximum number of retry attempts
        base_delay: Initial delay in seconds
        max_delay: Maximum delay between retries
        exceptions: Tuple of exceptions to catch and retry

    Example:
        @retry_with_backoff(max_retries=3, exceptions=(requests.RequestException,))
        def fetch_tcga_data(cohort_id):
            response = requests.get(f"https://api.gdc.cancer.gov/cases?filters={cohort_id}")
            response.raise_for_status()
            return response.json()
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs) -> Any:
            delay = base_delay
            last_exception = None

            for attempt in range(max_retries + 1):
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    last_exception = e

                    if attempt == max_retries:
                        logger.error(f"❌ {func.__name__} failed after {max_retries} retries: {e}")
                        raise

                    logger.warning(f"⚠️  {func.__name__} attempt {attempt + 1}/{max_retries} failed: {e}")
                    logger.info(f"⏳ Retrying in {delay:.1f}s...")

                    time.sleep(delay)
                    delay = min(delay * 2, max_delay)  # Exponential backoff

            # Should never reach here, but just in case
            raise last_exception

        return wrapper
    return decorator

def optional_api_call(
    fallback_value: Any,
    log_failure: bool = True
):
    """
    Decorator for optional API calls that should degrade gracefully.

    If API call fails, returns fallback value instead of raising exception.

    Example:
        @optional_api_call(fallback_value=None)
        def fetch_tcga_cohort(cancer_type):
            # If TCGA API is down, return None instead of crashing
            return tcga_api.get_cohort(cancer_type)
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs) -> Any:
            try:
                return func(*args, **kwargs)
            except Exception as e:
                if log_failure:
                    logger.warning(f"⚠️  Optional API call {func.__name__} failed: {e}")
                    logger.info(f"💡 Continuing with fallback value: {fallback_value}")
                return fallback_value

        return wrapper
    return decorator
```

**Example usage in mcp-tcga:**

```python
# servers/mcp-tcga/src/mcp_tcga/__main__.py

from shared.utils.api_retry import retry_with_backoff, optional_api_call
import requests

@app.tool()
@retry_with_backoff(max_retries=3, exceptions=(requests.RequestException,))
async def query_tcga_cohorts(cancer_type: str):
    """Query TCGA cohorts with automatic retry."""

    # This will retry up to 3 times with exponential backoff
    response = requests.get(
        "https://api.gdc.cancer.gov/cases",
        params={"filters": f'{{"op":"in","content":{{"field":"primary_site","value":["{cancer_type}"]}}}}'}
    )
    response.raise_for_status()
    return response.json()

@app.tool()
@optional_api_call(fallback_value=None)
async def fetch_survival_data(cohort_id: str):
    """Fetch survival data - optional, degrades gracefully if API down."""

    # If this fails, returns None instead of crashing entire analysis
    response = requests.get(f"https://api.gdc.cancer.gov/survival/{cohort_id}")
    response.raise_for_status()
    return response.json()
```

**Files to create/modify:**
1. Create `shared/utils/api_retry.py`
2. Modify `servers/mcp-tcga/src/mcp_tcga/__main__.py` - Add retry to all API calls
3. Modify `servers/mcp-huggingface/src/mcp_huggingface/__main__.py` - Add retry
4. Modify `servers/mcp-seqera/src/mcp_seqera/__main__.py` - Add retry

**Acceptance Criteria:**
- [ ] Retry utilities created and tested
- [ ] All external API calls use retry decorator
- [ ] Optional API calls degrade gracefully with fallback values
- [ ] Logs clearly indicate retry attempts
- [ ] Test with simulated API failures (mock server down)

**Estimated Time:** 4-6 hours

---

#### **WI-6: Cost Tracking & Monitoring**
**Addresses:** R2 (Cost Overruns) - Severity 8/10
**Effort:** 3-4 hours
**Impact:** 🟡 MEDIUM - Enables cost awareness

**Problem:**
- No way to track actual costs per patient analysis
- Can't identify cost drivers until bill arrives
- No warnings before expensive operations

**Solution - Add Cost Tracking:**

**File to create:** `shared/utils/cost_tracker.py`

```python
"""Cost tracking utilities for MCP servers."""

import time
import logging
from pathlib import Path
from datetime import datetime
from typing import Optional

logger = logging.getLogger(__name__)

class CostTracker:
    """Track costs for patient analysis."""

    def __init__(self, patient_id: str, output_dir: str = "./costs"):
        self.patient_id = patient_id
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        self.costs = {
            "patient_id": patient_id,
            "start_time": datetime.now().isoformat(),
            "compute_costs": [],
            "api_costs": [],
            "token_costs": [],
            "total_cost_usd": 0.0
        }

    def add_compute_cost(
        self,
        operation: str,
        cpu_hours: float,
        cost_per_hour: float = 0.50,
        notes: Optional[str] = None
    ):
        """Track computational costs (e.g., STAR alignment, DeepCell)."""
        cost = cpu_hours * cost_per_hour

        self.costs["compute_costs"].append({
            "operation": operation,
            "cpu_hours": cpu_hours,
            "cost_per_hour_usd": cost_per_hour,
            "cost_usd": cost,
            "notes": notes
        })

        self.costs["total_cost_usd"] += cost

        logger.info(f"💰 Compute cost: {operation} = ${cost:.2f} ({cpu_hours:.2f} CPU-hours @ ${cost_per_hour}/hr)")

        if cost > 10.0:
            logger.warning(f"⚠️  High compute cost detected: ${cost:.2f}")

    def add_api_cost(
        self,
        service: str,
        num_calls: int,
        cost_per_call: float,
        notes: Optional[str] = None
    ):
        """Track external API costs (TCGA, HuggingFace, Seqera)."""
        cost = num_calls * cost_per_call

        self.costs["api_costs"].append({
            "service": service,
            "num_calls": num_calls,
            "cost_per_call_usd": cost_per_call,
            "cost_usd": cost,
            "notes": notes
        })

        self.costs["total_cost_usd"] += cost

        logger.info(f"💰 API cost: {service} = ${cost:.2f} ({num_calls} calls @ ${cost_per_call}/call)")

    def add_token_cost(
        self,
        input_tokens: int,
        output_tokens: int,
        model: str = "claude-sonnet-4.5"
    ):
        """Track Claude token costs."""
        # Pricing as of Dec 2025
        prices = {
            "claude-sonnet-4.5": {"input": 0.000003, "output": 0.000015},  # $3/$15 per 1M tokens
            "claude-opus-4.5": {"input": 0.000015, "output": 0.000075}
        }

        pricing = prices.get(model, prices["claude-sonnet-4.5"])
        cost = (input_tokens * pricing["input"]) + (output_tokens * pricing["output"])

        self.costs["token_costs"].append({
            "model": model,
            "input_tokens": input_tokens,
            "output_tokens": output_tokens,
            "cost_usd": cost
        })

        self.costs["total_cost_usd"] += cost

        logger.info(f"💰 Token cost: {model} = ${cost:.2f} ({input_tokens} in, {output_tokens} out)")

    def save_report(self):
        """Save cost report to JSON file."""
        import json

        self.costs["end_time"] = datetime.now().isoformat()

        output_file = self.output_dir / f"{self.patient_id}_cost_report.json"
        with open(output_file, 'w') as f:
            json.dump(self.costs, f, indent=2)

        logger.info(f"💰 Total cost for {self.patient_id}: ${self.costs['total_cost_usd']:.2f}")
        logger.info(f"📊 Cost report saved: {output_file}")

        # Warn if over budget
        if self.costs['total_cost_usd'] > 50.0:
            logger.warning(f"⚠️  Cost exceeds $50 threshold: ${self.costs['total_cost_usd']:.2f}")

    def __enter__(self):
        """Context manager support."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Auto-save on exit."""
        self.save_report()
```

**Example usage:**

```python
# In PatientOne workflow

from shared.utils.cost_tracker import CostTracker

async def run_patient_one_analysis(patient_id: str):
    """Run complete PatientOne workflow with cost tracking."""

    with CostTracker(patient_id) as tracker:
        # TEST_1: Clinical + Genomic
        tracker.add_api_cost("TCGA", num_calls=5, cost_per_call=0.0)  # Free
        tracker.add_token_cost(input_tokens=1500, output_tokens=2000)

        # TEST_2: Multi-Omics
        tracker.add_compute_cost("Multi-omics preprocessing", cpu_hours=0.25, cost_per_hour=0.50)
        tracker.add_token_cost(input_tokens=2500, output_tokens=3500)

        # TEST_3: Spatial (expensive!)
        tracker.add_compute_cost("STAR alignment", cpu_hours=1.5, cost_per_hour=2.0, notes="16 cores")
        tracker.add_token_cost(input_tokens=1800, output_tokens=2200)

        # TEST_4: Imaging
        tracker.add_compute_cost("DeepCell segmentation", cpu_hours=0.5, cost_per_hour=4.0, notes="GPU")
        tracker.add_api_cost("HuggingFace", num_calls=10, cost_per_call=0.02)

        # TEST_5: Integration
        tracker.add_token_cost(input_tokens=3500, output_tokens=5000)

        # Auto-saves cost report on exit
        # Output: PAT001_cost_report.json with itemized costs
```

**Files to create/modify:**
1. Create `shared/utils/cost_tracker.py`
2. Add cost tracking to PatientOne test scripts
3. Update `COST_ANALYSIS.md` with "How to Track Actual Costs" section

**Acceptance Criteria:**
- [ ] CostTracker class created and tested
- [ ] Can track compute, API, and token costs separately
- [ ] Saves JSON report per patient
- [ ] Logs warnings for costs > $50
- [ ] Documentation explains how to use

**Estimated Time:** 3-4 hours

---

### 🟢 MEDIUM PRIORITY (Week 2 - Nice to Have)

These items improve documentation and governance but aren't blocking for research use.

---

#### **WI-7: Data Governance Documentation**
**Addresses:** R1 (Patient Privacy) - Severity 9/10 (but lower priority for research data)
**Effort:** 2-3 hours
**Impact:** 🟢 MEDIUM - Important for governance

**Problem:**
- No documented data handling practices
- Unclear who can access what data
- No guidance on de-identification

**Solution:**

**File to create:** `docs/DATA_GOVERNANCE.md`

```markdown
# Data Governance Guide

## Data Classification

### Research Data (Current Scope)
- **Type:** PDX models, cell lines, de-identified patient samples
- **Sensitivity:** Medium
- **Regulations:** IRB approval, institutional data policies
- **Access:** Authorized research personnel only

### Clinical Patient Data (Future Scope)
- **Type:** Identifiable patient genomic, clinical, imaging data
- **Sensitivity:** HIGH
- **Regulations:** HIPAA, GDPR, institutional compliance
- **Access:** Strictly controlled, BAA required

## Access Control Guidelines

### Recommended File Permissions

```bash
# Data directories (read-only for processing)
chmod 750 /path/to/MULTIOMICS_DATA_DIR
chmod 750 /path/to/SPATIAL_DATA_DIR
chmod 750 /path/to/IMAGE_DATA_DIR

# Output directories (write access for results)
chmod 770 /path/to/results

# Restrict to authorized group
chown -R :bioinformatics /path/to/data
```

### User Roles

| Role | Access | Permissions |
|------|--------|-------------|
| Bioinformatics Lead | All data | Read/write, can add users |
| Analyst | Project-specific data | Read-only |
| Clinician | De-identified results only | Read-only reports |
| IT Admin | No data access | System maintenance only |

## De-identification for Clinical Data

### HIPAA Safe Harbor Method (18 Identifiers to Remove)

Before processing ANY clinical data, remove:

1. Names (patient, relatives, employers)
2. Geographic subdivisions smaller than state
3. Dates (except year)
4. Phone/fax numbers
5. Email addresses
6. Social Security numbers
7. Medical record numbers
8. Health plan numbers
9. Account numbers
10. Certificate/license numbers
11. Vehicle identifiers
12. Device identifiers
13. URLs
14. IP addresses
15. Biometric identifiers
16. Full-face photos
17. Any other unique identifying number

### De-identification Script

```python
# tools/deidentify_patient_data.py

import pandas as pd
import hashlib
from datetime import datetime

def deidentify_patient_id(patient_id: str, salt: str) -> str:
    """Convert real patient ID to de-identified pseudonym."""
    return "PAT" + hashlib.sha256(f"{patient_id}{salt}".encode()).hexdigest()[:8]

def deidentify_dates(date_str: str) -> str:
    """Convert full dates to year only."""
    try:
        date_obj = datetime.strptime(date_str, "%Y-%m-%d")
        return str(date_obj.year)
    except:
        return "REDACTED"

def deidentify_clinical_file(input_file: str, output_file: str, salt: str):
    """De-identify clinical data file."""
    df = pd.read_csv(input_file)

    # Remove direct identifiers
    cols_to_remove = ['patient_name', 'mrn', 'ssn', 'dob', 'phone', 'email', 'address']
    df = df.drop(columns=[col for col in cols_to_remove if col in df.columns])

    # Pseudonymize patient ID
    if 'patient_id' in df.columns:
        df['patient_id'] = df['patient_id'].apply(lambda x: deidentify_patient_id(x, salt))

    # Convert dates to year only
    date_cols = [col for col in df.columns if 'date' in col.lower()]
    for col in date_cols:
        df[col] = df[col].apply(deidentify_dates)

    # Save de-identified version
    df.to_csv(output_file, index=False)
    print(f"✅ De-identified data saved: {output_file}")
    print(f"⚠️  VERIFY manually before using in production!")
```

## Logging & Audit Trail

### What to Log (Securely)
- ✅ User who initiated analysis
- ✅ Timestamp
- ✅ Input file paths (NOT contents)
- ✅ Analysis type and parameters
- ✅ Output file paths
- ✅ Errors and warnings

### What NOT to Log
- ❌ Patient identifiers (names, MRNs)
- ❌ Genomic sequences
- ❌ Clinical notes
- ❌ PHI/PII in any form

## Data Retention

- **Raw input data:** Retain per IRB/institutional policy (typically 3-7 years)
- **Analysis results:** Retain 3 years minimum
- **Logs:** Retain 1 year
- **Deletion:** Secure deletion (not just rm, use shred or equivalent)

## Claude API Considerations

### Current Setup (Claude Desktop)
- Data sent to Anthropic cloud
- Subject to Anthropic's data retention policy
- Review: https://www.anthropic.com/legal/privacy

### For Clinical Data (Future)
- ❌ Do NOT use cloud Claude API with identifiable patient data
- ✅ Consider self-hosted deployment
- ✅ Ensure Business Associate Agreement (BAA) in place

## Compliance Checklist

Before processing clinical patient data:

- [ ] IRB approval obtained
- [ ] Data de-identified per HIPAA Safe Harbor
- [ ] Access controls configured (file permissions, user roles)
- [ ] Secure logging implemented (no PHI in logs)
- [ ] Data retention policy documented
- [ ] BAA signed with Anthropic (if using cloud API)
- [ ] Security risk assessment completed
- [ ] Regular audits scheduled
```

**Acceptance Criteria:**
- [ ] DATA_GOVERNANCE.md created
- [ ] De-identification script provided
- [ ] Access control guidelines clear
- [ ] Link from main README

**Estimated Time:** 2-3 hours

---

## Summary & Timeline

### Week 1 (Priority: Critical items)

| Day | Work Item | Hours | Cumulative |
|-----|-----------|-------|------------|
| Mon | WI-1: Server Implementation Status Doc | 4-6h | 4-6h |
| Tue | WI-2: Runtime DRY_RUN Warnings | 3-4h | 7-10h |
| Wed | WI-3: Input Validation (part 1) | 4h | 11-14h |
| Thu | WI-3: Input Validation (part 2) | 2-4h | 13-18h |
| Fri | WI-4: Research Use Disclaimers | 2-3h | 15-21h |

**Week 1 Total:** 15-21 hours (~3-4 days solo)

### Week 2 (Priority: High + Medium items)

| Day | Work Item | Hours | Cumulative |
|-----|-----------|-------|------------|
| Mon | WI-5: Error Handling & Retry Logic | 4-6h | 4-6h |
| Tue | WI-6: Cost Tracking | 3-4h | 7-10h |
| Wed | WI-7: Data Governance Docs | 2-3h | 9-13h |
| Thu | Testing & validation of all changes | 4h | 13-17h |
| Fri | Documentation polish, merge to main | 2h | 15-19h |

**Week 2 Total:** 15-19 hours (~3-4 days solo)

**Grand Total:** 30-40 hours (1-2 weeks solo developer)

---

## Risk Reduction Matrix

**Before this work:**

| Risk | Severity | Mitigation Status |
|------|----------|-------------------|
| R3: Mocked servers in production | 🔴 9/10 | ❌ Not mitigated |
| R5: Poor data quality | 🟡 5/10 | ❌ Not mitigated |
| R7: Incorrect recommendations | 🟡 6/10 | ❌ Not mitigated |
| R4: API failures | 🟡 6/10 | ❌ Not mitigated |
| R2: Cost overruns | 🔴 8/10 | ❌ Not mitigated |

**After this work:**

| Risk | Severity | Mitigation Status | Reduction |
|------|----------|-------------------|-----------|
| R3: Mocked servers | 🔴 9/10 → 🟡 4/10 | ✅ Documented + Warnings | **55% reduction** |
| R5: Data quality | 🟡 5/10 → 🟢 2/10 | ✅ Input validation | **60% reduction** |
| R7: Incorrect recommendations | 🟡 6/10 → 🟡 3/10 | ✅ Disclaimers added | **50% reduction** |
| R4: API failures | 🟡 6/10 → 🟢 3/10 | ✅ Retry logic | **50% reduction** |
| R2: Cost overruns | 🔴 8/10 → 🟡 5/10 | ✅ Cost tracking | **37% reduction** |

**Overall Risk Reduction: ~50% across critical areas**

---

## Next Steps After This Work

Once these 7 work items are complete, the next recommended priorities are:

1. **Implement Real TCGA Integration** (1 week)
   - Addresses R3 by de-mocking a critical server
   - Enables real cohort comparisons

2. **Add Comprehensive Test Suite with Real Data** (1 week)
   - Test with 3-5 real PDX datasets
   - Validate against known results
   - Addresses R5 (data quality) further

3. **Pilot Deployment** (2-3 weeks)
   - Process 10 real PDX samples
   - Monitor costs, errors, user feedback
   - Iterate based on findings

---

**Questions or clarifications needed on any work item?**
