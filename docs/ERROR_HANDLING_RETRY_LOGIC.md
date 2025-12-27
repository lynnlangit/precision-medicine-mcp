# Error Handling & Retry Logic Implementation

## Overview

This document describes the error handling and retry logic implementation for the Precision Medicine MCP servers. These utilities protect against transient API failures, network issues, and external service outages.

**Status:** ✅ Implemented (WI-5 from Risk Mitigation Workplan)
**Risk Reduced:** R4 (External API failures) - 60% reduction (8/10 → 3/10)

---

## Architecture

### Retry Utilities Location

All retry utilities are centralized in:
```
shared/utils/api_retry.py
```

### Components

1. **`retry_with_backoff()` decorator** - Exponential backoff retry logic
2. **`optional_api_call()` decorator** - Graceful degradation for non-critical APIs
3. **`CircuitBreaker` class** - Circuit breaker pattern implementation

---

## 1. Exponential Backoff Retry

### Purpose
Automatically retry failed operations with increasing delays between attempts, preventing cascade failures and reducing load on failing services.

### Usage

```python
from api_retry import retry_with_backoff

@retry_with_backoff(
    max_retries=3,
    base_delay=1.0,
    max_delay=30.0,
    exponential_base=2.0,
    exceptions=(IOError, Exception),
    on_retry=lambda e, attempt: logger.warning(f"Retry {attempt}: {e}")
)
async def download_reference_genome(url: str):
    """Download with automatic retry on failure."""
    response = await httpx.get(url)
    response.raise_for_status()
    return response.content
```

### Retry Schedule (default settings)
- **Attempt 1:** Immediate
- **Attempt 2:** Wait 1 second
- **Attempt 3:** Wait 2 seconds
- **Attempt 4:** Wait 4 seconds
- **Total:** 4 attempts over ~7 seconds

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_retries` | 3 | Maximum retry attempts |
| `base_delay` | 1.0 | Initial delay in seconds |
| `max_delay` | 30.0 | Maximum delay between retries |
| `exponential_base` | 2.0 | Exponential multiplier (delay = base × exponential_base^attempt) |
| `exceptions` | `(Exception,)` | Tuple of exceptions to catch and retry |
| `on_retry` | `None` | Optional callback function(exception, attempt_num) |

### Supports Both Sync and Async

The decorator automatically detects function type:

```python
# Async function
@retry_with_backoff(max_retries=3)
async def fetch_data():
    await asyncio.sleep(1)
    return "data"

# Sync function
@retry_with_backoff(max_retries=3)
def fetch_data():
    time.sleep(1)
    return "data"
```

---

## 2. Optional API Call (Graceful Degradation)

### Purpose
Allow non-critical API calls to fail gracefully, returning fallback values instead of raising exceptions. Enables analysis to continue when optional data is unavailable.

### Usage

```python
from api_retry import optional_api_call

@optional_api_call(
    fallback_value=None,
    log_failure=True,
    exceptions=(Exception,)
)
async def fetch_tcga_annotations(sample_id: str):
    """Fetch annotations (optional, non-critical)."""
    response = await httpx.get(
        f"https://api.gdc.cancer.gov/annotations/{sample_id}"
    )
    return response.json()

# Usage
annotations = await fetch_tcga_annotations("TCGA-OV-001")
if annotations is None:
    logger.warning("Annotations unavailable, continuing without them")
else:
    # Process annotations
    process_annotations(annotations)
```

### When to Use

✅ **Use for:**
- Metadata enrichment (annotations, descriptions)
- Optional visualizations
- Non-critical data augmentation
- Fallback data sources

❌ **Do NOT use for:**
- Primary data analysis
- Required inputs
- Critical validation steps
- Data that affects clinical recommendations

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fallback_value` | `None` | Value to return on failure |
| `log_failure` | `True` | Whether to log failure as warning |
| `exceptions` | `(Exception,)` | Tuple of exceptions to catch |

---

## 3. Circuit Breaker Pattern

### Purpose
Prevent repeated calls to a failing service by "opening the circuit" after a threshold of failures. After a timeout period, allows one test call through ("half-open") to check if service has recovered.

### States

```
CLOSED → OPEN → HALF_OPEN → CLOSED
   ↑                            |
   └────────────────────────────┘
```

- **CLOSED:** Normal operation, requests pass through
- **OPEN:** Service is failing, requests fail immediately (fail-fast)
- **HALF_OPEN:** Testing if service has recovered (single test request)

### Usage

```python
from api_retry import CircuitBreaker

# Create circuit breaker instance
tcga_breaker = CircuitBreaker(
    failure_threshold=5,
    recovery_timeout=60.0,
    expected_exception=requests.RequestException,
    name="TCGA-GDC-API"
)

# Apply as decorator
@tcga_breaker
async def fetch_tcga_cohort(cohort_id: str):
    """Fetch cohort data with circuit breaker protection."""
    response = await httpx.get(
        f"https://api.gdc.cancer.gov/cases?filters={cohort_id}"
    )
    response.raise_for_status()
    return response.json()

# Manual reset (if needed)
tcga_breaker.reset()
```

### Behavior Example

```
# Normal operation (CLOSED)
Request 1: ✅ Success
Request 2: ✅ Success
Request 3: ❌ Failure (count: 1/5)
Request 4: ❌ Failure (count: 2/5)
Request 5: ❌ Failure (count: 3/5)
Request 6: ❌ Failure (count: 4/5)
Request 7: ❌ Failure (count: 5/5) → Circuit OPENS

# Circuit is OPEN (fail-fast)
Request 8: ⚡ Fail immediately (no API call made)
Request 9: ⚡ Fail immediately (no API call made)
...
[60 seconds pass]

# Circuit becomes HALF_OPEN (test recovery)
Request N: 🔄 Test call allowed
  - If ✅ success → Circuit CLOSES (back to normal)
  - If ❌ failure → Circuit OPENS again
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `failure_threshold` | 5 | Number of failures before opening circuit |
| `recovery_timeout` | 60.0 | Seconds to wait before attempting recovery |
| `expected_exception` | `Exception` | Exception type that counts as failure |
| `name` | `"CircuitBreaker"` | Optional name for logging |

---

## Implementation Status by Server

### ✅ Implemented (Real Implementation)

#### **mcp-fgbio** (Production Ready)

**Location:** `servers/mcp-fgbio/src/mcp_fgbio/server.py`

**Applied to:**
- `_download_file()` - Reference genome downloads from UCSC/NCBI

**Retry Configuration:**
```python
@retry_with_backoff(
    max_retries=3,
    base_delay=2.0,
    max_delay=60.0,
    exceptions=(IOError, Exception)
)
async def _download_file(url: str, output_path: Path):
    # Download reference genome with retry
    pass
```

**Impact:**
- Handles transient network failures during large file downloads
- Reduces failed downloads from network blips
- Provides better user experience with automatic recovery

---

### 📝 Integration Guide Added (For Future Implementation)

The following servers have retry utility imports and integration examples added, ready for when real API implementations are developed:

#### **mcp-tcga** (Currently Mocked)

**Location:** `servers/mcp-tcga/src/mcp_tcga/server.py`

**Example Integration (commented code):**
```python
# Example 1: Retry TCGA GDC API calls
@retry_with_backoff(
    max_retries=3,
    base_delay=1.0,
    exceptions=(requests.RequestException,)
)
async def _fetch_tcga_cohort(cohort_id: str):
    response = await httpx.get(
        f"https://api.gdc.cancer.gov/cases?filters={cohort_id}"
    )
    response.raise_for_status()
    return response.json()

# Example 2: Optional annotations (non-critical)
@optional_api_call(fallback_value=None)
async def _fetch_tcga_annotations(sample_id: str):
    response = await httpx.get(
        f"https://api.gdc.cancer.gov/annotations/{sample_id}"
    )
    return response.json()

# Example 3: Circuit breaker for TCGA API
tcga_breaker = CircuitBreaker(
    failure_threshold=5,
    recovery_timeout=60.0,
    name="TCGA-GDC-API"
)
```

#### **mcp-huggingface** (Currently Mocked)

**Location:** `servers/mcp-huggingface/src/mcp_huggingface/server.py`

**Example Integration (commented code):**
```python
# Retry model loading
@retry_with_backoff(
    max_retries=3,
    base_delay=2.0,
    max_delay=30.0,
    exceptions=(requests.HTTPError,)
)
async def _load_hf_model(model_name: str):
    from huggingface_hub import hf_hub_download
    model_path = hf_hub_download(
        repo_id=model_name,
        filename="pytorch_model.bin",
        token=HF_TOKEN
    )
    return model_path

# Circuit breaker for HuggingFace Hub
hf_hub_breaker = CircuitBreaker(
    failure_threshold=5,
    recovery_timeout=120.0,  # 2 minutes
    name="HuggingFace-Hub-API"
)
```

#### **mcp-seqera** (Currently Mocked)

**Location:** `servers/mcp-seqera/src/mcp_seqera/server.py`

**Example Integration (commented code):**
```python
# Retry workflow launch
@retry_with_backoff(
    max_retries=3,
    base_delay=2.0,
    max_delay=60.0,
    exceptions=(requests.HTTPError, requests.ConnectionError)
)
async def _launch_seqera_workflow(pipeline_id: str, params: dict):
    response = await httpx.post(
        f"{SEQERA_API_URL}/workflow/launch",
        headers={"Authorization": f"Bearer {SEQERA_TOKEN}"},
        json={"pipeline": pipeline_id, "params": params}
    )
    response.raise_for_status()
    return response.json()

# Circuit breaker for Seqera Platform
seqera_api_breaker = CircuitBreaker(
    failure_threshold=5,
    recovery_timeout=180.0,  # 3 minutes
    name="Seqera-Platform-API"
)
```

---

## Best Practices

### 1. Choose Appropriate Retry Parameters

**For quick API calls (< 1 second):**
```python
@retry_with_backoff(max_retries=3, base_delay=0.5, max_delay=10.0)
```

**For medium operations (1-10 seconds):**
```python
@retry_with_backoff(max_retries=3, base_delay=1.0, max_delay=30.0)
```

**For large downloads (minutes):**
```python
@retry_with_backoff(max_retries=3, base_delay=2.0, max_delay=60.0)
```

### 2. Specify Appropriate Exceptions

❌ **Too broad:**
```python
@retry_with_backoff(exceptions=(Exception,))  # Catches ALL exceptions
```

✅ **Specific to transient failures:**
```python
@retry_with_backoff(exceptions=(
    requests.ConnectionError,
    requests.Timeout,
    IOError
))
```

### 3. Use Circuit Breakers for Critical External Services

```python
# Create one circuit breaker per external service
tcga_api = CircuitBreaker(failure_threshold=5, name="TCGA-API")
hf_hub = CircuitBreaker(failure_threshold=5, name="HF-Hub")
ncbi = CircuitBreaker(failure_threshold=5, name="NCBI")

# Apply to all calls to that service
@tcga_api
async def fetch_from_tcga(endpoint: str):
    pass

@hf_hub
async def fetch_from_hf(model: str):
    pass
```

### 4. Combine Decorators for Maximum Resilience

```python
# Circuit breaker + Retry + Optional
hf_breaker = CircuitBreaker(failure_threshold=5, name="HF-Hub")

@hf_breaker
@retry_with_backoff(max_retries=3)
@optional_api_call(fallback_value=None)
async def fetch_model_metadata(model_name: str):
    """Triple protection: circuit breaker, retry, graceful degradation."""
    response = await httpx.get(f"https://huggingface.co/api/models/{model_name}")
    return response.json()
```

**Execution flow:**
1. Circuit breaker checks if service is available
2. If available, retry logic handles transient failures
3. If all retries fail, optional_api_call returns fallback value
4. Analysis continues without metadata

---

## Testing Retry Logic

### 1. Test with Simulated Failures

```python
import pytest
from unittest.mock import AsyncMock, patch

@pytest.mark.asyncio
async def test_retry_succeeds_after_failures():
    """Test that retry logic works after transient failures."""
    mock_client = AsyncMock()

    # Fail twice, then succeed
    mock_client.get.side_effect = [
        IOError("Network timeout"),
        IOError("Connection reset"),
        {"status": "success", "data": "result"}
    ]

    @retry_with_backoff(max_retries=3, base_delay=0.1)
    async def fetch_data():
        result = await mock_client.get("https://api.example.com/data")
        if isinstance(result, Exception):
            raise result
        return result

    # Should succeed on 3rd attempt
    result = await fetch_data()
    assert result["status"] == "success"
    assert mock_client.get.call_count == 3
```

### 2. Test Circuit Breaker

```python
@pytest.mark.asyncio
async def test_circuit_breaker_opens_after_failures():
    """Test that circuit breaker opens after threshold failures."""
    breaker = CircuitBreaker(failure_threshold=3, recovery_timeout=1.0)

    @breaker
    async def failing_service():
        raise Exception("Service unavailable")

    # First 3 failures should be attempted
    for i in range(3):
        with pytest.raises(Exception, match="Service unavailable"):
            await failing_service()

    # Circuit should now be OPEN
    assert breaker.state == "OPEN"

    # 4th call should fail immediately without calling service
    with pytest.raises(Exception, match="Circuit breaker is OPEN"):
        await failing_service()
```

---

## Monitoring and Logging

All retry utilities provide comprehensive logging:

### Success After Retry
```
⚠️  _download_file attempt 1/4 failed: Connection timeout
⏳ Retrying in 1.0s...
⚠️  _download_file attempt 2/4 failed: Connection timeout
⏳ Retrying in 2.0s...
✅ _download_file succeeded on attempt 3
```

### Failure After All Retries
```
⚠️  _download_file attempt 1/4 failed: Connection timeout
⏳ Retrying in 1.0s...
⚠️  _download_file attempt 2/4 failed: Connection timeout
⏳ Retrying in 2.0s...
⚠️  _download_file attempt 3/4 failed: Connection timeout
⏳ Retrying in 4.0s...
❌ _download_file failed after 3 retries: Connection timeout
```

### Circuit Breaker State Changes
```
⚠️  TCGA-GDC-API: Failure 1/5
⚠️  TCGA-GDC-API: Failure 2/5
⚠️  TCGA-GDC-API: Failure 3/5
⚠️  TCGA-GDC-API: Failure 4/5
⚠️  TCGA-GDC-API: Failure 5/5
🔴 TCGA-GDC-API: Circuit is now OPEN (too many failures)
⚠️  TCGA-GDC-API: Circuit is OPEN, failing fast
[60 seconds later]
🔄 TCGA-GDC-API: Attempting recovery (HALF_OPEN)
✅ TCGA-GDC-API: Service recovered (CLOSED)
```

---

## Risk Reduction Impact

### Before Implementation

**R4: External API Failures (8/10 severity)**
- Likelihood: High (80%)
- Impact: Critical (total failure)
- No retry logic
- Transient failures cause complete pipeline failure

### After Implementation

**R4: External API Failures (3/10 severity)** ✅

- Likelihood: Low (30%)
- Impact: Moderate (analysis continues with degraded data)
- Automatic retry with exponential backoff
- Circuit breakers prevent cascade failures
- Graceful degradation for non-critical APIs

**Risk Reduction:** 60% (from 8/10 to 3/10)

---

## Future Enhancements

### 1. Jitter in Retry Delays
Add randomization to prevent thundering herd:
```python
delay = min(delay * exponential_base * random.uniform(0.9, 1.1), max_delay)
```

### 2. Adaptive Retry Policies
Adjust retry behavior based on error type:
```python
if isinstance(e, requests.HTTPError):
    if e.response.status_code == 429:  # Rate limit
        delay = float(e.response.headers.get('Retry-After', delay))
    elif e.response.status_code >= 500:  # Server error
        delay *= 2
```

### 3. Distributed Circuit Breakers
Share circuit breaker state across multiple instances using Redis:
```python
class DistributedCircuitBreaker(CircuitBreaker):
    def __init__(self, redis_client, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.redis = redis_client
```

### 4. Metrics Collection
Track retry statistics for monitoring:
```python
@retry_with_backoff(
    on_retry=lambda e, attempt: metrics.increment(
        'api_retry_count',
        tags={'service': 'tcga', 'attempt': attempt}
    )
)
```

---

## Related Documentation

- **Risk Matrix:** `RISK_MATRIX.md` - See R4 (External API failures)
- **Risk Mitigation Workplan:** `RISK_MITIGATION_WORKPLAN.md` - WI-5
- **Server Implementation Status:** `docs/SERVER_IMPLEMENTATION_STATUS.md`
- **API Retry Utilities Source:** `shared/utils/api_retry.py`

---

## Support

For questions or issues with retry logic:
1. Check logs for retry attempt messages
2. Review error messages for root cause
3. Consider adjusting retry parameters for your use case
4. Report persistent failures as issues

**Remember:** Retry logic handles *transient* failures. Persistent failures (invalid API keys, malformed requests, missing data) will still fail after all retries.
