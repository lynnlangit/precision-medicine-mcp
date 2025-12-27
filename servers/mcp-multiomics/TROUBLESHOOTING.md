# Pydantic V2 Boolean Parsing Fix

## Problem

When `MULTIOMICS_DRY_RUN="false"` was set in the Claude Desktop config, the server was still running in DRY_RUN mode. This prevented the Stouffer's meta-analysis from using directionality features.

## Root Cause

**Pydantic V2 Settings** parses environment variable string `"false"` as boolean `True`. This is a known behavior where any non-empty string is considered truthy.

```python
# Before fix:
os.environ['MULTIOMICS_DRY_RUN'] = 'false'
config = MultiOmicsConfig()
print(config.dry_run)  # Prints: True (WRONG!)
```

Additionally, environment variables weren't being read at all because:
1. No `env_prefix` was configured in SettingsConfigDict
2. Pydantic was looking for `data_dir` instead of `MULTIOMICS_DATA_DIR`

## Solution

### 1. Added `env_prefix` to SettingsConfigDict

**File:** `src/mcp_multiomics/config.py`

```python
model_config = SettingsConfigDict(
    env_file=".env",
    env_file_encoding="utf-8",
    case_sensitive=False,
    env_prefix="MULTIOMICS_",  # Added this line
)
```

This tells Pydantic Settings to look for environment variables prefixed with `MULTIOMICS_`:
- `data_dir` field reads from `MULTIOMICS_DATA_DIR`
- `cache_dir` field reads from `MULTIOMICS_CACHE_DIR`
- `dry_run` field reads from `MULTIOMICS_DRY_RUN`

### 2. Added `model_validator` to Fix Boolean Parsing

```python
@model_validator(mode='after')
def parse_boolean_env_vars(self):
    """Fix boolean parsing from environment variables.

    Pydantic V2 Settings parses string 'false' as boolean True,
    so we need to explicitly check environment variable values.
    """
    env_val = os.getenv('MULTIOMICS_DRY_RUN')
    if env_val is not None:
        self.dry_run = env_val.lower() not in ('false', '0', 'no', 'off', '')
    return self
```

The validator runs **after** Pydantic Settings has loaded all values, then:
- Checks the raw environment variable value
- Converts string `"false"` to boolean `False`
- Handles common falsy values: `"false"`, `"0"`, `"no"`, `"off"`, `""`

### 3. Removed Directory Creation from `__init__`

The original `__init__` method created directories immediately:

```python
def __init__(self, **kwargs):
    super().__init__(**kwargs)
    if not self.dry_run:
        self.data_dir.mkdir(parents=True, exist_ok=True)  # Failed at import time!
```

This caused `OSError: [Errno 30] Read-only file system: '/workspace'` because:
- The config module is imported with `config = MultiOmicsConfig()` at the bottom
- This happens **before** Claude Desktop injects environment variables
- It tried to create `/workspace` (the default) instead of the actual path

**Solution:** Moved directory creation to a separate method:

```python
def ensure_directories(self):
    """Create directories if in real mode. Call this before file operations."""
    if not self.dry_run:
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
```

Tools that need to write files can call `config.ensure_directories()` before file operations.

## Verification

### Boolean Parsing Test

```bash
# Test: DRY_RUN=false → False
MULTIOMICS_DRY_RUN="false" python -c "from mcp_multiomics.config import config; print(config.dry_run)"
# Output: False ✅

# Test: DRY_RUN=true → True
MULTIOMICS_DRY_RUN="true" python -c "from mcp_multiomics.config import config; print(config.dry_run)"
# Output: True ✅

# Test: Not set → True (default)
python -c "from mcp_multiomics.config import config; print(config.dry_run)"
# Output: True ✅
```

### Environment Variable Parsing Test

```bash
MULTIOMICS_DATA_DIR="/path/to/data" \
MULTIOMICS_CACHE_DIR="/path/to/cache" \
MULTIOMICS_DRY_RUN="false" \
python -c "from mcp_multiomics.config import config; \
print(f'DRY_RUN: {config.dry_run}'); \
print(f'Data dir: {config.data_dir}'); \
print(f'Cache dir: {config.cache_dir}')"

# Output:
# DRY_RUN: False ✅
# Data dir: /path/to/data ✅
# Cache dir: /path/to/cache ✅
```

### Test Suite Results

All 29 tests pass:
- 14 integration tests (data loading, alignment, normalization)
- 15 Stouffer's meta-analysis tests (including directionality)
- 84% code coverage

```bash
MULTIOMICS_DRY_RUN="false" pytest tests/ -v
# Result: 29 passed in 1.86s ✅
```

### Directionality Tests

```bash
MULTIOMICS_DRY_RUN="false" pytest tests/test_stouffer.py -v -k "directionality"
# Result: 3 passed (test_p_to_z_with_directionality, test_meta_analyze_with_directionality, test_meta_analysis_with_directionality) ✅
```

## Claude Desktop Configuration

The Claude Desktop config file is correctly set:

**File:** `~/Library/Application Support/Claude/claude_desktop_config.json`

```json
"multiomics": {
  "command": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-multiomics/venv/bin/python",
  "args": ["-m", "mcp_multiomics"],
  "cwd": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-multiomics",
  "env": {
    "PYTHONPATH": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/servers/mcp-multiomics/src",
    "MULTIOMICS_DATA_DIR": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/data/multiomics",
    "MULTIOMICS_CACHE_DIR": "/Users/lynnlangit/Documents/GitHub/precision-medicine-mcp/data/cache/multiomics",
    "MULTIOMICS_DRY_RUN": "false"  ✅
  }
}
```

## Next Steps for User

**To activate the fix in Claude Desktop:**

1. **Restart Claude Desktop completely:**
   - Quit Claude Desktop (⌘Q)
   - Relaunch Claude Desktop

2. **Test the multiomics server:**
   ```
   I have p-values from RNA and Protein analyses for genes TP53, MYC, KRAS:

   RNA p-values: [0.001, 0.002, 0.05]
   RNA effect sizes (log2FC): [2.5, 1.8, 1.2]

   Protein p-values: [0.005, 0.01, 0.03]
   Protein effect sizes (log2FC): [2.0, 1.6, 1.1]

   Please combine these using Stouffer's method with directionality.
   ```

3. **Expected result:**
   - Status: Success (no DRY_RUN mode message)
   - Directionality used: **Yes** ✅
   - Combined Z-scores should be positive (matching effect size direction)

## Technical Details

### Why `mode='after'` for the Validator?

Using `mode='after'` ensures the validator runs **after** Pydantic Settings has:
1. Read the `.env` file
2. Read environment variables
3. Applied the `env_prefix`
4. Set all field values

If we used `mode='before'`, it would run before environment variables are loaded.

### Why Not Use `default_factory`?

The `default_factory` approach was attempted:

```python
dry_run: bool = Field(
    default_factory=lambda: _parse_bool_env("MULTIOMICS_DRY_RUN", default=True)
)
```

This failed because:
1. The factory function runs at field definition time
2. At that point, environment variables might not be set yet
3. It caused the `__init__` method to run at import time
4. Led to directory creation errors before paths were configured

### Alternative Solutions Not Used

1. **Custom Pydantic Type:** Could create a custom `BooleanEnv` type, but the validator is simpler
2. **Field Serializer:** Only affects output, not input parsing
3. **Config Preprocessor:** Would require overriding `__init__` before super()

The `model_validator(mode='after')` approach is the cleanest solution that:
- Works with Pydantic Settings architecture
- Doesn't interfere with environment variable loading
- Is easy to understand and maintain

## Files Modified

1. **src/mcp_multiomics/config.py** (lines 24-114)
   - Added `env_prefix="MULTIOMICS_"` to SettingsConfigDict
   - Simplified `dry_run` field (removed `default_factory`)
   - Added `parse_boolean_env_vars()` model validator
   - Replaced `__init__` directory creation with `ensure_directories()` method

2. **configs/claude_desktop_config.json** (line 89)
   - Changed `"MULTIOMICS_DRY_RUN": "true"` → `"false"`

## Status

✅ **FIXED** - All tests passing, ready for Claude Desktop testing

---

**Created:** November 11, 2025
**Author:** Claude (Sonnet 4.5)
**Issue:** Pydantic V2 Boolean parsing and environment variable configuration
**Resolution:** Added env_prefix and model_validator to properly parse boolean environment variables
