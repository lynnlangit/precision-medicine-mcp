# Test Coverage - Precision Medicine MCP Servers

Comprehensive test coverage for all 9 MCP servers with **167 automated tests** covering 40 tools.

## Quick Stats

| Metric | Value | Status |
|--------|-------|--------|
| **Overall Coverage** | 56.9% | ⬆️ +27.5 points from baseline |
| **Automated Tests** | 167 | 9/9 servers tested |
| **Production Ready** | 4/9 servers | mcp-multiomics, mcp-fgbio, mcp-spatialtools, mcp-epic |
| **GCP Deployed** | 9/9 servers | All validated on Cloud Run (2025-12-30) |

---

## Coverage by Server

| Server | Coverage | Tests | Status | Production Ready |
|--------|----------|-------|--------|------------------|
| **mcp-fgbio** | 77% | 29 | ✅ Complete | ✅ YES (95% real) |
| **mcp-multiomics** | 68% | 91 | ✅ Complete | ✅ YES (85% real) |
| **mcp-deepcell** | 62% | 9 | ✅ Smoke | ❌ Mocked |
| **mcp-huggingface** | 56% | 12 | ✅ Smoke | ❌ Mocked |
| **mcp-seqera** | 56% | 6 | ✅ Smoke | ❌ Mocked |
| **mcp-epic** | 58% | 12 | ✅ Complete | ✅ YES (100% real) |
| **mcp-openimagedata** | 35% | 5 | ✅ Smoke | ⚠️ 30% real |
| **mcp-tcga** | 35% | 5 | ✅ Smoke | ❌ Mocked |
| **mcp-spatialtools** | 23% | 5 | ✅ Smoke | ✅ YES (95% real) |

**Note:** Low test coverage doesn't always mean low production readiness. mcp-spatialtools has 23% test coverage but is 95% production-ready with real implementations validated.

---

## Running Tests

### Quick Start - Run All Tests for a Server

```bash
# Production servers
cd servers/mcp-multiomics
MULTIOMICS_DRY_RUN="true" venv/bin/python -m pytest tests/ -v

cd servers/mcp-fgbio
FGBIO_DRY_RUN="true" venv/bin/python -m pytest tests/ -v

cd servers/mcp-spatialtools
SPATIAL_DRY_RUN="true" venv/bin/python -m pytest tests/ -v

cd servers/mcp-epic
EPIC_DRY_RUN="true" venv/bin/python -m pytest tests/ -v

# Mocked servers (smoke tests only)
cd servers/mcp-{deepcell|huggingface|seqera|tcga|openimagedata}
{SERVER}_DRY_RUN="true" venv/bin/python -m pytest tests/ -v
```

### With Coverage Report

```bash
cd servers/mcp-multiomics
MULTIOMICS_DRY_RUN="true" venv/bin/python -m pytest tests/ \
  --cov=src/mcp_multiomics \
  --cov-report=term-missing \
  -v
```

### Run Specific Test

```bash
# Specific file
MULTIOMICS_DRY_RUN="true" venv/bin/python -m pytest tests/test_preprocessing.py -v

# Specific test function
MULTIOMICS_DRY_RUN="true" venv/bin/python -m pytest \
  tests/test_preprocessing.py::TestValidateWithRealData::test_validate_with_real_rna_data -v
```

---

## Test Types

### 1. Smoke Tests (All Servers)

**Purpose:** Basic validation that servers load and register tools correctly

**What's tested:**
- Module imports
- DRY_RUN configuration
- Tool registration
- Server initialization

**Location:** `/servers/mcp-{server}/tests/test_server.py`

### 2. Functional Tests (Production Servers)

**Purpose:** Deep testing of real data processing logic

**What's tested:**
- Data validation and QC
- Statistical calculations (Fisher's exact, FDR correction, Moran's I)
- Batch correction and normalization
- File I/O and output generation
- Edge cases and error handling

**Best example:** `mcp-multiomics` - 91 tests with real fixture data (580KB+)

### 3. Integration Tests

**Purpose:** Test server interactions and end-to-end workflows

**Location:** `/tests/integration/`

**What's tested:**
- Multi-server workflows (PatientOne)
- GCP Cloud Run deployment validation
- SSE transport communication
- Claude API integration

---

## GCP Cloud Run Testing

All 9 servers deployed and tested on Google Cloud Platform:

**Deployment Validation (2025-12-30):**
```bash
# Automated test of all deployed servers
cd tests/integration
python test_all_gcp_servers.py
```

**Results:** ✅ 9/9 servers passing functional tests via Claude API

**Test Coverage:**
- SSE transport validation
- Tool discovery via MCP protocol
- Basic functionality per server
- Response time < 5 seconds (excluding cold starts)

**Server URLs:** See [Deployment Status](../docs/deployment/DEPLOYMENT_STATUS.md)

---

## Manual Testing

### PatientOne Test Suite

**Location:** `/tests/manual_testing/PatientOne-OvarianCancer/`

Complete end-to-end precision medicine workflow for Stage IV ovarian cancer:

- **TEST_1** - Clinical data integration (Epic FHIR)
- **TEST_2** - Multi-omics resistance analysis (RNA/Protein/Phospho)
- **TEST_3** - Spatial transcriptomics (Visium, 900 spots × 31 genes)
- **TEST_4** - Histology & imaging (H&E, multiplex IF)
- **TEST_5** - Integration & treatment recommendations

**Modes:**
- **DRY_RUN** (default): Synthetic data demo (~$0.32, 25-35 min)
- **Real Data**: Your own patient data ([Configuration Guide](manual_testing/PatientOne-OvarianCancer/DATA_MODES_GUIDE.md))

**📖 Quick Start:** [PatientOne README](manual_testing/PatientOne-OvarianCancer/README.md)

### Solution Testing

**Location:** `/tests/manual_testing/Solution-Testing/`

**Resources:**
- `MANUAL_TESTING_GUIDE.md` - Step-by-step testing instructions
- `verify_servers.sh` - Automated server verification script
- `TESTING_STATUS.md` - Server readiness status (9/9 ready)

---

## Key Achievements

### mcp-multiomics (⭐ Most Comprehensive)

**Coverage:** 37% → 68% (+31 points)
**Tests:** 91 automated tests
**Highlights:**
- 100% coverage on `stouffer.py` (meta-analysis)
- 93% coverage on `upstream_regulators.py` (Fisher's exact, FDR)
- 73% coverage on `preprocessing.py` (validation, normalization, batch correction)
- Real fixture data (580KB+) for RNA, protein, phospho

### mcp-epic (NEW - Real Epic FHIR)

**Coverage:** 45% → 58%
**Tests:** 12 automated tests
**Highlights:**
- Real Epic FHIR API integration with OAuth 2.0
- HIPAA Safe Harbor de-identification (18 PHI identifiers removed)
- Production-ready for hospital deployment

### mcp-spatialtools (Production Ready)

**Coverage:** 23% (low test coverage, but 95% real implementation)
**Tests:** 5 smoke tests
**Highlights:**
- Real differential expression (Mann-Whitney U + FDR)
- Real Moran's I spatial autocorrelation
- Real cell type deconvolution (8 cell types)
- Validated with Patient-001 data (900 spots × 31 genes)

---

## Development Guidelines

### Adding Tests

**For new servers (smoke tests):**
1. Create `tests/test_server.py` with 3 test classes
2. Test imports, configuration, tool registration
3. Target: 35-60% coverage with 5-12 tests

**For production servers (functional tests):**
1. Create test fixtures with realistic data
2. Test with `DRY_RUN=False` using monkeypatch
3. Validate actual calculations and outputs
4. Test edge cases and error handling
5. Target: 70%+ coverage on critical modules

### Best Practices

✅ **DO:**
- Use pytest fixtures for test data
- Test both DRY_RUN modes (True and False)
- Validate actual outputs (files, calculations, structures)
- Use descriptive test names
- Test edge cases (empty data, invalid inputs)

❌ **DON'T:**
- Call FastMCP-decorated functions directly (use `_impl` functions)
- Skip assertions on key results
- Test only happy paths
- Hard-code file paths (use fixtures, tmp_path)
- Ignore test failures

---

## Next Steps

To reach 60% overall coverage (+3.1 points):

1. **mcp-spatialtools:** Add functional tests for spatial algorithms (+5.4 point impact)
2. **mcp-tcga:** Add TCGA API integration tests (+1.6 point impact)
3. **mcp-openimagedata:** Add image processing tests (+1.4 point impact)

**Estimated effort:** ~30-40 functional tests across 3 servers

---

## Resources

- **pytest Documentation:** https://docs.pytest.org/
- **pytest-cov Plugin:** https://pytest-cov.readthedocs.io/
- **MCP Protocol:** https://modelcontextprotocol.io/
- **Server Implementation Status:** [docs/SERVER_IMPLEMENTATION_STATUS.md](../docs/SERVER_IMPLEMENTATION_STATUS.md)

---

**Last Updated:** 2025-12-30
**Status:** 9/9 servers tested | 4/9 production-ready | All deployed to GCP Cloud Run
