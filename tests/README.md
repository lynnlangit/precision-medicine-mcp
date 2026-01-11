# 🧪 Test Coverage - Precision Medicine MCP Servers

**167 automated tests** covering 40 tools across 9 MCP servers.

---

## 📊 Quick Stats

| Metric | Value | Status |
|--------|-------|--------|
| **Overall Coverage** | 56.9% | ⬆️ +27.5 points from baseline |
| **Automated Tests** | 167 | ✅ 9/9 servers tested |
| **Production Ready** | 4/9 servers | mcp-multiomics, mcp-fgbio, mcp-spatialtools, mcp-epic |
| **GCP Deployed** | 9/9 servers | ✅ All validated on Cloud Run |

---

## 📋 Coverage by Server

📋 **[See Server Status →](../servers/README.md#-server-status)** - Complete server status and implementation details

| Server | Coverage | Tests | Status | Production Ready |
|--------|----------|-------|--------|------------------|
| 🧬 **mcp-fgbio** | 77% | 29 | ✅ Complete | ✅ YES (95% real) |
| 🔬 **mcp-multiomics** | 68% | 91 | ✅ Complete | ✅ YES (85% real) |
| 🖼️ **mcp-deepcell** | 62% | 9 | ✅ Smoke | ❌ Mocked |
| 🤖 **mcp-huggingface** | 56% | 12 | ✅ Smoke | ❌ Mocked |
| ⚙️ **mcp-seqera** | 56% | 6 | ✅ Smoke | ❌ Mocked |
| 🏥 **mcp-epic** | 58% | 12 | ✅ Complete | ✅ YES (100% real) |
| 🖼️ **mcp-openimagedata** | 35% | 5 | ✅ Smoke | ⚠️ 30% real |
| 🧪 **mcp-tcga** | 35% | 5 | ✅ Smoke | ❌ Mocked |
| 📍 **mcp-spatialtools** | 23% | 5 | ✅ Smoke | ✅ YES (95% real) |

**Note:** Low test coverage ≠ low production readiness. mcp-spatialtools has 23% coverage but is 95% production-ready with real implementations.

---

## 📂 Test Organization

```
/tests/
├── unit/                  # Unit tests for all servers (pytest)
│   ├── mcp-fgbio/        # 29 tests (77% coverage)
│   ├── mcp-multiomics/   # 91 tests (68% coverage) ⭐ Most comprehensive
│   ├── mcp-spatialtools/ # 5 tests (production-ready)
│   └── mcp-epic/         # 12 tests (real Epic FHIR)
├── integration/          # GCP deployment validation
├── manual_testing/       # End-to-end test suites
│   ├── PatientOne-OvarianCancer/  # Full workflow (TESTS 1-5)
│   └── Solution-Testing/          # Server verification
└── verification/         # Server utilities
```

---

## 🚀 Running Tests

### Run All Tests for a Server

**Production servers:**
```bash
# From repository root
MULTIOMICS_DRY_RUN="true" servers/mcp-multiomics/venv/bin/python -m pytest tests/unit/mcp-multiomics/ -v
FGBIO_DRY_RUN="true" servers/mcp-fgbio/venv/bin/python -m pytest tests/unit/mcp-fgbio/ -v
SPATIAL_DRY_RUN="true" servers/mcp-spatialtools/venv/bin/python -m pytest tests/unit/mcp-spatialtools/ -v
```

**Mocked servers:**
```bash
DEEPCELL_DRY_RUN="true" servers/mcp-deepcell/venv/bin/python -m pytest tests/unit/mcp-deepcell/ -v
TCGA_DRY_RUN="true" servers/mcp-tcga/venv/bin/python -m pytest tests/unit/mcp-tcga/ -v
```

### With Coverage Report

```bash
MULTIOMICS_DRY_RUN="true" servers/mcp-multiomics/venv/bin/python -m pytest tests/unit/mcp-multiomics/ \
  --cov=servers/mcp-multiomics/src/mcp_multiomics --cov-report=term-missing -v
```

### Run Specific Test

```bash
# Specific file
MULTIOMICS_DRY_RUN="true" servers/mcp-multiomics/venv/bin/python -m pytest tests/unit/mcp-multiomics/test_preprocessing.py -v

# Specific test function
MULTIOMICS_DRY_RUN="true" servers/mcp-multiomics/venv/bin/python -m pytest \
  tests/unit/mcp-multiomics/test_preprocessing.py::TestValidateWithRealData::test_validate_with_real_rna_data -v
```

---

## 🧩 Test Types

### ✅ Smoke Tests (All Servers)
- Module imports and DRY_RUN configuration
- Tool registration and server initialization
- Location: `/tests/unit/mcp-{server}/test_server.py`

### 🔬 Functional Tests (Production Servers)
- Real data processing and statistical calculations
- Batch correction, normalization, edge cases
- Best example: `mcp-multiomics` - 91 tests with 580KB+ fixture data
- Location: `/tests/unit/mcp-{server}/` (e.g., `test_preprocessing.py`, `test_integration.py`)

### 🔗 Integration Tests
- Multi-server workflows (PatientOne)
- GCP Cloud Run deployment validation
- SSE transport and Claude API integration
- Location: `/tests/integration/`

---

## ☁️ GCP Cloud Run Testing

**Validation:** ✅ 9/9 servers passing functional tests via Claude API (2025-12-30)

```bash
cd tests/integration
python test_all_gcp_servers.py
```

**Coverage:**
- SSE transport validation
- Tool discovery via MCP protocol
- Response time < 5 seconds (excluding cold starts)

📖 **Server URLs:** [Deployment Status](../docs/deployment/DEPLOYMENT_STATUS.md)

---

## 🏥 PatientOne Test Suite

**Location:** `/tests/manual_testing/PatientOne-OvarianCancer/`

Complete end-to-end precision medicine workflow for Stage IV ovarian cancer:

- 🧬 **TEST_1** - Clinical data (Epic FHIR)
- 🔬 **TEST_2** - Multi-omics resistance (RNA/Protein/Phospho)
- 📍 **TEST_3** - Spatial transcriptomics (900 spots × 31 genes)
- 🖼️ **TEST_4** - Histology & imaging (H&E, MxIF)
- 🎯 **TEST_5** - Integration & treatment recommendations

**Modes:**
- **DRY_RUN** (default): Synthetic data demo (~$0.32, 25-35 min)
- **Real Data**: Your own patient data ([Configuration Guide](manual_testing/PatientOne-OvarianCancer/DATA_MODES_GUIDE.md))

📖 **[PatientOne Quick Start →](manual_testing/PatientOne-OvarianCancer/README.md)**

---

## 📚 Development Guidelines

### Adding Tests

**For new servers (smoke tests):**
1. Create `/tests/unit/mcp-{server}/test_server.py`
2. Test imports, configuration, tool registration
3. Target: 35-60% coverage with 5-12 tests

**For production servers (functional tests):**
1. Create fixtures in `/tests/unit/mcp-{server}/fixtures/` with realistic data
2. Test with `DRY_RUN=False` using monkeypatch
3. Validate actual calculations and outputs
4. Target: 70%+ coverage on critical modules

### Best Practices

**✅ DO:**
- Use pytest fixtures for test data
- Test both DRY_RUN modes (True/False)
- Validate actual outputs (files, calculations)
- Use descriptive test names

**❌ DON'T:**
- Call FastMCP-decorated functions directly (use `_impl` functions)
- Skip assertions on key results
- Hard-code file paths (use fixtures, tmp_path)
- Ignore test failures

---

## 🎯 Next Steps to 60% Coverage

1. **mcp-spatialtools:** Add functional tests for spatial algorithms (+5.4 point impact)
2. **mcp-tcga:** Add TCGA API integration tests (+1.6 point impact)
3. **mcp-openimagedata:** Add image processing tests (+1.4 point impact)

**Estimated effort:** ~30-40 functional tests across 3 servers

---

## 🔗 Resources

- 📖 [pytest Documentation](https://docs.pytest.org/)
- 📊 [pytest-cov Plugin](https://pytest-cov.readthedocs.io/)
- 🔌 [MCP Protocol](https://modelcontextprotocol.io/)
- ✅ [Server Implementation Status](../docs/SERVER_IMPLEMENTATION_STATUS.md)

---

**Last Updated:** 2026-01-11
**Status:** 167 tests | 56.9% coverage | 9/9 servers on GCP ✅
