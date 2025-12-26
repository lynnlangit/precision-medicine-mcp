"""Pytest configuration and fixtures."""

import os
from pathlib import Path

import pytest


@pytest.fixture(autouse=True, scope="function")
def setup_test_env(tmp_path, monkeypatch):
    """Set up test environment with temp directories."""
    # Create temp directories
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir(exist_ok=True)

    data_dir = tmp_path / "data"
    data_dir.mkdir(exist_ok=True)

    # Directly patch the config object's directories
    from mcp_multiomics.config import config

    monkeypatch.setattr(config, "cache_dir", cache_dir)
    monkeypatch.setattr(config, "data_dir", data_dir)
    # Set dry_run based on environment variable (default to True for tests)
    dry_run_env = os.getenv('MULTIOMICS_DRY_RUN', 'true')
    dry_run_value = dry_run_env.lower() not in ('false', '0', 'no', 'off', '')
    monkeypatch.setattr(config, "dry_run", dry_run_value)

    yield {
        "cache_dir": cache_dir,
        "data_dir": data_dir,
    }


@pytest.fixture
def fixture_dir():
    """Return path to test fixtures directory."""
    return Path(__file__).parent / "fixtures"


@pytest.fixture
def rna_path(fixture_dir):
    """Return path to RNA test data."""
    return str(fixture_dir / "sample_rna.csv")


@pytest.fixture
def protein_path(fixture_dir):
    """Return path to protein test data."""
    return str(fixture_dir / "sample_protein.csv")


@pytest.fixture
def phospho_path(fixture_dir):
    """Return path to phospho test data."""
    return str(fixture_dir / "sample_phospho.csv")


@pytest.fixture
def metadata_path(fixture_dir):
    """Return path to metadata."""
    return str(fixture_dir / "sample_metadata.csv")


@pytest.fixture
def sample_p_values():
    """Return sample p-values for Stouffer's method testing."""
    return {
        "rna": [0.001, 0.05, 0.3, 0.8],
        "protein": [0.002, 0.04, 0.25, 0.75],
        "phospho": [0.01, 0.06, 0.28, 0.82],
    }


@pytest.fixture
def sample_effect_sizes():
    """Return sample effect sizes for directionality testing."""
    return {
        "rna": [2.5, 1.2, -0.3, 0.1],
        "protein": [1.8, 1.5, -0.2, 0.15],
        "phospho": [1.2, 0.8, -0.4, -0.05],
    }
