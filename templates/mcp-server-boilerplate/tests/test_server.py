"""Tests for mcp-{{SERVER_NAME}} server."""

import pytest
import os
import sys
from pathlib import Path

# Add src to path for testing
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


def test_imports():
    """Test that server module imports successfully."""
    try:
        from mcp_template import server
        assert server is not None
    except ImportError as e:
        pytest.fail(f"Failed to import server module: {e}")


def test_dry_run_mode():
    """Test DRY_RUN mode is enabled by default in test environment."""
    from mcp_template.server import DRY_RUN

    # In test/dev environment, DRY_RUN should default to true
    assert DRY_RUN is True, "DRY_RUN should be enabled by default"


def test_server_initialization():
    """Test FastMCP server initializes correctly."""
    from mcp_template.server import mcp

    assert mcp is not None
    assert mcp.name == "{{SERVER_NAME}}"


def test_environment_variables():
    """Test environment variables are read correctly."""
    from mcp_template.server import DATA_DIR, CACHE_DIR

    # Should have default values
    assert DATA_DIR is not None
    assert CACHE_DIR is not None
    assert "{{SERVER_NAME}}" in DATA_DIR.lower() or "{{SERVER_NAME}}" in CACHE_DIR.lower()


@pytest.mark.asyncio
async def test_example_tool_dry_run():
    """Test example tool in DRY_RUN mode."""
    from mcp_template.server import {{TOOL_EXAMPLE}}

    # Call the _impl method directly (FastMCP pattern for testing)
    result = await {{TOOL_EXAMPLE}}._impl(
        data_file="/fake/path/data.csv",
        param1="test_value",
        param2=10
    )

    # Verify DRY_RUN response
    assert result["status"] == "DRY_RUN"
    assert "message" in result
    assert result["result_key1"] is not None


@pytest.mark.asyncio
async def test_list_available_data_dry_run():
    """Test list_available_data tool in DRY_RUN mode."""
    from mcp_template.server import list_available_data

    result = await list_available_data._impl()

    assert result["status"] == "DRY_RUN"
    assert "files" in result
    assert result["file_count"] > 0
    assert isinstance(result["files"], list)


# Add more tests for each tool...
# Pattern:
# @pytest.mark.asyncio
# async def test_your_tool_name_dry_run():
#     """Test your_tool in DRY_RUN mode."""
#     from mcp_template.server import your_tool
#
#     result = await your_tool._impl(...)
#
#     assert result["status"] == "DRY_RUN"
#     # Add assertions for expected mock data


# Optional: Test with real data (requires fixtures)
# @pytest.mark.asyncio
# async def test_example_tool_with_fixtures(tmp_path):
#     """Test example tool with actual test data."""
#     # This test requires DRY_RUN=false and test fixtures
#     pytest.skip("Real data test - requires fixtures and DRY_RUN=false")
#
#     # Create test data
#     test_file = tmp_path / "test_data.csv"
#     # ... write test data ...
#
#     # Set DRY_RUN to false temporarily
#     import mcp_template.server
#     original_dry_run = mcp_template.server.DRY_RUN
#     mcp_template.server.DRY_RUN = False
#
#     try:
#         result = await {{TOOL_EXAMPLE}}._impl(
#             data_file=str(test_file),
#             param1="test",
#             param2=5
#         )
#
#         assert result["status"] == "success"
#         # Add assertions for real results
#
#     finally:
#         # Restore DRY_RUN
#         mcp_template.server.DRY_RUN = original_dry_run
