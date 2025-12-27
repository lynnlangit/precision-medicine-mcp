"""Smoke tests for mcp-deepcell server."""

import os
import sys
import pytest

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))


class TestServerImport:
    """Test server can be imported and initialized."""

    def test_import_server_module(self):
        """Test server module imports successfully."""
        from mcp_deepcell import server
        assert server is not None

    def test_mcp_server_exists(self):
        """Test MCP server object is created."""
        from mcp_deepcell import server
        assert hasattr(server, 'mcp')
        assert server.mcp is not None

    def test_main_function_exists(self):
        """Test main function is defined."""
        from mcp_deepcell import server
        assert hasattr(server, 'main')
        assert callable(server.main)


class TestConfiguration:
    """Test server configuration."""

    def test_dry_run_variable_exists(self):
        """Test DRY_RUN variable is defined."""
        from mcp_deepcell import server
        assert hasattr(server, 'DRY_RUN')
        assert isinstance(server.DRY_RUN, bool)

    def test_dry_run_default_true(self):
        """Test DRY_RUN defaults to True."""
        # Clear environment
        os.environ.pop('DEEPCELL_DRY_RUN', None)

        # Reload module
        import importlib
        from mcp_deepcell import server
        importlib.reload(server)

        assert server.DRY_RUN is True

    def test_dry_run_respects_env_false(self):
        """Test DRY_RUN respects environment variable."""
        os.environ['DEEPCELL_DRY_RUN'] = 'false'

        import importlib
        from mcp_deepcell import server
        importlib.reload(server)

        assert server.DRY_RUN is False

        # Cleanup
        os.environ['DEEPCELL_DRY_RUN'] = 'true'


class TestTools:
    """Test tool registration."""

    def test_segment_cells_exists(self):
        """Test segment_cells tool is registered."""
        from mcp_deepcell import server
        assert hasattr(server, 'segment_cells')

    def test_classify_cell_states_exists(self):
        """Test classify_cell_states tool is registered."""
        from mcp_deepcell import server
        assert hasattr(server, 'classify_cell_states')


class TestResources:
    """Test resource registration."""

    def test_get_membrane_model_info_exists(self):
        """Test get_membrane_model_info resource is registered."""
        from mcp_deepcell import server
        assert hasattr(server, 'get_membrane_model_info')
