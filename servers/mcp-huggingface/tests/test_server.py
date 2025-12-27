"""Smoke tests for mcp-huggingface server."""

import os
import sys
import pytest

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))


class TestServerImport:
    """Test server can be imported and initialized."""

    def test_import_server_module(self):
        """Test server module imports successfully."""
        from mcp_huggingface import server
        assert server is not None

    def test_mcp_server_exists(self):
        """Test MCP server object is created."""
        from mcp_huggingface import server
        assert hasattr(server, 'mcp')
        assert server.mcp is not None

    def test_main_function_exists(self):
        """Test main function is defined."""
        from mcp_huggingface import server
        assert hasattr(server, 'main')
        assert callable(server.main)


class TestConfiguration:
    """Test server configuration."""

    def test_dry_run_variable_exists(self):
        """Test DRY_RUN variable is defined."""
        from mcp_huggingface import server
        assert hasattr(server, 'DRY_RUN')
        assert isinstance(server.DRY_RUN, bool)

    def test_hf_token_variable_exists(self):
        """Test HF_TOKEN variable is defined."""
        from mcp_huggingface import server
        assert hasattr(server, 'HF_TOKEN')
        assert isinstance(server.HF_TOKEN, str)

    def test_dry_run_default_true(self):
        """Test DRY_RUN defaults to True."""
        os.environ.pop('HF_DRY_RUN', None)

        import importlib
        from mcp_huggingface import server
        importlib.reload(server)

        assert server.DRY_RUN is True

    def test_dry_run_respects_env_false(self):
        """Test DRY_RUN respects environment variable."""
        os.environ['HF_DRY_RUN'] = 'false'

        import importlib
        from mcp_huggingface import server
        importlib.reload(server)

        assert server.DRY_RUN is False

        # Cleanup
        os.environ['HF_DRY_RUN'] = 'true'


class TestTools:
    """Test tool registration."""

    def test_load_genomic_model_exists(self):
        """Test load_genomic_model tool is registered."""
        from mcp_huggingface import server
        assert hasattr(server, 'load_genomic_model')

    def test_predict_cell_type_exists(self):
        """Test predict_cell_type tool is registered."""
        from mcp_huggingface import server
        assert hasattr(server, 'predict_cell_type')

    def test_embed_sequences_exists(self):
        """Test embed_sequences tool is registered."""
        from mcp_huggingface import server
        assert hasattr(server, 'embed_sequences')


class TestResources:
    """Test resource registration."""

    def test_get_dnabert_info_exists(self):
        """Test get_dnabert_info resource is registered."""
        from mcp_huggingface import server
        assert hasattr(server, 'get_dnabert_info')

    def test_get_geneformer_info_exists(self):
        """Test get_geneformer_info resource is registered."""
        from mcp_huggingface import server
        assert hasattr(server, 'get_geneformer_info')
