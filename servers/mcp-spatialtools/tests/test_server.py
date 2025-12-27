"""Smoke tests for mcp-spatialtools server."""

import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))


class TestServerImport:
    def test_import_server_module(self):
        from mcp_spatialtools import server
        assert server is not None

    def test_mcp_server_exists(self):
        from mcp_spatialtools import server
        assert hasattr(server, 'mcp')

    def test_main_function_exists(self):
        from mcp_spatialtools import server
        assert callable(server.main)


class TestConfiguration:
    def test_dry_run_variable_exists(self):
        from mcp_spatialtools import server
        assert hasattr(server, 'DRY_RUN')
        assert isinstance(server.DRY_RUN, bool)


class TestTools:
    def test_tools_registered(self):
        from mcp_spatialtools import server
        # Spatialtools has 8 tools
        assert hasattr(server, 'filter_quality')
        assert hasattr(server, 'split_by_region')
        assert hasattr(server, 'align_spatial_data')
        assert hasattr(server, 'merge_tiles')
        assert hasattr(server, 'calculate_spatial_autocorrelation')
        assert hasattr(server, 'perform_differential_expression')
        assert hasattr(server, 'perform_batch_correction')
        assert hasattr(server, 'perform_pathway_enrichment')
