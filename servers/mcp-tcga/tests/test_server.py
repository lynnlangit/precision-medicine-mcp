"""Smoke tests for mcp-tcga server."""

import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))


class TestServerImport:
    def test_import_server_module(self):
        from mcp_tcga import server
        assert server is not None

    def test_mcp_server_exists(self):
        from mcp_tcga import server
        assert hasattr(server, 'mcp')

    def test_main_function_exists(self):
        from mcp_tcga import server
        assert callable(server.main)


class TestConfiguration:
    def test_dry_run_variable_exists(self):
        from mcp_tcga import server
        assert hasattr(server, 'DRY_RUN')
        assert isinstance(server.DRY_RUN, bool)


class TestTools:
    def test_tools_registered(self):
        from mcp_tcga import server
        # TCGA server has 5 tools
        assert hasattr(server, 'query_tcga_cohorts')
        assert hasattr(server, 'fetch_expression_data')
        assert hasattr(server, 'compare_to_cohort')
        assert hasattr(server, 'get_survival_data')
        assert hasattr(server, 'get_mutation_data')
