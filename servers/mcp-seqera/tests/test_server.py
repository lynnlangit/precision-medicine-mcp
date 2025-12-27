"""Smoke tests for mcp-seqera server."""

import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))


class TestServerImport:
    def test_import_server_module(self):
        from mcp_seqera import server
        assert server is not None

    def test_mcp_server_exists(self):
        from mcp_seqera import server
        assert hasattr(server, 'mcp')

    def test_main_function_exists(self):
        from mcp_seqera import server
        assert callable(server.main)


class TestConfiguration:
    def test_dry_run_variable_exists(self):
        from mcp_seqera import server
        assert hasattr(server, 'DRY_RUN')
        assert isinstance(server.DRY_RUN, bool)

    def test_seqera_token_exists(self):
        from mcp_seqera import server
        assert hasattr(server, 'SEQERA_TOKEN')


class TestTools:
    def test_tools_registered(self):
        from mcp_seqera import server
        # Server should have tools registered
        assert hasattr(server, 'launch_nextflow_pipeline')
        assert hasattr(server, 'monitor_workflow_status')
        assert hasattr(server, 'list_available_pipelines')
