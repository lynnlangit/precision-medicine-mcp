"""Smoke tests for mcp-openimagedata server."""

import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))


class TestServerImport:
    def test_import_server_module(self):
        from mcp_openimagedata import server
        assert server is not None

    def test_mcp_server_exists(self):
        from mcp_openimagedata import server
        assert hasattr(server, 'mcp')

    def test_main_function_exists(self):
        from mcp_openimagedata import server
        assert callable(server.main)


class TestConfiguration:
    def test_dry_run_variable_exists(self):
        from mcp_openimagedata import server
        assert hasattr(server, 'DRY_RUN')
        assert isinstance(server.DRY_RUN, bool)


class TestTools:
    def test_tools_registered(self):
        from mcp_openimagedata import server
        assert hasattr(server, 'fetch_histology_image')
        assert hasattr(server, 'register_image_to_spatial')
        assert hasattr(server, 'extract_image_features')
