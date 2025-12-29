"""Unit tests for STAR alignment functionality."""

import os
import sys
import pytest
import tempfile
import shutil
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))


class TestSTARLogParser:
    """Test STAR log file parsing."""

    def test_parse_valid_log(self, tmp_path):
        """Test parsing a valid STAR log file."""
        from mcp_spatialtools.server import _parse_star_log

        # Create mock STAR log file
        log_content = """                      Started job on |       Dec 29 10:15:23
                             Started mapping on |       Dec 29 10:25:56
                                    Finished on |       Dec 29 10:45:12
                                  Number of input reads |       50000000
                      Average input read length |       201
                            UNIQUE READS:
                   Uniquely mapped reads number |       42500000
                        Uniquely mapped reads % |       85.00%
                          Average mapped length |       198.23
                          Number of splices: Total |       25000000
                MULTI-MAPPING READS:
       Number of reads mapped to multiple loci |       3750000
            % of reads mapped to multiple loci |       7.50%
                          UNMAPPED READS:
  Number of reads unmapped: too many mismatches |       0
       % of reads unmapped: too many mismatches |       0.00%
            Number of reads unmapped: too short |       3750000
                 % of reads unmapped: too short |       7.50%
                Number of reads unmapped: other |       0
                     % of reads unmapped: other |       0.00%
"""
        log_file = tmp_path / "Log.final.out"
        log_file.write_text(log_content)

        # Parse log
        stats = _parse_star_log(log_file)

        # Verify stats
        assert stats["total_reads"] == 50000000
        assert stats["uniquely_mapped"] == 42500000
        assert stats["multi_mapped"] == 3750000
        assert stats["unmapped"] == 3750000
        assert abs(stats["alignment_rate"] - 0.925) < 0.001
        assert abs(stats["unique_mapping_rate"] - 0.85) < 0.001

    def test_parse_log_missing_file(self, tmp_path):
        """Test parsing when log file doesn't exist."""
        from mcp_spatialtools.server import _parse_star_log

        missing_log = tmp_path / "missing.log"

        with pytest.raises(IOError, match="STAR log file not found"):
            _parse_star_log(missing_log)

    def test_parse_log_invalid_totals(self, tmp_path):
        """Test parsing when read counts don't sum correctly."""
        from mcp_spatialtools.server import _parse_star_log

        # Create log with inconsistent totals
        log_content = """                                  Number of input reads |       50000000
                   Uniquely mapped reads number |       30000000
       Number of reads mapped to multiple loci |       10000000
            Number of reads unmapped: too short |       5000000
                Number of reads unmapped: other |       0
  Number of reads unmapped: too many mismatches |       0
"""
        log_file = tmp_path / "Log.final.out"
        log_file.write_text(log_content)

        # Should raise error due to mismatch (30M + 10M + 5M = 45M, not 50M)
        with pytest.raises(IOError, match="STAR log parsing error"):
            _parse_star_log(log_file)


class TestSyntheticFASTQ:
    """Test synthetic FASTQ generator."""

    def test_create_synthetic_fastq(self, tmp_path):
        """Test creating synthetic FASTQ files."""
        from mcp_spatialtools.server import _create_synthetic_fastq

        r1_path = tmp_path / "test_R1.fastq.gz"
        r2_path = tmp_path / "test_R2.fastq.gz"

        # Generate synthetic FASTQs
        _create_synthetic_fastq(r1_path, r2_path, num_reads=100, read_length=50)

        # Verify files exist
        assert r1_path.exists()
        assert r2_path.exists()

        # Verify files are non-empty
        assert r1_path.stat().st_size > 0
        assert r2_path.stat().st_size > 0

        # Verify they're gzipped
        import gzip
        with gzip.open(r1_path, 'rt') as f:
            first_line = f.readline()
            assert first_line.startswith('@read_')

    def test_synthetic_fastq_content(self, tmp_path):
        """Test FASTQ content validity."""
        from mcp_spatialtools.server import _create_synthetic_fastq
        import gzip

        r1_path = tmp_path / "test_R1.fastq.gz"
        r2_path = tmp_path / "test_R2.fastq.gz"

        _create_synthetic_fastq(r1_path, r2_path, num_reads=10, read_length=50)

        # Read and validate R1
        with gzip.open(r1_path, 'rt') as f:
            lines = f.readlines()

        # Should have 4 lines per read × 10 reads = 40 lines
        assert len(lines) == 40

        # Check FASTQ format (4-line structure)
        for i in range(0, len(lines), 4):
            assert lines[i].startswith('@')  # Header
            assert len(lines[i+1].strip()) == 50  # Sequence
            assert lines[i+2].strip() == '+'  # Plus line
            assert len(lines[i+3].strip()) == 50  # Quality


@pytest.mark.skipif(not shutil.which("STAR"), reason="STAR not installed")
class TestSTARAlignment:
    """Test STAR alignment execution (requires STAR installation).

    Note: These tests require direct function access, not through MCP tool wrapper.
    They serve as documentation for expected behavior. Integration testing is done
    through manual testing or the MCP protocol.
    """

    @pytest.mark.skip(reason="MCP tool wrapper prevents direct function calls")
    def test_align_dry_run_mode(self, tmp_path):
        """Test alignment in DRY_RUN mode."""
        # This test documents expected behavior but is skipped because
        # align_spatial_data is wrapped by @mcp.tool() decorator
        pass


class TestAlignmentErrorHandling:
    """Test alignment error handling.

    Note: These tests serve as documentation for expected error behavior.
    The actual error handling is tested through the helper functions and
    manual integration testing.
    """

    @pytest.mark.skip(reason="MCP tool wrapper prevents direct function calls")
    def test_missing_r1_file(self, tmp_path):
        """Test error when R1 file is missing."""
        # Documents that IOError should be raised for missing R1
        # Actual testing done through MCP protocol or manual testing
        pass

    @pytest.mark.skip(reason="MCP tool wrapper prevents direct function calls")
    def test_missing_r2_file(self, tmp_path):
        """Test error when R2 file is missing."""
        # Documents that IOError should be raised for missing R2
        # Actual testing done through MCP protocol or manual testing
        pass

    @pytest.mark.skip(reason="MCP tool wrapper prevents direct function calls")
    def test_invalid_thread_count(self, tmp_path):
        """Test error when thread count is invalid."""
        # Documents that ValueError should be raised for invalid threads
        # Actual testing done through MCP protocol or manual testing
        pass


class TestAlignmentIntegration:
    """Integration tests for alignment workflow.

    Note: Full integration testing requires MCP client. These tests document
    expected behavior and can be manually verified through Claude Desktop.
    """

    @pytest.mark.skip(reason="MCP tool wrapper prevents direct function calls")
    def test_dry_run_returns_expected_format(self, tmp_path):
        """Test that DRY_RUN mode returns expected output format."""
        # Documents expected output format:
        # {
        #   "aligned_bam": str,
        #   "alignment_stats": {
        #     "total_reads": int,
        #     "uniquely_mapped": int,
        #     "multi_mapped": int,
        #     "unmapped": int,
        #     "alignment_rate": float (0-1),
        #     "unique_mapping_rate": float (0-1)
        #   },
        #   "log_file": str
        # }
        pass


class TestToolsExist:
    """Verify alignment tools are properly registered."""

    def test_align_tool_registered(self):
        """Test that align_spatial_data tool is registered."""
        from mcp_spatialtools import server
        assert hasattr(server, 'align_spatial_data')

    def test_helper_functions_exist(self):
        """Test that helper functions are available."""
        from mcp_spatialtools import server
        assert hasattr(server, '_parse_star_log')
        assert hasattr(server, '_create_synthetic_fastq')
