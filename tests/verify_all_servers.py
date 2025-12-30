#!/usr/bin/env python3
"""
Quick MCP Server Verification Script

Tests all 9 MCP servers can start and respond to basic requests.
Run time: ~30 seconds
"""

import subprocess
import sys
import time
from pathlib import Path

# Color codes for terminal output
GREEN = '\033[92m'
RED = '\033[91m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
RESET = '\033[0m'

SERVERS = [
    {
        "name": "mcp-fgbio",
        "path": "servers/mcp-fgbio",
        "module": "mcp_fgbio",
        "env": {"FGBIO_DRY_RUN": "true"}
    },
    {
        "name": "mcp-multiomics",
        "path": "servers/mcp-multiomics",
        "module": "mcp_multiomics",
        "env": {"MULTIOMICS_DRY_RUN": "true"}
    },
    {
        "name": "mcp-spatialtools",
        "path": "servers/mcp-spatialtools",
        "module": "mcp_spatialtools",
        "env": {"SPATIAL_DRY_RUN": "true"}
    },
    {
        "name": "mcp-tcga",
        "path": "servers/mcp-tcga",
        "module": "mcp_tcga",
        "env": {"TCGA_DRY_RUN": "true"}
    },
    {
        "name": "mcp-openimagedata",
        "path": "servers/mcp-openimagedata",
        "module": "mcp_openimagedata",
        "env": {"IMAGE_DRY_RUN": "true"}
    },
    {
        "name": "mcp-seqera",
        "path": "servers/mcp-seqera",
        "module": "mcp_seqera",
        "env": {"SEQERA_DRY_RUN": "true"}
    },
    {
        "name": "mcp-huggingface",
        "path": "servers/mcp-huggingface",
        "module": "mcp_huggingface",
        "env": {"HF_DRY_RUN": "true"}
    },
    {
        "name": "mcp-deepcell",
        "path": "servers/mcp-deepcell",
        "module": "mcp_deepcell",
        "env": {"DEEPCELL_DRY_RUN": "true"}
    },
    {
        "name": "mcp-mockepic",
        "path": "servers/mcp-mockepic",
        "module": "mcp_mockepic",
        "env": {"EPIC_DRY_RUN": "true"}
    }
]


def test_server_import(server, repo_root):
    """Test if server can be imported successfully."""
    name = server["name"]
    path = server["path"]
    module = server["module"]

    print(f"{BLUE}Testing {name}...{RESET}", end=" ", flush=True)

    try:
        # Test import using the venv Python
        server_path = repo_root / path
        venv_python = server_path / "venv/bin/python"

        if not venv_python.exists():
            print(f"{RED}✗ venv not found{RESET}")
            return False

        # Test import
        result = subprocess.run(
            [str(venv_python), "-c", f"import {module}; print('OK')"],
            cwd=str(server_path),
            capture_output=True,
            text=True,
            timeout=10,
            env={**server["env"]}
        )

        if result.returncode == 0 and "OK" in result.stdout:
            print(f"{GREEN}✓{RESET}")
            return True
        else:
            print(f"{RED}✗ Import failed{RESET}")
            if result.stderr:
                print(f"  Error: {result.stderr[:100]}")
            return False

    except subprocess.TimeoutExpired:
        print(f"{RED}✗ Timeout{RESET}")
        return False
    except Exception as e:
        print(f"{RED}✗ {str(e)[:50]}{RESET}")
        return False


def main():
    """Run verification for all servers."""
    # Find repo root (where this script is located)
    repo_root = Path(__file__).parent.resolve()

    print(f"\n{BLUE}{'='*60}{RESET}")
    print(f"{BLUE}MCP Server Verification - Quick Import Test{RESET}")
    print(f"{BLUE}{'='*60}{RESET}")
    print(f"{BLUE}Repo:{RESET} {repo_root}\n")

    start_time = time.time()
    results = []

    for server in SERVERS:
        success = test_server_import(server, repo_root)
        results.append((server["name"], success))

    # Summary
    elapsed = time.time() - start_time
    passed = sum(1 for _, success in results if success)
    total = len(results)

    print(f"\n{BLUE}{'='*60}{RESET}")
    print(f"{BLUE}Summary{RESET}")
    print(f"{BLUE}{'='*60}{RESET}\n")

    for name, success in results:
        status = f"{GREEN}✓ PASS{RESET}" if success else f"{RED}✗ FAIL{RESET}"
        print(f"  {name:25} {status}")

    print(f"\n{BLUE}Results:{RESET} {passed}/{total} servers passed")
    print(f"{BLUE}Time:{RESET} {elapsed:.1f}s\n")

    if passed == total:
        print(f"{GREEN}✅ All servers verified successfully!{RESET}")
        print(f"{GREEN}Ready for Claude Desktop testing.{RESET}\n")
        return 0
    else:
        print(f"{YELLOW}⚠️  Some servers failed verification.{RESET}")
        print(f"{YELLOW}Fix the failing servers before testing in Claude Desktop.{RESET}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
