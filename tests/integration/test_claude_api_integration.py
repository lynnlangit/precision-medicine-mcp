#!/usr/bin/env python3
"""
Test GCP-deployed MCP servers using Claude API

Prerequisites:
    pip install anthropic
    export ANTHROPIC_API_KEY=your_key_here

Usage:
    python test_claude_api_integration.py
"""

import anthropic
import os
import sys
from pathlib import Path

# Color codes
GREEN = '\033[92m'
RED = '\033[91m'
BLUE = '\033[94m'
YELLOW = '\033[93m'
NC = '\033[0m'


def load_deployment_urls():
    """Load deployed server URLs from infrastructure/deployment_urls.txt"""
    urls_file = Path(__file__).parent.parent.parent / "infrastructure" / "deployment_urls.txt"

    if not urls_file.exists():
        print(f"{RED}Error: infrastructure/deployment_urls.txt not found{NC}")
        print("Run ./scripts/deployment/deploy_to_gcp.sh first")
        sys.exit(1)

    servers = {}
    with open(urls_file) as f:
        for line in f:
            if '=' in line:
                name, url = line.strip().split('=', 1)
                servers[name] = url

    return servers


def create_mcp_config(servers):
    """Create MCP server configuration for Claude API"""
    mcp_servers = []

    for name, url in servers.items():
        # Convert server name to friendly name (remove mcp- prefix)
        friendly_name = name.replace('mcp-', '')

        mcp_servers.append({
            "type": "url",
            "url": f"{url}/sse",
            "name": friendly_name,
            # Add authentication if needed:
            # "authorization_token": f"Bearer {get_gcp_token()}"
        })

    return mcp_servers


def test_server_connection(client, server_name):
    """Test a single MCP server via Claude API"""
    print(f"{BLUE}Testing {server_name}...{NC}", end=" ", flush=True)

    try:
        # Simple test query
        response = client.beta.messages.create(
            model="claude-opus-4-5",
            max_tokens=500,
            messages=[{
                "role": "user",
                "content": f"List the available tools from the {server_name} MCP server."
            }],
            mcp_servers=[{
                "type": "url",
                "url": load_deployment_urls()[f"mcp-{server_name}"] + "/sse",
                "name": server_name,
            }],
            tools=[{
                "type": "mcp_toolset",
                "mcp_server_name": server_name
            }],
            betas=["mcp-client-2025-11-20"]
        )

        if response.content:
            print(f"{GREEN}✓ PASS{NC}")
            return True
        else:
            print(f"{RED}✗ FAIL (no response){NC}")
            return False

    except Exception as e:
        print(f"{RED}✗ FAIL{NC}")
        print(f"  Error: {str(e)[:100]}")
        return False


def main():
    """Main test function"""
    print(f"\n{BLUE}{'='*60}{NC}")
    print(f"{BLUE}Claude API Integration Test - GCP MCP Servers{NC}")
    print(f"{BLUE}{'='*60}{NC}\n")

    # Check for API key
    if not os.getenv("ANTHROPIC_API_KEY"):
        print(f"{RED}Error: ANTHROPIC_API_KEY not set{NC}")
        print("Set it with: export ANTHROPIC_API_KEY=your_key_here")
        sys.exit(1)

    # Load deployed servers
    servers = load_deployment_urls()
    print(f"{BLUE}Found {len(servers)} deployed servers{NC}\n")

    # Initialize Claude API client
    client = anthropic.Anthropic()

    # Test each server
    results = {}
    for full_name in servers.keys():
        server_name = full_name.replace('mcp-', '')
        results[server_name] = test_server_connection(client, server_name)

    # Summary
    print(f"\n{BLUE}{'='*60}{NC}")
    print(f"{BLUE}Summary{NC}")
    print(f"{BLUE}{'='*60}{NC}\n")

    passed = sum(1 for v in results.values() if v)
    failed = sum(1 for v in results.values() if not v)

    for server, success in results.items():
        status = f"{GREEN}✓ PASS{NC}" if success else f"{RED}✗ FAIL{NC}"
        print(f"  {server:20} {status}")

    print(f"\n{BLUE}Results:{NC} {passed}/{len(results)} servers passed\n")

    if passed == len(results):
        print(f"{GREEN}✅ All servers working with Claude API!{NC}\n")
        return 0
    else:
        print(f"{YELLOW}⚠️  Some servers failed. Check logs above.{NC}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
