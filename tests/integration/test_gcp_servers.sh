#!/bin/bash

# Test all deployed GCP Cloud Run MCP servers
# Tests health endpoints and basic connectivity

set -e

# Color codes
GREEN='\033[0;32m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
URLS_FILE="${REPO_ROOT}/infrastructure/deployment_urls.txt"

print_header() {
    echo -e "${BLUE}============================================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}============================================================${NC}"
}

test_server() {
    local server_name=$1
    local url=$2

    echo -n "Testing ${server_name}... "

    # Test SSE endpoint - should return 405 Method Not Allowed for GET without SSE headers
    # but that still means the server is running
    http_code=$(curl -s -o /dev/null -w "%{http_code}" --max-time 10 "${url}/sse")

    if [ "$http_code" = "405" ] || [ "$http_code" = "200" ]; then
        echo -e "${GREEN}✓ PASS (HTTP $http_code)${NC}"
        return 0
    else
        echo -e "${RED}✗ FAIL (HTTP $http_code)${NC}"
        return 1
    fi
}

main() {
    print_header "GCP MCP Server Health Check"

    if [ ! -f "$URLS_FILE" ]; then
        echo -e "${RED}Error: deployment_urls.txt not found${NC}"
        echo "Run ./deploy_to_gcp.sh first"
        exit 1
    fi

    echo ""
    PASSED=0
    FAILED=0

    while IFS='=' read -r server url; do
        if test_server "$server" "$url"; then
            ((PASSED++))
        else
            ((FAILED++))
        fi
    done < "$URLS_FILE"

    echo ""
    print_header "Summary"
    echo -e "${GREEN}Passed: ${PASSED}${NC}"
    if [ $FAILED -gt 0 ]; then
        echo -e "${RED}Failed: ${FAILED}${NC}"
    fi
    echo ""

    if [ $FAILED -eq 0 ]; then
        echo -e "${GREEN}✅ All servers healthy!${NC}"
        return 0
    else
        echo -e "${RED}⚠️  Some servers failed health checks${NC}"
        return 1
    fi
}

main "$@"
