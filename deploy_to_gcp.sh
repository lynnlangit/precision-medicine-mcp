#!/usr/bin/env bash

# GCP Cloud Run Deployment Script for All 9 MCP Servers
# Prerequisites:
#   - gcloud CLI installed and authenticated
#   - Project configured with billing enabled
#   - Cloud Run API enabled

set -e  # Exit on error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
PROJECT_ID="${GCP_PROJECT_ID:-precision-medicine-poc}"
REGION="${GCP_REGION:-us-central1}"
REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"

# Server configurations (bash 3 compatible)
SERVERS=(
    "mcp-fgbio:3000"
    "mcp-multiomics:3001"
    "mcp-spatialtools:3002"
    "mcp-tcga:3003"
    "mcp-openimagedata:3004"
    "mcp-seqera:3005"
    "mcp-huggingface:3006"
    "mcp-deepcell:3007"
    "mcp-mockepic:3008"
)

# Functions
print_header() {
    echo -e "${BLUE}============================================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}============================================================${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

print_info() {
    echo -e "${YELLOW}→ $1${NC}"
}

# Check prerequisites
check_prerequisites() {
    print_header "Checking Prerequisites"

    # Check if gcloud is installed
    if ! command -v gcloud &> /dev/null; then
        print_error "gcloud CLI not found. Install from: https://cloud.google.com/sdk/docs/install"
        exit 1
    fi
    print_success "gcloud CLI found"

    # Check if authenticated
    if ! gcloud auth list --filter=status:ACTIVE --format="value(account)" &> /dev/null; then
        print_error "Not authenticated. Run: gcloud auth login"
        exit 1
    fi
    print_success "Authenticated to GCP"

    # Set project
    gcloud config set project "${PROJECT_ID}"
    print_success "Project set to: ${PROJECT_ID}"

    # Enable required APIs
    print_info "Enabling required APIs..."
    gcloud services enable \
        cloudbuild.googleapis.com \
        run.googleapis.com \
        artifactregistry.googleapis.com \
        --quiet
    print_success "APIs enabled"

    echo ""
}

# Deploy a single server
deploy_server() {
    local server_name=$1
    local port=$2
    local server_path="${REPO_ROOT}/servers/${server_name}"

    print_header "Deploying ${server_name}"

    # Check if Dockerfile exists
    if [ ! -f "${server_path}/Dockerfile" ]; then
        print_error "Dockerfile not found at ${server_path}/Dockerfile"
        return 1
    fi

    print_info "Building and deploying from ${server_path}..."

    # Deploy to Cloud Run
    # Note: Build from repo root to access shared/utils
    # Dockerfile is at servers/${server_name}/Dockerfile
    cd "${REPO_ROOT}"
    gcloud run deploy "${server_name}" \
        --source "." \
        --dockerfile "servers/${server_name}/Dockerfile" \
        --platform managed \
        --region "${REGION}" \
        --allow-unauthenticated \
        --memory 2Gi \
        --cpu 2 \
        --min-instances 0 \
        --max-instances 10 \
        --timeout 300 \
        --quiet
    cd - > /dev/null

    if [ $? -eq 0 ]; then
        # Get the service URL
        SERVICE_URL=$(gcloud run services describe "${server_name}" \
            --platform managed \
            --region "${REGION}" \
            --format 'value(status.url)')

        print_success "Deployed: ${SERVICE_URL}"
        echo "${server_name}=${SERVICE_URL}" >> "${REPO_ROOT}/deployment_urls.txt"
    else
        print_error "Failed to deploy ${server_name}"
        return 1
    fi

    echo ""
}

# Test deployed server
test_server() {
    local server_name=$1
    local service_url=$2

    print_info "Testing ${server_name} at ${service_url}..."

    # Test health endpoint
    if curl -s -f "${service_url}/health" > /dev/null 2>&1; then
        print_success "${server_name} is healthy"
        return 0
    else
        print_error "${server_name} health check failed"
        return 1
    fi
}

# Main deployment function
main() {
    print_header "GCP Cloud Run Deployment - 9 MCP Servers"
    echo -e "${BLUE}Project: ${PROJECT_ID}${NC}"
    echo -e "${BLUE}Region: ${REGION}${NC}"
    echo ""

    # Ask for confirmation
    read -p "Deploy all 9 servers to GCP Cloud Run? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_info "Deployment cancelled"
        exit 0
    fi

    # Check prerequisites
    check_prerequisites

    # Clear previous deployment URLs
    > "${REPO_ROOT}/deployment_urls.txt"

    # Deploy all servers
    print_header "Deploying All Servers"

    DEPLOYED=0
    FAILED=0

    for server_config in "${SERVERS[@]}"; do
        server_name="${server_config%%:*}"
        port="${server_config##*:}"
        if deploy_server "$server_name" "$port"; then
            ((DEPLOYED++))
        else
            ((FAILED++))
        fi
    done

    echo ""
    print_header "Deployment Summary"
    echo -e "${GREEN}Deployed: ${DEPLOYED}/9${NC}"
    if [ $FAILED -gt 0 ]; then
        echo -e "${RED}Failed: ${FAILED}/9${NC}"
    fi
    echo ""

    # Display all URLs
    if [ -f "${REPO_ROOT}/deployment_urls.txt" ]; then
        print_header "Deployed Server URLs"
        cat "${REPO_ROOT}/deployment_urls.txt"
        echo ""
    fi

    # Test deployed servers
    print_header "Health Check Tests"
    while IFS='=' read -r server url; do
        test_server "$server" "$url" || true
    done < "${REPO_ROOT}/deployment_urls.txt"

    echo ""
    print_header "Next Steps"
    echo -e "${YELLOW}1. Copy deployment_urls.txt URLs for your Claude API configuration${NC}"
    echo -e "${YELLOW}2. Update authentication tokens (if needed)${NC}"
    echo -e "${YELLOW}3. Test servers with: curl \${SERVICE_URL}/health${NC}"
    echo -e "${YELLOW}4. See GCP_TESTING_GUIDE.md for Claude API integration${NC}"
    echo ""

    print_success "Deployment complete!"
}

# Run main function
main "$@"
