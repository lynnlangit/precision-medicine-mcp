#!/bin/bash
# Quick deploy script without special character issues
set -e

cd "$(dirname "$0")"

# Load API key from .env
if [ ! -f .env ]; then
    echo "Error: .env file not found"
    exit 1
fi

# Source the .env file
set -a
source .env
set +a

# Verify API key is set
if [ -z "$ANTHROPIC_API_KEY" ]; then
    echo "Error: ANTHROPIC_API_KEY not set in .env"
    exit 1
fi

echo "=========================================="
echo "Deploying Streamlit to Cloud Run"
echo "=========================================="
echo "API Key length: ${#ANTHROPIC_API_KEY} characters"
echo ""

# Deploy
gcloud run deploy streamlit-mcp-chat \
    --source . \
    --platform managed \
    --region us-central1 \
    --project precision-medicine-poc \
    --allow-unauthenticated \
    --memory 1Gi \
    --cpu 1 \
    --min-instances 0 \
    --max-instances 5 \
    --timeout 300 \
    --set-env-vars "ANTHROPIC_API_KEY=${ANTHROPIC_API_KEY},ENVIRONMENT=development" \
    --port 8501 \
    --quiet

echo ""
echo "✅ Deployment complete!"

# Get URL
SERVICE_URL=$(gcloud run services describe streamlit-mcp-chat \
    --region us-central1 \
    --project precision-medicine-poc \
    --format 'value(status.url)')

echo "🌐 Streamlit URL: $SERVICE_URL"
