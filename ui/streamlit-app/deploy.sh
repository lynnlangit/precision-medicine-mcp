#!/bin/bash
# Deploy Streamlit MCP Chat UI to GCP Cloud Run

set -e

# Configuration
PROJECT_ID="precision-medicine-poc"
REGION="us-central1"
SERVICE_NAME="streamlit-mcp-chat"

# Check if ANTHROPIC_API_KEY is set
if [ -z "$ANTHROPIC_API_KEY" ]; then
    echo "Error: ANTHROPIC_API_KEY environment variable not set"
    echo "Usage: export ANTHROPIC_API_KEY=your_key_here && ./deploy.sh"
    exit 1
fi

# Check if GEMINI_API_KEY is set (optional)
if [ -n "$GEMINI_API_KEY" ]; then
    echo "Note: GEMINI_API_KEY is set - Gemini provider will be available"
    GEMINI_ENV="GEMINI_API_KEY=$GEMINI_API_KEY,"
else
    echo "Note: GEMINI_API_KEY not set - Only Claude provider will be available"
    GEMINI_ENV=""
fi

echo "=========================================="
echo "Deploying Streamlit MCP Chat to Cloud Run"
echo "=========================================="
echo "Project: $PROJECT_ID"
echo "Region: $REGION"
echo "Service: $SERVICE_NAME"
echo ""

# Deploy to Cloud Run
echo "Building and deploying..."
gcloud run deploy "$SERVICE_NAME" \
    --source . \
    --platform managed \
    --region "$REGION" \
    --project "$PROJECT_ID" \
    --allow-unauthenticated \
    --memory 1Gi \
    --cpu 1 \
    --min-instances 0 \
    --max-instances 5 \
    --timeout 300 \
    --set-env-vars ${GEMINI_ENV}ANTHROPIC_API_KEY="$ANTHROPIC_API_KEY",ENVIRONMENT=development \
    --port 8501 \
    --quiet

echo ""
echo "=========================================="
echo "✅ Deployment Complete!"
echo "=========================================="
echo ""

# Get the service URL
SERVICE_URL=$(gcloud run services describe "$SERVICE_NAME" \
    --region "$REGION" \
    --project "$PROJECT_ID" \
    --format 'value(status.url)')

echo "🌐 Your Streamlit app is live at:"
echo "   $SERVICE_URL"
echo ""
echo "📊 View logs:"
echo "   gcloud logging read \"resource.type=cloud_run_revision AND resource.labels.service_name=$SERVICE_NAME\" --limit=50 --project=$PROJECT_ID"
echo ""
echo "🔧 Manage service:"
echo "   https://console.cloud.google.com/run/detail/$REGION/$SERVICE_NAME?project=$PROJECT_ID"
echo ""
