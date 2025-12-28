# Phase 1B: Epic FHIR Integration

**Timeline:** Week 2-3
**Budget Focus:** Use Epic sandbox (free) for development
**Goal:** Replace mcp-mockepic with real Epic FHIR integration

---

## Overview

Build mcp-epic server to retrieve real patient clinical data from Epic EHR via FHIR API, with automatic de-identification for research use.

---

## Prerequisites

- [x] Phase 1A GCP infrastructure completed
- [ ] Epic on FHIR credentials obtained
- [ ] Access to Epic sandbox environment
- [ ] Python 3.11+ development environment

---

## Task 1: Epic FHIR Access Setup (60 mins)

### 1.1 Register Application in Epic

**Option A: Use Epic Sandbox (Recommended for Development)**

1. Go to https://fhir.epic.com/
2. Create account if needed
3. Register a new app:
   - **App Name:** Precision Medicine MCP
   - **Redirect URI:** `http://localhost:8080/callback`
   - **FHIR Version:** R4
   - **Scopes:**
     - `patient/*.read`
     - `Patient.read`
     - `Observation.read`
     - `Condition.read`
     - `MedicationStatement.read`
     - `Procedure.read`

4. Save credentials:
   - Client ID
   - Client Secret (if provided)
   - FHIR endpoint: `https://fhir.epic.com/interconnect-fhir-oauth/api/FHIR/R4/`

**Option B: Hospital Epic Instance (Production)**

Work with hospital IT to:
1. Register app in Epic App Orchard
2. Obtain production credentials
3. Configure OAuth flow
4. Get FHIR endpoint URL

### 1.2 Store Credentials in GCP Secret Manager

```bash
# Source GCP environment
source .env.gcp

# Update Epic credentials in Secret Manager
echo -n "YOUR_EPIC_CLIENT_ID" | \
  gcloud secrets versions add epic-client-id --data-file=- --project=$GCP_PROJECT_ID

echo -n "YOUR_EPIC_CLIENT_SECRET" | \
  gcloud secrets versions add epic-client-secret --data-file=- --project=$GCP_PROJECT_ID

echo -n "https://fhir.epic.com/interconnect-fhir-oauth/api/FHIR/R4/" | \
  gcloud secrets versions add epic-fhir-endpoint --data-file=- --project=$GCP_PROJECT_ID
```

---

## Task 2: Create mcp-epic Server (4-6 hours)

### 2.1 Server Structure

```bash
cd /Users/lynnlangit/Documents/GitHub/spatial-mcp/servers

# Create new server directory
mkdir -p mcp-epic/{src/mcp_epic,tests}
cd mcp-epic

# Create virtual environment
python3.11 -m venv venv
source venv/bin/activate

# Install dependencies
pip install fastmcp fhirclient google-cloud-secret-manager google-cloud-healthcare pydantic pydantic-settings
```

### 2.2 Create pyproject.toml

```toml
[build-system]
requires = ["setuptools>=61.0", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "mcp-epic"
version = "0.1.0"
description = "MCP server for Epic EHR FHIR integration"
authors = [
    {name = "Precision Medicine MCP Team"}
]
readme = "README.md"
requires-python = ">=3.11"
license = {text = "Apache-2.0"}
dependencies = [
    "fastmcp>=0.2.0",
    "fhirclient>=4.1.0",
    "pydantic>=2.0.0",
    "pydantic-settings>=2.0.0",
    "google-cloud-secret-manager>=2.16.0",
    "google-cloud-healthcare>=1.16.0",
    "httpx>=0.27.0",
]

[project.urls]
Homepage = "https://github.com/lynnlangit/precision-medicine-mcp"
Issues = "https://github.com/lynnlangit/precision-medicine-mcp/issues"
Documentation = "https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs"

[project.optional-dependencies]
dev = [
    "pytest>=8.0.0",
    "pytest-asyncio>=0.23.0",
    "pytest-cov>=4.1.0",
    "pytest-mock>=3.12.0",
    "black>=24.0.0",
    "ruff>=0.3.0",
    "mypy>=1.8.0",
]

[tool.setuptools.packages.find]
where = ["src"]

[tool.pytest.ini_options]
testpaths = ["tests"]
asyncio_mode = "auto"
addopts = "-v --tb=short --cov=src/mcp_epic --cov-report=term-missing --cov-report=html"

[tool.black]
line-length = 100
target-version = ["py311"]

[tool.ruff]
line-length = 100
target-version = "py311"
```

### 2.3 Create Server Implementation

`src/mcp_epic/server.py`:

```python
"""MCP Epic server - Epic EHR FHIR integration with de-identification."""

import json
import logging
import os
from datetime import datetime, timedelta
from typing import Any, Dict, List, Optional

from fastmcp import FastMCP
from fhirclient import client
from fhirclient.models import (
    patient,
    observation,
    condition,
    medicationstatement,
    procedure,
)
from google.cloud import secretmanager
from pydantic import BaseModel, Field

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize MCP server
mcp = FastMCP("epic")


# ============================================================================
# CONFIGURATION
# ============================================================================

class EpicConfig(BaseModel):
    """Epic FHIR configuration."""

    project_id: str = Field(default_factory=lambda: os.getenv("GCP_PROJECT_ID", ""))
    client_id_secret: str = Field(
        default_factory=lambda: os.getenv(
            "EPIC_CLIENT_ID_SECRET",
            "projects/precision-medicine-poc/secrets/epic-client-id",
        )
    )
    client_secret_secret: str = Field(
        default_factory=lambda: os.getenv(
            "EPIC_CLIENT_SECRET_SECRET",
            "projects/precision-medicine-poc/secrets/epic-client-secret",
        )
    )
    fhir_endpoint_secret: str = Field(
        default_factory=lambda: os.getenv(
            "EPIC_FHIR_ENDPOINT_SECRET",
            "projects/precision-medicine-poc/secrets/epic-fhir-endpoint",
        )
    )
    deidentify: bool = Field(
        default_factory=lambda: os.getenv("EPIC_DEIDENTIFY", "true").lower() == "true"
    )


config = EpicConfig()


# ============================================================================
# SECRET MANAGEMENT
# ============================================================================

def get_secret(secret_name: str) -> str:
    """Retrieve secret from GCP Secret Manager."""
    try:
        client = secretmanager.SecretManagerServiceClient()
        response = client.access_secret_version(request={"name": f"{secret_name}/versions/latest"})
        return response.payload.data.decode("UTF-8")
    except Exception as e:
        logger.error(f"Failed to retrieve secret {secret_name}: {e}")
        raise


def get_fhir_client() -> client.FHIRClient:
    """Create authenticated FHIR client."""
    settings = {
        "app_id": "precision-medicine-mcp",
        "api_base": get_secret(config.fhir_endpoint_secret),
    }

    # Create client
    smart = client.FHIRClient(settings=settings)

    # Note: In production, implement proper OAuth flow
    # For sandbox/development, may use client credentials
    return smart


# ============================================================================
# DE-IDENTIFICATION
# ============================================================================

def deidentify_patient_data(data: Dict[str, Any]) -> Dict[str, Any]:
    """Apply HIPAA Safe Harbor de-identification."""
    if not config.deidentify:
        return data

    deidentified = data.copy()

    # Remove HIPAA identifiers
    identifiers_to_remove = [
        "name",
        "telecom",
        "address",
        "photo",
        "contact",
        "managingOrganization",
    ]

    for identifier in identifiers_to_remove:
        if identifier in deidentified:
            del deidentified[identifier]

    # Hash patient ID
    if "id" in deidentified:
        import hashlib

        deidentified["id"] = hashlib.sha256(deidentified["id"].encode()).hexdigest()[:16]

    # Date shift birthDate (keep year, remove month/day)
    if "birthDate" in deidentified:
        try:
            birth_date = datetime.fromisoformat(deidentified["birthDate"])
            deidentified["birthDate"] = str(birth_date.year)
        except:
            deidentified["birthDate"] = "REDACTED"

    # Add de-identification notice
    deidentified["_deidentified"] = True
    deidentified["_deidentification_method"] = "HIPAA Safe Harbor"

    return deidentified


# ============================================================================
# TOOLS
# ============================================================================

@mcp.tool()
async def get_patient_demographics(patient_id: str) -> Dict[str, Any]:
    """Retrieve patient demographics from Epic.

    Args:
        patient_id: Epic patient identifier

    Returns:
        Patient demographics (de-identified if configured)

    Example:
        >>> result = await get_patient_demographics("eM5CWtq15N0WJeuCet5bJlQ3")
    """
    try:
        smart = get_fhir_client()
        patient_resource = patient.Patient.read(patient_id, smart.server)

        patient_data = patient_resource.as_json()

        # De-identify
        if config.deidentify:
            patient_data = deidentify_patient_data(patient_data)

        return {
            "status": "success",
            "patient": patient_data,
            "retrieved_at": datetime.now().isoformat(),
        }

    except Exception as e:
        logger.error(f"Failed to retrieve patient {patient_id}: {e}")
        return {"status": "error", "message": str(e)}


@mcp.tool()
async def get_patient_conditions(patient_id: str) -> Dict[str, Any]:
    """Retrieve patient conditions/diagnoses from Epic.

    Args:
        patient_id: Epic patient identifier

    Returns:
        List of conditions with codes and dates

    Example:
        >>> result = await get_patient_conditions("eM5CWtq15N0WJeuCet5bJlQ3")
    """
    try:
        smart = get_fhir_client()

        # Search for conditions
        search = condition.Condition.where(struct={"patient": patient_id})
        conditions = search.perform_resources(smart.server)

        condition_list = []
        for cond in conditions:
            cond_json = cond.as_json()

            # Extract key information
            condition_data = {
                "code": cond_json.get("code", {}).get("coding", [{}])[0].get("display", "Unknown"),
                "code_system": cond_json.get("code", {}).get("coding", [{}])[0].get("system", ""),
                "onset_date": cond_json.get("onsetDateTime", ""),
                "clinical_status": cond_json.get("clinicalStatus", {})
                .get("coding", [{}])[0]
                .get("code", ""),
            }

            if config.deidentify:
                # Hash condition ID
                import hashlib

                condition_data["id"] = hashlib.sha256(
                    cond_json.get("id", "").encode()
                ).hexdigest()[:16]

            condition_list.append(condition_data)

        return {
            "status": "success",
            "patient_id": patient_id if not config.deidentify else "REDACTED",
            "conditions": condition_list,
            "count": len(condition_list),
            "retrieved_at": datetime.now().isoformat(),
        }

    except Exception as e:
        logger.error(f"Failed to retrieve conditions for {patient_id}: {e}")
        return {"status": "error", "message": str(e)}


@mcp.tool()
async def get_patient_medications(patient_id: str) -> Dict[str, Any]:
    """Retrieve patient medications from Epic.

    Args:
        patient_id: Epic patient identifier

    Returns:
        List of medications with dosages and dates

    Example:
        >>> result = await get_patient_medications("eM5CWtq15N0WJeuCet5bJlQ3")
    """
    try:
        smart = get_fhir_client()

        # Search for medication statements
        search = medicationstatement.MedicationStatement.where(struct={"patient": patient_id})
        medications = search.perform_resources(smart.server)

        medication_list = []
        for med in medications:
            med_json = med.as_json()

            # Extract key information
            medication_data = {
                "medication": med_json.get("medicationCodeableConcept", {})
                .get("coding", [{}])[0]
                .get("display", "Unknown"),
                "status": med_json.get("status", ""),
                "effective_period": med_json.get("effectivePeriod", {}),
                "dosage": med_json.get("dosage", [{}])[0].get("text", "") if med_json.get("dosage") else "",
            }

            if config.deidentify:
                import hashlib

                medication_data["id"] = hashlib.sha256(
                    med_json.get("id", "").encode()
                ).hexdigest()[:16]

            medication_list.append(medication_data)

        return {
            "status": "success",
            "patient_id": patient_id if not config.deidentify else "REDACTED",
            "medications": medication_list,
            "count": len(medication_list),
            "retrieved_at": datetime.now().isoformat(),
        }

    except Exception as e:
        logger.error(f"Failed to retrieve medications for {patient_id}: {e}")
        return {"status": "error", "message": str(e)}


@mcp.tool()
async def get_patient_observations(
    patient_id: str, category: Optional[str] = None, limit: int = 100
) -> Dict[str, Any]:
    """Retrieve patient observations (labs, vitals, etc.) from Epic.

    Args:
        patient_id: Epic patient identifier
        category: Filter by category (e.g., "laboratory", "vital-signs")
        limit: Maximum number of observations to return

    Returns:
        List of observations with values and dates

    Example:
        >>> result = await get_patient_observations("eM5CWtq15N0WJeuCet5bJlQ3", category="laboratory")
    """
    try:
        smart = get_fhir_client()

        # Build search parameters
        search_params = {"patient": patient_id, "_count": limit}
        if category:
            search_params["category"] = category

        # Search for observations
        search = observation.Observation.where(struct=search_params)
        observations = search.perform_resources(smart.server)

        observation_list = []
        for obs in observations:
            obs_json = obs.as_json()

            # Extract key information
            observation_data = {
                "code": obs_json.get("code", {}).get("coding", [{}])[0].get("display", "Unknown"),
                "value": obs_json.get("valueQuantity", {}).get("value", ""),
                "unit": obs_json.get("valueQuantity", {}).get("unit", ""),
                "date": obs_json.get("effectiveDateTime", ""),
                "category": obs_json.get("category", [{}])[0]
                .get("coding", [{}])[0]
                .get("display", ""),
            }

            if config.deidentify:
                import hashlib

                observation_data["id"] = hashlib.sha256(
                    obs_json.get("id", "").encode()
                ).hexdigest()[:16]

            observation_list.append(observation_data)

        return {
            "status": "success",
            "patient_id": patient_id if not config.deidentify else "REDACTED",
            "observations": observation_list,
            "count": len(observation_list),
            "retrieved_at": datetime.now().isoformat(),
        }

    except Exception as e:
        logger.error(f"Failed to retrieve observations for {patient_id}: {e}")
        return {"status": "error", "message": str(e)}


# ============================================================================
# SERVER ENTRYPOINT
# ============================================================================

def main() -> None:
    """Run the MCP Epic server."""
    logger.info("Starting mcp-epic server...")

    if config.deidentify:
        logger.warning("=" * 80)
        logger.warning("⚠️  DE-IDENTIFICATION ENABLED")
        logger.warning("⚠️  All patient data will be de-identified per HIPAA Safe Harbor")
        logger.warning("=" * 80)
    else:
        logger.warning("=" * 80)
        logger.warning("⚠️  DE-IDENTIFICATION DISABLED")
        logger.warning("⚠️  Patient data will contain PHI - handle accordingly")
        logger.warning("=" * 80)

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
```

### 2.4 Create __init__.py and __main__.py

`src/mcp_epic/__init__.py`:
```python
"""MCP Epic server."""

from .server import main

__version__ = "0.1.0"
```

`src/mcp_epic/__main__.py`:
```python
if __name__ == "__main__":
    from .server import main
    main()
```

---

## Task 3: Testing (2-3 hours)

### 3.1 Create Test Files

`tests/test_epic_server.py`:

```python
"""Tests for mcp-epic server."""

import pytest
from unittest.mock import Mock, patch
from mcp_epic.server import (
    get_patient_demographics,
    get_patient_conditions,
    deidentify_patient_data,
)


@pytest.fixture
def mock_fhir_client():
    """Mock FHIR client."""
    with patch("mcp_epic.server.get_fhir_client") as mock:
        yield mock


@pytest.fixture
def sample_patient_data():
    """Sample patient FHIR data."""
    return {
        "resourceType": "Patient",
        "id": "test-patient-123",
        "name": [{"family": "Doe", "given": ["John"]}],
        "birthDate": "1970-01-15",
        "gender": "male",
        "address": [{"city": "Boston", "state": "MA"}],
        "telecom": [{"system": "phone", "value": "555-1234"}],
    }


class TestDeidentification:
    """Test de-identification functions."""

    def test_deidentify_removes_name(self, sample_patient_data):
        """Test that name is removed."""
        result = deidentify_patient_data(sample_patient_data)
        assert "name" not in result

    def test_deidentify_removes_address(self, sample_patient_data):
        """Test that address is removed."""
        result = deidentify_patient_data(sample_patient_data)
        assert "address" not in result

    def test_deidentify_removes_telecom(self, sample_patient_data):
        """Test that telecom is removed."""
        result = deidentify_patient_data(sample_patient_data)
        assert "telecom" not in result

    def test_deidentify_shifts_date(self, sample_patient_data):
        """Test that birthDate is shifted to year only."""
        result = deidentify_patient_data(sample_patient_data)
        assert result["birthDate"] == "1970"

    def test_deidentify_hashes_id(self, sample_patient_data):
        """Test that ID is hashed."""
        result = deidentify_patient_data(sample_patient_data)
        assert result["id"] != "test-patient-123"
        assert len(result["id"]) == 16  # SHA256 hash truncated

    def test_deidentify_adds_notice(self, sample_patient_data):
        """Test that de-identification notice is added."""
        result = deidentify_patient_data(sample_patient_data)
        assert result["_deidentified"] is True
        assert "Safe Harbor" in result["_deidentification_method"]


class TestPatientDemographics:
    """Test patient demographics retrieval."""

    @pytest.mark.asyncio
    async def test_get_patient_demographics_success(self, mock_fhir_client, sample_patient_data):
        """Test successful patient retrieval."""
        # Mock FHIR client response
        mock_patient = Mock()
        mock_patient.as_json.return_value = sample_patient_data

        with patch("mcp_epic.server.patient.Patient.read", return_value=mock_patient):
            result = await get_patient_demographics("test-patient-123")

            assert result["status"] == "success"
            assert result["patient"]["resourceType"] == "Patient"
            assert result["patient"]["_deidentified"] is True

    @pytest.mark.asyncio
    async def test_get_patient_demographics_error(self, mock_fhir_client):
        """Test error handling."""
        with patch("mcp_epic.server.patient.Patient.read", side_effect=Exception("FHIR error")):
            result = await get_patient_demographics("invalid-id")

            assert result["status"] == "error"
            assert "FHIR error" in result["message"]


# Run tests
if __name__ == "__main__":
    pytest.main([__file__, "-v", "--cov=mcp_epic", "--cov-report=term-missing"])
```

### 3.2 Run Tests

```bash
cd /Users/lynnlangit/Documents/GitHub/spatial-mcp/servers/mcp-epic

# Install package in editable mode
pip install -e ".[dev]"

# Run tests
pytest tests/ -v --cov=src/mcp_epic --cov-report=term-missing
```

---

## Task 4: Integration with Claude Desktop (30 mins)

### 4.1 Update Claude Desktop Config

Edit `~/Library/Application Support/Claude/claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "epic": {
      "command": "uv",
      "args": [
        "run",
        "--directory",
        "/Users/lynnlangit/Documents/GitHub/spatial-mcp/servers/mcp-epic",
        "python",
        "-m",
        "mcp_epic"
      ],
      "env": {
        "GCP_PROJECT_ID": "precision-medicine-poc",
        "EPIC_CLIENT_ID_SECRET": "projects/precision-medicine-poc/secrets/epic-client-id",
        "EPIC_CLIENT_SECRET_SECRET": "projects/precision-medicine-poc/secrets/epic-client-secret",
        "EPIC_FHIR_ENDPOINT_SECRET": "projects/precision-medicine-poc/secrets/epic-fhir-endpoint",
        "EPIC_DEIDENTIFY": "true",
        "GOOGLE_APPLICATION_CREDENTIALS": "/Users/lynnlangit/Documents/GitHub/spatial-mcp/mcp-server-key.json"
      }
    }
  }
}
```

### 4.2 Test in Claude Desktop

Restart Claude Desktop and test with:

```
Retrieve patient demographics for Epic test patient ID: eM5CWtq15N0WJeuCet5bJlQ3
```

---

## Task 5: Documentation (30 mins)

### 5.1 Create README

`servers/mcp-epic/README.md`:

```markdown
# mcp-epic: Epic EHR FHIR Integration

MCP server for retrieving patient clinical data from Epic EHR via FHIR API.

## Features

- ✅ Patient demographics retrieval
- ✅ Conditions/diagnoses
- ✅ Medications
- ✅ Observations (labs, vitals)
- ✅ Automatic HIPAA Safe Harbor de-identification
- ✅ GCP Secret Manager integration

## Configuration

Environment variables:
- `GCP_PROJECT_ID` - GCP project ID
- `EPIC_CLIENT_ID_SECRET` - Secret Manager path for Epic client ID
- `EPIC_CLIENT_SECRET_SECRET` - Secret Manager path for Epic client secret
- `EPIC_FHIR_ENDPOINT_SECRET` - Secret Manager path for FHIR endpoint
- `EPIC_DEIDENTIFY` - Enable de-identification (default: true)

## Usage

### Get Patient Demographics
```python
result = await get_patient_demographics("patient-id")
```

### Get Patient Conditions
```python
result = await get_patient_conditions("patient-id")
```

### Get Patient Medications
```python
result = await get_patient_medications("patient-id")
```

### Get Patient Observations
```python
result = await get_patient_observations("patient-id", category="laboratory")
```

## De-identification

When `EPIC_DEIDENTIFY=true` (default), the following HIPAA identifiers are removed:
- Names
- Addresses
- Phone numbers
- Email addresses
- Photos
- Medical record numbers

Dates are shifted to year only, and IDs are hashed.

## Testing

```bash
pytest tests/ -v --cov=src/mcp_epic
```

## Production Deployment

1. Obtain Epic production credentials
2. Update Secret Manager secrets
3. Set `EPIC_DEIDENTIFY=true`
4. Deploy to GCP Cloud Run or Compute Engine

## License

Apache 2.0
```

---

## Validation Checklist

- [ ] mcp-epic server created and installed
- [ ] Epic sandbox credentials obtained and stored
- [ ] De-identification tested and validated
- [ ] All tests passing (>70% coverage)
- [ ] Claude Desktop integration working
- [ ] Documentation complete
- [ ] Code committed to repository

---

## Next Steps

1. Test with Epic sandbox patient data
2. Validate de-identification meets HIPAA Safe Harbor
3. Obtain hospital Epic credentials for production
4. Proceed to [Phase 1C: Complete Spatial Analysis](PHASE1C_SPATIALTOOLS.md)

---

## Estimated Effort Summary

- Epic setup: 1 hour
- Server implementation: 4-6 hours
- Testing: 2-3 hours
- Integration: 30 mins
- Documentation: 30 mins

**Total: 8-11 hours**
