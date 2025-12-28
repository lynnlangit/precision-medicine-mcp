"""Tests for FHIR tools."""

import pytest
from unittest.mock import AsyncMock, patch, MagicMock
from mcp_epic.server import (
    get_patient_demographics,
    get_patient_conditions,
    get_patient_observations,
    get_patient_medications,
)


@pytest.fixture
def mock_patient_response():
    """Mock patient FHIR resource."""
    return {
        "resourceType": "Patient",
        "id": "patient-001",
        "name": [{"family": "TestPatient", "given": ["Jane"]}],
        "gender": "female",
        "birthDate": "1968-03-15",
        "telecom": [{"system": "phone", "value": "555-0123"}],
        "address": [{"city": "TestCity", "state": "IA"}],
    }


@pytest.fixture
def mock_conditions_bundle():
    """Mock conditions search bundle."""
    return {
        "resourceType": "Bundle",
        "type": "searchset",
        "total": 1,
        "entry": [
            {
                "resource": {
                    "resourceType": "Condition",
                    "id": "condition-001",
                    "code": {"text": "Stage IV Ovarian Cancer"},
                    "subject": {
                        "reference": "Patient/patient-001",
                        "display": "Jane TestPatient",
                    },
                    "clinicalStatus": {
                        "coding": [{"code": "active", "display": "Active"}]
                    },
                }
            }
        ],
    }


@pytest.fixture
def mock_observations_bundle():
    """Mock observations search bundle."""
    return {
        "resourceType": "Bundle",
        "type": "searchset",
        "total": 2,
        "entry": [
            {
                "resource": {
                    "resourceType": "Observation",
                    "id": "obs-ca125",
                    "code": {"text": "CA-125"},
                    "subject": {"reference": "Patient/patient-001"},
                    "valueQuantity": {"value": 487, "unit": "U/mL"},
                    "issued": "2024-12-01T14:30:00Z",
                }
            },
            {
                "resource": {
                    "resourceType": "Observation",
                    "id": "obs-brca",
                    "code": {"text": "BRCA1/2"},
                    "subject": {"reference": "Patient/patient-001"},
                    "valueCodeableConcept": {"text": "Negative"},
                }
            },
        ],
    }


@pytest.fixture
def mock_medications_bundle():
    """Mock medications search bundle."""
    return {
        "resourceType": "Bundle",
        "type": "searchset",
        "total": 1,
        "entry": [
            {
                "resource": {
                    "resourceType": "MedicationStatement",
                    "id": "med-carboplatin",
                    "medicationCodeableConcept": {"text": "Carboplatin"},
                    "subject": {"reference": "Patient/patient-001"},
                    "status": "completed",
                    "dateAsserted": "2022-12-05",
                }
            }
        ],
    }


class TestGetPatientDemographics:
    """Tests for get_patient_demographics tool."""

    @pytest.mark.asyncio
    async def test_success(self, mock_patient_response):
        """Test successful patient retrieval."""
        with patch("mcp_epic.server.fhir_client") as mock_client:
            mock_client.get = AsyncMock(return_value=mock_patient_response)

            result = await get_patient_demographics("patient-001")

            assert result["status"] == "success"
            assert "data" in result
            assert result["data"]["_deidentified"] is True
            assert "name" not in result["data"]  # PHI removed
            assert "gender" in result["data"]  # Non-PHI preserved
            mock_client.get.assert_called_once_with("Patient/patient-001")

    @pytest.mark.asyncio
    async def test_http_error(self):
        """Test handling of HTTP errors."""
        with patch("mcp_epic.server.fhir_client") as mock_client:
            error = MagicMock()
            error.response.status_code = 404
            mock_client.get = AsyncMock(side_effect=Exception("Not found"))

            result = await get_patient_demographics("patient-999")

            assert result["status"] == "error"
            assert "error" in result
            assert "message" in result

    @pytest.mark.asyncio
    async def test_deidentification_applied(self, mock_patient_response):
        """Test that de-identification is applied."""
        with patch("mcp_epic.server.fhir_client") as mock_client:
            mock_client.get = AsyncMock(return_value=mock_patient_response)

            result = await get_patient_demographics("patient-001")

            data = result["data"]
            # PHI should be removed
            assert "name" not in data
            assert "telecom" not in data
            assert "address" not in data
            # ID should be hashed
            assert data["id"].startswith("deidentified-")
            # Birth date should be year only
            assert data["birthDate"] == "1968"


class TestGetPatientConditions:
    """Tests for get_patient_conditions tool."""

    @pytest.mark.asyncio
    async def test_success(self, mock_conditions_bundle):
        """Test successful conditions retrieval."""
        with patch("mcp_epic.server.fhir_client") as mock_client:
            mock_client.search = AsyncMock(return_value=mock_conditions_bundle)

            result = await get_patient_conditions("patient-001")

            assert result["status"] == "success"
            assert result["count"] == 1
            assert "data" in result
            mock_client.search.assert_called_once_with(
                "Condition", {"patient": "patient-001"}
            )

    @pytest.mark.asyncio
    async def test_multiple_conditions(self):
        """Test retrieval of multiple conditions."""
        bundle = {
            "resourceType": "Bundle",
            "total": 3,
            "entry": [
                {"resource": {"resourceType": "Condition", "id": f"cond-{i}"}}
                for i in range(3)
            ],
        }

        with patch("mcp_epic.server.fhir_client") as mock_client:
            mock_client.search = AsyncMock(return_value=bundle)

            result = await get_patient_conditions("patient-001")

            assert result["count"] == 3
            assert len(result["data"]["entry"]) == 3

    @pytest.mark.asyncio
    async def test_no_conditions(self):
        """Test patient with no conditions."""
        empty_bundle = {
            "resourceType": "Bundle",
            "type": "searchset",
            "total": 0,
        }

        with patch("mcp_epic.server.fhir_client") as mock_client:
            mock_client.search = AsyncMock(return_value=empty_bundle)

            result = await get_patient_conditions("patient-001")

            assert result["status"] == "success"
            assert result["count"] == 0


class TestGetPatientObservations:
    """Tests for get_patient_observations tool."""

    @pytest.mark.asyncio
    async def test_success(self, mock_observations_bundle):
        """Test successful observations retrieval."""
        with patch("mcp_epic.server.fhir_client") as mock_client:
            mock_client.search = AsyncMock(return_value=mock_observations_bundle)

            result = await get_patient_observations("patient-001")

            assert result["status"] == "success"
            assert result["count"] == 2
            mock_client.search.assert_called_once_with(
                "Observation", {"patient": "patient-001"}
            )

    @pytest.mark.asyncio
    async def test_with_category_filter(self, mock_observations_bundle):
        """Test observations with category filter."""
        with patch("mcp_epic.server.fhir_client") as mock_client:
            mock_client.search = AsyncMock(return_value=mock_observations_bundle)

            result = await get_patient_observations("patient-001", category="laboratory")

            assert result["status"] == "success"
            mock_client.search.assert_called_once_with(
                "Observation", {"patient": "patient-001", "category": "laboratory"}
            )

    @pytest.mark.asyncio
    async def test_clinical_data_preserved(self, mock_observations_bundle):
        """Test that clinical data is preserved after de-identification."""
        with patch("mcp_epic.server.fhir_client") as mock_client:
            mock_client.search = AsyncMock(return_value=mock_observations_bundle)

            result = await get_patient_observations("patient-001")

            # Check first observation
            obs = result["data"]["entry"][0]["resource"]
            assert obs["code"]["text"] == "CA-125"
            assert obs["valueQuantity"]["value"] == 487
            assert obs["valueQuantity"]["unit"] == "U/mL"


class TestGetPatientMedications:
    """Tests for get_patient_medications tool."""

    @pytest.mark.asyncio
    async def test_success(self, mock_medications_bundle):
        """Test successful medications retrieval."""
        with patch("mcp_epic.server.fhir_client") as mock_client:
            mock_client.search = AsyncMock(return_value=mock_medications_bundle)

            result = await get_patient_medications("patient-001")

            assert result["status"] == "success"
            assert result["count"] == 1
            mock_client.search.assert_called_once_with(
                "MedicationStatement", {"patient": "patient-001"}
            )

    @pytest.mark.asyncio
    async def test_medication_data_preserved(self, mock_medications_bundle):
        """Test that medication data is preserved."""
        with patch("mcp_epic.server.fhir_client") as mock_client:
            mock_client.search = AsyncMock(return_value=mock_medications_bundle)

            result = await get_patient_medications("patient-001")

            med = result["data"]["entry"][0]["resource"]
            assert med["medicationCodeableConcept"]["text"] == "Carboplatin"
            assert med["status"] == "completed"

    @pytest.mark.asyncio
    async def test_active_medications(self):
        """Test retrieval of active medications."""
        bundle = {
            "resourceType": "Bundle",
            "total": 2,
            "entry": [
                {
                    "resource": {
                        "resourceType": "MedicationStatement",
                        "id": "med-1",
                        "status": "active",
                        "medicationCodeableConcept": {"text": "Bevacizumab"},
                    }
                },
                {
                    "resource": {
                        "resourceType": "MedicationStatement",
                        "id": "med-2",
                        "status": "completed",
                        "medicationCodeableConcept": {"text": "Carboplatin"},
                    }
                },
            ],
        }

        with patch("mcp_epic.server.fhir_client") as mock_client:
            mock_client.search = AsyncMock(return_value=bundle)

            result = await get_patient_medications("patient-001")

            assert result["count"] == 2
            statuses = [
                entry["resource"]["status"] for entry in result["data"]["entry"]
            ]
            assert "active" in statuses
            assert "completed" in statuses
