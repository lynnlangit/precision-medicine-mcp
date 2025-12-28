"""Tests for HIPAA de-identification functions."""

import pytest
from datetime import datetime
from mcp_epic.server import (
    deidentify_patient,
    deidentify_resource,
    deidentify_bundle,
)


class TestDeidentifyPatient:
    """Tests for patient de-identification."""

    def test_removes_name(self):
        """Test that patient names are removed."""
        patient = {
            "resourceType": "Patient",
            "id": "patient-001",
            "name": [{"family": "Doe", "given": ["John"]}],
            "gender": "male",
        }
        result = deidentify_patient(patient)
        assert "name" not in result
        assert result["_deidentified"] is True

    def test_removes_telecom(self):
        """Test that phone numbers and emails are removed."""
        patient = {
            "resourceType": "Patient",
            "id": "patient-001",
            "telecom": [
                {"system": "phone", "value": "555-1234"},
                {"system": "email", "value": "john@example.com"},
            ],
        }
        result = deidentify_patient(patient)
        assert "telecom" not in result

    def test_removes_address(self):
        """Test that addresses are removed."""
        patient = {
            "resourceType": "Patient",
            "id": "patient-001",
            "address": [
                {"line": ["123 Main St"], "city": "Boston", "state": "MA"}
            ],
        }
        result = deidentify_patient(patient)
        assert "address" not in result

    def test_hashes_patient_id(self):
        """Test that patient ID is hashed."""
        patient = {
            "resourceType": "Patient",
            "id": "patient-001",
        }
        result = deidentify_patient(patient)
        assert result["id"] != "patient-001"
        assert result["id"].startswith("deidentified-")

    def test_hashes_identifier_values(self):
        """Test that identifier values are hashed."""
        patient = {
            "resourceType": "Patient",
            "id": "patient-001",
            "identifier": [{"system": "http://hospital.org", "value": "MRN123"}],
        }
        result = deidentify_patient(patient)
        assert result["identifier"][0]["value"].startswith("HASH-")
        assert "MRN123" not in result["identifier"][0]["value"]

    def test_shifts_birthdate_to_year_only(self):
        """Test that birth dates are shifted to year only."""
        patient = {
            "resourceType": "Patient",
            "id": "patient-001",
            "birthDate": "1968-03-15",
        }
        result = deidentify_patient(patient)
        assert result["birthDate"] == "1968"

    def test_preserves_gender(self):
        """Test that gender is preserved (not PHI per HIPAA)."""
        patient = {
            "resourceType": "Patient",
            "id": "patient-001",
            "gender": "female",
        }
        result = deidentify_patient(patient)
        assert result["gender"] == "female"

    def test_adds_deidentification_metadata(self):
        """Test that de-identification metadata is added."""
        patient = {
            "resourceType": "Patient",
            "id": "patient-001",
        }
        result = deidentify_patient(patient)
        assert result["_deidentified"] is True
        assert result["_deidentification_method"] == "HIPAA Safe Harbor"
        assert "_deidentification_date" in result


class TestDeidentifyResource:
    """Tests for general resource de-identification."""

    def test_hashes_resource_id(self):
        """Test that resource IDs are hashed."""
        resource = {
            "resourceType": "Condition",
            "id": "condition-123",
        }
        result = deidentify_resource(resource)
        assert result["id"] != "condition-123"
        assert result["id"].startswith("deidentified-")

    def test_hashes_patient_reference(self):
        """Test that patient references are hashed."""
        resource = {
            "resourceType": "Condition",
            "id": "condition-123",
            "subject": {
                "reference": "Patient/patient-001",
                "display": "John Doe",
            },
        }
        result = deidentify_resource(resource)
        assert "patient-001" not in result["subject"]["reference"]
        assert result["subject"]["reference"].startswith("Patient/deidentified-")
        assert "display" not in result["subject"]

    def test_removes_dates(self):
        """Test that dates are removed or shifted to year only."""
        resource = {
            "resourceType": "Observation",
            "id": "obs-123",
            "issued": "2024-12-01T14:30:00Z",
            "recordedDate": "2024-12-01",
        }
        result = deidentify_resource(resource)
        assert result["issued"] == "2024"
        assert result["recordedDate"] == "2024"

    def test_preserves_clinical_data(self):
        """Test that clinical data (non-PHI) is preserved."""
        resource = {
            "resourceType": "Observation",
            "id": "obs-123",
            "code": {"text": "CA-125"},
            "valueQuantity": {"value": 487, "unit": "U/mL"},
        }
        result = deidentify_resource(resource)
        assert result["code"] == {"text": "CA-125"}
        assert result["valueQuantity"] == {"value": 487, "unit": "U/mL"}


class TestDeidentifyBundle:
    """Tests for FHIR Bundle de-identification."""

    def test_deidentifies_all_entries(self):
        """Test that all bundle entries are de-identified."""
        bundle = {
            "resourceType": "Bundle",
            "type": "searchset",
            "total": 2,
            "entry": [
                {
                    "resource": {
                        "resourceType": "Condition",
                        "id": "condition-1",
                        "subject": {"reference": "Patient/patient-001"},
                    }
                },
                {
                    "resource": {
                        "resourceType": "Condition",
                        "id": "condition-2",
                        "subject": {"reference": "Patient/patient-001"},
                    }
                },
            ],
        }
        result = deidentify_bundle(bundle)

        # Check both entries are de-identified
        for entry in result["entry"]:
            resource = entry["resource"]
            assert resource["_deidentified"] is True
            assert "patient-001" not in resource["subject"]["reference"]
            assert resource["id"].startswith("deidentified-")

    def test_preserves_bundle_structure(self):
        """Test that bundle structure is preserved."""
        bundle = {
            "resourceType": "Bundle",
            "type": "searchset",
            "total": 1,
            "entry": [
                {
                    "resource": {
                        "resourceType": "Observation",
                        "id": "obs-1",
                    }
                }
            ],
        }
        result = deidentify_bundle(bundle)
        assert result["resourceType"] == "Bundle"
        assert result["type"] == "searchset"
        assert result["total"] == 1

    def test_handles_empty_bundle(self):
        """Test that empty bundles are handled correctly."""
        bundle = {
            "resourceType": "Bundle",
            "type": "searchset",
            "total": 0,
        }
        result = deidentify_bundle(bundle)
        assert result["total"] == 0
        assert "entry" not in result or result["entry"] == []


class TestDeidentificationConsistency:
    """Tests for de-identification consistency."""

    def test_same_id_produces_same_hash(self):
        """Test that the same ID always produces the same hash."""
        patient1 = {"resourceType": "Patient", "id": "patient-001"}
        patient2 = {"resourceType": "Patient", "id": "patient-001"}

        result1 = deidentify_patient(patient1)
        result2 = deidentify_patient(patient2)

        assert result1["id"] == result2["id"]

    def test_different_ids_produce_different_hashes(self):
        """Test that different IDs produce different hashes."""
        patient1 = {"resourceType": "Patient", "id": "patient-001"}
        patient2 = {"resourceType": "Patient", "id": "patient-002"}

        result1 = deidentify_patient(patient1)
        result2 = deidentify_patient(patient2)

        assert result1["id"] != result2["id"]
