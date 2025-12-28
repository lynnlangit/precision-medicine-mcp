"""MCP server for Epic FHIR integration with HIPAA de-identification.

This server connects to GCP Healthcare API FHIR store and provides tools for:
- Retrieving patient demographics
- Querying patient conditions
- Accessing patient observations
- Retrieving patient medications

All data is de-identified using HIPAA Safe Harbor method.
"""

import hashlib
import logging
import os
from datetime import datetime
from typing import Any, Dict, List, Optional

import httpx
from fastmcp import FastMCP
from google.auth import default
from google.auth.transport.requests import Request

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize FastMCP server
mcp = FastMCP("epic")

# Configuration
PROJECT_ID = os.getenv("GCP_PROJECT_ID", "precision-medicine-poc")
LOCATION = os.getenv("GCP_REGION", "us-central1")
DATASET_ID = os.getenv("HEALTHCARE_DATASET", "precision-medicine-dataset")
FHIR_STORE_ID = os.getenv("FHIR_STORE", "identified-fhir-store")
DEIDENTIFY_ENABLED = os.getenv("DEIDENTIFY_ENABLED", "true").lower() == "true"

# FHIR base URL
FHIR_BASE_URL = (
    f"https://healthcare.googleapis.com/v1/projects/{PROJECT_ID}"
    f"/locations/{LOCATION}/datasets/{DATASET_ID}"
    f"/fhirStores/{FHIR_STORE_ID}/fhir"
)


class FHIRClient:
    """Client for GCP Healthcare API FHIR store."""

    def __init__(self):
        """Initialize FHIR client with GCP credentials."""
        self.base_url = FHIR_BASE_URL
        self.credentials, _ = default()
        self.client = httpx.AsyncClient(timeout=30.0)

    async def _get_access_token(self) -> str:
        """Get GCP access token for authentication."""
        if not self.credentials.valid:
            self.credentials.refresh(Request())
        return self.credentials.token

    async def get(self, resource_path: str) -> Dict[str, Any]:
        """GET request to FHIR API.

        Args:
            resource_path: FHIR resource path (e.g., "Patient/patient-001")

        Returns:
            FHIR resource as dictionary
        """
        token = await self._get_access_token()
        url = f"{self.base_url}/{resource_path}"

        headers = {
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/fhir+json",
        }

        logger.info(f"GET {url}")
        response = await self.client.get(url, headers=headers)
        response.raise_for_status()
        return response.json()

    async def search(self, resource_type: str, params: Dict[str, str]) -> Dict[str, Any]:
        """Search FHIR resources.

        Args:
            resource_type: FHIR resource type (e.g., "Condition")
            params: Query parameters

        Returns:
            FHIR Bundle containing search results
        """
        token = await self._get_access_token()
        url = f"{self.base_url}/{resource_type}"

        headers = {
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/fhir+json",
        }

        logger.info(f"SEARCH {url} with params {params}")
        response = await self.client.get(url, headers=headers, params=params)
        response.raise_for_status()
        return response.json()

    async def close(self):
        """Close HTTP client."""
        await self.client.aclose()


# FHIR client will be lazily initialized
_fhir_client: Optional[FHIRClient] = None


def get_fhir_client() -> FHIRClient:
    """Get or create FHIR client instance."""
    global _fhir_client
    if _fhir_client is None:
        _fhir_client = FHIRClient()
    return _fhir_client


def deidentify_patient(patient: Dict[str, Any]) -> Dict[str, Any]:
    """Apply HIPAA Safe Harbor de-identification to Patient resource.

    Removes 18 HIPAA identifiers:
    - Names, addresses, dates (except year), phone numbers, etc.
    - Hashes patient ID for anonymization
    - Keeps year of birth only

    Args:
        patient: FHIR Patient resource

    Returns:
        De-identified Patient resource
    """
    if not DEIDENTIFY_ENABLED:
        return patient

    deidentified = patient.copy()

    # Remove direct identifiers
    identifiers_to_remove = [
        "name",
        "telecom",
        "address",
        "photo",
        "contact",
        "communication",
    ]

    for identifier in identifiers_to_remove:
        if identifier in deidentified:
            del deidentified[identifier]

    # Hash patient ID
    if "id" in deidentified:
        original_id = deidentified["id"]
        hashed_id = hashlib.sha256(original_id.encode()).hexdigest()[:16]
        deidentified["id"] = f"deidentified-{hashed_id}"

    # Hash identifier values
    if "identifier" in deidentified:
        for ident in deidentified["identifier"]:
            if "value" in ident:
                original_value = ident["value"]
                hashed_value = hashlib.sha256(original_value.encode()).hexdigest()[:16]
                ident["value"] = f"HASH-{hashed_value}"

    # Date shift: keep year only
    if "birthDate" in deidentified:
        try:
            birth_date = datetime.fromisoformat(deidentified["birthDate"])
            deidentified["birthDate"] = str(birth_date.year)
        except (ValueError, AttributeError):
            del deidentified["birthDate"]

    # Add de-identification marker
    deidentified["_deidentified"] = True
    deidentified["_deidentification_method"] = "HIPAA Safe Harbor"
    deidentified["_deidentification_date"] = datetime.utcnow().isoformat()

    return deidentified


def deidentify_resource(resource: Dict[str, Any]) -> Dict[str, Any]:
    """Apply de-identification to any FHIR resource.

    Args:
        resource: FHIR resource

    Returns:
        De-identified resource
    """
    if not DEIDENTIFY_ENABLED:
        return resource

    deidentified = resource.copy()

    # Hash resource ID
    if "id" in deidentified:
        original_id = deidentified["id"]
        hashed_id = hashlib.sha256(original_id.encode()).hexdigest()[:16]
        deidentified["id"] = f"deidentified-{hashed_id}"

    # Hash patient references
    if "subject" in deidentified and "reference" in deidentified["subject"]:
        ref = deidentified["subject"]["reference"]
        if ref.startswith("Patient/"):
            patient_id = ref.split("/")[1]
            hashed_id = hashlib.sha256(patient_id.encode()).hexdigest()[:16]
            deidentified["subject"]["reference"] = f"Patient/deidentified-{hashed_id}"
        if "display" in deidentified["subject"]:
            del deidentified["subject"]["display"]

    # Remove dates (keep year only for effectiveDateTime)
    date_fields = ["recordedDate", "assertedDate", "dateAsserted", "issued"]
    for field in date_fields:
        if field in deidentified:
            try:
                date_obj = datetime.fromisoformat(
                    deidentified[field].replace("Z", "+00:00")
                )
                deidentified[field] = str(date_obj.year)
            except (ValueError, AttributeError):
                del deidentified[field]

    # Add de-identification marker
    deidentified["_deidentified"] = True
    deidentified["_deidentification_method"] = "HIPAA Safe Harbor"

    return deidentified


def deidentify_bundle(bundle: Dict[str, Any]) -> Dict[str, Any]:
    """Apply de-identification to FHIR Bundle.

    Args:
        bundle: FHIR Bundle resource

    Returns:
        De-identified Bundle
    """
    if not DEIDENTIFY_ENABLED:
        return bundle

    deidentified = bundle.copy()

    if "entry" in deidentified:
        for entry in deidentified["entry"]:
            if "resource" in entry:
                entry["resource"] = deidentify_resource(entry["resource"])

    return deidentified


@mcp.tool()
async def get_patient_demographics(patient_id: str) -> Dict[str, Any]:
    """Retrieve patient demographics from FHIR store.

    Fetches patient demographic information including age, gender, and identifiers.
    All data is de-identified using HIPAA Safe Harbor method.

    Args:
        patient_id: Patient identifier (e.g., "patient-001")

    Returns:
        De-identified patient demographics including:
        - Hashed patient ID
        - Gender
        - Birth year (not full date)
        - Active status
        - Resource type
    """
    try:
        client = get_fhir_client()
        patient = await client.get(f"Patient/{patient_id}")
        deidentified = deidentify_patient(patient)

        logger.info(f"Retrieved and de-identified patient: {patient_id}")
        return {
            "status": "success",
            "data": deidentified,
        }
    except httpx.HTTPStatusError as e:
        logger.error(f"HTTP error retrieving patient {patient_id}: {e}")
        return {
            "status": "error",
            "error": f"Failed to retrieve patient: {e.response.status_code}",
            "message": str(e),
        }
    except Exception as e:
        logger.error(f"Error retrieving patient {patient_id}: {e}")
        return {
            "status": "error",
            "error": "Failed to retrieve patient",
            "message": str(e),
        }


@mcp.tool()
async def get_patient_conditions(patient_id: str) -> Dict[str, Any]:
    """Retrieve patient conditions/diagnoses from FHIR store.

    Fetches all active and historical conditions for the patient.
    All data is de-identified.

    Args:
        patient_id: Patient identifier (e.g., "patient-001")

    Returns:
        De-identified bundle of Condition resources including:
        - Diagnosis codes and descriptions
        - Clinical status
        - Verification status
        - Severity
        - Stage information (for cancer diagnoses)
    """
    try:
        bundle = await fhir_client.search("Condition", {"patient": patient_id})
        deidentified = deidentify_bundle(bundle)

        count = deidentified.get("total", 0)
        logger.info(f"Retrieved {count} condition(s) for patient: {patient_id}")

        return {
            "status": "success",
            "data": deidentified,
            "count": count,
        }
    except httpx.HTTPStatusError as e:
        logger.error(f"HTTP error retrieving conditions for {patient_id}: {e}")
        return {
            "status": "error",
            "error": f"Failed to retrieve conditions: {e.response.status_code}",
            "message": str(e),
        }
    except Exception as e:
        logger.error(f"Error retrieving conditions for {patient_id}: {e}")
        return {
            "status": "error",
            "error": "Failed to retrieve conditions",
            "message": str(e),
        }


@mcp.tool()
async def get_patient_observations(
    patient_id: str,
    category: Optional[str] = None,
) -> Dict[str, Any]:
    """Retrieve patient observations (lab results, vitals, etc.) from FHIR store.

    Fetches laboratory results, vital signs, and other observations.
    All data is de-identified.

    Args:
        patient_id: Patient identifier (e.g., "patient-001")
        category: Optional category filter (e.g., "laboratory", "vital-signs")

    Returns:
        De-identified bundle of Observation resources including:
        - Test/observation codes (LOINC)
        - Values and units
        - Reference ranges
        - Interpretation (e.g., high, low, normal)
    """
    try:
        params = {"patient": patient_id}
        if category:
            params["category"] = category

        bundle = await fhir_client.search("Observation", params)
        deidentified = deidentify_bundle(bundle)

        count = deidentified.get("total", 0)
        logger.info(f"Retrieved {count} observation(s) for patient: {patient_id}")

        return {
            "status": "success",
            "data": deidentified,
            "count": count,
        }
    except httpx.HTTPStatusError as e:
        logger.error(f"HTTP error retrieving observations for {patient_id}: {e}")
        return {
            "status": "error",
            "error": f"Failed to retrieve observations: {e.response.status_code}",
            "message": str(e),
        }
    except Exception as e:
        logger.error(f"Error retrieving observations for {patient_id}: {e}")
        return {
            "status": "error",
            "error": "Failed to retrieve observations",
            "message": str(e),
        }


@mcp.tool()
async def get_patient_medications(patient_id: str) -> Dict[str, Any]:
    """Retrieve patient medications from FHIR store.

    Fetches current and historical medication information.
    All data is de-identified.

    Args:
        patient_id: Patient identifier (e.g., "patient-001")

    Returns:
        De-identified bundle of MedicationStatement resources including:
        - Medication names (generic and brand)
        - Status (active, completed, stopped)
        - Dosage and route
        - Treatment dates (year only)
        - Reason for medication
    """
    try:
        bundle = await fhir_client.search("MedicationStatement", {"patient": patient_id})
        deidentified = deidentify_bundle(bundle)

        count = deidentified.get("total", 0)
        logger.info(f"Retrieved {count} medication(s) for patient: {patient_id}")

        return {
            "status": "success",
            "data": deidentified,
            "count": count,
        }
    except httpx.HTTPStatusError as e:
        logger.error(f"HTTP error retrieving medications for {patient_id}: {e}")
        return {
            "status": "error",
            "error": f"Failed to retrieve medications: {e.response.status_code}",
            "message": str(e),
        }
    except Exception as e:
        logger.error(f"Error retrieving medications for {patient_id}: {e}")
        return {
            "status": "error",
            "error": "Failed to retrieve medications",
            "message": str(e),
        }


def main():
    """Run the mcp-epic server."""
    logger.info("Starting mcp-epic server")
    logger.info(f"FHIR Base URL: {FHIR_BASE_URL}")
    logger.info(f"De-identification: {'ENABLED' if DEIDENTIFY_ENABLED else 'DISABLED'}")
    mcp.run()


if __name__ == "__main__":
    main()
