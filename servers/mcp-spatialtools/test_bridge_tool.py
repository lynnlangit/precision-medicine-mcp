#!/usr/bin/env python3
"""Test the clinical-spatial bridge tool."""

import asyncio
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from mcp_spatialtools.server import get_spatial_data_for_patient


async def test_bridge_tool():
    """Test the clinical-spatial bridge with Patient 001."""
    print("=" * 80)
    print("Testing Clinical-Spatial Bridge Tool")
    print("=" * 80)
    print()

    # Test 1: Basic lookup without clinical context
    print("Test 1: Basic patient lookup")
    result = await get_spatial_data_for_patient.fn(
        patient_id="patient-001",
        include_clinical_context=False
    )
    print(f"Status: {result.get('status')}")
    print(f"Spatial dataset: {result.get('spatial_dataset')}")
    print(f"Data directory: {result.get('data_directory')}")
    print(f"Available files: {result.get('available_files')}")
    print(f"Data ready: {result.get('data_ready')}")
    print()

    # Test 2: With clinical context (ovarian cancer patient)
    print("Test 2: With clinical context - Ovarian Cancer Patient")
    result = await get_spatial_data_for_patient.fn(
        patient_id="patient-001",
        conditions=["Stage IV High-Grade Serous Ovarian Carcinoma", "platinum-resistant"],
        medications=["Bevacizumab", "Carboplatin", "Paclitaxel"],
        biomarkers={"CA-125": 487, "BRCA": "negative"},
        tissue_type="tumor",
        include_clinical_context=True
    )

    print(f"Status: {result.get('status')}")
    print()

    print("Clinical Summary:")
    if "clinical_summary" in result:
        cs = result["clinical_summary"]
        print(f"  Patient ID: {cs['patient_id']}")
        print(f"  Conditions: {cs['conditions']}")
        print(f"  Medications: {cs['medications']}")
        print(f"  Biomarkers: {cs['biomarkers']}")
        print(f"  Tissue type: {cs['tissue_type']}")
    print()

    print(f"Genes of Interest ({result.get('num_genes_of_interest')}):")
    for gene in result.get("genes_of_interest", [])[:15]:  # Show first 15
        print(f"  - {gene}")
    print()

    print("Suggested Analyses:")
    for analysis in result.get("suggested_analyses", []):
        print(f"  • {analysis}")
    print()

    print("Data Files:")
    for file_type, file_path in result.get("files", {}).items():
        status = "✅ EXISTS" if file_path and os.path.exists(file_path) else "❌ MISSING"
        print(f"  {file_type}: {status}")
        if file_path:
            print(f"    Path: {file_path}")
    print()

    # Test 3: Non-existent patient
    print("Test 3: Non-existent patient (error handling)")
    result = await get_spatial_data_for_patient.fn(
        patient_id="patient-999"
    )
    print(f"Status: {result.get('status')}")
    if result.get('status') == 'error':
        print(f"Error: {result.get('error')}")
        print(f"Searched path: {result.get('searched_path')}")
    print()

    print("=" * 80)
    print("Bridge Tool Testing Complete!")
    print("=" * 80)


if __name__ == "__main__":
    # Set environment
    os.environ["SPATIAL_DATA_DIR"] = "/Users/lynnlangit/Documents/GitHub/spatial-mcp/data"
    os.environ["SPATIAL_DRY_RUN"] = "false"

    asyncio.run(test_bridge_tool())
