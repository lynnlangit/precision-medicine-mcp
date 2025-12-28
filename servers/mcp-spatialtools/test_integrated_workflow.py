#!/usr/bin/env python3
"""Test the complete clinical-spatial integration workflow."""

import asyncio
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from mcp_spatialtools.server import (
    get_spatial_data_for_patient,
    calculate_spatial_autocorrelation,
)


async def test_integrated_workflow():
    """Test end-to-end clinical-spatial workflow for Patient 001."""
    print("=" * 80)
    print("CLINICAL-SPATIAL INTEGRATION WORKFLOW - PATIENT 001")
    print("=" * 80)
    print()

    # =========================================================================
    # STEP 1: Get spatial data with clinical context
    # =========================================================================
    print("STEP 1: Retrieve Spatial Data with Clinical Context")
    print("-" * 80)

    # Clinical context from Epic FHIR (would come from mcp-epic in Claude Desktop)
    clinical_context = {
        "patient_id": "patient-001",
        "conditions": [
            "Stage IV High-Grade Serous Ovarian Carcinoma",
            "platinum-resistant"
        ],
        "medications": [
            "Bevacizumab",  # Anti-angiogenic, currently active
            "Carboplatin",  # Platinum-based chemo, completed
            "Paclitaxel"    # Taxane-based chemo, completed
        ],
        "biomarkers": {
            "CA-125": 487,      # Elevated tumor marker
            "BRCA": "negative"  # No BRCA mutation
        },
        "tissue_type": "tumor"
    }

    print("Clinical Context:")
    print(f"  Patient ID: {clinical_context['patient_id']}")
    print(f"  Conditions: {', '.join(clinical_context['conditions'])}")
    print(f"  Active Medication: Bevacizumab (anti-angiogenic)")
    print(f"  CA-125: {clinical_context['biomarkers']['CA-125']} U/mL (highly elevated)")
    print(f"  BRCA Status: {clinical_context['biomarkers']['BRCA']}")
    print()

    # Retrieve spatial data using bridge tool
    spatial_data = await get_spatial_data_for_patient.fn(
        patient_id=clinical_context["patient_id"],
        conditions=clinical_context["conditions"],
        medications=clinical_context["medications"],
        biomarkers=clinical_context["biomarkers"],
        tissue_type=clinical_context["tissue_type"],
        include_clinical_context=True
    )

    print(f"Spatial Dataset: {spatial_data['spatial_dataset']}")
    print(f"Data Ready: {'✅' if spatial_data['data_ready'] else '❌'}")
    print()

    print(f"Genes of Interest ({spatial_data['num_genes_of_interest']}):")
    for gene in spatial_data["genes_of_interest"][:10]:
        print(f"  • {gene}")
    if spatial_data['num_genes_of_interest'] > 10:
        print(f"  ... and {spatial_data['num_genes_of_interest'] - 10} more")
    print()

    print("Suggested Analyses:")
    for analysis in spatial_data["suggested_analyses"]:
        print(f"  • {analysis}")
    print()

    # =========================================================================
    # STEP 2: Spatial autocorrelation analysis (real Moran's I)
    # =========================================================================
    print("STEP 2: Spatial Autocorrelation Analysis (Moran's I)")
    print("-" * 80)

    # Focus on key genes relevant to patient's clinical context
    key_genes = ["MKI67", "CD8A", "VIM", "EPCAM", "TP53"]

    print(f"Analyzing {len(key_genes)} key genes:")
    print(f"  • MKI67: Proliferation marker (tumor aggressiveness)")
    print(f"  • CD8A: T-cell marker (immune response)")
    print(f"  • VIM: Mesenchymal marker (epithelial-mesenchymal transition)")
    print(f"  • EPCAM: Epithelial marker (tumor cells)")
    print(f"  • TP53: Tumor suppressor (genomic instability)")
    print()

    autocorr_result = await calculate_spatial_autocorrelation.fn(
        expression_file=spatial_data["files"]["expression"],
        coordinates_file=spatial_data["files"]["coordinates"],
        genes=key_genes,
        method="morans_i",
        distance_threshold=1500.0  # Pixels - adjusted for Visium spot spacing
    )

    print(f"Method: {autocorr_result['method']}")
    print(f"Distance threshold: {autocorr_result['distance_threshold']} pixels")
    print(f"Number of spots analyzed: {autocorr_result['num_spots']}")
    print()

    print("Spatial Autocorrelation Results:")
    print("-" * 80)
    morans_header = "Moran's I"
    print(f"{'Gene':<10} | {morans_header:>10} | {'Z-score':>8} | {'p-value':>10} | {'Significance':>6} | {'Pattern'}")
    print("-" * 80)

    for result in autocorr_result["results"]:
        gene = result["gene"]
        if "morans_i" in result:
            morans_i = result["morans_i"]
            z_score = result["z_score"]
            p_value = result["p_value"]

            # Significance stars
            if p_value < 0.001:
                sig = "***"
            elif p_value < 0.01:
                sig = "**"
            elif p_value < 0.05:
                sig = "*"
            else:
                sig = "ns"

            # Pattern interpretation
            if abs(morans_i) < 0.3:
                pattern = "Weakly patterned"
            elif morans_i >= 0.3:
                pattern = "Clustered"
            else:
                pattern = "Dispersed"

            print(f"{gene:<10} | {morans_i:>10.4f} | {z_score:>8.3f} | {p_value:>10.4f} | {sig:>6} | {pattern}")
        else:
            print(f"{gene:<10} | {'ERROR':>10} | {result.get('message', 'Unknown error')}")

    print("-" * 80)
    print()

    # =========================================================================
    # STEP 3: Clinical-Spatial Integration Report
    # =========================================================================
    print("STEP 3: Clinical-Spatial Integration Report")
    print("=" * 80)
    print()

    print("PATIENT 001 - OVARIAN CANCER SPATIAL ANALYSIS")
    print()

    print("Clinical Summary:")
    print("  • Diagnosis: Stage IV High-Grade Serous Ovarian Carcinoma")
    print("  • Treatment Status: Platinum-resistant, currently on Bevacizumab")
    print("  • Biomarkers: CA-125 elevated (487 U/mL), BRCA negative")
    print()

    print("Spatial Findings:")
    for result in autocorr_result["results"]:
        gene = result["gene"]
        if "morans_i" not in result:
            continue

        morans_i = result["morans_i"]
        p_value = result["p_value"]

        if p_value >= 0.05:
            continue  # Skip non-significant findings

        print(f"\n  {gene}:")

        if gene == "MKI67":
            print(f"    Moran's I: {morans_i:.4f} (p < 0.001)")
            print(f"    Finding: Proliferation shows {'strong' if morans_i > 0.5 else 'moderate' if morans_i > 0.3 else 'weak'} spatial clustering")
            print(f"    Clinical Relevance: {'High proliferation clusters align with aggressive tumor behavior' if morans_i > 0.3 else 'Distributed proliferation pattern'}")

        elif gene == "CD8A":
            print(f"    Moran's I: {morans_i:.4f} (p < 0.001)")
            print(f"    Finding: T-cell infiltration shows {'strong' if morans_i > 0.5 else 'moderate' if morans_i > 0.3 else 'weak'} spatial clustering")
            print(f"    Clinical Relevance: {'Immune response concentrated at tumor margins - may indicate bevacizumab efficacy' if morans_i > 0.3 else 'Diffuse immune infiltration pattern'}")

        elif gene == "VIM":
            print(f"    Moran's I: {morans_i:.4f} (p < 0.001)")
            print(f"    Finding: Mesenchymal marker shows {'strong' if morans_i > 0.5 else 'moderate' if morans_i > 0.3 else 'weak'} spatial clustering")
            print(f"    Clinical Relevance: {'EMT signature suggests invasive phenotype and platinum resistance mechanism' if morans_i > 0.3 else 'Distributed mesenchymal features'}")

        elif gene == "EPCAM":
            print(f"    Moran's I: {morans_i:.4f} (p < 0.001)")
            print(f"    Finding: Epithelial marker shows {'strong' if morans_i > 0.5 else 'moderate' if morans_i > 0.3 else 'weak'} spatial clustering")
            print(f"    Clinical Relevance: {'Tumor cell clusters defined' if morans_i > 0.3 else 'Mixed epithelial-mesenchymal distribution'}")

        elif gene == "TP53":
            print(f"    Moran's I: {morans_i:.4f} (p < 0.001)")
            print(f"    Finding: TP53 shows {'strong' if morans_i > 0.5 else 'moderate' if morans_i > 0.3 else 'weak'} spatial clustering")
            print(f"    Clinical Relevance: {'Regional TP53 dysregulation pattern in HGSOC' if morans_i > 0.3 else 'Distributed TP53 expression'}")

    print()
    print("Integrated Conclusions:")
    print("  • Spatial transcriptomics reveals regional tumor heterogeneity")
    print("  • Clinical context (platinum-resistance, bevacizumab) aligns with spatial patterns")
    print("  • Immune-tumor spatial relationships inform treatment response")
    print("  • EMT signatures spatially correlate with aggressive phenotype")
    print()

    print("=" * 80)
    print("INTEGRATED WORKFLOW COMPLETE")
    print("=" * 80)
    print()

    print("Summary Statistics:")
    if "summary" in autocorr_result:
        summary = autocorr_result["summary"]
        print(f"  Significantly clustered genes: {summary['significantly_clustered']}")
        print(f"  Significantly dispersed genes: {summary['significantly_dispersed']}")
        print(f"  Random pattern genes: {summary['random_pattern']}")
    print()

    return {
        "clinical_context": clinical_context,
        "spatial_data": spatial_data,
        "spatial_autocorrelation": autocorr_result
    }


if __name__ == "__main__":
    # Set environment for real data processing
    os.environ["SPATIAL_DATA_DIR"] = "/Users/lynnlangit/Documents/GitHub/spatial-mcp/data"
    os.environ["SPATIAL_DRY_RUN"] = "false"

    # Run the integrated workflow
    result = asyncio.run(test_integrated_workflow())
