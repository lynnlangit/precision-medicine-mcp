#!/usr/bin/env python3
"""Test real spatial autocorrelation implementation."""

import asyncio
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from mcp_spatialtools.server import calculate_spatial_autocorrelation


async def test_spatial_autocorrelation():
    """Test Moran's I with real Patient 001 data."""
    print("=" * 80)
    print("Testing Real Spatial Autocorrelation (Moran's I)")
    print("=" * 80)
    print()

    # Paths to Patient 001 spatial data
    expression_file = "/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial/visium_gene_expression.csv"
    coordinates_file = "/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial/visium_spatial_coordinates.csv"

    # Genes of interest for ovarian cancer patient (use actual gene names from data)
    genes_to_test = ["MKI67", "CD8A", "VIM", "EPCAM", "TP53"]

    print(f"Expression file: {expression_file}")
    print(f"Coordinates file: {coordinates_file}")
    print(f"Genes to analyze: {genes_to_test}")
    print()

    # Run spatial autocorrelation analysis
    result = await calculate_spatial_autocorrelation.fn(
        expression_file=expression_file,
        coordinates_file=coordinates_file,
        genes=genes_to_test,
        method="morans_i",
        distance_threshold=1500.0  # Pixels - spots are ~1000 pixels apart
    )

    print(f"Status: {result.get('status')}")
    print(f"Method: {result.get('method')}")
    print(f"Genes analyzed: {result.get('genes_analyzed')}")
    print(f"Number of spots: {result.get('num_spots')}")
    print(f"Distance threshold: {result.get('distance_threshold')}")
    print()

    print("Results by Gene:")
    print("-" * 80)
    for gene_result in result.get("results", []):
        gene = gene_result["gene"]
        if "morans_i" in gene_result:
            morans_i = gene_result["morans_i"]
            z_score = gene_result["z_score"]
            p_value = gene_result["p_value"]
            interpretation = gene_result["interpretation"]

            sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"

            print(f"{gene:12} | Moran's I: {morans_i:7.4f} | Z: {z_score:6.3f} | p: {p_value:.4f} {sig}")
            print(f"{'':12} | {interpretation}")
            print()
        else:
            print(f"{gene:12} | {gene_result.get('message', 'Error')}")
            print()

    print("-" * 80)
    print()

    if "summary" in result:
        summary = result["summary"]
        print("Summary:")
        print(f"  Significantly clustered: {summary['significantly_clustered']}")
        print(f"  Significantly dispersed: {summary['significantly_dispersed']}")
        print(f"  Random pattern: {summary['random_pattern']}")
    print()

    print("=" * 80)
    print("Spatial Autocorrelation Testing Complete!")
    print("=" * 80)


if __name__ == "__main__":
    # Set environment
    os.environ["SPATIAL_DATA_DIR"] = "/Users/lynnlangit/Documents/GitHub/spatial-mcp/data"
    os.environ["SPATIAL_DRY_RUN"] = "false"

    asyncio.run(test_spatial_autocorrelation())
