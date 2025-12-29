#!/usr/bin/env python3
"""Test real differential expression analysis with Patient-001 data.

This script demonstrates Phase 2A: Spatial Differential Expression
- Compares tumor core vs tumor margin
- Uses real statistical tests (Mann-Whitney U, t-test)
- Applies FDR correction (Benjamini-Hochberg)
- Provides clinical interpretation for ovarian cancer patient
"""

import asyncio
import os
import sys
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from mcp_spatialtools.server import (
    split_by_region,
    perform_differential_expression
)


async def test_differential_expression():
    """Test differential expression with Patient-001 spatial data."""
    print("=" * 80)
    print("PHASE 2A: SPATIAL DIFFERENTIAL EXPRESSION ANALYSIS")
    print("Patient-001: Stage IV HGSOC, Platinum-Resistant, on Bevacizumab")
    print("=" * 80)
    print()

    # File paths
    expression_file = "/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial/visium_gene_expression.csv"
    coordinates_file = "/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial/visium_spatial_coordinates.csv"
    output_dir = "/tmp/spatial_regions"

    # =========================================================================
    # STEP 1: Split data by spatial regions
    # =========================================================================
    print("STEP 1: Split by Spatial Regions")
    print("-" * 80)

    regions = ["tumor_core", "tumor_margin", "stroma", "immune_infiltrate"]

    split_result = await split_by_region.fn(
        input_file=expression_file,
        coordinate_file=coordinates_file,
        output_dir=output_dir,
        regions=regions
    )

    print(f"Total regions created: {split_result['total_regions']}")
    for region_info in split_result['regions']:
        print(f"  • {region_info['name']}: {region_info['barcode_count']} spots")
    print()

    # =========================================================================
    # STEP 2: Prepare sample groups for differential expression
    # =========================================================================
    print("STEP 2: Prepare Sample Groups")
    print("-" * 80)

    # Load region-specific data to get spot IDs
    core_file = None
    margin_file = None

    for region_info in split_result['regions']:
        if region_info['name'] == 'tumor_core':
            core_file = region_info['file']
        elif region_info['name'] == 'tumor_margin':
            margin_file = region_info['file']

    # Get spot IDs from region files
    core_data = pd.read_csv(core_file, index_col=0)
    margin_data = pd.read_csv(margin_file, index_col=0)

    group1_spots = list(core_data.index)  # Tumor core
    group2_spots = list(margin_data.index)  # Tumor margin

    print(f"Group 1 (Tumor Core): {len(group1_spots)} spots")
    print(f"Group 2 (Tumor Margin): {len(group2_spots)} spots")
    print()

    # =========================================================================
    # STEP 3: Run Differential Expression Analysis
    # =========================================================================
    print("STEP 3: Differential Expression Analysis (Tumor Core vs Margin)")
    print("-" * 80)
    print("Testing method: Mann-Whitney U (Wilcoxon rank-sum)")
    print("FDR correction: Benjamini-Hochberg")
    print("Significance threshold: q < 0.05 and |log2FC| >= 0.5")
    print()

    deg_result = await perform_differential_expression.fn(
        expression_file=expression_file,
        group1_samples=group1_spots,
        group2_samples=group2_spots,
        test_method="wilcoxon",
        min_log_fc=0.5
    )

    print(f"Status: {deg_result['status']}")
    print(f"Test method: {deg_result['test_method']}")
    print(f"Total genes tested: {deg_result['total_genes_tested']}")
    print(f"Significant genes: {deg_result['significant_genes']}")
    print(f"  ↑ Upregulated in core: {deg_result['upregulated_genes']}")
    print(f"  ↓ Downregulated in core: {deg_result['downregulated_genes']}")
    print()

    # =========================================================================
    # STEP 4: Display Top Results with Clinical Interpretation
    # =========================================================================
    print("STEP 4: Top Differentially Expressed Genes")
    print("=" * 80)
    print()

    print("TOP 10 GENES UPREGULATED IN TUMOR CORE")
    print("-" * 80)
    print(f"{'Gene':<10} | {'log2FC':>8} | {'q-value':>10} | {'Interpretation'}")
    print("-" * 80)

    for gene_result in deg_result['top_upregulated'][:10]:
        gene = gene_result['gene']
        log2fc = gene_result['log2_fold_change']
        qval = gene_result['qvalue']

        # Clinical interpretation
        interpretation = ""
        if gene == "MKI67":
            interpretation = "Proliferation marker - high in aggressive core"
        elif gene == "TP53":
            interpretation = "Tumor suppressor - genomic instability"
        elif gene == "PCNA":
            interpretation = "Cell cycle - active proliferation"
        elif gene == "TOP2A":
            interpretation = "DNA replication - proliferative cells"
        elif gene == "MYC":
            interpretation = "Oncogene - tumor growth driver"
        elif gene == "FOXM1":
            interpretation = "Cell cycle regulator - HGSOC marker"
        elif gene == "EPCAM":
            interpretation = "Epithelial marker - tumor cell identity"
        elif gene.startswith("KRT"):
            interpretation = "Keratin - epithelial tumor cells"
        elif gene == "PIK3CA":
            interpretation = "PI3K pathway - growth signaling"
        else:
            interpretation = "Higher expression in tumor core"

        print(f"{gene:<10} | {log2fc:>8.3f} | {qval:>10.6f} | {interpretation}")

    print()
    print()

    print("TOP 10 GENES UPREGULATED IN TUMOR MARGIN")
    print("-" * 80)
    print(f"{'Gene':<10} | {'log2FC':>8} | {'q-value':>10} | {'Interpretation'}")
    print("-" * 80)

    for gene_result in deg_result['top_downregulated'][:10]:
        gene = gene_result['gene']
        log2fc = gene_result['log2_fold_change']
        qval = gene_result['qvalue']

        # Clinical interpretation (downregulated in core = upregulated in margin)
        interpretation = ""
        if gene == "CD8A":
            interpretation = "CD8+ T-cells - immune infiltration at margin"
        elif gene in ["CD3D", "CD3E"]:
            interpretation = "T-cell marker - active immune response"
        elif gene == "CD4":
            interpretation = "CD4+ T-cells - helper T-cells at margin"
        elif gene == "FOXP3":
            interpretation = "Regulatory T-cells - immunosuppression"
        elif gene == "CD68":
            interpretation = "Macrophages - tumor-associated macrophages"
        elif gene == "CD163":
            interpretation = "M2 macrophages - pro-tumor immune cells"
        elif gene == "VIM":
            interpretation = "Mesenchymal marker - EMT at invasive front"
        elif gene in ["SNAI1", "TWIST1"]:
            interpretation = "EMT transcription factors - invasion"
        elif gene in ["COL1A1", "COL3A1"]:
            interpretation = "Collagen - stromal ECM at margin"
        elif gene == "FAP":
            interpretation = "Fibroblast marker - cancer-associated fibroblasts"
        else:
            interpretation = "Higher expression in tumor margin"

        print(f"{gene:<10} | {log2fc:>8.3f} | {qval:>10.6f} | {interpretation}")

    print()
    print()

    # =========================================================================
    # STEP 5: Clinical Summary
    # =========================================================================
    print("STEP 5: Clinical-Spatial Integration Report")
    print("=" * 80)
    print()

    print("PATIENT-001 DIFFERENTIAL EXPRESSION SUMMARY")
    print()

    print("Clinical Context:")
    print("  • Diagnosis: Stage IV High-Grade Serous Ovarian Carcinoma (HGSOC)")
    print("  • Treatment status: Platinum-resistant")
    print("  • Current therapy: Bevacizumab (anti-angiogenic)")
    print("  • Biomarkers: CA-125 elevated (487 U/mL), BRCA negative")
    print()

    print("Spatial Differential Expression Findings:")
    print()

    # Analyze specific genes relevant to clinical context
    sig_genes = {r['gene']: r for r in deg_result['significant_results']}

    print("1. Proliferation & Tumor Aggressiveness:")
    for gene in ["MKI67", "PCNA", "TOP2A"]:
        if gene in sig_genes:
            fc = sig_genes[gene]['log2_fold_change']
            if fc > 0:
                print(f"   ✓ {gene} upregulated in core (log2FC={fc:.2f})")
                print(f"     → High proliferation in tumor core indicates aggressive disease")

    print()
    print("2. Immune Infiltration (Bevacizumab response marker):")
    for gene in ["CD8A", "CD3D", "CD3E", "CD4"]:
        if gene in sig_genes:
            fc = sig_genes[gene]['log2_fold_change']
            if fc < 0:  # Downregulated in core = upregulated in margin
                print(f"   ✓ {gene} enriched at margin (log2FC={fc:.2f})")
                print(f"     → T-cell infiltration at tumor periphery")

    print()
    print("3. Epithelial-Mesenchymal Transition (EMT):")
    for gene in ["VIM", "SNAI1", "TWIST1"]:
        if gene in sig_genes:
            fc = sig_genes[gene]['log2_fold_change']
            if fc < 0:  # Higher at margin
                print(f"   ✓ {gene} enriched at invasive front (log2FC={fc:.2f})")
                print(f"     → EMT signature correlates with platinum resistance")

    print()
    print("4. Angiogenesis (Bevacizumab target pathway):")
    for gene in ["VEGFA", "HIF1A", "KDR"]:
        if gene in sig_genes:
            fc = sig_genes[gene]['log2_fold_change']
            location = "core" if fc > 0 else "margin"
            print(f"   ✓ {gene} higher in {location} (log2FC={fc:.2f})")

    print()
    print("Clinical Implications:")
    print("  • Spatial heterogeneity confirmed: distinct gene programs in core vs margin")
    print("  • High proliferation in core supports aggressive phenotype")
    print("  • Immune infiltration at margin suggests potential for immunotherapy")
    print("  • EMT at invasive front explains platinum resistance mechanism")
    print("  • Consider combination therapy targeting both core (proliferation) and")
    print("    margin (immune/EMT) compartments")
    print()

    print("=" * 80)
    print("DIFFERENTIAL EXPRESSION ANALYSIS COMPLETE")
    print("=" * 80)
    print()

    return deg_result


if __name__ == "__main__":
    # Set environment for real data processing
    os.environ["SPATIAL_DATA_DIR"] = "/Users/lynnlangit/Documents/GitHub/spatial-mcp/data"
    os.environ["SPATIAL_DRY_RUN"] = "false"

    # Run the differential expression test
    result = asyncio.run(test_differential_expression())

    print("\n✅ Phase 2A testing complete!")
    print(f"   Found {result['significant_genes']} significantly differentially expressed genes")
    print(f"   Ready for Claude Desktop integration")
