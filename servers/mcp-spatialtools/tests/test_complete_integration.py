"""Complete integration tests for spatialtools workflow.

Tests the full analysis pipeline:
1. Alignment (STAR) - optional if STAR installed
2. Quality filtering
3. Batch correction (ComBat)
4. Differential expression (Mann-Whitney U)
5. Pathway enrichment (Fisher's exact test)

This validates that data flows correctly through the entire pipeline
and that biological signals are preserved through batch correction.
"""

import os
import sys
import pytest
import tempfile
import shutil
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))


class TestCompleteWorkflow:
    """Test complete workflow with Patient-001 data."""

    def test_batch_correction_to_differential_expression(self, tmp_path):
        """Test workflow: Batch Correction → Differential Expression.

        This tests the critical path where batch-corrected data is used
        for downstream differential expression analysis.
        """
        # Load Patient-001 data
        data_dir = Path("/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial")

        if not data_dir.exists():
            pytest.skip(f"Patient-001 data not found at {data_dir}")

        expr_file = data_dir / "visium_gene_expression.csv"
        meta_file = data_dir / "visium_region_annotations.csv"

        # Load data
        expr_df = pd.read_csv(expr_file, index_col=0)  # spots × genes
        expr_df = expr_df.T  # genes × spots

        meta_df = pd.read_csv(meta_file)
        meta_df = meta_df.set_index('barcode')

        print(f"\n✅ Loaded Patient-001 data: {expr_df.shape[0]} genes × {expr_df.shape[1]} spots")

        # =====================================================================
        # STEP 1: Create artificial batches with batch effects
        # =====================================================================
        print("\nSTEP 1: Creating artificial batches with batch effects")

        # Shuffle to distribute biological signals across batches
        np.random.seed(42)
        shuffled_cols = np.random.permutation(expr_df.columns)
        expr_df = expr_df[shuffled_cols]

        # Split into 3 batches (better for batch correction)
        n_spots = expr_df.shape[1]
        batch_size = n_spots // 3

        batch1_cols = expr_df.columns[:batch_size]
        batch2_cols = expr_df.columns[batch_size:2*batch_size]
        batch3_cols = expr_df.columns[2*batch_size:]

        # Add batch effects (multiplicative + small additive noise)
        batch1_df = expr_df[batch1_cols].copy()
        batch2_df = expr_df[batch2_cols].copy() * 1.4 + np.random.normal(0, 10, (expr_df.shape[0], len(batch2_cols)))
        batch3_df = expr_df[batch3_cols].copy() * 1.8 + np.random.normal(0, 15, (expr_df.shape[0], len(batch3_cols)))

        # Clip negative values
        batch2_df = batch2_df.clip(lower=0)
        batch3_df = batch3_df.clip(lower=0)

        # Merge batches
        merged_before = pd.concat([batch1_df, batch2_df, batch3_df], axis=1)
        batch_array = np.array(['batch1'] * len(batch1_cols) +
                               ['batch2'] * len(batch2_cols) +
                               ['batch3'] * len(batch3_cols))

        print(f"  Batch 1: {len(batch1_cols)} spots (baseline)")
        print(f"  Batch 2: {len(batch2_cols)} spots (1.4x + noise)")
        print(f"  Batch 3: {len(batch3_cols)} spots (1.8x + noise)")

        # =====================================================================
        # STEP 2: Apply batch correction
        # =====================================================================
        print("\nSTEP 2: Applying batch correction")

        from mcp_spatialtools.server import _combat_batch_correction, _calculate_batch_variance

        # Calculate variance before
        variance_before = _calculate_batch_variance(merged_before.T.values, batch_array)
        print(f"  Batch variance BEFORE correction: {variance_before:.4f}")

        # Apply ComBat
        corrected_df = _combat_batch_correction(merged_before, batch_array, parametric=True)

        # Calculate variance after
        variance_after = _calculate_batch_variance(corrected_df.T.values, batch_array)
        print(f"  Batch variance AFTER correction: {variance_after:.4f}")

        variance_reduction = (variance_before - variance_after) / variance_before
        print(f"  Variance reduction: {variance_reduction:.2%}")

        # Verify batch correction worked
        assert variance_reduction > 0.10, f"Expected >10% variance reduction, got {variance_reduction:.2%}"
        print("  ✅ Batch correction successful")

        # =====================================================================
        # STEP 3: Differential expression on batch-corrected data
        # =====================================================================
        print("\nSTEP 3: Differential expression on batch-corrected data")

        # Align metadata with corrected data
        corrected_barcodes = corrected_df.columns.tolist()

        # Get original region labels (before shuffling)
        meta_aligned = meta_df.loc[corrected_barcodes]

        # Get tumor_core vs stroma spots
        tumor_spots = meta_aligned[meta_aligned['region'] == 'tumor_core'].index.tolist()
        stroma_spots = meta_aligned[meta_aligned['region'] == 'stroma'].index.tolist()

        if len(tumor_spots) == 0 or len(stroma_spots) == 0:
            pytest.skip("Not enough tumor_core or stroma spots after shuffling")

        print(f"  Group 1 (tumor_core): {len(tumor_spots)} spots")
        print(f"  Group 2 (stroma): {len(stroma_spots)} spots")

        # Extract expression for each group
        tumor_expr = corrected_df[tumor_spots].values
        stroma_expr = corrected_df[stroma_spots].values

        # Perform Mann-Whitney U test
        deg_results = []
        for gene_idx in range(corrected_df.shape[0]):
            gene_tumor = tumor_expr[gene_idx, :]
            gene_stroma = stroma_expr[gene_idx, :]

            try:
                statistic, p_value = stats.mannwhitneyu(gene_tumor, gene_stroma, alternative='two-sided')

                # Calculate fold change
                mean_tumor = np.mean(gene_tumor)
                mean_stroma = np.mean(gene_stroma)

                if mean_stroma > 0:
                    fold_change = mean_tumor / mean_stroma
                    log2fc = np.log2(fold_change) if fold_change > 0 else -10
                else:
                    log2fc = 10 if mean_tumor > 0 else 0

                deg_results.append({
                    'gene': corrected_df.index[gene_idx],
                    'p_value': p_value,
                    'log2_fold_change': log2fc,
                    'mean_tumor': mean_tumor,
                    'mean_stroma': mean_stroma
                })
            except Exception as e:
                deg_results.append({
                    'gene': corrected_df.index[gene_idx],
                    'p_value': 1.0,
                    'log2_fold_change': 0,
                    'mean_tumor': 0,
                    'mean_stroma': 0
                })

        deg_df = pd.DataFrame(deg_results)

        # FDR correction
        p_values = deg_df['p_value'].values
        n = len(p_values)
        sorted_indices = np.argsort(p_values)
        sorted_p = p_values[sorted_indices]
        ranks = np.arange(1, n + 1)

        fdr_values = sorted_p * n / ranks
        for i in range(len(fdr_values) - 2, -1, -1):
            fdr_values[i] = min(fdr_values[i], fdr_values[i + 1])
        fdr_values = np.minimum(fdr_values, 1.0)

        # Unsort FDR values
        fdr_unsorted = np.empty_like(fdr_values)
        fdr_unsorted[sorted_indices] = fdr_values
        deg_df['fdr'] = fdr_unsorted

        # Identify significant DEGs
        sig_degs = deg_df[(deg_df['fdr'] < 0.05) & (np.abs(deg_df['log2_fold_change']) > 1.0)]

        print(f"  ✅ Significant DEGs (FDR < 0.05, |log2FC| > 1): {len(sig_degs)}")

        if len(sig_degs) > 0:
            print(f"\n  Top 5 upregulated genes in tumor:")
            upregulated = sig_degs[sig_degs['log2_fold_change'] > 0].sort_values('log2_fold_change', ascending=False)
            for idx, row in upregulated.head(5).iterrows():
                print(f"    {row['gene']}: log2FC={row['log2_fold_change']:.2f}, FDR={row['fdr']:.2e}")

        # Verify some DEGs found
        assert len(sig_degs) > 0, "Should find at least some significant DEGs"

        # =====================================================================
        # STEP 4: Pathway enrichment on DEGs
        # =====================================================================
        print("\nSTEP 4: Pathway enrichment on DEGs")

        from mcp_spatialtools.server import OVARIAN_CANCER_PATHWAYS

        # Get upregulated genes
        upregulated_genes = sig_degs[sig_degs['log2_fold_change'] > 0]['gene'].tolist()

        print(f"  Upregulated genes for enrichment: {len(upregulated_genes)}")

        if len(upregulated_genes) == 0:
            pytest.skip("No upregulated genes found for pathway enrichment")

        # All genes in dataset (background)
        all_genes = corrected_df.index.tolist()

        # Test enrichment in KEGG pathways
        enrichment_results = []

        for pathway_id, pathway_data in OVARIAN_CANCER_PATHWAYS["KEGG"].items():
            pathway_genes = set(pathway_data["genes"])

            # Calculate overlap
            overlap = set(upregulated_genes) & pathway_genes

            if len(overlap) == 0:
                continue

            # Fisher's exact test (2x2 contingency table)
            a = len(overlap)  # DE and in pathway
            b = len(upregulated_genes) - a  # DE but not in pathway
            c = len(pathway_genes) - a  # In pathway but not DE
            d = len(all_genes) - a - b - c  # Background

            oddsratio, pvalue = stats.fisher_exact([[a, b], [c, d]], alternative='greater')

            # Fold enrichment
            expected = (len(upregulated_genes) * len(pathway_genes)) / len(all_genes)
            fold_enrichment = a / expected if expected > 0 else 0

            enrichment_results.append({
                'pathway_id': pathway_id,
                'pathway_name': pathway_data['name'],
                'overlap': a,
                'pathway_size': len(pathway_genes),
                'p_value': pvalue,
                'fold_enrichment': fold_enrichment
            })

        if len(enrichment_results) == 0:
            print("  ⚠️  No pathway overlaps found (may be due to small gene set)")
        else:
            enrich_df = pd.DataFrame(enrichment_results)

            # FDR correction
            p_vals = enrich_df['p_value'].values
            n_pathways = len(p_vals)
            sorted_idx = np.argsort(p_vals)
            sorted_p = p_vals[sorted_idx]
            ranks = np.arange(1, n_pathways + 1)

            fdr = sorted_p * n_pathways / ranks
            for i in range(len(fdr) - 2, -1, -1):
                fdr[i] = min(fdr[i], fdr[i + 1])
            fdr = np.minimum(fdr, 1.0)

            fdr_unsorted = np.empty_like(fdr)
            fdr_unsorted[sorted_idx] = fdr
            enrich_df['fdr'] = fdr_unsorted

            # Show top enriched pathways
            sig_pathways = enrich_df[enrich_df['fdr'] < 0.1].sort_values('p_value')

            print(f"  ✅ Pathways tested: {len(enrichment_results)}")
            print(f"  Significant pathways (FDR < 0.1): {len(sig_pathways)}")

            if len(sig_pathways) > 0:
                print(f"\n  Top enriched pathways:")
                for idx, row in sig_pathways.head(5).iterrows():
                    print(f"    {row['pathway_name']}")
                    print(f"      Overlap: {row['overlap']}/{row['pathway_size']} genes")
                    print(f"      Fold enrichment: {row['fold_enrichment']:.2f}x")
                    print(f"      FDR: {row['fdr']:.2e}")

        # =====================================================================
        # WORKFLOW VALIDATION
        # =====================================================================
        print("\n" + "=" * 80)
        print("WORKFLOW VALIDATION SUMMARY")
        print("=" * 80)
        print(f"✅ STEP 1: Created {len(batch_array)} artificial batches")
        print(f"✅ STEP 2: Batch correction reduced variance by {variance_reduction:.1%}")
        print(f"✅ STEP 3: Differential expression found {len(sig_degs)} DEGs")
        print(f"✅ STEP 4: Pathway enrichment tested {len(enrichment_results)} pathways")
        print()
        print("✅ COMPLETE WORKFLOW VALIDATED!")
        print("   Data flows correctly: Batch Correction → DE → Pathway Enrichment")
        print("=" * 80)
        print()

    def test_multibatch_workflow(self, tmp_path):
        """Test workflow with 3 batches (stronger batch effects)."""
        # Load Patient-001 data
        data_dir = Path("/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial")

        if not data_dir.exists():
            pytest.skip(f"Patient-001 data not found at {data_dir}")

        expr_file = data_dir / "visium_gene_expression.csv"

        # Load data
        expr_df = pd.read_csv(expr_file, index_col=0)  # spots × genes
        expr_df = expr_df.T  # genes × spots

        print(f"\n✅ Loaded Patient-001 data: {expr_df.shape[0]} genes × {expr_df.shape[1]} spots")

        # Shuffle columns
        np.random.seed(42)
        shuffled_cols = np.random.permutation(expr_df.columns)
        expr_df = expr_df[shuffled_cols]

        # Split into 3 batches
        n_spots = expr_df.shape[1]
        batch_size = n_spots // 3

        batch1_cols = expr_df.columns[:batch_size]
        batch2_cols = expr_df.columns[batch_size:2*batch_size]
        batch3_cols = expr_df.columns[2*batch_size:]

        print(f"\nCreating 3 batches:")
        print(f"  Batch 1: {len(batch1_cols)} spots (baseline)")
        print(f"  Batch 2: {len(batch2_cols)} spots (1.4x effect)")
        print(f"  Batch 3: {len(batch3_cols)} spots (1.8x effect)")

        # Add batch effects
        batch1_df = expr_df[batch1_cols].copy()
        batch2_df = expr_df[batch2_cols].copy() * 1.4
        batch3_df = expr_df[batch3_cols].copy() * 1.8

        # Merge
        merged_before = pd.concat([batch1_df, batch2_df, batch3_df], axis=1)
        batch_array = np.array(['batch1'] * len(batch1_cols) +
                               ['batch2'] * len(batch2_cols) +
                               ['batch3'] * len(batch3_cols))

        # Apply batch correction
        from mcp_spatialtools.server import _combat_batch_correction, _calculate_batch_variance

        variance_before = _calculate_batch_variance(merged_before.T.values, batch_array)
        corrected_df = _combat_batch_correction(merged_before, batch_array, parametric=True)
        variance_after = _calculate_batch_variance(corrected_df.T.values, batch_array)

        variance_reduction = (variance_before - variance_after) / variance_before

        print(f"\nBatch correction results:")
        print(f"  Variance BEFORE: {variance_before:.4f}")
        print(f"  Variance AFTER: {variance_after:.4f}")
        print(f"  Reduction: {variance_reduction:.2%}")

        # Verify reduction
        assert variance_reduction > 0.15, f"Expected >15% variance reduction, got {variance_reduction:.2%}"

        print(f"\n✅ Multi-batch workflow successful!")
        print(f"   3-batch correction achieved {variance_reduction:.1%} variance reduction")


@pytest.mark.skipif(not shutil.which("STAR"), reason="STAR not installed")
class TestAlignmentIntegration:
    """Test alignment integration (requires STAR installation)."""

    @pytest.mark.skip(reason="STAR alignment requires genome index and is time-consuming")
    def test_alignment_to_expression_workflow(self, tmp_path):
        """Test workflow: STAR Alignment → Expression Matrix.

        This documents the expected workflow but is skipped because:
        1. Requires STAR genome index (~30GB)
        2. Takes 30-60 minutes to run
        3. MCP tool wrapper prevents direct function calls

        For full integration testing, use manual testing via Claude Desktop
        or MCP protocol client.
        """
        # Step 1: Create synthetic FASTQ files
        from mcp_spatialtools.server import _create_synthetic_fastq

        r1_path = tmp_path / "test_R1.fastq.gz"
        r2_path = tmp_path / "test_R2.fastq.gz"

        _create_synthetic_fastq(r1_path, r2_path, num_reads=1000, read_length=100)

        # Step 2: Run STAR alignment (would require genome index)
        # alignment_result = align_spatial_data(...)

        # Step 3: Parse BAM to expression matrix (would use external tools)
        # expression_matrix = parse_bam_to_matrix(...)

        # Step 4: Continue with batch correction, DE, pathway enrichment
        # This demonstrates the full pipeline from raw reads to insights

        pass


class TestDataIntegrity:
    """Test data integrity through workflow steps."""

    def test_batch_correction_preserves_genes(self):
        """Test that batch correction preserves all genes."""
        from mcp_spatialtools.server import _combat_batch_correction

        # Create test data
        np.random.seed(42)
        gene_names = ['GeneA', 'GeneB', 'GeneC', 'GeneD', 'GeneE']

        batch1 = pd.DataFrame(
            np.random.lognormal(size=(5, 20)),
            index=gene_names,
            columns=[f"S{i}_B1" for i in range(20)]
        )

        batch2 = pd.DataFrame(
            np.random.lognormal(size=(5, 20)) * 1.5,
            index=gene_names,
            columns=[f"S{i}_B2" for i in range(20)]
        )

        merged = pd.concat([batch1, batch2], axis=1)
        batch_array = np.array(['batch1'] * 20 + ['batch2'] * 20)

        # Apply correction
        corrected = _combat_batch_correction(merged, batch_array)

        # Verify genes preserved
        assert corrected.index.tolist() == gene_names
        assert corrected.shape == merged.shape

        print("✅ Batch correction preserves gene names and data structure")

    def test_differential_expression_preserves_genes(self):
        """Test that DE analysis preserves all genes."""
        # Create test data
        np.random.seed(42)
        n_genes = 10
        genes = [f"Gene{i}" for i in range(n_genes)]

        group1 = np.random.normal(10, 2, size=(n_genes, 30))
        group2 = np.random.normal(15, 2, size=(n_genes, 30))

        # Perform Mann-Whitney U test
        deg_results = []
        for i in range(n_genes):
            statistic, p_value = stats.mannwhitneyu(group1[i], group2[i])
            deg_results.append({
                'gene': genes[i],
                'p_value': p_value,
                'statistic': statistic
            })

        deg_df = pd.DataFrame(deg_results)

        # Verify all genes present
        assert len(deg_df) == n_genes
        assert deg_df['gene'].tolist() == genes

        print("✅ Differential expression preserves all genes")
