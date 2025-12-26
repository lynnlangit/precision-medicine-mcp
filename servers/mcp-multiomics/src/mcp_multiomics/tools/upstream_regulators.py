"""Upstream regulator analysis for multi-omics data.

Identifies potential upstream regulators (kinases, transcription factors, drugs)
that may explain observed expression changes. Based on enrichment analysis of
known regulator-target relationships.

Key insights per bioinformatician feedback:
- Kinases: Identify activated/inhibited protein kinases
- Transcription Factors: Find TFs regulating gene expression changes
- Drug Responses: Predict compounds that could modulate observed changes

This provides similar insights to IPA's Upstream Regulator Analysis but using
open-source databases and enrichment methods.
"""

import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from ..config import config

logger = logging.getLogger(__name__)


# Mock regulator databases (in production, these would be loaded from files or APIs)
# These are simplified examples - real implementation would use comprehensive databases

KINASE_TARGETS = {
    # Kinase -> list of known phosphorylation targets
    "AKT1": ["GSK3B", "FOXO1", "FOXO3", "TSC2", "MDM2", "BAD", "CDKN1B"],
    "MAPK1": ["ELK1", "MYC", "JUN", "FOS", "TP53", "CREB1", "ATF2"],
    "MTOR": ["RPS6KB1", "EIF4EBP1", "ULK1", "TFEB", "SREBF1"],
    "GSK3B": ["CTNNB1", "JUN", "MYC", "NFATC1", "TP53", "SNAI1"],
    "CDK1": ["RB1", "BRCA1", "TP53", "CDKN1A", "CDKN1B"],
    "PRKCA": ["MARCKS", "RAF1", "MAPK3", "NFKB1", "TP53"],
    "SRC": ["STAT3", "PTK2", "EGFR", "CTNNB1", "YES1"],
    "JAK2": ["STAT3", "STAT5A", "STAT5B", "STAT1", "SHP2"],
    "PI3K": ["AKT1", "PDK1", "PTEN", "FOXO1", "FOXO3"],
    "ERK1": ["ELK1", "RSK1", "MYC", "CREB1", "FOS"],
}

TF_TARGETS = {
    # Transcription Factor -> list of known target genes
    "TP53": ["CDKN1A", "BAX", "BBC3", "MDM2", "GADD45A", "FAS", "PUMA", "NOXA"],
    "MYC": ["CDK4", "CCND1", "TERT", "CAD", "ODC1", "LDHA", "PKM"],
    "NFKB1": ["IL6", "IL8", "TNF", "CXCL10", "CCL2", "ICAM1", "VCAM1"],
    "STAT3": ["BCL2", "MYC", "CCND1", "VEGF", "HIF1A", "SOCS3", "MCL1"],
    "HIF1A": ["VEGFA", "LDHA", "PDK1", "SLC2A1", "PGK1", "ENO1", "BNIP3"],
    "FOXO1": ["SOD2", "GADD45A", "G6PC", "PCK1", "CDKN1B", "BCL6"],
    "CREB1": ["FOS", "JUN", "BDNF", "NR4A1", "EGR1", "BCL2"],
    "JUN": ["MMP1", "MMP3", "IL8", "CCND1", "VEGFA", "FOS"],
    "FOS": ["MMP1", "IL8", "CCND1", "JUNB", "EGR1"],
    "E2F1": ["CCNE1", "CDC6", "MCM2", "PCNA", "RB1", "TP53"],
}

DRUG_TARGETS = {
    # Drug -> list of genes/proteins it affects
    "Alpelisib": ["PIK3CA", "AKT1", "MTOR", "FOXO1", "FOXO3"],  # PI3K inhibitor
    "Trametinib": ["MAP2K1", "MAPK1", "MAPK3", "ELK1", "MYC"],  # MEK inhibitor
    "Everolimus": ["MTOR", "RPS6KB1", "EIF4EBP1", "HIF1A"],  # mTOR inhibitor
    "Trastuzumab": ["ERBB2", "PI3K", "AKT1", "MAPK1"],  # HER2 inhibitor
    "Pembrolizumab": ["PDCD1", "CD274", "CTLA4"],  # PD-1 inhibitor
    "Olaparib": ["PARP1", "BRCA1", "BRCA2", "ATM", "RAD51"],  # PARP inhibitor
    "Venetoclax": ["BCL2", "BAX", "BAK1", "BID"],  # BCL2 inhibitor
    "Ibrutinib": ["BTK", "NFKB1", "BCL2", "MYC"],  # BTK inhibitor
}


def predict_upstream_regulators_impl(
    differential_genes: Dict[str, Dict[str, float]],
    regulator_types: Optional[List[str]] = None,
    fdr_threshold: float = 0.05,
    activation_zscore_threshold: float = 2.0,
) -> Dict[str, Any]:
    """Predict upstream regulators from differential expression data.

    Uses enrichment analysis to identify kinases, transcription factors, and drugs
    that may regulate the observed expression changes. Similar to IPA's Upstream
    Regulator Analysis.

    Method:
    1. For each regulator, test if its known targets are enriched in DEGs
    2. Calculate activation Z-score based on target expression direction
    3. Apply FDR correction to enrichment p-values
    4. Return ranked list of predicted regulators

    Args:
        differential_genes: Dict of gene -> {"log2fc": float, "p_value": float}
                           Genes should be differentially expressed (significant)
        regulator_types: List of regulator types to analyze:
                        ["kinase", "transcription_factor", "drug"]
                        (default: all types)
        fdr_threshold: FDR threshold for significant enrichment (default: 0.05)
        activation_zscore_threshold: |Z-score| threshold for activation/inhibition
                                     (default: 2.0, roughly p < 0.05)

    Returns:
        Dictionary with:
        - kinases: List of predicted kinases with activation state
        - transcription_factors: List of predicted TFs
        - drugs: List of predicted drug responses
        - statistics: Summary statistics
        - method: Enrichment method details
    """
    logger.info("Starting upstream regulator analysis")
    logger.info(f"Input: {len(differential_genes)} differential genes")

    # Return mock data in DRY_RUN mode
    if config.dry_run:
        logger.info("DRY_RUN mode detected - returning mock upstream regulator results")
        return {
            "kinases": [
                {
                    "name": "AKT1",
                    "activation_state": "Activated",
                    "z_score": 3.2,
                    "p_value": 0.0001,
                    "q_value": 0.001,
                    "targets_in_dataset": 5,
                    "targets_consistent": 4,
                    "predicted_effect": "Promotes cell survival and growth",
                },
                {
                    "name": "MTOR",
                    "activation_state": "Activated",
                    "z_score": 2.8,
                    "p_value": 0.0005,
                    "q_value": 0.003,
                    "targets_in_dataset": 4,
                    "targets_consistent": 3,
                    "predicted_effect": "Enhances protein synthesis",
                },
                {
                    "name": "GSK3B",
                    "activation_state": "Inhibited",
                    "z_score": -2.5,
                    "p_value": 0.001,
                    "q_value": 0.005,
                    "targets_in_dataset": 4,
                    "targets_consistent": 3,
                    "predicted_effect": "Reduced tumor suppression",
                },
            ],
            "transcription_factors": [
                {
                    "name": "TP53",
                    "activation_state": "Inhibited",
                    "z_score": -3.5,
                    "p_value": 0.00005,
                    "q_value": 0.0005,
                    "targets_in_dataset": 6,
                    "targets_consistent": 5,
                    "predicted_effect": "Loss of apoptosis and cell cycle control",
                },
                {
                    "name": "MYC",
                    "activation_state": "Activated",
                    "z_score": 3.1,
                    "p_value": 0.0002,
                    "q_value": 0.002,
                    "targets_in_dataset": 5,
                    "targets_consistent": 4,
                    "predicted_effect": "Drives proliferation and metabolism",
                },
            ],
            "drugs": [
                {
                    "name": "Alpelisib",
                    "prediction": "Inhibits pathway",
                    "z_score": -2.9,
                    "p_value": 0.0003,
                    "q_value": 0.002,
                    "targets_in_dataset": 4,
                    "mechanism": "PI3K inhibitor",
                    "clinical_indication": "PI3K pathway activation",
                },
                {
                    "name": "Everolimus",
                    "prediction": "Inhibits pathway",
                    "z_score": -2.6,
                    "p_value": 0.0008,
                    "q_value": 0.004,
                    "targets_in_dataset": 3,
                    "mechanism": "mTOR inhibitor",
                    "clinical_indication": "mTOR pathway activation",
                },
            ],
            "statistics": {
                "total_genes_analyzed": len(differential_genes),
                "kinases_tested": 10,
                "tfs_tested": 10,
                "drugs_tested": 8,
                "significant_kinases": 3,
                "significant_tfs": 2,
                "significant_drugs": 2,
                "method": "Fisher's exact test + activation Z-score",
                "fdr_method": "Benjamini-Hochberg",
            },
            "method": {
                "enrichment_test": "Fisher's exact test (one-sided)",
                "activation_score": "Z-score based on target expression direction",
                "interpretation": "Positive Z-score = Activated, Negative = Inhibited",
                "threshold": f"|Z| > {activation_zscore_threshold} considered significant",
            },
            "recommendation": "Focus on Alpelisib (PI3K inhibitor) - targets activated AKT/mTOR pathway",
            "status": "success (DRY_RUN mode)",
        }

    # Set default regulator types
    if regulator_types is None:
        regulator_types = ["kinase", "transcription_factor", "drug"]

    # Prepare differential gene data
    deg_genes = set(differential_genes.keys())
    deg_log2fc = {gene: data["log2fc"] for gene, data in differential_genes.items()}

    logger.info(f"Analyzing regulator types: {regulator_types}")

    results = {
        "kinases": [],
        "transcription_factors": [],
        "drugs": [],
    }

    # Analyze kinases
    if "kinase" in regulator_types:
        logger.info("Analyzing kinase regulators...")
        kinase_results = _analyze_regulators(
            KINASE_TARGETS,
            deg_genes,
            deg_log2fc,
            "kinase",
            activation_zscore_threshold,
        )
        results["kinases"] = kinase_results

    # Analyze transcription factors
    if "transcription_factor" in regulator_types:
        logger.info("Analyzing transcription factor regulators...")
        tf_results = _analyze_regulators(
            TF_TARGETS,
            deg_genes,
            deg_log2fc,
            "transcription_factor",
            activation_zscore_threshold,
        )
        results["transcription_factors"] = tf_results

    # Analyze drug responses
    if "drug" in regulator_types:
        logger.info("Analyzing drug response predictions...")
        drug_results = _analyze_regulators(
            DRUG_TARGETS,
            deg_genes,
            deg_log2fc,
            "drug",
            activation_zscore_threshold,
        )
        results["drugs"] = drug_results

    # Apply FDR correction across all regulators
    all_p_values = []
    all_regulators = []

    for reg_type, regs in results.items():
        for reg in regs:
            all_p_values.append(reg["p_value"])
            all_regulators.append((reg_type, reg))

    if len(all_p_values) > 0:
        _, q_values, _, _ = multipletests(all_p_values, method="fdr_bh")

        # Update q-values in results
        idx = 0
        for reg_type, regs in results.items():
            for reg in regs:
                reg["q_value"] = float(q_values[idx])
                idx += 1

        # Filter by FDR threshold
        results["kinases"] = [k for k in results["kinases"] if k["q_value"] <= fdr_threshold]
        results["transcription_factors"] = [
            tf for tf in results["transcription_factors"] if tf["q_value"] <= fdr_threshold
        ]
        results["drugs"] = [d for d in results["drugs"] if d["q_value"] <= fdr_threshold]

    # Generate statistics
    results["statistics"] = {
        "total_genes_analyzed": len(differential_genes),
        "kinases_tested": len(KINASE_TARGETS),
        "tfs_tested": len(TF_TARGETS),
        "drugs_tested": len(DRUG_TARGETS),
        "significant_kinases": len(results["kinases"]),
        "significant_tfs": len(results["transcription_factors"]),
        "significant_drugs": len(results["drugs"]),
        "method": "Fisher's exact test + activation Z-score",
        "fdr_method": "Benjamini-Hochberg",
        "fdr_threshold": fdr_threshold,
    }

    results["method"] = {
        "enrichment_test": "Fisher's exact test (one-sided)",
        "activation_score": "Z-score based on target expression direction",
        "interpretation": "Positive Z-score = Activated, Negative = Inhibited",
        "threshold": f"|Z| > {activation_zscore_threshold} considered significant",
    }

    # Generate recommendation
    if len(results["drugs"]) > 0:
        top_drug = results["drugs"][0]
        results["recommendation"] = (
            f"Consider {top_drug['name']} - targets pathway with "
            f"Z-score = {top_drug['z_score']:.2f}"
        )
    elif len(results["kinases"]) > 0:
        top_kinase = results["kinases"][0]
        results["recommendation"] = (
            f"Top predicted regulator: {top_kinase['name']} "
            f"({top_kinase['activation_state']})"
        )
    else:
        results["recommendation"] = "No significant upstream regulators found at current threshold"

    results["status"] = "success"

    logger.info(f"Upstream regulator analysis complete")
    logger.info(f"Found: {len(results['kinases'])} kinases, "
                f"{len(results['transcription_factors'])} TFs, "
                f"{len(results['drugs'])} drugs")

    return results


def _analyze_regulators(
    regulator_db: Dict[str, List[str]],
    deg_genes: set,
    deg_log2fc: Dict[str, float],
    regulator_type: str,
    zscore_threshold: float,
) -> List[Dict[str, Any]]:
    """Analyze regulators using enrichment and activation scoring.

    Args:
        regulator_db: Dict of regulator -> list of targets
        deg_genes: Set of differential genes
        deg_log2fc: Dict of gene -> log2 fold change
        regulator_type: "kinase", "transcription_factor", or "drug"
        zscore_threshold: Threshold for significant activation/inhibition

    Returns:
        List of regulator predictions with statistics
    """
    results = []

    total_genes_in_universe = 20000  # Approximate human protein-coding genes

    for regulator, targets in regulator_db.items():
        # Count overlapping targets
        targets_set = set(targets)
        overlap = targets_set & deg_genes

        if len(overlap) < 2:
            # Need at least 2 targets for meaningful analysis
            continue

        # Fisher's exact test for enrichment
        # 2x2 contingency table:
        # | In DEG | Not in DEG |
        # |--------|------------|
        # | a      | b          | <- Regulator targets
        # | c      | d          | <- Not regulator targets

        a = len(overlap)  # Targets in DEG
        b = len(targets_set) - a  # Targets not in DEG
        c = len(deg_genes) - a  # DEG not targets
        d = total_genes_in_universe - a - b - c  # Neither

        _, p_value = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

        # Calculate activation Z-score
        # Positive if targets are upregulated, negative if downregulated
        target_fc = [deg_log2fc[gene] for gene in overlap if gene in deg_log2fc]

        if len(target_fc) > 0:
            z_score = np.mean(target_fc) / (np.std(target_fc) + 1e-10) * np.sqrt(len(target_fc))
        else:
            z_score = 0.0

        # Determine activation state
        if abs(z_score) >= zscore_threshold:
            if z_score > 0:
                activation_state = "Activated"
            else:
                activation_state = "Inhibited"
        else:
            activation_state = "Ambiguous"

        # Calculate consistency (how many targets change in expected direction)
        if regulator_type == "drug":
            # For drugs, negative Z-score means drug would inhibit the pathway
            prediction = "Inhibits pathway" if z_score < 0 else "Mimics pathway"
        else:
            prediction = None

        targets_consistent = sum(1 for fc in target_fc if (fc > 0) == (z_score > 0))

        result = {
            "name": regulator,
            "activation_state": activation_state if regulator_type != "drug" else None,
            "prediction": prediction if regulator_type == "drug" else None,
            "z_score": float(z_score),
            "p_value": float(p_value),
            "q_value": None,  # Will be filled in after FDR correction
            "targets_in_dataset": len(overlap),
            "targets_consistent": int(targets_consistent),
        }

        results.append(result)

    # Sort by p-value
    results.sort(key=lambda x: x["p_value"])

    return results
