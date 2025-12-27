"""Multi-omics data integration functionality."""

import logging
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd

from ..config import config
from ..validation import (
    ValidationError,
    validate_multiomics_file,
    validate_metadata_file,
    validate_data_integration,
    format_validation_error,
)
from .utils import (
    align_samples,
    calculate_qc_metrics,
    filter_missing_features,
    load_omics_data,
    normalize_zscore,
    save_integrated_data,
)

logger = logging.getLogger(__name__)


def integrate_omics_data_impl(
    rna_path: str,
    protein_path: Optional[str] = None,
    phospho_path: Optional[str] = None,
    metadata_path: Optional[str] = None,
    normalize: bool = True,
    filter_missing: float = 0.5,
) -> Dict[str, Any]:
    """Implementation of multi-omics data integration.

    Args:
        rna_path: Path to RNA expression data
        protein_path: Path to protein abundance data (optional)
        phospho_path: Path to phosphorylation data (optional)
        metadata_path: Path to sample metadata (optional)
        normalize: Apply Z-score normalization
        filter_missing: Remove features with >X fraction missing

    Returns:
        Dictionary with integration results
    """
    logger.info("Starting multi-omics data integration")

    # =========================================================================
    # VALIDATION: Check input files before processing
    # =========================================================================

    validation_errors = []

    # Validate RNA file (required)
    logger.info(f"Validating RNA file: {rna_path}")
    rna_valid, rna_errors, rna_df = validate_multiomics_file(
        rna_path,
        required_columns=["gene_id"],  # Minimal requirement
        file_type="RNA",
        allow_missing_columns=True  # We'll check sample columns later
    )
    if not rna_valid:
        validation_errors.extend(rna_errors)
    elif rna_errors:  # Warnings
        for warning in rna_errors:
            logger.warning(warning)

    # Validate Protein file (optional)
    if protein_path:
        logger.info(f"Validating Protein file: {protein_path}")
        protein_valid, protein_errors, protein_df = validate_multiomics_file(
            protein_path,
            required_columns=["gene_id"],
            file_type="Protein",
            allow_missing_columns=True
        )
        if not protein_valid:
            validation_errors.extend(protein_errors)
        elif protein_errors:
            for warning in protein_errors:
                logger.warning(warning)

    # Validate Phospho file (optional)
    if phospho_path:
        logger.info(f"Validating Phospho file: {phospho_path}")
        phospho_valid, phospho_errors, phospho_df = validate_multiomics_file(
            phospho_path,
            required_columns=["gene_id"],
            file_type="Phospho",
            allow_missing_columns=True
        )
        if not phospho_valid:
            validation_errors.extend(phospho_errors)
        elif phospho_errors:
            for warning in phospho_errors:
                logger.warning(warning)

    # Validate Metadata file (optional)
    if metadata_path:
        logger.info(f"Validating Metadata file: {metadata_path}")
        meta_valid, meta_errors, meta_df = validate_metadata_file(
            metadata_path,
            required_columns=["Sample"]
        )
        if not meta_valid:
            validation_errors.extend(meta_errors)
        elif meta_errors:
            for warning in meta_errors:
                logger.warning(warning)

    # If any validation errors, raise exception with formatted message
    if validation_errors:
        error_msg = format_validation_error(validation_errors, rna_path)
        logger.error("Input validation failed")
        raise ValidationError(error_msg)

    logger.info("✅ All input files validated successfully")

    # =========================================================================
    # Load data files
    # =========================================================================

    dataframes = {}
    original_counts = {}

    logger.info(f"Loading RNA data from {rna_path}")
    dataframes["rna"] = load_omics_data(rna_path)
    original_counts["rna"] = dataframes["rna"].shape[0]

    if protein_path:
        logger.info(f"Loading protein data from {protein_path}")
        dataframes["protein"] = load_omics_data(protein_path)
        original_counts["protein"] = dataframes["protein"].shape[0]

    if phospho_path:
        logger.info(f"Loading phosphorylation data from {phospho_path}")
        dataframes["phospho"] = load_omics_data(phospho_path)
        original_counts["phospho"] = dataframes["phospho"].shape[0]

    # Load metadata if provided
    metadata = None
    if metadata_path:
        logger.info(f"Loading metadata from {metadata_path}")
        metadata = pd.read_csv(metadata_path)
        if "Sample" in metadata.columns:
            metadata = metadata.set_index("Sample")

    # Align samples across modalities
    logger.info("Aligning samples across modalities")
    aligned_data, common_samples = align_samples(dataframes)
    logger.info(f"Found {len(common_samples)} common samples")

    # Filter features with too much missing data
    if filter_missing < 1.0:
        logger.info(f"Filtering features with >{filter_missing*100}% missing data")
        aligned_data = {
            modality: filter_missing_features(df, filter_missing)
            for modality, df in aligned_data.items()
        }

    # Apply normalization
    if normalize:
        logger.info("Applying Z-score normalization")
        aligned_data = {
            modality: normalize_zscore(df)
            for modality, df in aligned_data.items()
        }

    # Calculate QC metrics
    qc_metrics = calculate_qc_metrics(aligned_data, original_counts)
    qc_metrics["normalization"] = "z-score" if normalize else "none"
    qc_metrics["missing_threshold"] = filter_missing

    # Save integrated data to cache
    cache_path = config.cache_dir / "integrated_data.pkl"
    config.cache_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Saving integrated data to {cache_path}")
    save_integrated_data(aligned_data, metadata, str(cache_path))

    # Prepare metadata summary
    metadata_summary = None
    if metadata is not None:
        # Filter metadata to common samples
        metadata_aligned = metadata.loc[
            metadata.index.intersection(common_samples)
        ]

        metadata_summary = {
            "samples": len(metadata_aligned),
        }

        # Count treatment groups if available
        if "Response" in metadata_aligned.columns:
            response_counts = metadata_aligned["Response"].value_counts().to_dict()
            for response, count in response_counts.items():
                metadata_summary[f"treatment_{response.lower()}"] = int(count)
        elif "Treatment" in metadata_aligned.columns:
            treatment_counts = metadata_aligned["Treatment"].value_counts().to_dict()
            for treatment, count in treatment_counts.items():
                metadata_summary[f"treatment_{treatment.lower()}"] = int(count)

    # Prepare response
    result = {
        "integrated_data": {
            modality: {
                "shape": list(df.shape),
                "features": list(df.index[:5]) + ["..."],  # Show first 5 features
                "samples": list(df.columns[:5]) + ["..."],  # Show first 5 samples
            }
            for modality, df in aligned_data.items()
        },
        "common_samples": common_samples,
        "feature_counts": {
            modality: df.shape[0] for modality, df in aligned_data.items()
        },
        "metadata": metadata_summary,
        "qc_metrics": qc_metrics,
        "cache_path": str(cache_path),
        "status": "success",
    }

    logger.info("Multi-omics integration complete")
    return result
