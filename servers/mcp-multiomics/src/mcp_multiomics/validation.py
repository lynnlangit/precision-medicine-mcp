"""Input validation utilities for multi-omics data."""

import pandas as pd
from pathlib import Path
from typing import Tuple, List, Dict, Any, Optional


class ValidationError(Exception):
    """Raised when input data fails validation."""
    pass


def validate_file_exists(file_path: str, file_type: str = "file") -> Tuple[bool, List[str]]:
    """
    Validate that a file exists and is accessible.

    Args:
        file_path: Path to file
        file_type: Description of file type (e.g., "RNA data", "metadata")

    Returns:
        (is_valid, error_messages)
    """
    errors = []

    if not file_path:
        errors.append(f"❌ No {file_type} path provided")
        errors.append(f"💡 Please specify a path to your {file_type} file")
        return False, errors

    path = Path(file_path)

    if not path.exists():
        errors.append(f"❌ {file_type.capitalize()} file not found: {file_path}")
        errors.append(f"💡 Check that the path is correct and the file exists")
        errors.append(f"💡 If using relative path, make sure you're in the correct directory")
        errors.append(f"💡 Current working directory: {Path.cwd()}")
        return False, errors

    if not path.is_file():
        errors.append(f"❌ Path is not a file: {file_path}")
        errors.append(f"💡 The path points to a directory, not a file")
        return False, errors

    if path.stat().st_size == 0:
        errors.append(f"❌ {file_type.capitalize()} file is empty: {file_path}")
        errors.append(f"💡 The file exists but contains no data")
        return False, errors

    return True, []


def validate_multiomics_file(
    file_path: str,
    required_columns: List[str],
    file_type: str = "multi-omics",
    allow_missing_columns: bool = False
) -> Tuple[bool, List[str], Optional[pd.DataFrame]]:
    """
    Validate multi-omics input file format.

    Args:
        file_path: Path to data file
        required_columns: List of required column names
        file_type: Description (e.g., "RNA", "protein", "phospho")
        allow_missing_columns: If True, warn but don't fail on missing columns

    Returns:
        (is_valid, error_messages, dataframe or None)
    """
    errors = []

    # Check file exists
    exists, exist_errors = validate_file_exists(file_path, file_type)
    if not exists:
        return False, exist_errors, None

    # Try to read as TSV/CSV
    df = None
    for sep, sep_name in [('\t', 'tab-separated'), (',', 'comma-separated')]:
        try:
            df = pd.read_csv(file_path, sep=sep, nrows=5)
            if len(df.columns) > 1:  # Successfully parsed
                break
        except Exception as e:
            continue

    if df is None or len(df.columns) <= 1:
        errors.append(f"❌ Cannot parse {file_type} file as tab-separated or comma-separated: {file_path}")
        errors.append(f"💡 File must be in TSV (tab-separated) or CSV (comma-separated) format")
        errors.append(f"💡 From Excel: File → Save As → Tab delimited text (.txt) or CSV (.csv)")
        errors.append(f"💡 Check that the file is not corrupted")
        return False, errors, None

    # Check required columns exist
    missing_cols = set(required_columns) - set(df.columns)
    if missing_cols:
        if allow_missing_columns:
            errors.append(f"⚠️  Warning: Missing optional columns: {list(missing_cols)}")
        else:
            errors.append(f"❌ Missing required columns in {file_type} file: {list(missing_cols)}")
            errors.append(f"💡 Found columns: {list(df.columns)[:10]}")
            if len(df.columns) > 10:
                errors.append(f"   ... and {len(df.columns) - 10} more")
            errors.append(f"💡 Required columns must include: {required_columns}")
            errors.append(f"💡 Column names are case-sensitive")
            return False, errors, None

    # Read full file for additional checks
    try:
        df_full = pd.read_csv(file_path, sep='\t' if '\t' in open(file_path, 'r').read(1000) else ',')
    except Exception as e:
        errors.append(f"⚠️  Warning: Could not read full file for validation: {str(e)}")
        return True, errors, df  # Still return success with partial data

    # Check for common issues
    warnings = []

    # Check for all-missing columns
    all_missing_cols = df_full.columns[df_full.isnull().all()].tolist()
    if all_missing_cols:
        warnings.append(f"⚠️  Columns with all missing values: {all_missing_cols[:5]}")
        warnings.append(f"💡 These columns will not contribute to analysis")

    # Check missing value percentage
    missing_pct = (df_full.isnull().sum() / len(df_full) * 100)
    high_missing = missing_pct[missing_pct > 50].index.tolist()
    if high_missing:
        warnings.append(f"⚠️  Columns with >50% missing values: {high_missing[:5]}")
        warnings.append(f"💡 High missing rates may affect analysis quality")
        warnings.append(f"💡 Consider removing these columns or using imputation")

    # Check for duplicate rows
    n_duplicates = df_full.duplicated().sum()
    if n_duplicates > 0:
        warnings.append(f"⚠️  Found {n_duplicates} duplicate rows")
        warnings.append(f"💡 Duplicates will be kept but may affect statistical tests")

    # Success with warnings
    if warnings:
        return True, warnings, df_full

    return True, [], df_full


def validate_metadata_file(
    file_path: str,
    required_columns: List[str] = ["Sample"],
) -> Tuple[bool, List[str], Optional[pd.DataFrame]]:
    """
    Validate metadata file format.

    Args:
        file_path: Path to metadata file
        required_columns: Required column names (default: ["Sample"])

    Returns:
        (is_valid, error_messages, dataframe or None)
    """
    errors = []

    # Check file exists
    exists, exist_errors = validate_file_exists(file_path, "metadata")
    if not exists:
        return False, exist_errors, None

    # Try to read file
    try:
        df = pd.read_csv(file_path, sep='\t' if '\t' in open(file_path, 'r').read(1000) else ',')
    except Exception as e:
        errors.append(f"❌ Cannot parse metadata file: {file_path}")
        errors.append(f"💡 Error: {str(e)}")
        errors.append(f"💡 Ensure file is in TSV or CSV format")
        return False, errors, None

    # Check required columns
    missing_cols = set(required_columns) - set(df.columns)
    if missing_cols:
        errors.append(f"❌ Missing required columns in metadata: {list(missing_cols)}")
        errors.append(f"💡 Found columns: {list(df.columns)}")
        errors.append(f"💡 Metadata must contain at minimum: {required_columns}")
        return False, errors, None

    # Check for common batch effect column
    if "Batch" not in df.columns:
        errors.append(f"⚠️  Warning: No 'Batch' column found in metadata")
        errors.append(f"💡 If your data has batch effects, add a 'Batch' column")
        errors.append(f"💡 This enables batch effect detection and correction")

    # Check sample names are unique
    if df["Sample"].duplicated().any():
        duplicates = df[df["Sample"].duplicated(keep=False)]["Sample"].tolist()
        errors.append(f"❌ Duplicate sample names in metadata: {duplicates[:5]}")
        errors.append(f"💡 Each sample must have a unique name")
        return False, errors, None

    return True, errors, df


def validate_data_integration(
    rna_samples: List[str],
    protein_samples: Optional[List[str]] = None,
    phospho_samples: Optional[List[str]] = None,
    metadata_samples: Optional[List[str]] = None
) -> Tuple[bool, List[str], Dict[str, Any]]:
    """
    Validate sample consistency across multiple modalities.

    Args:
        rna_samples: List of sample names from RNA data
        protein_samples: List of sample names from protein data (optional)
        phospho_samples: List of sample names from phospho data (optional)
        metadata_samples: List of sample names from metadata (optional)

    Returns:
        (is_valid, messages, info_dict)
    """
    messages = []
    info = {
        "rna_samples": len(rna_samples),
        "protein_samples": len(protein_samples) if protein_samples else 0,
        "phospho_samples": len(phospho_samples) if phospho_samples else 0,
        "metadata_samples": len(metadata_samples) if metadata_samples else 0,
        "common_samples": [],
        "rna_only": [],
        "protein_only": [],
        "phospho_only": []
    }

    # Find common samples
    all_modalities = [rna_samples]
    if protein_samples:
        all_modalities.append(protein_samples)
    if phospho_samples:
        all_modalities.append(phospho_samples)

    common = set(rna_samples)
    for samples in all_modalities[1:]:
        common = common.intersection(set(samples))

    info["common_samples"] = sorted(list(common))

    if len(common) == 0:
        messages.append(f"❌ No common samples found across modalities")
        messages.append(f"💡 RNA samples: {rna_samples[:5]}")
        if protein_samples:
            messages.append(f"💡 Protein samples: {protein_samples[:5]}")
        if phospho_samples:
            messages.append(f"💡 Phospho samples: {phospho_samples[:5]}")
        messages.append(f"💡 Check that sample names match exactly (case-sensitive)")
        messages.append(f"💡 Common issues: 'Sample_01' vs 'Sample-01' vs 'sample_01'")
        return False, messages, info

    if len(common) < 3:
        messages.append(f"⚠️  Warning: Very few common samples ({len(common)})")
        messages.append(f"💡 Multi-omics analysis requires at least 3 samples")
        messages.append(f"💡 More samples (10+) recommended for statistical power")

    # Check for modality-specific samples
    rna_only = set(rna_samples) - common
    if rna_only and len(rna_only) <= 5:
        info["rna_only"] = list(rna_only)
        messages.append(f"ℹ️  Samples only in RNA: {list(rna_only)}")
    elif rna_only:
        info["rna_only"] = list(rna_only)
        messages.append(f"ℹ️  {len(rna_only)} samples only in RNA")

    if protein_samples:
        protein_only = set(protein_samples) - common
        if protein_only and len(protein_only) <= 5:
            info["protein_only"] = list(protein_only)
            messages.append(f"ℹ️  Samples only in Protein: {list(protein_only)}")
        elif protein_only:
            info["protein_only"] = list(protein_only)
            messages.append(f"ℹ️  {len(protein_only)} samples only in Protein")

    # Check metadata consistency
    if metadata_samples:
        missing_in_metadata = common - set(metadata_samples)
        if missing_in_metadata:
            messages.append(f"⚠️  Samples in data but not in metadata: {list(missing_in_metadata)[:5]}")
            messages.append(f"💡 Add these samples to metadata or remove from data")

    # Success message
    messages.append(f"✅ Found {len(common)} common samples across all modalities")

    return True, messages, info


def format_validation_error(errors: List[str], file_path: str = None) -> str:
    """
    Format validation errors into a user-friendly message.

    Args:
        errors: List of error messages
        file_path: Optional file path to include in header

    Returns:
        Formatted error message string
    """
    header = "DATA VALIDATION FAILED"
    if file_path:
        header += f" - {Path(file_path).name}"

    formatted = f"""
╔═══════════════════════════════════════════════════════════════════════════╗
║  {header:^73}  ║
╚═══════════════════════════════════════════════════════════════════════════╝

"""

    for error in errors:
        formatted += error + "\n"

    formatted += """
═══════════════════════════════════════════════════════════════════════════

Need help? Check the data format guide:
https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/DATA_FORMATS.md

Or see example data files:
https://github.com/lynnlangit/precision-medicine-mcp/tree/main/data/multiomics
"""

    return formatted
