"""Input validation utilities for genomic data (FASTQ, VCF, etc.)."""

import gzip
from pathlib import Path
from typing import Tuple, List, Optional


class ValidationError(Exception):
    """Raised when input data fails validation."""
    pass


def validate_file_exists(file_path: str, file_type: str = "file") -> Tuple[bool, List[str]]:
    """
    Validate that a file exists and is accessible.

    Args:
        file_path: Path to file
        file_type: Description of file type

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
        errors.append(f"💡 If using relative path, ensure you're in the correct directory")
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


def validate_fastq_file(file_path: str) -> Tuple[bool, List[str], dict]:
    """
    Validate FASTQ file format and quality encoding.

    Args:
        file_path: Path to FASTQ file (.fastq, .fq, .fastq.gz, .fq.gz)

    Returns:
        (is_valid, messages, info_dict)
    """
    errors = []
    warnings = []
    info = {
        "is_gzipped": False,
        "quality_encoding": None,
        "first_read_id": None,
        "estimated_reads": None
    }

    # Check file exists
    exists, exist_errors = validate_file_exists(file_path, "FASTQ")
    if not exists:
        return False, exist_errors, info

    path = Path(file_path)

    # Check file extension
    valid_extensions = [".fastq", ".fq", ".fastq.gz", ".fq.gz"]
    if not any(str(path).endswith(ext) for ext in valid_extensions):
        warnings.append(f"⚠️  Unusual file extension: {path.suffix}")
        warnings.append(f"💡 Expected extensions: {', '.join(valid_extensions)}")

    # Determine if gzipped
    is_gzipped = str(path).endswith('.gz')
    info["is_gzipped"] = is_gzipped

    # Read first few lines
    try:
        if is_gzipped:
            with gzip.open(path, 'rt') as f:
                lines = [f.readline() for _ in range(8)]
        else:
            with open(path, 'r') as f:
                lines = [f.readline() for _ in range(8)]
    except Exception as e:
        errors.append(f"❌ Cannot read FASTQ file: {str(e)}")
        errors.append(f"💡 File may be corrupted or not in FASTQ format")
        if is_gzipped:
            errors.append(f"💡 If file is gzipped, ensure it's not corrupted: gunzip -t {file_path}")
        return False, errors, info

    # Validate FASTQ format (4 lines per read)
    if not lines[0].startswith('@'):
        errors.append(f"❌ Invalid FASTQ format: First line should start with '@'")
        errors.append(f"💡 Found: {lines[0][:50]}")
        errors.append(f"💡 FASTQ format: @read_id, sequence, +, quality_scores")
        return False, errors, info

    if len(lines) < 4:
        errors.append(f"❌ File too short to be valid FASTQ (< 4 lines)")
        return False, errors, info

    if not lines[2].startswith('+'):
        errors.append(f"❌ Invalid FASTQ format: Third line should start with '+'")
        errors.append(f"💡 Found: {lines[2][:50]}")
        return False, errors, info

    # Extract first read ID
    info["first_read_id"] = lines[0].strip()[1:]  # Remove '@'

    # Validate sequence and quality length match
    seq = lines[1].strip()
    qual = lines[3].strip()

    if len(seq) != len(qual):
        errors.append(f"❌ Sequence and quality length mismatch in first read")
        errors.append(f"💡 Sequence length: {len(seq)}, Quality length: {len(qual)}")
        errors.append(f"💡 Each base in sequence must have a quality score")
        return False, errors, info

    # Detect quality score encoding
    if qual:
        qual_min = ord(min(qual))
        qual_max = ord(max(qual))

        # Phred+33 (Sanger, standard): 33-73 (typical), max 126
        # Phred+64 (old Illumina): 64-104 (typical), max 126

        if qual_min >= 33 and qual_max <= 73:
            info["quality_encoding"] = "Phred+33 (Sanger/Illumina 1.8+)"
        elif qual_min >= 64 and qual_max <= 104:
            info["quality_encoding"] = "Phred+64 (old Illumina)"
            warnings.append(f"⚠️  Old Illumina quality encoding detected (Phred+64)")
            warnings.append(f"💡 Modern tools expect Phred+33")
            warnings.append(f"💡 Consider converting with: seqtk seq -Q64 -V input.fq > output.fq")
        elif qual_min < 33:
            errors.append(f"❌ Invalid quality scores (ASCII < 33): min={qual_min}")
            errors.append(f"💡 Quality scores must be ASCII 33-126")
            return False, errors, info
        elif qual_max > 126:
            errors.append(f"❌ Invalid quality scores (ASCII > 126): max={qual_max}")
            return False, errors, info
        else:
            info["quality_encoding"] = "Phred+33 (extended range)"

    # Estimate number of reads (rough estimate from file size)
    file_size_bytes = path.stat().st_size
    if is_gzipped:
        # Assume ~5x compression for FASTQ
        estimated_size = file_size_bytes * 5
    else:
        estimated_size = file_size_bytes

    # Assume average read length ~150bp, 4 lines per read, ~300 bytes per read
    estimated_reads = estimated_size // 300
    info["estimated_reads"] = estimated_reads

    if estimated_reads < 100000:  # < 100K reads
        warnings.append(f"⚠️  Very small FASTQ file (~{estimated_reads:,} reads)")
        warnings.append(f"💡 Typical sequencing runs have millions of reads")
        warnings.append(f"💡 Is this a test file or incomplete download?")

    # Success messages
    messages = []
    if warnings:
        messages.extend(warnings)

    messages.append(f"✅ FASTQ file format is valid")
    messages.append(f"ℹ️  Quality encoding: {info['quality_encoding']}")
    messages.append(f"ℹ️  Estimated reads: ~{estimated_reads:,}")
    messages.append(f"ℹ️  Compressed: {'Yes' if is_gzipped else 'No'}")

    return True, messages, info


def validate_vcf_file(file_path: str) -> Tuple[bool, List[str], dict]:
    """
    Validate VCF (Variant Call Format) file.

    Args:
        file_path: Path to VCF file (.vcf or .vcf.gz)

    Returns:
        (is_valid, messages, info_dict)
    """
    errors = []
    warnings = []
    info = {
        "is_gzipped": False,
        "vcf_version": None,
        "has_header": False,
        "sample_count": 0
    }

    # Check file exists
    exists, exist_errors = validate_file_exists(file_path, "VCF")
    if not exists:
        return False, exist_errors, info

    path = Path(file_path)

    # Check file extension
    if not (str(path).endswith('.vcf') or str(path).endswith('.vcf.gz')):
        warnings.append(f"⚠️  Unusual file extension: {path.suffix}")
        warnings.append(f"💡 Expected: .vcf or .vcf.gz")

    # Determine if gzipped
    is_gzipped = str(path).endswith('.gz')
    info["is_gzipped"] = is_gzipped

    # Read file
    try:
        if is_gzipped:
            with gzip.open(path, 'rt') as f:
                lines = [f.readline() for _ in range(100)]  # Read first 100 lines
        else:
            with open(path, 'r') as f:
                lines = [f.readline() for _ in range(100)]
    except Exception as e:
        errors.append(f"❌ Cannot read VCF file: {str(e)}")
        return False, errors, info

    # Check for VCF header
    if not lines[0].startswith('##fileformat=VCF'):
        errors.append(f"❌ Invalid VCF format: Missing ##fileformat=VCF header")
        errors.append(f"💡 Found: {lines[0][:50]}")
        errors.append(f"💡 VCF files must start with ##fileformat=VCFv4.x")
        return False, errors, info

    info["has_header"] = True

    # Extract VCF version
    version_line = lines[0].strip()
    if '=' in version_line:
        info["vcf_version"] = version_line.split('=')[1]

    # Find column header line (#CHROM...)
    column_line = None
    for line in lines:
        if line.startswith('#CHROM'):
            column_line = line.strip()
            break

    if not column_line:
        errors.append(f"❌ Missing column header line (#CHROM...)")
        errors.append(f"💡 VCF must have header line starting with #CHROM")
        return False, errors, info

    # Count samples (columns after FORMAT)
    columns = column_line.split('\t')
    required_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    if len(columns) < 8:
        errors.append(f"❌ Insufficient columns in VCF header")
        errors.append(f"💡 Found {len(columns)} columns, need at least 8")
        errors.append(f"💡 Required: {', '.join(required_cols)}")
        return False, errors, info

    # Count samples (columns after FORMAT)
    if len(columns) > 9:  # Has FORMAT + samples
        info["sample_count"] = len(columns) - 9
    else:
        info["sample_count"] = 0
        warnings.append(f"⚠️  No sample columns found (sites-only VCF)")

    messages = []
    if warnings:
        messages.extend(warnings)

    messages.append(f"✅ VCF file format is valid")
    messages.append(f"ℹ️  VCF version: {info['vcf_version']}")
    messages.append(f"ℹ️  Samples: {info['sample_count']}")

    return True, messages, info


def format_validation_error(errors: List[str], file_path: str = None) -> str:
    """
    Format validation errors into user-friendly message.

    Args:
        errors: List of error messages
        file_path: Optional file path

    Returns:
        Formatted error message string
    """
    header = "FILE VALIDATION FAILED"
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

Need help? Check the documentation:
- FASTQ format: https://en.wikipedia.org/wiki/FASTQ_format
- VCF format: https://samtools.github.io/hts-specs/VCFv4.2.pdf

Or see example data:
https://github.com/lynnlangit/precision-medicine-mcp/tree/main/data
"""

    return formatted
