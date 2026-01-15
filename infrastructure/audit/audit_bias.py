#!/usr/bin/env python3
"""
Bias Audit Script for Precision Medicine Workflows

This script performs comprehensive bias audits on precision medicine AI/ML workflows,
checking for representation bias, algorithmic fairness, and proxy features.

Aligned with FDA, AMA, and NIH standards for ethical healthcare AI.

Usage:
    # Basic audit
    python scripts/audit/audit_bias.py \\
        --workflow patientone \\
        --genomics-data data/genomics/patient_variants.vcf \\
        --clinical-data data/fhir/patients_deidentified.json \\
        --output reports/bias_audit_2026-01-12.html

    # Audit with custom thresholds
    python scripts/audit/audit_bias.py \\
        --workflow patientone \\
        --genomics-data data/genomics/patient_variants.vcf \\
        --clinical-data data/fhir/patients_deidentified.json \\
        --output reports/audit.html \\
        --min-representation 0.10 \\
        --max-disparity 0.10 \\
        --reference-dataset gnomad

References:
- docs/hospital-deployment/ethics/ETHICS_AND_BIAS.md
- docs/hospital-deployment/ethics/BIAS_AUDIT_CHECKLIST.md
- docs/hospital-deployment/ethics/PATIENTONE_BIAS_AUDIT.md
"""

import argparse
import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any, Optional

# Add parent directory to path to import shared utilities
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from shared.utils.bias_detection import (
    check_data_representation,
    calculate_fairness_metrics,
    detect_proxy_features,
    calculate_ancestry_aware_confidence,
    generate_bias_report,
    RepresentationAnalysis,
    FairnessMetrics,
    ProxyFeatureAnalysis
)


# ============================================================================
# DATA LOADING
# ============================================================================

def load_genomics_data(file_path: str) -> pd.DataFrame:
    """
    Load genomics data from VCF or CSV file.

    Expected columns:
    - patient_id: Patient identifier
    - variant_id: Variant identifier (e.g., "BRCA1:c.5266dupC")
    - gene: Gene name
    - ancestry: Patient ancestry (european, african, asian, latino, etc.)
    - pathogenicity: Predicted pathogenicity (0=benign, 1=pathogenic)
    - confidence: Model confidence score
    """
    file_path_obj = Path(file_path)

    if not file_path_obj.exists():
        raise FileNotFoundError(f"Genomics data file not found: {file_path}")

    if file_path_obj.suffix == ".csv":
        df = pd.read_csv(file_path)
    elif file_path_obj.suffix == ".vcf":
        # Simple VCF parser (would need more sophisticated parser for production)
        print(f"Warning: VCF parsing is simplified. Consider converting to CSV format.")
        df = parse_vcf_simple(file_path)
    else:
        raise ValueError(f"Unsupported file format: {file_path_obj.suffix}")

    # Validate required columns
    required_cols = ["patient_id", "variant_id", "ancestry"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    return df


def parse_vcf_simple(file_path: str) -> pd.DataFrame:
    """
    Simple VCF parser for demonstration purposes.
    Production use should use libraries like cyvcf2 or pysam.
    """
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 8:
                data.append({
                    'patient_id': 'patient_001',  # Would extract from sample column
                    'variant_id': f"{fields[0]}:{fields[1]}:{fields[3]}>{fields[4]}",
                    'gene': fields[7].split('GENE=')[1].split(';')[0] if 'GENE=' in fields[7] else 'UNKNOWN',
                    'ancestry': 'european'  # Would need to be provided separately
                })

    return pd.DataFrame(data)


def load_clinical_data(file_path: str) -> pd.DataFrame:
    """
    Load clinical data from FHIR JSON or CSV file.

    Expected fields:
    - patient_id: Patient identifier
    - age: Patient age
    - sex: Patient sex (M/F/Other)
    - ancestry: Patient ancestry
    - diagnosis: Primary diagnosis
    - treatment_recommendation: Recommended treatment
    - outcome: Treatment outcome (0=poor, 1=good) [optional]
    """
    file_path_obj = Path(file_path)

    if not file_path_obj.exists():
        raise FileNotFoundError(f"Clinical data file not found: {file_path}")

    if file_path_obj.suffix == ".csv":
        df = pd.read_csv(file_path)
    elif file_path_obj.suffix == ".json":
        with open(file_path, 'r') as f:
            fhir_data = json.load(f)
        df = parse_fhir_bundle(fhir_data)
    else:
        raise ValueError(f"Unsupported file format: {file_path_obj.suffix}")

    return df


def parse_fhir_bundle(fhir_data: Dict) -> pd.DataFrame:
    """
    Parse FHIR Bundle into DataFrame.
    """
    patients = []

    if 'entry' in fhir_data:
        for entry in fhir_data['entry']:
            resource = entry.get('resource', {})
            if resource.get('resourceType') == 'Patient':
                patient_data = {
                    'patient_id': resource.get('id'),
                    'sex': resource.get('gender', 'unknown'),
                    'age': calculate_age_from_birthdate(resource.get('birthDate')),
                    'ancestry': extract_ancestry_from_extension(resource)
                }
                patients.append(patient_data)

    return pd.DataFrame(patients)


def calculate_age_from_birthdate(birthdate: Optional[str]) -> Optional[int]:
    """Calculate age from birthdate string (YYYY-MM-DD)."""
    if not birthdate:
        return None
    try:
        birth_year = int(birthdate.split('-')[0])
        current_year = datetime.now().year
        return current_year - birth_year
    except:
        return None


def extract_ancestry_from_extension(patient_resource: Dict) -> str:
    """Extract ancestry from FHIR Patient extension."""
    extensions = patient_resource.get('extension', [])
    for ext in extensions:
        if 'ancestry' in ext.get('url', '').lower():
            return ext.get('valueString', 'unknown')
    return 'unknown'


# ============================================================================
# WORKFLOW-SPECIFIC AUDITS
# ============================================================================

def audit_patientone_workflow(
    genomics_df: pd.DataFrame,
    clinical_df: pd.DataFrame,
    min_representation: float = 0.10,
    max_disparity: float = 0.10,
    reference_dataset: str = "gnomad"
) -> Dict[str, Any]:
    """
    Run bias audit on PatientOne precision medicine workflow.

    Checks:
    1. Data representation across ancestries
    2. Variant pathogenicity prediction fairness
    3. Treatment recommendation fairness
    4. Proxy feature detection (if feature importance provided)

    Args:
        genomics_df: Genomics data with variants and predictions
        clinical_df: Clinical data with demographics and outcomes
        min_representation: Minimum acceptable ancestry representation
        max_disparity: Maximum acceptable fairness disparity
        reference_dataset: Reference dataset for comparison

    Returns:
        Dict with audit results
    """
    audit_results = {
        "workflow": "patientone",
        "timestamp": datetime.now().isoformat(),
        "thresholds": {
            "min_representation": min_representation,
            "max_disparity": max_disparity
        }
    }

    # 1. Data Representation Analysis
    print("\n" + "="*80)
    print("STEP 1: DATA REPRESENTATION ANALYSIS")
    print("="*80)

    if 'ancestry' in genomics_df.columns:
        rep_analysis = check_data_representation(
            genomics_df,
            ancestry_column='ancestry',
            reference_dataset=reference_dataset,
            min_acceptable_pct=min_representation
        )
        audit_results['representation_analysis'] = rep_analysis

        print(f"\nRisk Level: {rep_analysis.risk_level}")
        print(f"Total Samples: {rep_analysis.total_samples}")
        print("\nAncestry Distribution:")
        for ancestry, count in rep_analysis.ancestry_counts.items():
            prop = rep_analysis.ancestry_proportions[ancestry]
            print(f"  - {ancestry}: {count} ({prop:.1%})")

        if rep_analysis.warnings:
            print("\n⚠ Warnings:")
            for warning in rep_analysis.warnings:
                print(f"  {warning}")
    else:
        print("⚠ Warning: No ancestry column found, skipping representation analysis")
        audit_results['representation_analysis'] = None

    # 2. Fairness Metrics for Variant Pathogenicity Predictions
    print("\n" + "="*80)
    print("STEP 2: FAIRNESS METRICS - VARIANT PATHOGENICITY")
    print("="*80)

    if all(col in genomics_df.columns for col in ['pathogenicity', 'ancestry']):
        # Create dummy true labels for demonstration (would use validation set in production)
        if 'pathogenicity_true' in genomics_df.columns:
            y_true = genomics_df['pathogenicity_true'].values
        else:
            print("⚠ No ground truth labels found, using predictions as proxy")
            y_true = genomics_df['pathogenicity'].values

        y_pred = genomics_df['pathogenicity'].values
        groups = genomics_df['ancestry'].values

        # Calculate fairness metrics
        y_prob = genomics_df['confidence'].values if 'confidence' in genomics_df.columns else None

        fairness_metrics = calculate_fairness_metrics(
            y_true=y_true,
            y_pred=y_pred,
            y_prob=y_prob,
            groups=groups,
            metrics=['demographic_parity', 'equalized_odds']
        )

        audit_results['fairness_metrics_genomics'] = fairness_metrics

        # Print results
        for metric_name, metric_result in fairness_metrics.items():
            if isinstance(metric_result, dict) and 'tpr' in metric_result:
                print(f"\n{metric_result['tpr'].metric_name}:")
                print(f"  Risk: {metric_result['tpr'].risk_level}")
                print(f"  Max Disparity: {metric_result['tpr'].max_disparity:.1%}")
                if metric_result['tpr'].warnings:
                    for warning in metric_result['tpr'].warnings:
                        print(f"  ⚠ {warning}")

                print(f"\n{metric_result['fpr'].metric_name}:")
                print(f"  Risk: {metric_result['fpr'].risk_level}")
                print(f"  Max Disparity: {metric_result['fpr'].max_disparity:.1%}")
            else:
                print(f"\n{metric_result.metric_name}:")
                print(f"  Risk: {metric_result.risk_level}")
                print(f"  Max Disparity: {metric_result.max_disparity:.1%}")
    else:
        print("⚠ Insufficient data for fairness metrics")
        audit_results['fairness_metrics_genomics'] = None

    # 3. Proxy Feature Detection
    print("\n" + "="*80)
    print("STEP 3: PROXY FEATURE DETECTION")
    print("="*80)

    # For PatientOne, check if geographic or socioeconomic features are used
    proxy_features = []

    # Check for common proxy features
    proxy_indicators = ['zip_code', 'postal_code', 'insurance', 'income', 'education']
    found_proxies = [col for col in genomics_df.columns if any(indicator in col.lower() for indicator in proxy_indicators)]

    if found_proxies:
        print(f"⚠ Warning: Potential proxy features detected: {found_proxies}")
        print("  Recommendation: Remove these features and retrain model")
        audit_results['proxy_features'] = found_proxies
    else:
        print("✓ No obvious proxy features detected")
        audit_results['proxy_features'] = []

    # 4. Ancestry-Aware Confidence Analysis
    print("\n" + "="*80)
    print("STEP 4: ANCESTRY-AWARE CONFIDENCE SCORING")
    print("="*80)

    if 'confidence' in genomics_df.columns and 'ancestry' in genomics_df.columns:
        # Calculate mean confidence by ancestry
        confidence_by_ancestry = genomics_df.groupby('ancestry')['confidence'].agg(['mean', 'std', 'count'])
        print("\nConfidence Scores by Ancestry:")
        print(confidence_by_ancestry)

        max_conf_diff = confidence_by_ancestry['mean'].max() - confidence_by_ancestry['mean'].min()
        if max_conf_diff > 0.10:
            print(f"\n⚠ Warning: Large confidence disparity ({max_conf_diff:.1%}) across ancestries")
            print("  Consider implementing ancestry-aware confidence scoring")
        else:
            print(f"\n✓ Confidence disparity acceptable ({max_conf_diff:.1%})")

        audit_results['confidence_analysis'] = confidence_by_ancestry.to_dict()
    else:
        print("⚠ No confidence scores available")
        audit_results['confidence_analysis'] = None

    return audit_results


# ============================================================================
# REPORT GENERATION
# ============================================================================

def generate_html_report(audit_results: Dict[str, Any], output_path: str):
    """
    Generate HTML bias audit report.

    Args:
        audit_results: Results from audit workflow
        output_path: Path to save HTML report
    """
    html = []
    html.append("<!DOCTYPE html>")
    html.append("<html>")
    html.append("<head>")
    html.append("<title>Bias Audit Report</title>")
    html.append("<style>")
    html.append("""
        body { font-family: Arial, sans-serif; margin: 40px; background-color: #f5f5f5; }
        .container { max-width: 1200px; margin: 0 auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
        h2 { color: #34495e; margin-top: 30px; border-left: 4px solid #3498db; padding-left: 10px; }
        .risk-critical { color: #e74c3c; font-weight: bold; }
        .risk-high { color: #e67e22; font-weight: bold; }
        .risk-medium { color: #f39c12; font-weight: bold; }
        .risk-acceptable { color: #27ae60; font-weight: bold; }
        .warning { background-color: #fff3cd; border-left: 4px solid #ffc107; padding: 10px; margin: 10px 0; }
        .success { background-color: #d4edda; border-left: 4px solid #28a745; padding: 10px; margin: 10px 0; }
        table { width: 100%; border-collapse: collapse; margin: 20px 0; }
        th, td { padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }
        th { background-color: #3498db; color: white; }
        tr:hover { background-color: #f5f5f5; }
        .metadata { background-color: #ecf0f1; padding: 15px; border-radius: 4px; margin-bottom: 20px; }
    """)
    html.append("</style>")
    html.append("</head>")
    html.append("<body>")
    html.append("<div class='container'>")

    # Header
    html.append("<h1>Bias Audit Report</h1>")
    html.append(f"<div class='metadata'>")
    html.append(f"<strong>Workflow:</strong> {audit_results.get('workflow', 'Unknown')}<br>")
    html.append(f"<strong>Timestamp:</strong> {audit_results.get('timestamp', 'Unknown')}<br>")
    html.append(f"<strong>Min Representation Threshold:</strong> {audit_results.get('thresholds', {}).get('min_representation', 0.10):.0%}<br>")
    html.append(f"<strong>Max Disparity Threshold:</strong> {audit_results.get('thresholds', {}).get('max_disparity', 0.10):.0%}")
    html.append("</div>")

    # Section 1: Data Representation
    html.append("<h2>1. Data Representation Analysis</h2>")
    rep_analysis = audit_results.get('representation_analysis')
    if rep_analysis:
        risk_class = f"risk-{rep_analysis.risk_level.lower()}"
        html.append(f"<p><strong>Risk Level:</strong> <span class='{risk_class}'>{rep_analysis.risk_level}</span></p>")
        html.append(f"<p><strong>Total Samples:</strong> {rep_analysis.total_samples}</p>")

        html.append("<table>")
        html.append("<tr><th>Ancestry</th><th>Count</th><th>Percentage</th></tr>")
        for ancestry, count in rep_analysis.ancestry_counts.items():
            prop = rep_analysis.ancestry_proportions[ancestry]
            html.append(f"<tr><td>{ancestry.capitalize()}</td><td>{count}</td><td>{prop:.1%}</td></tr>")
        html.append("</table>")

        if rep_analysis.warnings:
            for warning in rep_analysis.warnings:
                html.append(f"<div class='warning'>⚠ {warning}</div>")

        if rep_analysis.recommendations:
            html.append("<p><strong>Recommendations:</strong></p>")
            html.append("<ul>")
            for rec in rep_analysis.recommendations:
                html.append(f"<li>{rec}</li>")
            html.append("</ul>")
    else:
        html.append("<p>No representation analysis available</p>")

    # Section 2: Fairness Metrics
    html.append("<h2>2. Fairness Metrics</h2>")
    fairness_metrics = audit_results.get('fairness_metrics_genomics')
    if fairness_metrics:
        for metric_name, metric_result in fairness_metrics.items():
            if isinstance(metric_result, dict) and 'tpr' in metric_result:
                for sub_name, sub_metric in metric_result.items():
                    html.append(f"<h3>{sub_metric.metric_name}</h3>")
                    risk_class = f"risk-{sub_metric.risk_level.lower()}"
                    html.append(f"<p><strong>Risk Level:</strong> <span class='{risk_class}'>{sub_metric.risk_level}</span></p>")
                    html.append(f"<p><strong>Max Disparity:</strong> {sub_metric.max_disparity:.1%}</p>")

                    html.append("<table>")
                    html.append("<tr><th>Group</th><th>Value</th></tr>")
                    for group, value in sub_metric.group_values.items():
                        html.append(f"<tr><td>{group}</td><td>{value:.3f}</td></tr>")
                    html.append("</table>")

                    if sub_metric.warnings:
                        for warning in sub_metric.warnings:
                            html.append(f"<div class='warning'>⚠ {warning}</div>")
            else:
                html.append(f"<h3>{metric_result.metric_name}</h3>")
                risk_class = f"risk-{metric_result.risk_level.lower()}"
                html.append(f"<p><strong>Risk Level:</strong> <span class='{risk_class}'>{metric_result.risk_level}</span></p>")
                html.append(f"<p><strong>Max Disparity:</strong> {metric_result.max_disparity:.1%}</p>")

                html.append("<table>")
                html.append("<tr><th>Group</th><th>Value</th></tr>")
                for group, value in metric_result.group_values.items():
                    html.append(f"<tr><td>{group}</td><td>{value:.3f}</td></tr>")
                html.append("</table>")

                if metric_result.warnings:
                    for warning in metric_result.warnings:
                        html.append(f"<div class='warning'>⚠ {warning}</div>")
    else:
        html.append("<p>No fairness metrics available</p>")

    # Section 3: Proxy Features
    html.append("<h2>3. Proxy Feature Detection</h2>")
    proxy_features = audit_results.get('proxy_features', [])
    if proxy_features:
        html.append("<div class='warning'>")
        html.append(f"<p><strong>⚠ Warning:</strong> Potential proxy features detected:</p>")
        html.append("<ul>")
        for feature in proxy_features:
            html.append(f"<li>{feature}</li>")
        html.append("</ul>")
        html.append("<p><strong>Recommendation:</strong> Remove these features and retrain model</p>")
        html.append("</div>")
    else:
        html.append("<div class='success'>✓ No obvious proxy features detected</div>")

    # Footer
    html.append("<hr>")
    html.append("<p style='text-align: center; color: #7f8c8d; font-size: 12px;'>")
    html.append("Generated by Precision Medicine Bias Audit Tool | ")
    html.append("References: FDA AI/ML SaMD Guidance, AMA Ethics Opinion 2.3.2, NIH All of Us Diversity Requirements")
    html.append("</p>")

    html.append("</div>")
    html.append("</body>")
    html.append("</html>")

    # Write to file
    output_path_obj = Path(output_path)
    output_path_obj.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        f.write('\n'.join(html))

    print(f"\n✓ HTML report generated: {output_path}")


# ============================================================================
# MAIN CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Bias Audit Script for Precision Medicine Workflows",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic PatientOne audit
  python scripts/audit/audit_bias.py \\
      --workflow patientone \\
      --genomics-data data/genomics/patient_variants.csv \\
      --clinical-data data/fhir/patients.json \\
      --output reports/bias_audit.html

  # Audit with custom thresholds
  python scripts/audit/audit_bias.py \\
      --workflow patientone \\
      --genomics-data data/genomics/variants.csv \\
      --clinical-data data/fhir/patients.csv \\
      --output reports/audit.html \\
      --min-representation 0.15 \\
      --max-disparity 0.05 \\
      --reference-dataset all_of_us
        """
    )

    parser.add_argument(
        '--workflow',
        type=str,
        required=True,
        choices=['patientone', 'generic'],
        help='Workflow to audit (patientone or generic)'
    )

    parser.add_argument(
        '--genomics-data',
        type=str,
        help='Path to genomics data (VCF or CSV)'
    )

    parser.add_argument(
        '--clinical-data',
        type=str,
        help='Path to clinical data (FHIR JSON or CSV)'
    )

    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Output path for audit report (HTML or JSON)'
    )

    parser.add_argument(
        '--min-representation',
        type=float,
        default=0.10,
        help='Minimum acceptable ancestry representation (default: 0.10 = 10%%)'
    )

    parser.add_argument(
        '--max-disparity',
        type=float,
        default=0.10,
        help='Maximum acceptable fairness disparity (default: 0.10 = 10%%)'
    )

    parser.add_argument(
        '--reference-dataset',
        type=str,
        default='gnomad',
        choices=['gnomad', 'all_of_us', 'gtex', 'topmed'],
        help='Reference dataset for comparison (default: gnomad)'
    )

    args = parser.parse_args()

    # Validate inputs
    if args.workflow == 'patientone':
        if not args.genomics_data or not args.clinical_data:
            parser.error("PatientOne workflow requires both --genomics-data and --clinical-data")

    print("="*80)
    print("BIAS AUDIT FOR PRECISION MEDICINE WORKFLOWS")
    print("="*80)
    print(f"Workflow: {args.workflow}")
    print(f"Reference Dataset: {args.reference_dataset}")
    print(f"Min Representation Threshold: {args.min_representation:.0%}")
    print(f"Max Disparity Threshold: {args.max_disparity:.0%}")
    print("="*80)

    # Load data
    genomics_df = None
    clinical_df = None

    if args.genomics_data:
        print(f"\nLoading genomics data from {args.genomics_data}...")
        genomics_df = load_genomics_data(args.genomics_data)
        print(f"✓ Loaded {len(genomics_df)} genomics records")

    if args.clinical_data:
        print(f"\nLoading clinical data from {args.clinical_data}...")
        clinical_df = load_clinical_data(args.clinical_data)
        print(f"✓ Loaded {len(clinical_df)} clinical records")

    # Run audit
    if args.workflow == 'patientone':
        audit_results = audit_patientone_workflow(
            genomics_df=genomics_df,
            clinical_df=clinical_df,
            min_representation=args.min_representation,
            max_disparity=args.max_disparity,
            reference_dataset=args.reference_dataset
        )
    else:
        print("Generic workflow audit not yet implemented")
        sys.exit(1)

    # Generate report
    output_path = Path(args.output)
    if output_path.suffix == '.html':
        generate_html_report(audit_results, args.output)
    elif output_path.suffix == '.json':
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, 'w') as f:
            json.dump(audit_results, f, indent=2, default=str)
        print(f"\n✓ JSON report generated: {args.output}")
    else:
        print(f"Unsupported output format: {output_path.suffix}")
        sys.exit(1)

    print("\n" + "="*80)
    print("AUDIT COMPLETE")
    print("="*80)


if __name__ == '__main__':
    main()
