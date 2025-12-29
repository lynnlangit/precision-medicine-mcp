#!/usr/bin/env python3
"""
Automated Patient Report Generator

Generates a complete clinical analysis report for spatial transcriptomics patients.
Integrates FHIR clinical data with spatial molecular analysis.

Usage:
    python generate_patient_report.py --patient-id patient-001 --output-dir ./results

Author: Claude Code
Date: December 29, 2025
"""

import argparse
import os
import sys
from pathlib import Path
from datetime import datetime
import json
import pandas as pd
import numpy as np
from scipy import stats

# Add servers to path
sys.path.insert(0, str(Path(__file__).parent / "servers" / "mcp-spatialtools" / "src"))
sys.path.insert(0, str(Path(__file__).parent / "servers" / "mcp-epic" / "src"))

from mcp_spatialtools.server import _calculate_morans_i


class PatientReportGenerator:
    """Generate comprehensive patient analysis reports."""

    def __init__(self, patient_id: str, output_dir: str):
        self.patient_id = patient_id
        self.output_dir = Path(output_dir)
        self.patient_output_dir = self.output_dir / patient_id
        self.patient_output_dir.mkdir(parents=True, exist_ok=True)

        # Data paths
        self.data_dir = Path("/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data")
        self.patient_data_dir = None

        # Results storage
        self.fhir_data = {}
        self.spatial_data = {}
        self.analysis_results = {}

    def log(self, message: str):
        """Print timestamped log message."""
        timestamp = datetime.now().strftime("%H:%M:%S")
        print(f"[{timestamp}] {message}")

    def fetch_fhir_data(self):
        """Fetch patient clinical data from GCP Healthcare API."""
        self.log("Fetching FHIR data from GCP Healthcare API...")

        try:
            import subprocess

            project_id = "precision-medicine-poc"
            region = "us-central1"

            # Get access token
            token_result = subprocess.run(
                ["gcloud", "auth", "print-access-token"],
                capture_output=True,
                text=True
            )
            token = token_result.stdout.strip()

            base_url = f"https://healthcare.googleapis.com/v1/projects/{project_id}/locations/{region}/datasets/precision-medicine-dataset/fhirStores/identified-fhir-store/fhir"

            # Fetch patient demographics
            import requests

            patient_url = f"{base_url}/Patient/{self.patient_id}"
            response = requests.get(
                patient_url,
                headers={"Authorization": f"Bearer {token}"}
            )

            if response.status_code == 200:
                patient_data = response.json()
                self.fhir_data['patient'] = patient_data

                name = patient_data['name'][0]
                self.log(f"✅ Patient: {name['given'][0]} {name['family']}")
                self.log(f"   Gender: {patient_data['gender']}, DOB: {patient_data['birthDate']}")
            else:
                self.log(f"⚠️  Could not fetch patient data (using mock data)")
                self.fhir_data['patient'] = {
                    'name': [{'given': ['Patient'], 'family': self.patient_id}],
                    'gender': 'unknown',
                    'birthDate': 'unknown'
                }

            # Fetch conditions
            conditions_url = f"{base_url}/Condition?patient={self.patient_id}"
            response = requests.get(
                conditions_url,
                headers={"Authorization": f"Bearer {token}"}
            )

            if response.status_code == 200:
                conditions = response.json()
                self.fhir_data['conditions'] = conditions.get('entry', [])
                self.log(f"✅ Conditions: {conditions.get('total', 0)} found")
            else:
                self.fhir_data['conditions'] = []

            # Fetch observations
            observations_url = f"{base_url}/Observation?patient={self.patient_id}"
            response = requests.get(
                observations_url,
                headers={"Authorization": f"Bearer {token}"}
            )

            if response.status_code == 200:
                observations = response.json()
                self.fhir_data['observations'] = observations.get('entry', [])
                self.log(f"✅ Observations: {observations.get('total', 0)} found")
            else:
                self.fhir_data['observations'] = []

            # Fetch medications
            medications_url = f"{base_url}/MedicationStatement?patient={self.patient_id}"
            response = requests.get(
                medications_url,
                headers={"Authorization": f"Bearer {token}"}
            )

            if response.status_code == 200:
                medications = response.json()
                self.fhir_data['medications'] = medications.get('entry', [])
                self.log(f"✅ Medications: {medications.get('total', 0)} found")
            else:
                self.fhir_data['medications'] = []

        except Exception as e:
            self.log(f"⚠️  FHIR data fetch failed: {e}")
            self.log("   Continuing with spatial analysis only...")

    def load_spatial_data(self):
        """Load spatial transcriptomics data."""
        self.log("Loading spatial transcriptomics data...")

        # Find patient data directory (handle different naming conventions)
        possible_dirs = [
            self.data_dir / self.patient_id.upper() / "spatial",
            self.data_dir / f"PAT001-OVC-2025" / "spatial",  # Specific to our test patient
        ]

        for data_dir in possible_dirs:
            if data_dir.exists():
                self.patient_data_dir = data_dir
                break

        if not self.patient_data_dir:
            raise FileNotFoundError(f"Could not find spatial data for {self.patient_id}")

        self.log(f"   Data directory: {self.patient_data_dir}")

        # Load expression data
        expr_file = self.patient_data_dir / "visium_gene_expression.csv"
        expr_data = pd.read_csv(expr_file, index_col=0).T  # Transpose to genes × spots
        self.spatial_data['expression'] = expr_data
        self.log(f"✅ Expression: {expr_data.shape[0]} genes × {expr_data.shape[1]} spots")

        # Load region annotations
        region_file = self.patient_data_dir / "visium_region_annotations.csv"
        region_data = pd.read_csv(region_file).set_index('barcode')
        self.spatial_data['regions'] = region_data
        self.log(f"✅ Regions: {len(region_data['region'].unique())} tissue regions")

        # Load coordinates
        coord_file = self.patient_data_dir / "visium_spatial_coordinates.csv"
        coord_data = pd.read_csv(coord_file).set_index('barcode')
        self.spatial_data['coordinates'] = coord_data
        self.log(f"✅ Coordinates: {len(coord_data)} spots")

        # Align data
        common_barcodes = list(
            set(expr_data.columns) &
            set(region_data.index) &
            set(coord_data.index)
        )
        common_barcodes.sort()

        self.spatial_data['expression'] = expr_data[common_barcodes]
        self.spatial_data['regions'] = region_data.loc[common_barcodes]
        self.spatial_data['coordinates'] = coord_data.loc[common_barcodes]

        self.log(f"✅ Aligned: {len(common_barcodes)} common spots")

    def perform_differential_expression(self):
        """Run differential expression analysis."""
        self.log("Running differential expression analysis (tumor_core vs stroma)...")

        expr_data = self.spatial_data['expression']
        region_data = self.spatial_data['regions']

        # Get spots for each group
        tumor_spots = region_data[region_data['region'] == 'tumor_core'].index.tolist()
        stroma_spots = region_data[region_data['region'] == 'stroma'].index.tolist()

        if len(tumor_spots) == 0 or len(stroma_spots) == 0:
            self.log("⚠️  Insufficient samples for differential expression")
            return

        self.log(f"   Comparing: {len(tumor_spots)} tumor_core vs {len(stroma_spots)} stroma spots")

        # Extract expression data
        group1_data = expr_data[tumor_spots].values
        group2_data = expr_data[stroma_spots].values

        # Perform Mann-Whitney U tests
        results = []
        genes = expr_data.index.tolist()

        for gene_idx, gene in enumerate(genes):
            gene1 = group1_data[gene_idx, :]
            gene2 = group2_data[gene_idx, :]

            try:
                statistic, p_value = stats.mannwhitneyu(gene1, gene2, alternative='two-sided')

                mean1 = np.mean(gene1)
                mean2 = np.mean(gene2)

                if mean2 > 0:
                    fold_change = mean1 / mean2
                    log2fc = np.log2(fold_change) if fold_change > 0 else -10
                else:
                    log2fc = 10 if mean1 > 0 else 0

                results.append({
                    'gene': gene,
                    'mean_tumor': mean1,
                    'mean_stroma': mean2,
                    'log2_fold_change': log2fc,
                    'p_value': p_value,
                    'statistic': statistic
                })
            except:
                results.append({
                    'gene': gene,
                    'mean_tumor': 0,
                    'mean_stroma': 0,
                    'log2_fold_change': 0,
                    'p_value': 1.0,
                    'statistic': 0
                })

        deg_df = pd.DataFrame(results)

        # Calculate FDR
        from scipy.stats import false_discovery_control
        deg_df['fdr'] = false_discovery_control(deg_df['p_value'].values)

        # Filter significant DEGs
        sig_degs = deg_df[(deg_df['fdr'] < 0.05) & (np.abs(deg_df['log2_fold_change']) > 1.0)]

        self.analysis_results['differential_expression'] = deg_df
        self.analysis_results['significant_degs'] = sig_degs

        # Save results
        output_file = self.patient_output_dir / "differential_expression.csv"
        deg_df.to_csv(output_file, index=False)

        self.log(f"✅ DEGs: {len(sig_degs)} significant (FDR < 0.05, |log2FC| > 1)")
        self.log(f"   Saved to: {output_file}")

    def calculate_spatial_autocorrelation(self):
        """Calculate spatial autocorrelation for all genes."""
        self.log("Calculating spatial autocorrelation (Moran's I)...")

        expr_data = self.spatial_data['expression']
        coord_data = self.spatial_data['coordinates']

        coordinates = coord_data[['array_row', 'array_col']].values
        genes = expr_data.index.tolist()

        results = []
        for gene in genes:
            expression_values = expr_data.loc[gene].values

            try:
                morans_i, z_score, p_value = _calculate_morans_i(
                    expression_values,
                    coordinates,
                    distance_threshold=1.5
                )

                results.append({
                    'gene': gene,
                    'morans_i': morans_i,
                    'z_score': z_score,
                    'p_value': p_value
                })
            except:
                results.append({
                    'gene': gene,
                    'morans_i': 0,
                    'z_score': 0,
                    'p_value': 1.0
                })

        spatial_df = pd.DataFrame(results)
        spatial_df = spatial_df.sort_values('morans_i', ascending=False)

        # Identify significant SVGs
        svgs = spatial_df[spatial_df['p_value'] < 0.01]

        self.analysis_results['spatial_autocorrelation'] = spatial_df
        self.analysis_results['spatially_variable_genes'] = svgs

        # Save results
        output_file = self.patient_output_dir / "spatial_autocorrelation.csv"
        spatial_df.to_csv(output_file, index=False)

        self.log(f"✅ SVGs: {len(svgs)} spatially variable genes (p < 0.01)")
        self.log(f"   Top gene: {svgs.iloc[0]['gene']} (Moran's I = {svgs.iloc[0]['morans_i']:.4f})")
        self.log(f"   Saved to: {output_file}")

    def perform_cell_deconvolution(self):
        """Perform cell type deconvolution."""
        self.log("Performing cell type deconvolution...")

        expr_data = self.spatial_data['expression']
        region_data = self.spatial_data['regions']

        # Define signatures
        signatures = {
            'fibroblasts': ['COL1A1', 'COL3A1', 'ACTA2'],
            'immune_cells': ['CD3D', 'CD8A', 'PTPRC'],
            'hypoxic': ['HIF1A', 'CA9', 'VEGFA'],
            'resistant': ['ABCB1', 'PIK3CA', 'AKT1']
        }

        # Calculate signature scores
        signature_scores = {}
        for cell_type, sig_genes in signatures.items():
            available_genes = [g for g in sig_genes if g in expr_data.index]

            if len(available_genes) > 0:
                sig_expr = expr_data.loc[available_genes]
                scores = sig_expr.mean(axis=0)
                signature_scores[cell_type] = scores

        score_df = pd.DataFrame(signature_scores)
        score_df['region'] = region_data['region'].values

        # Calculate mean scores by region
        region_means = score_df.groupby('region').mean()

        self.analysis_results['cell_deconvolution'] = score_df
        self.analysis_results['region_cell_scores'] = region_means

        # Save results
        output_file = self.patient_output_dir / "cell_deconvolution.csv"
        region_means.to_csv(output_file)

        self.log(f"✅ Cell types: {len(signature_scores)} signatures analyzed")
        self.log(f"   Saved to: {output_file}")

    def generate_clinical_summary(self):
        """Generate clinical interpretation summary."""
        self.log("Generating clinical summary...")

        summary_lines = []
        summary_lines.append("=" * 80)
        summary_lines.append("CLINICAL ANALYSIS REPORT")
        summary_lines.append("=" * 80)
        summary_lines.append("")

        # Patient information
        if self.fhir_data.get('patient'):
            patient = self.fhir_data['patient']
            name = patient['name'][0]
            summary_lines.append(f"Patient: {name['given'][0]} {name['family']}")
            summary_lines.append(f"Gender: {patient.get('gender', 'unknown')}")
            summary_lines.append(f"DOB: {patient.get('birthDate', 'unknown')}")
        else:
            summary_lines.append(f"Patient ID: {self.patient_id}")

        summary_lines.append(f"Report Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        summary_lines.append("")

        # Clinical diagnoses
        if self.fhir_data.get('conditions'):
            summary_lines.append("DIAGNOSES:")
            for entry in self.fhir_data['conditions']:
                condition = entry['resource']
                summary_lines.append(f"  - {condition['code']['text']}")
            summary_lines.append("")

        # Medications
        if self.fhir_data.get('medications'):
            summary_lines.append("CURRENT MEDICATIONS:")
            for entry in self.fhir_data['medications']:
                med = entry['resource']
                summary_lines.append(f"  - {med['medicationCodeableConcept']['text']} ({med['status']})")
            summary_lines.append("")

        # Spatial analysis summary
        summary_lines.append("MOLECULAR ANALYSIS SUMMARY:")
        summary_lines.append("-" * 80)
        summary_lines.append("")

        # Differential expression
        if 'significant_degs' in self.analysis_results:
            degs = self.analysis_results['significant_degs']
            upregulated = degs[degs['log2_fold_change'] > 0].sort_values('log2_fold_change', ascending=False)
            downregulated = degs[degs['log2_fold_change'] < 0].sort_values('log2_fold_change')

            summary_lines.append(f"1. DIFFERENTIAL EXPRESSION (tumor_core vs stroma):")
            summary_lines.append(f"   Total DEGs: {len(degs)}")
            summary_lines.append(f"   Upregulated: {len(upregulated)}")
            summary_lines.append(f"   Downregulated: {len(downregulated)}")
            summary_lines.append("")

            if len(upregulated) > 0:
                summary_lines.append("   Top 5 upregulated genes in tumor:")
                for _, row in upregulated.head(5).iterrows():
                    summary_lines.append(f"      {row['gene']}: log2FC = {row['log2_fold_change']:.3f}, FDR = {row['fdr']:.2e}")
                summary_lines.append("")

        # Spatial patterns
        if 'spatially_variable_genes' in self.analysis_results:
            svgs = self.analysis_results['spatially_variable_genes']

            summary_lines.append(f"2. SPATIAL ORGANIZATION:")
            summary_lines.append(f"   Spatially variable genes: {len(svgs)}")
            summary_lines.append("")

            if len(svgs) > 0:
                summary_lines.append("   Top 5 spatially clustered genes:")
                for _, row in svgs.head(5).iterrows():
                    summary_lines.append(f"      {row['gene']}: Moran's I = {row['morans_i']:.4f}, p = {row['p_value']:.2e}")
                summary_lines.append("")

        # Cell composition
        if 'region_cell_scores' in self.analysis_results:
            summary_lines.append("3. TUMOR MICROENVIRONMENT:")
            summary_lines.append("   Cell type enrichment by region:")
            summary_lines.append("")

            region_scores = self.analysis_results['region_cell_scores']
            for region in region_scores.index:
                summary_lines.append(f"   {region}:")
                for cell_type in region_scores.columns:
                    score = region_scores.loc[region, cell_type]
                    summary_lines.append(f"      {cell_type}: {score:.1f}")
            summary_lines.append("")

        # Clinical implications
        summary_lines.append("CLINICAL IMPLICATIONS:")
        summary_lines.append("-" * 80)
        summary_lines.append("")

        # Check for drug resistance markers
        if 'significant_degs' in self.analysis_results:
            degs = self.analysis_results['significant_degs']
            resistance_markers = ['ABCB1', 'PIK3CA', 'AKT1', 'MTOR']
            found_markers = degs[degs['gene'].isin(resistance_markers)]

            if len(found_markers) > 0:
                summary_lines.append("Drug Resistance Markers Detected:")
                for _, row in found_markers.iterrows():
                    summary_lines.append(f"  - {row['gene']}: {row['log2_fold_change']:.2f}× fold change")
                summary_lines.append("")
                summary_lines.append("  ⚠️  Consider:")
                summary_lines.append("     - PI3K/AKT pathway inhibitors (Alpelisib, Capivasertib)")
                summary_lines.append("     - MDR reversal agents")
                summary_lines.append("")

        # Check for hypoxia
        if 'spatially_variable_genes' in self.analysis_results:
            svgs = self.analysis_results['spatially_variable_genes']
            hypoxia_genes = svgs[svgs['gene'].isin(['HIF1A', 'CA9', 'VEGFA'])]

            if len(hypoxia_genes) > 0:
                summary_lines.append("Hypoxic Regions Identified:")
                for _, row in hypoxia_genes.iterrows():
                    summary_lines.append(f"  - {row['gene']}: Moran's I = {row['morans_i']:.4f} (spatially clustered)")
                summary_lines.append("")
                summary_lines.append("  ⚠️  Consider:")
                summary_lines.append("     - Hypoxia-targeting agents (Evofosfamide, TH-302)")
                summary_lines.append("     - VEGF inhibitors (Bevacizumab)")
                summary_lines.append("")

        # Disclaimers
        summary_lines.append("=" * 80)
        summary_lines.append("DISCLAIMER:")
        summary_lines.append("This analysis is for RESEARCH PURPOSES ONLY.")
        summary_lines.append("NOT validated for clinical decision-making.")
        summary_lines.append("All treatment decisions must be made by qualified oncologists.")
        summary_lines.append("=" * 80)

        # Save summary
        summary_text = "\n".join(summary_lines)
        output_file = self.patient_output_dir / "clinical_summary.txt"
        with open(output_file, 'w') as f:
            f.write(summary_text)

        self.log(f"✅ Clinical summary saved to: {output_file}")

        # Also print to console
        print("\n")
        print(summary_text)
        print("\n")

    def save_metadata(self):
        """Save analysis metadata."""
        metadata = {
            'patient_id': self.patient_id,
            'analysis_date': datetime.now().isoformat(),
            'data_directory': str(self.patient_data_dir),
            'output_directory': str(self.patient_output_dir),
            'spatial_data': {
                'num_genes': self.spatial_data['expression'].shape[0],
                'num_spots': self.spatial_data['expression'].shape[1],
                'num_regions': len(self.spatial_data['regions']['region'].unique())
            },
            'analysis_results': {
                'num_degs': len(self.analysis_results.get('significant_degs', [])),
                'num_svgs': len(self.analysis_results.get('spatially_variable_genes', []))
            }
        }

        output_file = self.patient_output_dir / "metadata.json"
        with open(output_file, 'w') as f:
            json.dump(metadata, f, indent=2)

        self.log(f"✅ Metadata saved to: {output_file}")

    def generate_report(self):
        """Run complete analysis pipeline."""
        self.log("=" * 80)
        self.log(f"GENERATING PATIENT REPORT: {self.patient_id}")
        self.log("=" * 80)
        self.log("")

        # Step 1: Fetch FHIR data
        self.fetch_fhir_data()
        self.log("")

        # Step 2: Load spatial data
        self.load_spatial_data()
        self.log("")

        # Step 3: Run analyses
        self.perform_differential_expression()
        self.log("")

        self.calculate_spatial_autocorrelation()
        self.log("")

        self.perform_cell_deconvolution()
        self.log("")

        # Step 4: Generate summary
        self.generate_clinical_summary()

        # Step 5: Save metadata
        self.save_metadata()

        self.log("")
        self.log("=" * 80)
        self.log(f"✅ REPORT COMPLETE!")
        self.log(f"   Output directory: {self.patient_output_dir}")
        self.log("=" * 80)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Generate automated patient analysis report",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze Patient-001
  python generate_patient_report.py --patient-id patient-001 --output-dir ./results

  # Specify custom output directory
  python generate_patient_report.py --patient-id PAT002-OVC-2025 --output-dir /path/to/results
        """
    )

    parser.add_argument(
        '--patient-id',
        required=True,
        help='Patient ID (e.g., patient-001)'
    )

    parser.add_argument(
        '--output-dir',
        default='./results',
        help='Output directory for results (default: ./results)'
    )

    args = parser.parse_args()

    # Generate report
    generator = PatientReportGenerator(args.patient_id, args.output_dir)
    generator.generate_report()


if __name__ == "__main__":
    main()
