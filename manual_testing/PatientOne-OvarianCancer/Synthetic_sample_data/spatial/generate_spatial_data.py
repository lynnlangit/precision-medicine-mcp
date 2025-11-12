#!/usr/bin/env python3
"""
Generate synthetic Visium spatial transcriptomics data for Stage IV HGSOC patient.
Simulates tumor heterogeneity with distinct spatial regions:
- Tumor core (high proliferation, TP53/PIK3CA signature)
- Tumor-stroma interface (EMT markers, CAFs)
- Immune infiltrated regions (T cells, macrophages)
- Necrotic regions (low viability, hypoxia markers)
"""

import numpy as np
import pandas as pd
from pathlib import Path

# Set random seed for reproducibility
np.random.seed(42)

# Define Visium array dimensions (standard Visium slide has hexagonal grid)
# We'll create a 30x30 grid for ~900 spots
n_rows, n_cols = 30, 30

# Key genes for ovarian cancer spatial analysis
genes = {
    # Tumor proliferation markers
    'MKI67': ('proliferation', 'high_in_tumor'),
    'PCNA': ('proliferation', 'high_in_tumor'),
    'TOP2A': ('proliferation', 'high_in_tumor'),

    # Cancer driver genes
    'TP53': ('driver', 'high_in_tumor'),
    'PIK3CA': ('driver', 'high_in_tumor'),
    'MYC': ('driver', 'high_in_tumor'),

    # Epithelial markers
    'EPCAM': ('epithelial', 'high_in_tumor'),
    'KRT8': ('epithelial', 'high_in_tumor'),
    'KRT18': ('epithelial', 'high_in_tumor'),

    # EMT markers (tumor-stroma interface)
    'VIM': ('EMT', 'high_at_interface'),
    'SNAI1': ('EMT', 'high_at_interface'),
    'TWIST1': ('EMT', 'high_at_interface'),
    'CDH2': ('EMT', 'high_at_interface'),

    # Stromal markers
    'COL1A1': ('stroma', 'high_in_stroma'),
    'COL3A1': ('stroma', 'high_in_stroma'),
    'ACTA2': ('CAF', 'high_in_stroma'),
    'FAP': ('CAF', 'high_in_stroma'),

    # Immune markers
    'CD3D': ('T_cell', 'high_in_immune'),
    'CD3E': ('T_cell', 'high_in_immune'),
    'CD8A': ('cytotoxic_T', 'high_in_immune'),
    'CD4': ('helper_T', 'high_in_immune'),
    'FOXP3': ('Treg', 'high_in_immune'),
    'CD68': ('macrophage', 'high_in_immune'),
    'CD163': ('M2_macrophage', 'high_in_immune'),

    # Hypoxia/necrosis markers
    'HIF1A': ('hypoxia', 'high_in_necrotic'),
    'CA9': ('hypoxia', 'high_in_necrotic'),
    'VEGFA': ('hypoxia', 'high_in_necrotic'),

    # Resistance pathway genes
    'AKT1': ('PI3K_pathway', 'high_in_tumor'),
    'MTOR': ('PI3K_pathway', 'high_in_tumor'),
    'ABCB1': ('drug_resistance', 'high_in_tumor'),
    'BCL2L1': ('anti_apoptotic', 'high_in_tumor'),
}

# Create spatial coordinates (hexagonal grid pattern)
coordinates = []
spot_ids = []

for row in range(n_rows):
    for col in range(n_cols):
        # Hexagonal spacing: offset every other row
        x = col * 100 + (50 if row % 2 == 1 else 0)
        y = row * 87  # 87 = 100 * sin(60°) for hexagonal packing

        spot_id = f"SPOT_{row:02d}_{col:02d}"
        spot_ids.append(spot_id)
        coordinates.append({
            'barcode': spot_id,
            'in_tissue': 1,
            'array_row': row,
            'array_col': col,
            'pxl_row_in_fullres': int(y * 10),  # Scale for full resolution image
            'pxl_col_in_fullres': int(x * 10)
        })

# Define spatial regions (based on distance from center)
center_row, center_col = n_rows // 2, n_cols // 2

def get_region(row, col):
    """Assign spots to regions based on spatial location."""
    dist_from_center = np.sqrt((row - center_row)**2 + (col - center_col)**2)

    if dist_from_center < 5:
        return 'tumor_core'
    elif dist_from_center < 8:
        return 'tumor_proliferative'
    elif dist_from_center < 10:
        return 'tumor_interface'
    elif dist_from_center < 13:
        return 'stroma_immune'
    elif dist_from_center < 15:
        return 'stroma'
    else:
        return 'necrotic_hypoxic'

# Generate expression data
expression_data = {}

for row in range(n_rows):
    for col in range(n_cols):
        spot_id = f"SPOT_{row:02d}_{col:02d}"
        region = get_region(row, col)

        spot_expr = {}

        for gene, (category, pattern) in genes.items():
            # Base expression with noise
            base = np.random.lognormal(4, 1.5)

            # Modulate by region
            if region == 'tumor_core' and pattern == 'high_in_tumor':
                expr = base * np.random.uniform(3, 5)
            elif region == 'tumor_proliferative' and pattern == 'high_in_tumor':
                expr = base * np.random.uniform(4, 6)
            elif region == 'tumor_interface' and pattern in ['high_in_tumor', 'high_at_interface']:
                expr = base * np.random.uniform(2, 4)
            elif region == 'stroma_immune':
                if pattern == 'high_in_immune':
                    expr = base * np.random.uniform(3, 5)
                elif pattern == 'high_in_stroma':
                    expr = base * np.random.uniform(2, 3)
                else:
                    expr = base * np.random.uniform(0.3, 0.8)
            elif region == 'stroma' and pattern == 'high_in_stroma':
                expr = base * np.random.uniform(4, 6)
            elif region == 'necrotic_hypoxic' and pattern == 'high_in_necrotic':
                expr = base * np.random.uniform(3, 5)
            else:
                expr = base * np.random.uniform(0.2, 0.6)

            spot_expr[gene] = expr

        expression_data[spot_id] = spot_expr

# Create DataFrames
coords_df = pd.DataFrame(coordinates)
expr_df = pd.DataFrame(expression_data).T

# Add region annotations
region_annotations = []
for row in range(n_rows):
    for col in range(n_cols):
        spot_id = f"SPOT_{row:02d}_{col:02d}"
        region = get_region(row, col)
        region_annotations.append({
            'barcode': spot_id,
            'region': region
        })

regions_df = pd.DataFrame(region_annotations)

# Save files
output_dir = Path(__file__).parent
coords_df.to_csv(output_dir / 'visium_spatial_coordinates.csv', index=False)
expr_df.to_csv(output_dir / 'visium_gene_expression.csv')
regions_df.to_csv(output_dir / 'visium_region_annotations.csv', index=False)

print("Spatial transcriptomics data generated successfully!")
print(f"\nFiles created:")
print(f"  1. visium_spatial_coordinates.csv - Spot coordinates ({len(coords_df)} spots)")
print(f"  2. visium_gene_expression.csv - Expression matrix ({len(genes)} genes x {len(expr_df)} spots)")
print(f"  3. visium_region_annotations.csv - Region annotations")
print(f"\nRegion summary:")
for region in regions_df['region'].value_counts().items():
    print(f"  {region[0]}: {region[1]} spots")
print(f"\nExpression range: {expr_df.min().min():.2f} - {expr_df.max().max():.2f}")
