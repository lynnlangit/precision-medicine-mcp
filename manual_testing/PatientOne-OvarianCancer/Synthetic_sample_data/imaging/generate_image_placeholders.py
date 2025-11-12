#!/usr/bin/env python3
"""
Generate placeholder TIFF imaging files for ovarian cancer patient.
Creates small synthetic images with appropriate metadata for testing.
"""

import numpy as np
from PIL import Image
from pathlib import Path

# Set random seed
np.random.seed(42)

def create_he_image(output_path, size=(512, 512)):
    """Create H&E histology placeholder with pink/purple staining pattern."""
    # Create RGB image with H&E-like appearance
    img = np.random.randint(180, 255, size=(*size, 3), dtype=np.uint8)

    # Add some structure (nuclei-like dark spots)
    num_nuclei = 200
    for _ in range(num_nuclei):
        x = np.random.randint(0, size[0])
        y = np.random.randint(0, size[1])
        radius = np.random.randint(3, 8)

        # Create circular nucleus (purple/blue)
        for i in range(max(0, x-radius), min(size[0], x+radius)):
            for j in range(max(0, y-radius), min(size[1], y+radius)):
                if (i-x)**2 + (j-y)**2 < radius**2:
                    img[i, j] = [80, 60, 120]  # Purple-ish

    # Add cytoplasm regions (pink)
    img[:, :, 0] = np.clip(img[:, :, 0] * 1.1, 0, 255)  # More red
    img[:, :, 1] = np.clip(img[:, :, 1] * 0.9, 0, 255)  # Less green

    pil_img = Image.fromarray(img, mode='RGB')
    pil_img.save(output_path, format='TIFF', compression='tiff_lzw')
    return output_path

def create_if_image(output_path, size=(512, 512), marker_name='DAPI'):
    """Create immunofluorescence placeholder."""
    # Create grayscale image
    if marker_name == 'DAPI':
        # Nuclear stain - discrete bright spots
        img = np.random.randint(10, 30, size=size, dtype=np.uint8)
        num_nuclei = 200
        for _ in range(num_nuclei):
            x = np.random.randint(0, size[0])
            y = np.random.randint(0, size[1])
            radius = np.random.randint(4, 7)
            for i in range(max(0, x-radius), min(size[0], x+radius)):
                for j in range(max(0, y-radius), min(size[1], y+radius)):
                    if (i-x)**2 + (j-y)**2 < radius**2:
                        img[i, j] = np.random.randint(180, 255)
    else:
        # Protein marker - more diffuse
        img = np.random.randint(20, 80, size=size, dtype=np.uint8)
        # Add some bright regions (positive cells)
        num_positive = 100
        for _ in range(num_positive):
            x = np.random.randint(0, size[0])
            y = np.random.randint(0, size[1])
            radius = np.random.randint(5, 12)
            for i in range(max(0, x-radius), min(size[0], x+radius)):
                for j in range(max(0, y-radius), min(size[1], y+radius)):
                    dist = np.sqrt((i-x)**2 + (j-y)**2)
                    if dist < radius:
                        intensity = int(200 - (dist/radius) * 100)
                        img[i, j] = max(img[i, j], intensity)

    pil_img = Image.fromarray(img, mode='L')
    pil_img.save(output_path, format='TIFF', compression='tiff_lzw')
    return output_path

def create_multiplex_if_image(output_path, size=(512, 512)):
    """Create multiplex IF placeholder with multiple channels."""
    # Create RGB image representing 3-channel multiplex IF
    # Channel 1 (Red): TP53
    # Channel 2 (Green): Ki67
    # Channel 3 (Blue): DAPI

    img = np.zeros((*size, 3), dtype=np.uint8)

    # DAPI channel (blue) - nuclei
    num_nuclei = 200
    for _ in range(num_nuclei):
        x = np.random.randint(0, size[0])
        y = np.random.randint(0, size[1])
        radius = np.random.randint(4, 7)
        for i in range(max(0, x-radius), min(size[0], x+radius)):
            for j in range(max(0, y-radius), min(size[1], y+radius)):
                if (i-x)**2 + (j-y)**2 < radius**2:
                    img[i, j, 2] = np.random.randint(150, 255)

    # Ki67 channel (green) - proliferation marker
    num_ki67 = 80
    for _ in range(num_ki67):
        x = np.random.randint(0, size[0])
        y = np.random.randint(0, size[1])
        radius = np.random.randint(5, 9)
        for i in range(max(0, x-radius), min(size[0], x+radius)):
            for j in range(max(0, y-radius), min(size[1], y+radius)):
                if (i-x)**2 + (j-y)**2 < radius**2:
                    img[i, j, 1] = np.random.randint(100, 200)

    # TP53 channel (red) - nuclear protein
    num_tp53 = 150
    for _ in range(num_tp53):
        x = np.random.randint(0, size[0])
        y = np.random.randint(0, size[1])
        radius = np.random.randint(4, 7)
        for i in range(max(0, x-radius), min(size[0], x+radius)):
            for j in range(max(0, y-radius), min(size[1], y+radius)):
                if (i-x)**2 + (j-y)**2 < radius**2:
                    img[i, j, 0] = np.random.randint(120, 220)

    pil_img = Image.fromarray(img, mode='RGB')
    pil_img.save(output_path, format='TIFF', compression='tiff_lzw')
    return output_path

# Generate all placeholder images
output_dir = Path(__file__).parent

print("Generating placeholder imaging files...\n")

# 1. H&E histology
he_path = output_dir / 'PAT001_tumor_HE_20x.tiff'
create_he_image(he_path)
print(f"✓ Created: {he_path.name} (H&E histology)")

# 2. Immunofluorescence images
if_markers = ['DAPI', 'CD8', 'CD3', 'KI67', 'PanCK']
for marker in if_markers:
    if_path = output_dir / f'PAT001_tumor_IF_{marker}.tiff'
    create_if_image(if_path, marker_name=marker)
    print(f"✓ Created: {if_path.name} (IF: {marker})")

# 3. Multiplex IF
multiplex_path = output_dir / 'PAT001_tumor_multiplex_IF_TP53_KI67_DAPI.tiff'
create_multiplex_if_image(multiplex_path)
print(f"✓ Created: {multiplex_path.name} (Multiplex IF)")

print(f"\nTotal files created: {len(if_markers) + 2}")
print("All placeholder imaging files generated successfully!")
