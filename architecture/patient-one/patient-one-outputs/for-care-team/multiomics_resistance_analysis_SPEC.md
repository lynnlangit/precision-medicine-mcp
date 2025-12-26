# Multi-Omics Resistance Analysis Visualization Specification

**File:** multiomics_resistance_analysis.png
**Version:** 2.0 (Enhanced with Preprocessing & Upstream Regulators)
**Date:** December 26, 2025
**Purpose:** Clinical decision support visualization for care team
**Format:** PNG, 300 DPI, 16" × 12" (4800 × 3600 pixels)

---

## Overview

This multi-panel figure visualizes the complete enhanced multi-omics workflow including:
- ⭐ NEW: Preprocessing QC (batch correction before/after)
- ENHANCED: Gene-level results with preprocessing context
- ⭐ NEW: Upstream regulator predictions (therapeutic targets)
- ENHANCED: Pathway summary with drug recommendations

---

## Layout Design

```
┌─────────────────────────────────────────────────────────────────┐
│                                                                   │
│  Multi-Omics Resistance Analysis: PAT001-OVC-2025                │
│  Stage IV HGSOC, Platinum-Resistant                              │
│  PDX Models (n=13 after QC)                                      │
│                                                                   │
└─────────────────────────────────────────────────────────────────┘

┌──────────────────────────┬──────────────────────────────────────┐
│                          │                                        │
│   Panel A: QC & Batch    │   Panel B: Gene-Level Results         │
│   Correction (NEW)       │   (Stouffer's Meta-Analysis)          │
│                          │                                        │
│   [PCA Before/After]     │   [Results Table + Heatmap]           │
│                          │                                        │
├──────────────────────────┴──────────────────────────────────────┤
│                                                                   │
│   Panel C: Upstream Regulator Predictions (NEW)                  │
│   [Kinases/TFs Network + Drug Targets]                           │
│                                                                   │
├───────────────────────────────────────────────────────────────── │
│                                                                   │
│   Panel D: Pathway Summary & Therapeutic Strategy                │
│   [PI3K/AKT/mTOR Pathway Diagram with Drug Targets]              │
│                                                                   │
└───────────────────────────────────────────────────────────────── ┘
```

**Panel Dimensions:**
- Panel A: 7" × 4" (top left)
- Panel B: 9" × 4" (top right)
- Panel C: 16" × 3.5" (middle)
- Panel D: 16" × 4" (bottom)

---

## Panel A: Data Quality & Batch Correction ⭐ NEW

**Purpose:** Show preprocessing QC and batch correction effectiveness

### Left Side: PCA Before Preprocessing

**Title:** "Before Batch Correction"

**Visualization:**
- 2D PCA plot (PC1 vs PC2)
- Points: 15 samples (includes outliers)
- Colors: By batch
  - Batch 1: Blue circles
  - Batch 2: Red circles
- Shape: By treatment response
  - Resistant: Filled circles
  - Sensitive: Open circles

**Annotations:**
- PC1 variance: 67%
- PC2 variance: 15%
- **PC1-Batch correlation: r = 0.82** (in red box)
- ⚠️ Warning label: "PC1 driven by batch (technical artifact)"

**Key Finding:**
Samples cluster by BATCH, not by treatment response → batch effects dominate

### Right Side: PCA After Preprocessing

**Title:** "After Batch Correction"

**Visualization:**
- 2D PCA plot (PC1 vs PC2)
- Points: 13 samples (outliers removed)
- Colors: By batch
  - Batch 1: Light blue circles
  - Batch 2: Light red circles
- Shape: By treatment response
  - Resistant: Filled circles (dark green)
  - Sensitive: Open circles (purple)

**Annotations:**
- PC1 variance: 42%
- PC2 variance: 23%
- **PC1-Batch correlation: r = 0.12** (in green box)
- ✅ Success label: "PC1 now biological (resistant vs sensitive)"

**Key Finding:**
Samples NOW cluster by TREATMENT RESPONSE, not batch → biological signal revealed

### Center: Arrow with Preprocessing Steps

**Arrow from left → right with text:**
```
ComBat Batch Correction
↓
PC1-batch: 0.82 → 0.12 (85% reduction)
KNN Imputation: 2000 protein values
Outliers removed: 2 samples
```

### Bottom Caption

"**Quality Control:** Batch correction successfully removed technical variation.
PC1 now reflects biological differences (resistant vs sensitive), not batch artifacts.
This preprocessing step is CRITICAL for valid multi-omics analysis."

---

## Panel B: Gene-Level Results (Stouffer's Meta-Analysis)

**Purpose:** Show resistance genes across all 3 modalities

### Top Section: Results Table

**Title:** "Resistance Genes Across RNA, Protein, and Phosphorylation"

**Table Columns:**
| Gene | RNA FC | Prot FC | Phos FC | Z-score | q-value | Direction | Modalities |
|------|--------|---------|---------|---------|---------|-----------|------------|
| AKT1 | +2.1 | +1.9 | +2.3 | 4.5 | <0.0001 | ↑ UP | ●●● |
| PIK3CA | +2.3 | +2.0 | +1.8 | 4.2 | 0.0001 | ↑ UP | ●●● |
| ABCB1 | +2.5 | +2.2 | +1.9 | 4.1 | 0.0001 | ↑ UP | ●●● |
| PTEN | -2.1 | -1.9 | -1.7 | -3.9 | 0.0002 | ↓ DOWN | ●●● |
| MTOR | +1.9 | +1.7 | +1.5 | 3.8 | 0.0003 | ↑ UP | ●●● |
| BCL2L1 | +1.8 | +1.6 | +1.4 | 3.2 | 0.002 | ↑ UP | ●●● |
| TP53 | -1.5 | -1.3 | -1.1 | -2.8 | 0.005 | ↓ DOWN | ●●● |

**Column Formatting:**
- FC columns: Color scale (red positive, blue negative)
- Z-score: Bold font, size by magnitude
- q-value: Scientific notation, bold if < 0.001
- Direction: Large arrows (↑ red, ↓ blue)
- Modalities: Filled circles (●) for each modality supporting

**Note at bottom:**
"*Results from 13 PDX samples after quality control and batch correction.
Stouffer's meta-analysis combines evidence across RNA, protein, and phosphorylation.
All genes show concordant dysregulation across all 3 modalities (●●●).*"

### Bottom Section: Expression Heatmap

**Heatmap showing 7 genes × 13 samples × 3 modalities**

**Rows:** 7 resistance genes (grouped by direction)
- Upregulated: AKT1, PIK3CA, ABCB1, MTOR, BCL2L1
- Downregulated: PTEN, TP53

**Columns:** 13 samples (7 resistant, 6 sensitive)
- Annotation bar: Resistant (red) vs Sensitive (blue)

**Color Scale:**
- Red (+3): Highly upregulated
- White (0): No change
- Blue (-3): Highly downregulated

**Side Panel (3 mini-heatmaps):**
- RNA (left)
- Protein (center)
- Phospho (right)

**Key Finding Visual:**
- Clear separation: Resistant samples (left) vs Sensitive (right)
- Concordance: All 3 modalities show same pattern

---

## Panel C: Upstream Regulator Predictions ⭐ NEW

**Purpose:** Identify therapeutic targets and drug recommendations

### Left Section: Kinase Activation (40% width)

**Title:** "Activated Kinases (Therapeutic Targets)"

**Visualization:** Horizontal bar chart

**Bars (sorted by Z-score):**
1. **AKT1** | Z = 3.2 | q = 0.001 | ████████████████ (dark green)
   - Targets: GSK3B, FOXO1, MDM2, TSC2, mTOR (5 genes)
   - Drug: Capivasertib (AKT inhibitor)

2. **PI3K** | Z = 3.0 | q = 0.002 | ███████████████ (green)
   - Targets: AKT1, PDK1, PIK3R1, PTEN, mTOR (6 genes)
   - Drug: Alpelisib (PI3K inhibitor)

3. **MTOR** | Z = 2.8 | q = 0.003 | ██████████████ (light green)
   - Targets: RPS6KB1, EIF4EBP1, ULK1, TFEB (4 genes)
   - Drug: Everolimus (mTOR inhibitor)

**X-axis:** Activation Z-score (0 to 3.5)
**Y-axis:** Kinase names

**Annotations:**
- Green bars = ACTIVATED (therapeutic targets)
- Number of targets shown in parentheses
- FDA-approved drugs listed

### Center Section: Transcription Factors (30% width)

**Title:** "Transcription Factors"

**Visualization:** Diverging bar chart (positive right, negative left)

**Bars:**
1. **MYC** | Z = +2.9 | q = 0.002 | ██████████████ (right, red)
   - State: ACTIVATED
   - Role: Proliferation, metabolism

2. **TP53** | Z = -3.5 | q = 0.0001 | ████████████████ (left, blue)
   - State: INHIBITED
   - Role: Tumor suppression LOST

**Annotations:**
- Red (right) = Activated TFs
- Blue (left) = Inhibited TFs
- ⚠️ Highlight on TP53: "Loss of tumor suppression"

### Right Section: Drug Recommendations (30% width)

**Title:** "FDA-Approved Drug Targets"

**Visualization:** Drug cards with icons

**Card 1: Alpelisib (Piqray®)**
```
┌─────────────────────────┐
│ 📍 Alpelisib             │
│                          │
│ Target: PI3K alpha       │
│ Mechanism: Selective     │
│   PI3K inhibitor         │
│                          │
│ FDA Status: ✅ Approved  │
│ (Breast cancer, 2019)    │
│                          │
│ Evidence: PI3K pathway   │
│ activated (Z=3.0)        │
└─────────────────────────┘
```

**Card 2: Capivasertib (AZD5363)**
```
┌─────────────────────────┐
│ 📍 Capivasertib          │
│                          │
│ Target: AKT (pan-AKT)    │
│ Mechanism: ATP-          │
│   competitive inhibitor  │
│                          │
│ Status: ⏳ Phase III     │
│ Clinical Trials          │
│                          │
│ Evidence: AKT1 activated │
│ (Z=3.2)                  │
└─────────────────────────┘
```

**Card 3: Everolimus (Afinitor®)**
```
┌─────────────────────────┐
│ 📍 Everolimus            │
│                          │
│ Target: mTOR             │
│ Mechanism: mTORC1        │
│   inhibitor              │
│                          │
│ FDA Status: ✅ Approved  │
│ (RCC, breast, NET)       │
│                          │
│ Evidence: mTOR activated │
│ (Z=2.8)                  │
└─────────────────────────┘
```

### Bottom: Clinical Trial Recommendation

**Box spanning full width:**
```
┌─────────────────────────────────────────────────────────────────┐
│ 🏥 Recommended Clinical Trial: NCT03602859                      │
│                                                                   │
│ "Alpelisib + Capivasertib in PTEN-deficient Solid Tumors"        │
│ Phase II | Eligibility: Platinum-resistant HGSOC with PI3K       │
│ pathway activation | Primary endpoint: Objective response rate   │
│                                                                   │
│ Rationale: Patient shows PI3K (Z=3.0) + AKT1 (Z=3.2) activation, │
│ PTEN loss, consistent with trial eligibility criteria.           │
└─────────────────────────────────────────────────────────────────┘
```

---

## Panel D: Pathway Summary & Therapeutic Strategy

**Purpose:** Mechanistic explanation and treatment strategy

### Pathway Diagram

**Title:** "PI3K/AKT/mTOR Pathway Activation in Platinum Resistance"

**Visualization:** Signaling cascade with drug targets

```
                Growth Factors
                      ↓
┌──────────────────────────────────────────────────────┐
│                 Cell Membrane                         │
└──────────────────────────────────────────────────────┘
                      ↓
                 Receptor
                      ↓
               ┌─────┴─────┐
               │  PI3K ↑   │ ← 📍 Alpelisib
               │ (Z=3.0)   │
               └─────┬─────┘
                      ↓
                   PIP3
                      ↓
               ┌─────┴─────┐
       PTEN ←─│  PDK1      │
       (LOST) │             │
       ↓ ↓    └─────┬─────┘
               ┌─────┴─────┐
               │  AKT1 ↑   │ ← 📍 Capivasertib
               │ (Z=3.2)   │
               └─────┬─────┘
                      ↓
          ┌───────────┼───────────┐
          ↓           ↓           ↓
       GSK3B       FOXO1        MDM2
     (inactive)  (inactive)       ↓
                                 TP53 ↓
                              (INHIBITED)
                                       ↓
                          Loss of Tumor Suppression

                      ↓
               ┌─────┴─────┐
               │  MTOR ↑   │ ← 📍 Everolimus
               │ (Z=2.8)   │
               └─────┬─────┘
                      ↓
          ┌───────────┼───────────┐
          ↓           ↓           ↓
       Protein    Autophagy    Lipid
     Synthesis   (inhibited) Synthesis
          ↓
    ┌─────┴─────┐
    │  ABCB1 ↑  │ Drug Efflux Pump
    │ (Z=4.1)   │ (Platinum Resistance)
    └───────────┘
```

**Color Coding:**
- Green boxes with ↑ = Activated (upregulated)
- Red boxes with ↓ = Inhibited/Lost (downregulated)
- Blue 📍 = Drug target sites
- Dashed lines = Inhibitory interactions
- Solid arrows = Activating interactions

**Annotations:**
- Z-scores shown for each activated node
- Drug names at target sites
- PTEN loss highlighted in red box

### Bottom Text Boxes

**Left Box: Mechanism**
```
Resistance Mechanism:
1. PTEN loss removes brake on PI3K
2. PI3K hyperactivation → PIP3 accumulation
3. AKT1 activation → multiple downstream effects:
   - GSK3B inhibition (anti-apoptotic)
   - FOXO1 inhibition (stress response blocked)
   - MDM2 activation → TP53 degradation
4. mTOR activation → protein synthesis, ABCB1 upregulation
5. ABCB1 (MDR1) pumps platinum drugs out of cells

Result: Multi-layered platinum resistance
```

**Right Box: Therapeutic Strategy**
```
Recommended Treatment Strategy:

Primary: Dual PI3K/AKT Inhibition
- Alpelisib 300 mg PO daily
- Capivasertib 400 mg PO BID (4 days on, 3 off)

Rationale:
- Single-agent PI3K inhibitor → compensatory AKT activation
- Dual inhibition prevents pathway escape
- Synergistic effects in PTEN-deficient models

Monitoring:
- Baseline: Phospho-AKT (S473), Phospho-S6 (S235/236)
- Week 4: Tumor markers (CA-125), imaging
- Week 8: Response assessment (RECIST)

Toxicity Management:
- Hyperglycemia (main toxicity): Monitor glucose closely
- Diarrhea: Loperamide as needed
- Consider metformin for glucose control
```

---

## Color Palette

**Main Colors:**
- Upregulated genes: #D32F2F (red)
- Downregulated genes: #1976D2 (blue)
- Activated pathways: #388E3C (green)
- Inhibited pathways: #F57C00 (orange)
- Batch 1: #64B5F6 (light blue)
- Batch 2: #EF5350 (light red)
- Resistant samples: #C62828 (dark red)
- Sensitive samples: #7B1FA2 (purple)
- Drug targets: #FFB300 (gold)
- Warnings: #FF5722 (deep orange)
- Success: #4CAF50 (green)

**Background:**
- Panel backgrounds: #FAFAFA (light gray)
- Figure background: white

---

## Typography

**Fonts:**
- Title: Arial Bold, 24pt
- Panel titles: Arial Bold, 18pt
- Section headers: Arial Bold, 14pt
- Body text: Arial, 11pt
- Table text: Arial, 10pt
- Annotations: Arial, 9pt
- Gene names: Arial Bold Italic, 12pt

**Special Formatting:**
- Z-scores: Bold
- q-values < 0.001: Bold + red
- Drug names: Bold + gold background
- Warnings: Bold + orange background
- Success markers: Bold + green background

---

## Technical Implementation

### Python Code Skeleton (matplotlib/seaborn)

```python
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

# Figure setup
fig = plt.figure(figsize=(16, 12), dpi=300)
gs = fig.add_gridspec(3, 2, height_ratios=[1.2, 1, 1.3], hspace=0.3, wspace=0.2)

# Panel A: QC & Batch Correction
ax_pca_before = fig.add_subplot(gs[0, 0])
ax_pca_after = fig.add_subplot(gs[0, 0])
# ... PCA plotting code

# Panel B: Gene Results
ax_table = fig.add_subplot(gs[0, 1])
ax_heatmap = fig.add_subplot(gs[0, 1])
# ... table and heatmap code

# Panel C: Upstream Regulators
ax_kinases = fig.add_subplot(gs[1, :])
# ... bar chart code

# Panel D: Pathway
ax_pathway = fig.add_subplot(gs[2, :])
# ... pathway diagram code (networkx or manual)

plt.savefig('multiomics_resistance_analysis.png', dpi=300, bbox_inches='tight')
```

---

## Data Sources

All data comes from TEST_2_MULTIOMICS_ENHANCED.txt execution:

1. **Panel A Data:**
   - Source: `validate_multiomics_data` + `preprocess_multiomics_data` + `visualize_data_quality`
   - PC1-batch before: 0.82
   - PC1-batch after: 0.12

2. **Panel B Data:**
   - Source: `calculate_stouffer_meta`
   - 7 resistance genes with Z-scores and q-values

3. **Panel C Data:**
   - Source: `predict_upstream_regulators`
   - Kinases: AKT1, PI3K, MTOR
   - TFs: TP53, MYC
   - Drugs: Alpelisib, Capivasertib, Everolimus

4. **Panel D Data:**
   - Derived from pathway analysis of differential genes
   - Drug target annotations from upstream regulators

---

## Quality Checklist

Before finalizing, verify:

- [ ] All 13 samples shown in PCA plots
- [ ] PC1-batch correlation: 0.82 → 0.12 clearly labeled
- [ ] All 7 resistance genes in results table
- [ ] Z-scores and q-values match TEST_2 expected results
- [ ] All 3 kinases (AKT1, PI3K, MTOR) shown with correct Z-scores
- [ ] TP53 inhibition clearly marked (negative Z-score)
- [ ] All 3 drugs (Alpelisib, Capivasertib, Everolimus) shown
- [ ] Clinical trial NCT03602859 referenced
- [ ] Pathway diagram shows PTEN loss and drug targets
- [ ] Color-blind friendly palette used
- [ ] 300 DPI resolution for print quality
- [ ] All text readable at 100% zoom
- [ ] Legend included for all symbols/colors

---

## Clinical Use Case

**Target Audience:** Medical oncologists, clinical researchers
**Use In:** Tumor board presentations, clinical trial enrollment decisions
**Key Messages:**
1. Preprocessing enabled valid multi-omics analysis
2. PI3K/AKT/mTOR pathway is activated (therapeutic target)
3. Drug recommendations: Alpelisib + Capivasertib
4. Clinical trial: NCT03602859 matches patient profile

---

**Specification Status:** ✅ Complete and ready for implementation
**Next Step:** Generate visualization using Python (matplotlib/seaborn) or R (ggplot2)
**Validation:** Compare generated figure against this specification before clinical use
