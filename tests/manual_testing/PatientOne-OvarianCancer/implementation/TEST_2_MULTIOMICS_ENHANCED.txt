TEST 2: Multi-Omics Resistance Analysis (ENHANCED WORKFLOW)
================================================================

Patient ID: PAT001-OVC-2025

⚠️ UPDATED: Now includes preprocessing pipeline and upstream regulator prediction

## Multi-Omics PDX Analysis (use mcp-multiomics - 9 tools)

For patient PAT001-OVC-2025, analyze platinum resistance in PDX models using the COMPLETE workflow:

### Data Files Location:
Files are in: `patient-data/PAT001-OVC-2025/multiomics/`
- sample_metadata.csv (includes Batch column: Batch 1 or Batch 2)
- pdx_rna_seq.csv (raw, not preprocessed)
- pdx_proteomics.csv (raw TMT data, ~18 samples/batch, BATCH EFFECTS EXPECTED)
- pdx_phosphoproteomics.csv (raw)

### Analysis Steps (Enhanced Workflow):

## STEP 0: DATA PREPROCESSING ⭐ NEW (CRITICAL for real proteomics data)

1. **Validate Data Quality:**
   ```
   Use mcp-multiomics tool: validate_multiomics_data

   Inputs:
   - rna_path: pdx_rna_seq.csv
   - protein_path: pdx_proteomics.csv
   - phospho_path: pdx_phosphoproteomics.csv
   - metadata_path: sample_metadata.csv (must have 'Batch' column)

   Check for:
   - Batch effects (expect PC1-batch correlation > 0.7 in protein data)
   - Missing value patterns
   - Sample naming consistency
   - Outlier samples
   ```

   **Expected Validation Results:**
   - ⚠️ CRITICAL: Batch effects detected in protein data (PC1-batch r=0.82)
   - ⚠️ WARNING: Missing values ~30-40% in phospho data
   - ⚠️ WARNING: Sample_07 identified as outlier
   - ✅ Sample names consistent across modalities

2. **Preprocess Data:**
   ```
   Use mcp-multiomics tool: preprocess_multiomics_data

   Apply preprocessing pipeline:
   - Batch correction: ComBat (CRITICAL for proteomics)
   - Imputation: KNN (k=5) for missing values
   - Normalization: Quantile normalization
   - Outlier removal: MAD threshold 3.0
   - Output: Save preprocessed data to /preprocessed/
   ```

   **Expected Preprocessing Results:**
   - Batch effects REDUCED: PC1-batch correlation 0.82 → 0.15 ✅
   - Missing values FILLED: ~2000 protein values imputed
   - Outliers REMOVED: Sample_07 excluded
   - Final: 14 samples (7 resistant, 7 sensitive)

3. **Visualize QC:**
   ```
   Use mcp-multiomics tool: visualize_data_quality

   Generate before/after plots:
   - PCA colored by batch (verify batch correction worked)
   - Correlation heatmaps
   - Missing value patterns

   Verify: PC1-batch correlation < 0.3 after preprocessing
   ```

## STEP 1: DATA INTEGRATION

4. **Integrate Preprocessed Data:**
   ```
   Use mcp-multiomics tool: integrate_omics_data

   Load PREPROCESSED files (not raw):
   - RNA: /preprocessed/pdx_rna_seq_preprocessed.csv
   - Protein: /preprocessed/pdx_proteomics_preprocessed.csv
   - Phospho: /preprocessed/pdx_phosphoproteomics_preprocessed.csv
   - Metadata: sample_metadata.csv

   Samples: 14 (7 resistant, 7 sensitive)
   ```

## STEP 2: ASSOCIATION TESTING & META-ANALYSIS

5. **Focus on KEY RESISTANCE GENES** (to reduce context):
   - **PI3K/AKT pathway:** PIK3CA, AKT1, MTOR, PTEN
   - **Drug resistance:** ABCB1 (MDR1)
   - **Anti-apoptotic:** BCL2L1
   - **Tumor suppressor:** TP53

   Extract data for these 7 genes from all 3 modalities

6. **Run Stouffer's Meta-Analysis:**
   ```
   Use mcp-multiomics tool: calculate_stouffer_meta

   For each of the 7 genes:
   - Input: NOMINAL p-values from differential expression (RNA, Protein, Phospho)
   - Input: log2 fold changes (for directionality)
   - Method: Stouffer's Z-score method
   - FDR correction: Applied AFTER combination (α = 0.05)

   Output:
   - meta_p_values (combined, still nominal)
   - q_values (FDR-corrected - USE THESE for significance)
   - meta_z_scores (with directionality)
   ```

## STEP 3: UPSTREAM REGULATOR PREDICTION ⭐ NEW

7. **Predict Therapeutic Targets:**
   ```
   Use mcp-multiomics tool: predict_upstream_regulators

   Input: Significant genes from Stouffer's (q < 0.05)
   - PIK3CA: log2FC=2.3, p=0.0001
   - AKT1: log2FC=2.1, p=0.0003
   - MTOR: log2FC=1.9, p=0.0005
   - ABCB1: log2FC=2.5, p=0.0002
   - BCL2L1: log2FC=1.8, p=0.001
   - PTEN: log2FC=-2.1, p=0.0001 (downregulated)
   - TP53: log2FC=-1.5, p=0.002 (inhibited)

   Analyze for:
   - Kinases (activation state)
   - Transcription factors
   - Drug responses
   ```

   **Expected Upstream Regulator Results:**

   **Kinases (Activated):**
   - AKT1: Z-score=3.2, q=0.001, Targets: 5 genes activated
   - MTOR: Z-score=2.8, q=0.003, Targets: 4 genes activated
   - PI3K: Z-score=3.0, q=0.002, Targets: 6 genes activated

   **Kinases (Inhibited):**
   - GSK3B: Z-score=-2.5, q=0.005 (tumor suppression reduced)

   **Transcription Factors:**
   - TP53: Z-score=-3.5, q=0.0001 (INHIBITED - loss of tumor suppression)
   - MYC: Z-score=2.9, q=0.002 (activated - proliferation)

   **Drug Targets Identified:**
   - Alpelisib (PI3K inhibitor): Targets activated PI3K pathway
   - Capivasertib (AKT inhibitor): Targets activated AKT1
   - Everolimus (mTOR inhibitor): Targets activated MTOR
   - Clinical Indication: PI3K/AKT/mTOR pathway activation

## Expected Results:

**Sample Summary (After Preprocessing):**
- Total: 14 PDX samples (outlier removed)
- Resistant: 7 samples (PDX_R001-R007)
- Sensitive: 7 samples (PDX_S001-S007)
- Batch correction: APPLIED (PC1-batch 0.82 → 0.15)

**Gene-Level Results (Stouffer's Meta-Analysis):**
| Gene   | RNA FC | Prot FC | Phos FC | Z-score | q-value | Direction |
|--------|--------|---------|---------|---------|---------|-----------|
| PIK3CA | +2.3   | +2.0    | +1.8    | 4.2     | 0.0001  | UP ↑      |
| AKT1   | +2.1   | +1.9    | +2.3    | 4.5     | <0.0001 | UP ↑      |
| MTOR   | +1.9   | +1.7    | +1.5    | 3.8     | 0.0003  | UP ↑      |
| ABCB1  | +2.5   | +2.2    | +1.9    | 4.1     | 0.0001  | UP ↑      |
| BCL2L1 | +1.8   | +1.6    | +1.4    | 3.2     | 0.002   | UP ↑      |
| PTEN   | -2.1   | -1.9    | -1.7    | -3.9    | 0.0002  | DOWN ↓    |
| TP53   | -1.5   | -1.3    | -1.1    | -2.8    | 0.005   | DOWN ↓    |

**Upstream Regulators (NEW):**
- **Activated Kinases:** AKT1, MTOR, PI3K (Z > 2.5, q < 0.005)
- **Inhibited TFs:** TP53 (Z = -3.5, loss of tumor suppression)
- **Drug Recommendations:** PI3K inhibitors (alpelisib), AKT inhibitors (capivasertib), mTOR inhibitors (everolimus)

**Pathway Analysis:**
- ✅ PI3K/AKT/mTOR pathway: **ACTIVATED** in resistant samples
- ✅ Evidence: PIK3CA, AKT1, MTOR upregulated; PTEN downregulated
- ✅ Mechanism: Loss of PTEN → PI3K hyperactivation → AKT/mTOR signaling
- ✅ Clinical significance: Platinum resistance via survival signaling
- ✅ Therapeutic strategy: PI3K/AKT/mTOR inhibitor combinations

## Output Format:

Please provide:

1. **Preprocessing Summary:**
   - Validation results (batch effects detected)
   - Batch correction applied (before/after PC1-batch correlation)
   - Imputation stats (values filled)
   - QC verification (plots generated)

2. **Sample Summary:**
   - Total samples: 14 (post-QC)
   - Resistant: 7, Sensitive: 7
   - Confirm data loaded from all 3 PREPROCESSED modalities

3. **Gene-Level Results:**
   For each of the 7 genes:
   - Expression in resistant vs sensitive (RNA/Protein/Phospho)
   - Stouffer's combined Z-score
   - q-value (FDR-corrected)
   - Direction (up/down in resistant)

4. **Upstream Regulator Predictions (NEW):**
   - Top activated kinases with Z-scores
   - Inhibited tumor suppressors
   - Drug recommendations with mechanisms
   - Clinical trial suggestions

5. **Pathway Interpretation:**
   - Is PI3K/AKT pathway activated? (Yes)
   - Clinical significance for platinum resistance
   - Therapeutic targets with evidence levels
   - Monitoring strategy

## Validation Checkpoints:

**Preprocessing:**
✅ Validation: Batch effects detected (PC1-batch r=0.82)
✅ Batch correction: PC1-batch reduced to r=0.15
✅ Imputation: ~2000 missing values filled
✅ QC plots: Generated and verified

**Analysis:**
✅ Data loaded: 14 samples from 3 PREPROCESSED modalities
✅ Resistant genes: PIK3CA, AKT1, MTOR, ABCB1, BCL2L1 upregulated
✅ Tumor suppressors: PTEN, TP53 downregulated
✅ Stouffer's Z-scores: >3 for top genes
✅ FDR correction: Applied AFTER combination (α=0.05)
✅ Pathway: PI3K/AKT/mTOR activated

**Upstream Regulators (NEW):**
✅ Kinases identified: AKT1, MTOR, PI3K activated
✅ TFs identified: TP53 inhibited
✅ Drug targets: Alpelisib, Capivasertib, Everolimus recommended
✅ Activation Z-scores: Calculated with directionality

## Technical Notes:

**Why Preprocessing Matters:**
- Real proteomics data has batch effects (~18 samples/MS run)
- Without ComBat correction, PC1 = batch artifact, NOT biology
- Missing values are systematic (different proteins per batch)
- KNN imputation preserves biological structure

**Why Correct FDR Timing Matters:**
- HAllA returns NOMINAL p-values (not FDR-corrected)
- Stouffer's combines NOMINAL p-values across modalities
- FDR applied AFTER combination (more statistical power)
- Using pre-corrected q-values would be WRONG (over-conservative)

**Why Upstream Regulators Matter:**
- Identifies druggable targets (kinases)
- Predicts pathway activation state
- Provides clinical trial recommendations
- IPA-like analysis without expensive software
