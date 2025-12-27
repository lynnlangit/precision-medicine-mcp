# Output Disclaimers & Safety Warnings

**Purpose:** Standardized disclaimers and safety warnings to prevent misuse of analysis results for clinical decision-making.

---

## Research Use Only Disclaimer

**Add to ALL analysis outputs (especially mcp-multiomics, PatientOne integration):**

```
⚠️  RESEARCH USE ONLY - NOT FOR CLINICAL DECISION-MAKING ⚠️

This analysis is provided for RESEARCH PURPOSES ONLY and has not been
clinically validated. Results should NOT be used for:
  • Patient diagnosis
  • Treatment selection
  • Clinical decision-making
  • Publication without independent validation

All findings must be:
  1. Reviewed by qualified bioinformatics personnel
  2. Validated with orthogonal methods (e.g., qPCR, Western blot, IHC)
  3. Interpreted by board-certified oncologist before any clinical consideration

The developers assume NO liability for clinical decisions based on this output.
```

---

## Uncertainty Quantification

**For statistical results (p-values, fold changes, pathway predictions):**

### Confidence Level Indicators

Use these standard confidence levels in all statistical outputs:

```
📊 Statistical Confidence:

⭐⭐⭐ HIGH CONFIDENCE
  - p < 0.001 (after FDR correction)
  - Concordant across multiple modalities (RNA + Protein + Phospho)
  - Effect size log2FC > 2.0
  - Validated in independent datasets
  → Findings likely robust, still require experimental validation

⭐⭐ MODERATE CONFIDENCE
  - p < 0.05 (after FDR correction)
  - Some modality disagreement
  - Effect size log2FC 1.0-2.0
  → Findings plausible, require careful validation

⭐ LOW CONFIDENCE
  - p < 0.1 (exploratory threshold)
  - Single modality only
  - Effect size log2FC < 1.0
  → Exploratory finding only, hypothesis generation
```

### Always Report

For every statistical result, include:
- **p-values** (both raw and FDR-corrected)
- **Effect sizes** (log2FC, correlation coefficients, not just significance)
- **Sample size** and statistical power
- **Missing data percentage** per modality
- **Batch effects** detected and correction applied (if any)
- **Concordance** across modalities (for multi-omics)

### Example Usage

```python
result = {
    "gene": "PIK3CA",
    "meta_p_value": 0.0001,
    "meta_q_value": 0.002,  # FDR-corrected
    "log2FC": 2.3,
    "confidence": "HIGH",  # ⭐⭐⭐
    "evidence": {
        "rna_p": 0.001,
        "protein_p": 0.003,
        "phospho_p": 0.002,
        "concordance": "All modalities agree on upregulation"
    },
    "sample_size": 15,
    "missing_data_pct": 5,
    "recommendation": "High-confidence finding. Validate with qPCR and Western blot before clinical consideration."
}
```

---

## Data Quality Warnings

**When data quality is suboptimal, add explicit warnings:**

### Template: Data Quality Issues

```
⚠️  DATA QUALITY WARNING

The following data quality issues were detected:

Issue 1: Missing Values
  - RNA: 12% missing (imputed using KNN, k=5)
  - Protein: 18% missing (imputed using KNN, k=5)
  - Phospho: 25% missing (imputed using KNN, k=5)
  → Impact: High missing rates may affect downstream analysis accuracy

Issue 2: Batch Effects
  - Detected: YES (PC1 correlates with batch, r=0.68, p<0.001)
  - Corrected: YES (using limma removeBatchEffect)
  - Residual correlation after correction: r=0.15, p=0.08
  → Impact: Batch correction applied, but some residual effects remain

Issue 3: Low Sequencing Depth
  - Samples with <10M reads: 3/15 (20%)
  - Median depth: 12M reads (recommended: >20M)
  → Impact: Low-depth samples may have reduced power to detect differential expression

Issue 4: Outlier Samples
  - Detected: 2 samples (Sample_07, Sample_12)
  - MAD > 3.0 from median
  - Action: Retained in analysis, flagged for review
  → Impact: Outliers may skew results, consider sensitivity analysis

RECOMMENDATION:
Interpret results with caution due to data quality limitations.
Consider re-sequencing low-depth samples and validating top findings
with higher-quality data before publication.
```

### Template: Good Quality Data

```
✅ DATA QUALITY: GOOD

Quality checks passed:
  • Missing values: <10% across all modalities
  • No significant batch effects detected (PC1 batch correlation r=0.08, p=0.52)
  • Sequencing depth: All samples >15M reads
  • No outlier samples detected (all within MAD < 2.5)

Data is suitable for downstream analysis with standard interpretation.
```

---

## Upstream Regulator Predictions

**For kinase/TF/drug target predictions (mcp-multiomics tool):**

```
⚠️  UPSTREAM REGULATOR PREDICTIONS - COMPUTATIONAL INFERENCE ⚠️

These predictions are based on:
  • Database annotations (PhosphoSitePlus, TRRUST, TRANSFAC)
  • Enrichment analysis (Fisher's exact test)
  • Activation state inference (Z-score method)

Limitations:
  • Predictions are COMPUTATIONAL, not experimental
  • Database annotations may be incomplete or context-specific
  • Activation states are inferred, not directly measured
  • Drug targets are predicted based on pathway activation, not drug efficacy

Required Validation:
  1. Kinase predictions: Validate with kinase assays or phospho-specific antibodies
  2. TF predictions: Validate with ChIP-seq or reporter assays
  3. Drug predictions: Validate in appropriate cell/animal models
  4. Clinical use: Requires clinical trial evidence

Do NOT use these predictions for:
  • Direct patient treatment decisions
  • Drug selection without clinical evidence
  • Publication without experimental validation
```

---

## Pathway Enrichment Results

**For pathway analysis outputs:**

```
⚠️  PATHWAY ENRICHMENT - HYPOTHESIS GENERATION TOOL ⚠️

Pathway enrichment identifies biological processes potentially affected
by your gene/protein changes. These are HYPOTHESES, not proven mechanisms.

Interpretation Guidelines:
  • q-value < 0.05: Statistically enriched pathway
  • q-value < 0.01: Strongly enriched pathway
  • Overlap: Higher gene overlap = stronger evidence

Limitations:
  • Based on annotation databases (may be incomplete)
  • Does not prove causal mechanism
  • Context-dependent (cell type, species, condition)
  • Biased toward well-studied pathways

Follow-up Validation:
  1. Review individual genes in pathway
  2. Check literature for pathway relevance to your biological question
  3. Validate key pathway components with targeted experiments
  4. Consider pathway crosstalk and indirect effects
```

---

## Multi-Omics Integration

**For integrated RNA/Protein/Phospho analysis:**

```
ℹ️  MULTI-OMICS INTEGRATION NOTES

This analysis integrates RNA, protein, and phosphorylation data using
Stouffer's meta-analysis. Understanding the relationships:

RNA → Protein → Phosphorylation

Expected Concordance:
  • RNA↑ + Protein↑ = High confidence (transcriptional regulation)
  • RNA↔ + Protein↑ = Moderate confidence (post-transcriptional regulation)
  • RNA↑ + Protein↔ + Phospho↑ = High confidence (post-translational activation)

Discordance May Indicate:
  • Post-transcriptional regulation (miRNA, RNA binding proteins)
  • Post-translational modifications (ubiquitination, acetylation)
  • Protein stability changes
  • Technical artifacts or batch effects

Interpretation:
  • Concordant signals (all modalities agree) = HIGH confidence
  • Mixed signals = Requires mechanistic investigation
  • Discordant signals = May indicate complex regulation OR technical issues
```

---

## Spatial Transcriptomics

**For spatial analysis results:**

```
⚠️  SPATIAL TRANSCRIPTOMICS - RESOLUTION LIMITATIONS ⚠️

10x Visium spatial transcriptomics provides tissue-level resolution,
not single-cell resolution.

Technical Specifications:
  • Spot diameter: 55 μm
  • Center-to-center distance: 100 μm
  • Cells per spot: 1-10 (typically 3-5)
  • Genes detected: ~5,000-10,000 per spot

Limitations:
  • Cannot resolve individual cells within spots
  • Gene expression is AVERAGE across cells in each spot
  • Cell type deconvolution is computational inference
  • Spatial resolution lower than microscopy

Appropriate Uses:
  ✅ Tissue region identification
  ✅ Regional gene expression patterns
  ✅ Tumor microenvironment architecture
  ✅ Immune infiltration patterns

Inappropriate Uses:
  ❌ Single-cell resolution claims
  ❌ Precise cell-cell interaction mapping
  ❌ Exact cell type proportions (use as estimates only)
```

---

## Treatment Recommendations

**For any therapeutic suggestions:**

```
🔴 CRITICAL: TREATMENT RECOMMENDATIONS - EXPERIMENTAL ONLY 🔴

ANY treatment recommendations in this analysis are:
  • Based on COMPUTATIONAL predictions from genomic/transcriptomic data
  • NOT based on clinical trial evidence
  • NOT FDA-approved for this indication
  • NOT a substitute for clinical judgment

Required Before Clinical Consideration:
  1. ✅ Review by board-certified oncologist
  2. ✅ Tumor board discussion
  3. ✅ Published clinical evidence for indication
  4. ✅ FDA approval or clinical trial availability
  5. ✅ Patient eligibility assessment
  6. ✅ Risk/benefit analysis

Examples of Appropriate Use:
  ✅ Identifying potential drug targets for research
  ✅ Prioritizing genes for functional studies
  ✅ Generating hypotheses for clinical trial design
  ✅ Identifying resistance mechanisms for investigation

Examples of INAPPROPRIATE Use:
  ❌ Prescribing drugs directly based on predictions
  ❌ Making treatment decisions without clinical evidence
  ❌ Bypassing standard-of-care therapies
  ❌ Clinical use without IRB approval

REMEMBER: Computational predictions ≠ Clinical recommendations
```

---

## PatientOne Workflow Disclaimer

**For comprehensive PatientOne analysis outputs:**

```
╔═══════════════════════════════════════════════════════════════════════════╗
║                     PATIENTONE WORKFLOW DISCLAIMER                         ║
╚═══════════════════════════════════════════════════════════════════════════╝

This PatientOne analysis integrates:
  • Clinical data (demographics, treatment history)
  • Genomic variants (somatic mutations, CNVs)
  • Multi-omics data (RNA, Protein, Phosphorylation)
  • Spatial transcriptomics (tumor microenvironment)
  • Histology imaging (H&E, multiplex IF)

⚠️  CRITICAL LIMITATIONS:

1. RESEARCH USE ONLY
   This workflow has NOT been clinically validated.
   Do NOT use for patient care decisions.

2. SYNTHETIC/DE-IDENTIFIED DATA
   If using real patient data, ensure:
   ✅ IRB approval obtained
   ✅ Data de-identified per HIPAA guidelines
   ✅ Informed consent for research use

3. AI-GENERATED INSIGHTS
   Claude AI orchestrates the analysis but:
   • May misinterpret complex data
   • May hallucinate connections not supported by data
   • Cannot replace expert human interpretation
   • Requires validation by qualified bioinformaticians

4. VALIDATION REQUIRED
   ALL findings must be validated:
   • Statistical findings → Independent cohort
   • Pathway predictions → Mechanistic studies
   • Drug targets → Experimental models
   • Clinical recommendations → Clinical trials

5. NO LIABILITY
   The developers and AI provider assume NO liability for:
   • Treatment decisions based on this analysis
   • Patient outcomes
   • Misinterpretation of results
   • Technical errors or bugs

╔═══════════════════════════════════════════════════════════════════════════╗
║  CONSULT A QUALIFIED ONCOLOGIST BEFORE ANY CLINICAL DECISION               ║
╚═══════════════════════════════════════════════════════════════════════════╝
```

---

## Example: Complete Analysis Output with Disclaimers

```python
def generate_multiomics_report():
    """Example of properly disclaimed analysis output."""

    report = """
⚠️  RESEARCH USE ONLY - NOT FOR CLINICAL DECISION-MAKING ⚠️

═══════════════════════════════════════════════════════════════════════════
MULTI-OMICS TREATMENT RESISTANCE ANALYSIS
═══════════════════════════════════════════════════════════════════════════

Patient: PAT001-OVC-2025 (DE-IDENTIFIED)
Analysis Date: 2025-12-27
Modalities: RNA-seq, Proteomics, Phosphoproteomics (15 PDX samples)

───────────────────────────────────────────────────────────────────────────
KEY FINDINGS
───────────────────────────────────────────────────────────────────────────

1. PI3K/AKT/mTOR Pathway Activation ⭐⭐⭐ HIGH CONFIDENCE
   - PIK3CA: log2FC=2.3, q=0.0001 (RNA+Protein concordant)
   - AKT1 phospho(S473): log2FC=1.8, q=0.002
   - mTOR: log2FC=1.5, q=0.01

   Validation Required:
   → Western blot for p-AKT(S473), p-mTOR(S2448)
   → PI3K kinase activity assay

2. TP53 Loss of Function ⭐⭐⭐ HIGH CONFIDENCE
   - Mutation: TP53 c.524G>A (p.R175H) - known hotspot
   - RNA: log2FC=-2.1, q<0.0001
   - Protein: log2FC=-1.8, q=0.001

   Clinical Context:
   → TP53 mutations present in >96% of HGSOC
   → Validates expected HGSOC molecular profile

3. Immune Exclusion Phenotype ⭐⭐ MODERATE CONFIDENCE
   - CD8A: log2FC=-1.2, q=0.02
   - CD8B: log2FC=-1.1, q=0.03
   - Based on spatial transcriptomics (tumor core vs invasive margin)

   Limitations:
   → Spatial resolution: 55μm spots (3-5 cells per spot)
   → Requires validation with IHC for CD8+ T cells

───────────────────────────────────────────────────────────────────────────
PREDICTED THERAPEUTIC TARGETS (COMPUTATIONAL INFERENCE)
───────────────────────────────────────────────────────────────────────────

⚠️  These are COMPUTATIONAL PREDICTIONS, not clinical recommendations

1. PIK3CA Inhibitors (e.g., Alpelisib)
   Rationale: PIK3CA upregulation and pathway activation
   Evidence Level: ⭐⭐ MODERATE
   Clinical Status: FDA-approved for HR+ breast cancer, not ovarian

   ⚠️  OFF-LABEL for ovarian cancer
   ⚠️  Requires tumor board review and clinical trial availability

2. mTOR Inhibitors (e.g., Everolimus)
   Rationale: mTOR pathway activation
   Evidence Level: ⭐⭐ MODERATE
   Clinical Status: FDA-approved for other indications

   ⚠️  Limited clinical evidence in platinum-resistant ovarian cancer
   ⚠️  Consider in context of clinical trial

3. Immune Checkpoint Inhibitors (e.g., Pembrolizumab + combinations)
   Rationale: Immune exclusion suggests need for immune activation
   Evidence Level: ⭐ LOW (exploratory)

   ⚠️  HGSOC generally has low response to ICI monotherapy
   ⚠️  Consider combination strategies in clinical trial setting

───────────────────────────────────────────────────────────────────────────
DATA QUALITY ASSESSMENT
───────────────────────────────────────────────────────────────────────────

✅ GOOD DATA QUALITY

  • Missing values: RNA 8%, Protein 15%, Phospho 22%
  • Batch effects: Corrected (PC1 batch correlation: 0.68 → 0.12)
  • Sample size: 15 (7 resistant, 8 sensitive) - adequate power
  • Outliers: None detected (all within MAD < 2.5)

───────────────────────────────────────────────────────────────────────────
REQUIRED VALIDATION STEPS
───────────────────────────────────────────────────────────────────────────

Before clinical consideration:

1. ✅ Orthogonal Validation
   - Western blot: p-AKT(S473), p-mTOR(S2448), TP53
   - IHC: CD8+ T cell infiltration
   - qPCR: PIK3CA, AKT1, mTOR expression

2. ✅ Functional Studies
   - PI3K inhibitor sensitivity assays in patient-derived cells
   - Combination therapy testing (PI3K + immune checkpoint)

3. ✅ Clinical Review
   - Tumor board presentation
   - Review of standard-of-care options
   - Clinical trial matching
   - Patient eligibility and preferences

4. ✅ Regulatory Compliance
   - IRB approval for off-label use (if applicable)
   - Informed consent
   - Insurance pre-authorization

═══════════════════════════════════════════════════════════════════════════

⚠️  FINAL REMINDER:

This analysis is for RESEARCH PURPOSES ONLY.
All findings require experimental validation and clinical review.
Do NOT use for patient treatment decisions without proper validation.

Consult a board-certified oncologist before any clinical consideration.

═══════════════════════════════════════════════════════════════════════════

Analysis generated by: Precision Medicine MCP (mcp-multiomics v2.0)
Disclaimer version: 1.0 (2025-12-27)

"""
    return report
```

---

## Quick Reference: When to Use Each Disclaimer

| Output Type | Required Disclaimers |
|-------------|---------------------|
| **Statistical Results** | Research Use Only + Confidence Levels + Data Quality |
| **Pathway Enrichment** | Research Use Only + Pathway Disclaimer |
| **Upstream Regulators** | Research Use Only + Upstream Regulator Disclaimer |
| **Drug Predictions** | Research Use Only + Treatment Recommendations (CRITICAL) |
| **Spatial Analysis** | Research Use Only + Spatial Transcriptomics Limitations |
| **Multi-Omics Integration** | Research Use Only + Multi-Omics Notes + Data Quality |
| **PatientOne Complete** | PatientOne Workflow Disclaimer (comprehensive) |
| **Any Clinical Context** | Treatment Recommendations (CRITICAL) |

---

## Implementation Checklist

- [ ] Add `docs/DISCLAIMERS.md` to repository
- [ ] Update mcp-multiomics tools to include appropriate disclaimers
- [ ] Add PatientOne disclaimer to README and test documentation
- [ ] Include confidence levels in all statistical outputs
- [ ] Add data quality reporting to preprocessing tools
- [ ] Train users on disclaimer interpretation
- [ ] Review disclaimers with legal counsel (if deploying in clinical setting)

---

**Last Updated:** December 27, 2025
**Maintained by:** Precision Medicine MCP Team
**Review Schedule:** Update before any production deployment or change in use case
