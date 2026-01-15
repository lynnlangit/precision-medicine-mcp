# Clinician-in-the-Loop (CitL) Workflow Guide

**Target Audience:** Oncologists, Pathologists, Radiologists
**Purpose:** Step-by-step instructions for reviewing and validating AI-generated precision medicine reports
**Time Required:** 20-30 minutes per report
**Version:** 1.0
**Last Updated:** 2026-01-13

---

## Table of Contents

1. [Overview](#overview)
2. [Workflow Diagram](#workflow-diagram)
3. [Step-by-Step Instructions](#step-by-step-instructions)
4. [Section-by-Section Guidance](#section-by-section-guidance)
5. [Decision Scenarios](#decision-scenarios)
6. [Troubleshooting](#troubleshooting)
7. [Best Practices](#best-practices)
8. [FAQs](#faqs)

---

## Overview

### What is CitL Review?

The **Clinician-in-the-Loop (CitL) validation workflow** is a formal review and sign-off process where you, as the treating oncologist or specialist, validate AI-generated precision medicine analysis results before they are incorporated into clinical decision-making.

**Think of it as:** The AI system is your research assistant who has analyzed all the patient's multi-omics data (genomics, spatial transcriptomics, imaging) and prepared a draft report. Your job is to review the draft for accuracy, clinical validity, and guideline compliance before approving it for use.

### Why CitL Matters

**AI as Co-Pilot, Not Autopilot:**
- High-stakes oncology decisions require human expertise and judgment
- AI can miss clinical context, guideline nuances, or patient-specific factors
- Regulatory and ethical standards require human accountability for clinical decisions
- CitL creates an audit trail demonstrating appropriate oversight

**What You'll Validate:**
1. **Molecular findings accuracy** - Are the reported mutations, gene expression changes, and pathway alterations clinically valid?
2. **Clinical guideline compliance** - Do treatment recommendations follow NCCN and institutional protocols?
3. **Data quality** - Are sample sizes adequate? Are statistical thresholds met?
4. **Treatment appropriateness** - Are recommended therapies suitable for this specific patient?

### Your Role & Responsibility

**As the reviewing clinician, you attest that:**
- You have carefully reviewed all findings in the draft report
- You have assessed alignment with NCCN and institutional clinical guidelines
- You have evaluated automated quality flags and determined whether they are acceptable
- Your approval decision reflects your independent clinical judgment
- You understand this review becomes part of the permanent medical record (10-year retention per HIPAA)

**Decision Authority:**
- **APPROVE** - Report is accurate and ready for clinical use (proceed to tumor board, treatment planning)
- **REVISE** - Report requires modifications (specify issues, bioinformatician will re-analyze)
- **REJECT** - Report has critical errors requiring complete re-analysis (escalate to PI/bioinformatics team)

---

## Workflow Diagram

```
┌─────────────────────────────────────────────────────────────┐
│  AI Analysis Pipeline (Automated)                           │
│  ├─ TEST_1: Clinical data extraction (EHR/FHIR)            │
│  ├─ TEST_2: Multi-omics integration                        │
│  ├─ TEST_3: Spatial transcriptomics analysis               │
│  ├─ TEST_4: Imaging analysis                               │
│  └─ TEST_5: Cross-modal integration & synthesis            │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│  STEP 1: Generate Draft Report (Automated, ~30 sec)        │
│  ├─ Run quality checks (sample size, FDR, completeness)    │
│  ├─ Generate draft_report.json                             │
│  └─ Generate clinical_summary.txt                          │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│  STEP 2: Clinician Review (Manual, 20-30 min) ⭐ YOU       │
│  ├─ Open draft report files                                │
│  ├─ Complete CitL Review Template (7 sections)             │
│  ├─ Validate molecular findings (10 top findings)          │
│  ├─ Check NCCN guideline compliance                        │
│  ├─ Assess quality flags                                   │
│  └─ Make decision: APPROVE / REVISE / REJECT               │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│  STEP 3: Submit Review (Automated, ~5 sec)                 │
│  ├─ Validate review against JSON schema                    │
│  ├─ Generate digital signature (SHA-256)                   │
│  └─ Log to Cloud Logging audit trail                       │
└─────────────────────────────────────────────────────────────┘
                            ↓
                    ┌───────┴───────┐
                    │   Decision?   │
                    └───────┬───────┘
                            ↓
        ┌───────────────────┼───────────────────┐
        ↓                   ↓                   ↓
┌───────────────┐  ┌────────────────┐  ┌──────────────┐
│  APPROVE      │  │  REVISE        │  │  REJECT      │
├───────────────┤  ├────────────────┤  ├──────────────┤
│ Finalize      │  │ Re-analyze     │  │ Escalate to  │
│ report        │  │ with adjusted  │  │ PI/team      │
│               │  │ parameters     │  │              │
│ Ready for     │  │                │  │ Schedule     │
│ tumor board   │  │ Resubmit for   │  │ review       │
│               │  │ review         │  │ meeting      │
└───────────────┘  └────────────────┘  └──────────────┘
```

---

## Step-by-Step Instructions

### STEP 1: Access Draft Report (Automated - Already Complete)

**What Happened Before You Arrived:**
The bioinformatics team or automated pipeline has:
1. ✅ Run AI analysis pipeline (TEST_1 through TEST_5)
2. ✅ Generated draft report with quality checks
3. ✅ Notified you that a report is ready for review

**Files You'll Receive:**
```
results/PAT001-OVC-2025/
├── draft_report.json          # Structured findings (JSON format)
├── quality_checks.json         # Automated QC results
├── clinical_summary.txt        # Human-readable summary
└── citl_review_template.md     # Form to complete (copy provided)
```

**Your Action:** Open these files in your preferred editor or viewer.

---

### STEP 2: Review Clinical Summary (5-10 minutes)

**File to Open:** `clinical_summary.txt`

**What to Look For:**

1. **Patient Demographics & Context**
   - Verify patient ID matches your expectations
   - Review cancer type, stage, prior treatments
   - Note any contraindications or special considerations

2. **Quality Check Status**
   ```
   ✅ PASS: Sample sizes adequate (≥30 spots per region)
   ✅ PASS: FDR thresholds met (FDR < 0.05 for significant findings)
   ✅ PASS: Data completeness (>95% complete)
   ⚠️  WARNING: One region has minimum size (35 spots)
   ```

   **Interpretation:**
   - **All PASS:** Proceed with review, data quality is good
   - **Warnings present:** Note the warning, assess whether it affects clinical validity (you'll address this in Section 4)
   - **Critical failures:** Contact bioinformatics team before proceeding

3. **Top Molecular Findings**
   - Review list of top 10 genes/mutations/pathways identified
   - Note confidence scores (FDR values)
   - Identify any findings that seem inconsistent with clinical presentation

4. **Treatment Recommendations**
   - Review suggested therapies
   - Note molecular targets and expected efficacy
   - Consider patient-specific factors (comorbidities, prior treatments, performance status)

**Red Flags to Watch For:**
- ❌ Findings contradict known clinical presentation (e.g., ER+ markers in triple-negative breast cancer)
- ❌ Treatment recommendations contraindicated by patient history
- ❌ Suspicious FDR values (too many findings with FDR exactly = 0.05)
- ❌ Sample sizes below 30 spots in multiple regions

**If Red Flags Present:** Make note in your review form and consider REVISE or REJECT decision.

---

### STEP 3: Complete Review Form (15-20 minutes)

**File to Complete:** `citl_review_template.md` (copy to `citl_review_completed.json`)

**Overview of 7 Sections:**
1. **High-Level Decision** (REQUIRED) - APPROVE / REVISE / REJECT
2. **Per-Finding Validation** (REQUIRED) - Validate top 10 molecular findings
3. **Clinical Guideline Compliance** (REQUIRED) - NCCN + institutional protocols
4. **Quality Flags Assessment** - Assess automated QC flags
5. **Treatment Recommendations Review** - Evaluate treatment suggestions
6. **Attestation & Signature** (REQUIRED) - Digital signature
7. **Revision Instructions** - Complete if REVISE or REJECT

---

## Section-by-Section Guidance

### SECTION 1: High-Level Decision ⭐ REQUIRED

**What to Do:**
1. Select ONE decision: `[x] APPROVE` or `[x] REVISE` or `[x] REJECT`
2. Write 2-3 sentence rationale explaining your decision

**Decision Criteria:**

**APPROVE when:**
- All findings are consistent with clinical presentation
- No critical quality flags
- Treatment recommendations follow NCCN guidelines
- Data quality is adequate for clinical decision-making
- You are confident presenting these results at tumor board

**REVISE when:**
- Findings are mostly valid but require minor corrections
- Quality flags need addressing (e.g., exclude small regions)
- Treatment recommendations need refinement
- Additional tests or analyses would strengthen conclusions
- You can provide specific instructions for re-analysis

**REJECT when:**
- Critical errors in findings (e.g., wrong patient data, major inconsistencies)
- Data quality is inadequate (multiple failed QC checks)
- Results are clinically implausible
- Complete re-analysis required (not just parameter adjustments)

**Example Rationales:**

✅ **Good APPROVE rationale:**
> "All findings are consistent with clinical presentation and imaging. Molecular results align with observed platinum resistance. Treatment recommendations follow NCCN guidelines for recurrent HGSOC. Quality checks passed with one minor warning that does not affect clinical validity."

✅ **Good REVISE rationale:**
> "Findings are valid but tumor_necrotic region (35 spots) should be excluded due to insufficient sample size. FDR values for DEGs in this region may be unreliable. Request re-analysis excluding this region and increasing minimum spots per region to 50."

✅ **Good REJECT rationale:**
> "Critical inconsistency: Report identifies BRCA1 germline mutation, but patient's prior genetic testing (2024-03-15) showed BRCA wild-type. Suspect sample mix-up or data error. Recommend verifying sample identity before re-analysis."

❌ **Poor rationales:**
- "Looks good" (too vague, no clinical reasoning)
- "I disagree with the AI" (lacks specifics)
- "Not sure about this" (uncertain - use REVISE instead)

---

### SECTION 2: Per-Finding Validation ⭐ REQUIRED

**What to Do:**
For each of the 10 top molecular findings listed in `draft_report.json`, indicate validation status:
- `[x] CONFIRMED` - Finding is clinically valid and supported by data
- `[x] UNCERTAIN` - Requires additional validation or expert consultation
- `[x] INCORRECT` - Finding appears to be erroneous

**How to Validate Each Finding:**

**Finding Structure:**
```
Finding 1: TP53 (R175H) - Somatic Mutation
- Finding ID: MUT_1
- Confidence Score: FDR = 5.04e-20
- Clinical Significance: Tumor suppressor loss, TP53 pathway disrupted
- Validation Status: [CONFIRM / UNCERTAIN / INCORRECT]
- Comments: [Your assessment]
```

**Validation Questions:**

1. **Does this finding align with clinical presentation?**
   - Example: TP53 mutations are expected in 96% of HGSOC cases → CONFIRMED
   - Example: HER2 amplification in histologically confirmed TNBC → INCORRECT

2. **Is the confidence score appropriate?**
   - FDR < 0.01 → High confidence
   - FDR 0.01-0.05 → Moderate confidence (acceptable)
   - FDR > 0.05 → Low confidence → Mark UNCERTAIN, request re-analysis

3. **Does the clinical significance make sense?**
   - Check if gene function matches reported effect
   - Example: "BCL2L1 overexpression → apoptosis resistance" → Biologically plausible → CONFIRMED

4. **Can you verify with other data?**
   - Cross-reference with imaging (e.g., high Ki67 expression should correlate with high proliferation on imaging)
   - Cross-reference with prior genomic reports if available
   - Check protein expression (IHC) if available

**Examples:**

✅ **CONFIRMED Example:**
```
Finding 1: TP53 (R175H) - Somatic Mutation
- Validation Status: [x] CONFIRMED
- Comments: "TP53 R175H is a hotspot mutation in HGSOC (present in >95% of cases).
  Consistent with p53 staining pattern on IHC (diffuse nuclear positivity).
  Very high confidence (FDR = 5.04e-20). CONFIRMED."
```

⚠️ **UNCERTAIN Example:**
```
Finding 5: NRAS (Q61K) - Somatic Mutation
- Validation Status: [x] UNCERTAIN
- Comments: "NRAS mutations uncommon in HGSOC (<5% prevalence). Confidence score
  is moderate (FDR = 0.03). Recommend orthogonal validation with targeted
  sequencing (Sanger or ddPCR) before incorporating into treatment decisions."
```

❌ **INCORRECT Example:**
```
Finding 8: EGFR (Amplification) - Copy Number Variant
- Validation Status: [x] INCORRECT
- Comments: "EGFR amplification is rare in HGSOC (<2% prevalence) and not supported
  by IHC staining (EGFR 0/3+). Suspected false positive. Recommend excluding
  this finding from final report."
```

**Pro Tip:** If you have 2+ UNCERTAIN or 1+ INCORRECT findings, consider REVISE decision in Section 1.

---

### SECTION 3: Clinical Guideline Compliance ⭐ REQUIRED

**What to Do:**
Assess whether treatment recommendations align with:
1. **NCCN Guidelines** (National Comprehensive Cancer Network)
2. **Institutional Protocols** (your hospital's formulary and treatment pathways)

**For Each Guideline Set:**

Select alignment level:
- `[x] ALIGNED` - All recommendations follow guidelines
- `[x] PARTIAL` - Some recommendations deviate (explain below)
- `[x] NOT_ALIGNED` - Significant deviations from guidelines (explain below)

**How to Assess:**

1. **Open NCCN Guidelines for Patient's Cancer Type**
   - Example: NCCN Guidelines for Ovarian Cancer (current version)
   - Navigate to relevant section (e.g., "Recurrent or Persistent Disease")

2. **Compare Recommended Therapies**
   - Check if suggested drugs are in NCCN-recommended regimens
   - Verify sequencing (e.g., platinum-based before PARP inhibitors)
   - Check Category ratings (Category 1 preferred, Category 2A acceptable, Category 2B conditional)

3. **Document Deviations**
   - If recommendation is Category 2B or off-label, explain justification
   - If recommendation deviates from institutional formulary, note alternative options

**Examples:**

✅ **ALIGNED Example:**
```
NCCN Guidelines Alignment: [x] ALIGNED
NCCN Deviations: None

Institutional Protocol Alignment: [x] ALIGNED
Institutional Deviations: None

Comment: All recommended therapies (carboplatin/paclitaxel, olaparib, bevacizumab)
are NCCN Category 1 or 2A recommendations for platinum-sensitive recurrent HGSOC.
All drugs are on institutional formulary. No deviations.
```

⚠️ **PARTIAL Example:**
```
NCCN Guidelines Alignment: [x] PARTIAL
NCCN Deviations:
- "Alpelisib (PI3K inhibitor) is Category 2B (limited data) for PIK3CA-mutant
  ovarian cancer. NCCN prefers clinical trial enrollment. However, patient has
  PIK3CA E545K mutation with strong rationale for PI3K targeting. Recommend
  discussing with patient: clinical trial (preferred) vs. off-label alpelisib."

Institutional Protocol Alignment: [x] ALIGNED
Institutional Deviations: None
```

❌ **NOT_ALIGNED Example (should trigger REVISE):**
```
NCCN Guidelines Alignment: [x] NOT_ALIGNED
NCCN Deviations:
- "Recommendation for imatinib (PDGFR inhibitor) is not supported by NCCN
  guidelines for HGSOC. Imatinib has no established efficacy in this disease.
  Recommend removing this recommendation."

Institutional Protocol Alignment: [x] NOT_ALIGNED
Institutional Deviations:
- "Bevacizumab continuation beyond 15 months is not per institutional protocol
  (protocol limits to 12-15 months due to GI perforation risk). Recommend
  discontinuing bevacizumab at 15 months per protocol."
```

**Pro Tip:** PARTIAL alignment is common and acceptable if deviations are justified. NOT_ALIGNED typically requires REVISE decision.

---

### SECTION 4: Quality Flags Assessment

**What to Do:**
Review automated quality flags raised during analysis. For each flag, assess:
- `[x] ACCEPTABLE` - Flag noted but not concerning, does not affect clinical validity
- `[x] REQUIRES_ACTION` - Flag must be addressed before approval

**Common Quality Flags:**

**1. Sample Size Warning**
```
Flag: sample_size_warning
Severity: warning
Message: "Minimum region size: 35 spots (tumor_necrotic region)"
Recommendation: "Prefer ≥50 spots for robust statistical power. Consider excluding
                small regions or interpreting results with caution."

Your Assessment: [ACCEPTABLE / REQUIRES_ACTION]
```

**When to mark ACCEPTABLE:**
- Region is 35-49 spots AND FDR values are very low (< 1e-10) → Strong signal compensates for smaller sample
- Region is biologically unimportant (e.g., necrotic tissue) → Can be excluded without affecting conclusions

**When to mark REQUIRES_ACTION:**
- Region is <35 spots → Too small for reliable statistics
- FDR values are marginal (0.01-0.05) → Insufficient power → Request re-analysis excluding small regions

**2. Missing Data Warning**
```
Flag: missing_data_warning
Severity: warning
Message: "5% of spots have missing gene expression values (technical dropout)"
Recommendation: "Acceptable if <10%. Values were imputed using k-NN method."

Your Assessment: [ACCEPTABLE / REQUIRES_ACTION]
```

**When ACCEPTABLE:** <10% missing, imputation method is appropriate
**When REQUIRES_ACTION:** >10% missing, indicates poor data quality

**3. Batch Effect Warning**
```
Flag: batch_effect_detected
Severity: critical
Message: "Significant batch effect detected between sample batches (p < 0.001)"
Recommendation: "Apply batch correction (ComBat or Harmony) before analysis."

Your Assessment: [ACCEPTABLE / REQUIRES_ACTION]
```

**When ACCEPTABLE:** Batch correction was already applied, residual effect is minimal
**When REQUIRES_ACTION:** No batch correction applied → REVISE and request correction

**Example Assessment:**
```
Flag 1: sample_size_warning (tumor_necrotic region, 35 spots)
Reviewer Assessment: [x] ACCEPTABLE
Comments: "Minimum region size of 35 spots is acceptable. While below ideal threshold
of 50, the strong effect sizes (log2FC > 3) and very low FDR values (< 1e-10) provide
sufficient confidence. Additionally, tumor_necrotic region findings are corroborated
by imaging (necrotic areas on CT). ACCEPTABLE."
```

**Pro Tip:** If ANY flag is marked REQUIRES_ACTION, you must select REVISE in Section 1 and provide re-analysis instructions in Section 7.

---

### SECTION 5: Treatment Recommendations Review

**What to Do:**
For each treatment recommendation, indicate:
- `[x] AGREE` - Recommendation is appropriate for this patient
- `[x] DISAGREE` - Recommendation is not appropriate (explain below)

**Evaluation Criteria:**

1. **Molecular Target Match**
   - Does the patient have the molecular alteration targeted by this therapy?
   - Example: PIK3CA inhibitor requires PIK3CA mutation → Check if present

2. **Expected Efficacy**
   - Is the reported response rate realistic based on literature?
   - Is the patient likely to benefit given disease stage and prior treatments?

3. **Patient-Specific Factors**
   - Contraindications (e.g., bevacizumab contraindicated if recent GI perforation)
   - Performance status (ECOG 0-2 required for most trials)
   - Organ function (adequate renal/hepatic function)
   - Prior treatment history (e.g., platinum-free interval)

4. **Availability**
   - Is this therapy available at your institution?
   - Is it covered by insurance?
   - Are there clinical trials available?

**Examples:**

✅ **AGREE Example:**
```
Recommendation 1: PI3K inhibitor (Alpelisib) + PARP inhibitor (Olaparib)
Molecular Target: PIK3CA E545K mutation + BRCA1 germline mutation
Expected Efficacy: 40-50% ORR in PIK3CA-mutant, PARP-sensitive tumors

Agreement: [x] AGREE
Comments: "Patient has both PIK3CA E545K (confirmed by NGS) and BRCA1 germline
mutation (confirmed by prior germline testing). Platinum-free interval is 8 months
(platinum-sensitive). Combination is rational based on synthetic lethality. Patient
meets eligibility criteria for clinical trial NCT12345678. AGREE - recommend
discussing clinical trial enrollment."
```

❌ **DISAGREE Example:**
```
Recommendation 3: Bevacizumab continuation (cycle 20+)
Molecular Target: VEGFA overexpression
Expected Efficacy: Modest PFS benefit (2-3 months)

Agreement: [x] DISAGREE
Comments: "Patient has already received 19 cycles of bevacizumab (>18 months).
Institutional protocol limits bevacizumab to 12-15 months due to cumulative GI
perforation risk. Patient has history of diverticulosis (relative contraindication).
Risk exceeds benefit at this point. DISAGREE - recommend discontinuing bevacizumab
and considering alternative anti-angiogenic approaches if indicated."
```

**Pro Tip:** It's okay to DISAGREE with some recommendations - this demonstrates appropriate clinical oversight. Document your reasoning clearly.

---

### SECTION 6: Attestation & Signature ⭐ REQUIRED

**What to Do:**
1. Check ALL attestation boxes (must be checked to submit):
   - `[x]` I have reviewed all findings in this report
   - `[x]` I have assessed clinical guideline compliance
   - `[x]` I have evaluated quality flags
   - `[x]` My decision reflects my clinical judgment
   - `[x]` I understand this review is part of the medical record

2. Fill in reviewer information:
   - Full Name: `Dr. Sarah Johnson`
   - Credentials: `MD, Gynecologic Oncology`
   - Email: `sarah.johnson@hospital.org`
   - Role: `oncologist`
   - Date: `2026-01-13` (YYYY-MM-DD)
   - Time: `14:30 PST` (HH:MM timezone)

3. Leave Digital Signature blank (auto-generated by submission script)

**Legal & Regulatory Notes:**

⚠️ **This attestation is legally binding:**
- Your review becomes part of the permanent medical record
- 10-year retention per HIPAA requirements
- Audit trail includes your name, email, timestamp, and decision
- Digital signature (SHA-256 hash) ensures immutability

⚠️ **You are attesting that:**
- You personally reviewed the report (not delegated to resident/fellow)
- Your decision reflects your independent clinical judgment (not influenced by administrative pressure)
- You are qualified to evaluate this cancer type and molecular findings

---

### SECTION 7: Revision Instructions (Required if REVISE or REJECT)

**When to Complete:** ONLY if you selected REVISE or REJECT in Section 1

**What to Do:**

1. **List Specific Issues to Address** (minimum 1 issue)
   ```
   Issue 1: Exclude tumor_necrotic region (35 spots) due to insufficient sample size.
            FDR values for DEGs in this region are unreliable.

   Issue 2: Re-run differential expression analysis with stricter FDR threshold
            (FDR < 0.01 instead of FDR < 0.05) to reduce false positives.

   Issue 3: Validate TP53 mutation with orthogonal sequencing method (Sanger or ddPCR)
            to confirm hotspot mutation.
   ```

2. **Specify Re-analysis Parameters**
   ```
   Regions to exclude:
   - tumor_necrotic (n=35 spots)

   FDR threshold adjustment:
   - Change from FDR < 0.05 to FDR < 0.01 for differential expression

   Minimum spots per region:
   - Increase from 30 to 50 spots minimum per region

   Additional tests required:
   - Sanger sequencing confirmation for TP53 R175H
   ```

3. **Set Resubmission Timeline**
   ```
   Expected resubmission date: 2026-01-20
   Estimated time for corrections: 3-5 days
   ```

**Pro Tip:** Be specific! Vague instructions like "improve data quality" are not actionable. Provide exact parameters to adjust.

---

## Decision Scenarios

### Scenario 1: Clean Approval (Most Common)

**Situation:**
- All quality checks passed
- Findings consistent with clinical presentation
- NCCN-aligned treatment recommendations
- No critical quality flags

**Your Actions:**
1. Section 1: Select **APPROVE**
2. Section 2: Mark 8-10 findings as **CONFIRMED**, 0-2 as **UNCERTAIN**
3. Section 3: Mark **ALIGNED** for both NCCN and institutional
4. Section 4: Mark all flags as **ACCEPTABLE** (if any)
5. Section 5: **AGREE** with most/all treatment recommendations
6. Section 6: Complete attestation
7. Section 7: Leave blank (not required for APPROVE)

**Estimated Time:** 20-25 minutes

**Next Steps:** Submit review → Finalize report → Present at tumor board

---

### Scenario 2: Revise Due to Quality Flag

**Situation:**
- Findings are valid but one region has inadequate sample size
- Quality flag raised: sample_size_warning
- You determine sample size is too small for reliable statistics

**Your Actions:**
1. Section 1: Select **REVISE**
   - Rationale: "Quality flag for tumor_necrotic region (35 spots) requires action. Sample size below threshold for reliable differential expression. Request re-analysis excluding this region."

2. Section 2: Mark findings from adequate regions as **CONFIRMED**, findings from small region as **UNCERTAIN**

3. Section 3: Mark **ALIGNED** (no guideline issues)

4. Section 4: Mark sample_size_warning as **REQUIRES_ACTION**
   - Comment: "35 spots insufficient for robust statistics. FDR values marginal (0.02-0.04). Exclude tumor_necrotic region and re-run analysis."

5. Section 5: **AGREE** with treatments (assuming they're not based on unreliable findings)

6. Section 6: Complete attestation

7. Section 7: **MUST COMPLETE**
   ```
   Issues to Address:
   1. Exclude tumor_necrotic region (35 spots) from differential expression analysis

   Re-analysis Parameters:
   - Regions to exclude: tumor_necrotic
   - Minimum spots per region: Increase to 50 spots

   Resubmission Timeline: 2026-01-18 (3 days)
   ```

**Estimated Time:** 25-35 minutes (additional time for Section 7)

**Next Steps:** Submit review → Bioinformatician re-runs analysis → You review version 2

---

### Scenario 3: Reject Due to Critical Error

**Situation:**
- Findings contradict known clinical data
- Example: Report shows BRCA1 mutation but prior genetic testing was wild-type
- Suspect sample mix-up or data error

**Your Actions:**
1. Section 1: Select **REJECT**
   - Rationale: "Critical inconsistency: Report identifies BRCA1 germline mutation, but patient's prior genetic testing (2024-03-15, GeneDx report #12345) showed BRCA wild-type. Suspect sample mix-up or data error. Cannot proceed with clinical decision-making until discrepancy is resolved. Recommend verifying sample identity (fingerprinting) and re-running from raw data."

2. Section 2: Mark contradictory findings as **INCORRECT**, others as **UNCERTAIN**

3. Section 3: Cannot assess (report validity in question)

4. Section 4: Mark critical flags as **REQUIRES_ACTION**

5. Section 5: Cannot assess (recommendations based on invalid data)

6. Section 6: Complete attestation

7. Section 7: **MUST COMPLETE**
   ```
   Issues to Address:
   1. Verify sample identity - confirm correct patient specimen used
   2. Reconcile BRCA1 discrepancy (report shows BRCA1 mutation, prior test wild-type)
   3. If sample mix-up confirmed, re-run entire analysis pipeline with correct sample
   4. If data processing error, re-run from raw sequencing data

   Re-analysis Parameters:
   - Verify sample barcode matches patient ID
   - Perform sample fingerprinting (SNP concordance check)
   - If sample correct: Investigate BRCA1 call (tumor vs. germline, technical artifact)

   Resubmission Timeline: 2026-01-27 (2 weeks for investigation + re-analysis)
   ```

**Estimated Time:** 30-40 minutes (requires investigation and documentation)

**Next Steps:** Submit review → Escalate to PI/bioinformatics → Schedule team meeting → Corrective action → Complete re-analysis

---

## Troubleshooting

### Problem: Draft report files missing or incomplete

**Symptoms:**
- `draft_report.json` not found
- Files are empty or corrupted

**Solutions:**
1. Contact bioinformatics team - analysis may still be running
2. Check file paths - ensure you're looking in correct patient directory
3. Re-generate draft report:
   ```bash
   python tools/reports/generate_patient_report.py --patient-id PAT001-OVC-2025 --output-dir ./results --generate-draft
   ```

---

### Problem: I don't understand a molecular finding

**Example:** "What does 'PIK3CA E545K' mean?"

**Solutions:**
1. **Look up gene function:**
   - GeneCards: https://www.genecards.org/
   - OMIM: https://www.omim.org/
   - Example: PIK3CA = phosphoinositide 3-kinase, E545K = hotspot mutation in helical domain

2. **Check clinical significance:**
   - ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/
   - OncoKB: https://www.oncokb.org/
   - Example: PIK3CA E545K = Oncogenic, targetable with PI3K inhibitors

3. **Mark as UNCERTAIN and request expert consultation:**
   ```
   Finding: PIK3CA E545K
   Validation Status: [x] UNCERTAIN
   Comments: "Requesting molecular pathology consultation to confirm clinical
   significance and treatment implications before finalizing approval."
   ```

---

### Problem: Treatment recommendation seems off-label or experimental

**Example:** Report suggests drug not in standard NCCN guidelines

**Solutions:**
1. **Check if supported by molecular rationale:**
   - Is there a strong biological justification?
   - Is there clinical trial data or case reports supporting use?

2. **Consult NCCN guidelines carefully:**
   - Drug may be Category 2B (acceptable with lower evidence)
   - May be recommended for specific molecular subtype

3. **Document in Section 3 (Guidelines) and Section 5 (Treatments):**
   ```
   NCCN Alignment: [x] PARTIAL
   Deviation: "Alpelisib is Category 2B for PIK3CA-mutant ovarian cancer.
   Prefer clinical trial enrollment per NCCN. However, patient has exhausted
   standard options and has strong molecular rationale (PIK3CA E545K).
   Acceptable as off-label use after discussing trial availability."
   ```

4. **Mark AGREE with caveat:**
   ```
   Recommendation: Alpelisib (PI3K inhibitor)
   Agreement: [x] AGREE
   Comments: "Off-label use acceptable given PIK3CA E545K mutation and exhausted
   standard therapies. Recommend prior authorization and informed consent
   discussion with patient."
   ```

---

### Problem: Quality flag is unclear or I'm unsure if it's acceptable

**Example:** "Batch effect detected - is this a problem?"

**Solutions:**
1. **Read the flag details carefully:**
   - What is the severity? (critical / warning / info)
   - What does the recommendation say?
   - Has correction been applied?

2. **Consult with bioinformatics team:**
   - Ask: "Was batch correction applied?"
   - Ask: "Does this affect findings validity?"

3. **If unsure, mark REQUIRES_ACTION and request clarification:**
   ```
   Flag: batch_effect_detected
   Reviewer Assessment: [x] REQUIRES_ACTION
   Comments: "Unclear whether batch correction was applied. Request confirmation
   from bioinformatics team. If not applied, request re-analysis with ComBat or
   Harmony batch correction before approval."
   ```

4. **Select REVISE in Section 1** and document in Section 7

---

### Problem: Schema validation errors when submitting review

**Symptoms:**
```
❌ Schema validation failed:
   'status' is a required property
   Path: decision → status
```

**Solutions:**
1. **Check required fields:**
   - Section 1: Decision status and rationale
   - Section 2: At least one finding validated
   - Section 3: NCCN and institutional alignment selected
   - Section 6: All attestation checkboxes checked

2. **Check data types:**
   - Dates must be YYYY-MM-DD format
   - Email must be valid email format
   - Enums must exactly match (e.g., "APPROVE" not "Approve")

3. **Check conditional requirements:**
   - If status = REVISE or REJECT → Section 7 is required
   - If quality flags exist → Section 4 is required

4. **Validate JSON syntax:**
   - Use online JSON validator: https://jsonlint.com/
   - Check for missing commas, brackets, quotes

---

## Best Practices

### Time Management

**Efficient Review (20-25 minutes):**
1. **Skim clinical summary (5 min)** - Get high-level overview
2. **Section 1 decision (2 min)** - Quick initial assessment (APPROVE/REVISE/REJECT)
3. **Section 2 findings (8 min)** - Validate top 10 findings (~1 min each)
4. **Section 3 guidelines (3 min)** - Check NCCN alignment
5. **Section 4-5 (5 min)** - Quick flag and treatment review
6. **Section 6 attestation (2 min)** - Complete signature
7. **Buffer (5 min)** - Double-check, finalize

**When to Take More Time:**
- Complex molecular findings you're unfamiliar with (consult references)
- Multiple quality flags requiring assessment
- Off-label or experimental treatment recommendations
- REVISE or REJECT decisions (Section 7 requires detailed instructions)

---

### Communication with Bioinformatics Team

**When to Reach Out:**
1. **Before starting review:** If files are missing or analysis incomplete
2. **During review:** If you encounter unclear findings or flags
3. **After REVISE decision:** Confirm they understand re-analysis instructions
4. **After REJECT decision:** Schedule meeting to discuss corrective action

**How to Communicate:**
- Be specific: "Exclude tumor_necrotic region" not "improve data quality"
- Provide rationale: Explain why change is needed
- Set timeline: Expected resubmission date
- Use patient ID consistently: Avoid confusion with multiple patients

---

### Documentation Standards

**Good Comments:**
✅ Specific, clinical reasoning
✅ References to supporting data (IHC, imaging, prior reports)
✅ Actionable instructions (for REVISE)

**Example:**
> "TP53 R175H confirmed. Hotspot mutation present in 96% of HGSOC cases (TCGA data). Consistent with p53 IHC staining pattern (diffuse nuclear positivity, consistent with missense mutation). Very high confidence (FDR = 5.04e-20). CONFIRMED."

**Poor Comments:**
❌ Vague, no reasoning
❌ No clinical context
❌ Non-actionable

**Example:**
> "Looks good"

---

### Quality Assurance

**Self-Check Before Submitting:**
- [ ] All REQUIRED sections completed (1, 2, 3, 6)
- [ ] Decision (APPROVE/REVISE/REJECT) matches findings
- [ ] If REVISE/REJECT: Section 7 completed with specific instructions
- [ ] If APPROVE: You're confident presenting at tumor board
- [ ] Attestation checkboxes all checked
- [ ] Reviewer information complete and accurate
- [ ] No typos in patient ID or critical fields

---

## FAQs

**Q: How long does a typical review take?**
A: 20-30 minutes for straightforward APPROVE, 25-35 minutes for REVISE, 30-40 minutes for REJECT.

**Q: Can I delegate review to a fellow or resident?**
A: No. The attestation requires that YOU personally reviewed the report. Fellows/residents can assist with data review, but final validation and signature must be yours.

**Q: What if I'm not confident in molecular findings?**
A: Mark findings as UNCERTAIN and request expert consultation (molecular pathologist, geneticist). You can still APPROVE with UNCERTAIN findings if they don't affect treatment decisions, or select REVISE and request additional validation.

**Q: What happens after I submit an APPROVE review?**
A: The finalize_patient_report.py script generates a final approved report with status "clinically_approved". This report is ready for tumor board presentation and clinical decision-making.

**Q: What happens after I submit a REVISE review?**
A: Bioinformatics team receives your revision instructions (Section 7), re-runs analysis with adjusted parameters, and generates a new draft report (version 2). You'll then review version 2 using the same CitL workflow.

**Q: What happens after I submit a REJECT review?**
A: Your REJECT decision triggers escalation to PI/bioinformatics team lead. A team meeting is scheduled to investigate the critical error, determine corrective action, and create a plan for re-analysis.

**Q: Can I change my decision after submitting?**
A: No. Once submitted, the review is immutable (digital signature ensures tamper-evidence). If you realize you made an error, contact bioinformatics team and PI to discuss next steps. A new review version may be required.

**Q: What if I disagree with NCCN guidelines?**
A: Document your reasoning in Section 3. If you have strong clinical justification for deviating from guidelines (e.g., patient-specific factors, recent literature), you can mark PARTIAL and explain the deviation. This is acceptable and demonstrates appropriate clinical judgment.

**Q: What if quality checks failed but findings still seem valid?**
A: Assess each flag individually in Section 4. Some flags are informational and don't affect validity (mark ACCEPTABLE). Others require action (mark REQUIRES_ACTION and select REVISE). Use your clinical judgment.

**Q: Can I approve a report with UNCERTAIN findings?**
A: Yes, IF the uncertain findings don't materially affect treatment decisions. Document your reasoning: "Finding X marked UNCERTAIN pending additional validation, but primary treatment recommendation (Y) is supported by other CONFIRMED findings." If uncertain findings are critical to treatment decisions, select REVISE and request additional validation.

**Q: How do I convert the markdown template to JSON?**
A: Currently manual (copy values from .md to .json following schema structure). A web-based form interface is planned for future release to automate this conversion.

**Q: Where can I find example completed reviews?**
A: See `docs/hospital-deployment/citl-workflows/CITL_EXAMPLES.md` for three complete example reviews (APPROVE, REVISE, REJECT scenarios).

**Q: Who can I contact for help?**
A:
- Technical questions (files, scripts, schema): Bioinformatics team
- Clinical questions (molecular findings): Molecular pathology or genetics
- Workflow questions (process, timeline): Clinic administrator or PI

---

## Summary Checklist

Before submitting your review, confirm:

**Required Sections Complete:**
- [ ] Section 1: Decision (APPROVE/REVISE/REJECT) and rationale
- [ ] Section 2: All 10 findings validated (CONFIRMED/UNCERTAIN/INCORRECT)
- [ ] Section 3: NCCN and institutional guideline compliance assessed
- [ ] Section 6: All attestation boxes checked, reviewer info complete

**Conditional Sections:**
- [ ] Section 4: Quality flags assessed (if flags exist)
- [ ] Section 5: Treatment recommendations reviewed (if recommendations exist)
- [ ] Section 7: Revision instructions (if REVISE or REJECT)

**Quality Checks:**
- [ ] Decision is consistent with findings (e.g., APPROVE but all findings INCORRECT → inconsistent)
- [ ] Comments are specific and clinically justified
- [ ] Patient ID is correct
- [ ] You're confident in your decision

**Ready to Submit:**
```bash
python tools/reports/citl_submit_review.py \
  --patient-id PAT001-OVC-2025 \
  --review-file ./results/PAT001-OVC-2025/citl_review_completed.json
```

---

**Document Version:** 1.0
**Last Updated:** 2026-01-13
**Part of:** Precision Medicine MCP - Clinician-in-the-Loop Validation Workflow
**Contact:** bioinformatics@hospital.org
