# Clinician-in-the-Loop (CitL) Review Form

**Patient ID:** {{PATIENT_ID}}
**Report Date:** {{REPORT_DATE}}
**Reviewer:** {{REVIEWER_NAME}}
**Review Date:** {{REVIEW_DATE}}

---

## Instructions

This form captures your clinical review and validation of AI-generated precision medicine analysis results. Complete all required sections below. This review becomes part of the permanent medical record and audit trail (10-year retention per HIPAA).

**Time estimate:** 20-30 minutes
**Required:** All sections except comments fields
**Format:** Convert to JSON using provided template or automation tool

---

## SECTION 1: HIGH-LEVEL DECISION ⭐ REQUIRED

**Overall Assessment:**
- [ ] **APPROVE** - Report is accurate and ready for clinical use
- [ ] **REVISE** - Report requires modifications (complete Section 7)
- [ ] **REJECT** - Report has critical errors, requires re-analysis (complete Section 7)

**Rationale for Decision:** (2-3 sentences explaining your decision)

```
[Your rationale here - must be 10-1000 characters]
```

---

## SECTION 2: PER-FINDING VALIDATION ⭐ REQUIRED

For each major molecular finding in the draft report, indicate validation status. Review the top 10 findings listed in `draft_report.json` under `key_molecular_findings`.

### Finding 1: {{GENE_NAME_1}} - {{FINDING_TYPE_1}}
- **Finding ID:** {{FINDING_ID_1}}
- **Confidence Score:** FDR = {{FDR_VALUE_1}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_1}}
- **Validation Status:**
  - [ ] **CONFIRMED** - Finding is clinically valid and supported by data
  - [ ] **UNCERTAIN** - Requires additional validation or expert consultation
  - [ ] **INCORRECT** - Finding appears to be erroneous
- **Comments:** (Optional - explain your validation decision)
  ```
  [Your comments here]
  ```

### Finding 2: {{GENE_NAME_2}} - {{FINDING_TYPE_2}}
- **Finding ID:** {{FINDING_ID_2}}
- **Confidence Score:** FDR = {{FDR_VALUE_2}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_2}}
- **Validation Status:**
  - [ ] CONFIRMED
  - [ ] UNCERTAIN
  - [ ] INCORRECT
- **Comments:**
  ```
  [Your comments here]
  ```

### Finding 3: {{GENE_NAME_3}} - {{FINDING_TYPE_3}}
- **Finding ID:** {{FINDING_ID_3}}
- **Confidence Score:** FDR = {{FDR_VALUE_3}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_3}}
- **Validation Status:**
  - [ ] CONFIRMED
  - [ ] UNCERTAIN
  - [ ] INCORRECT
- **Comments:**
  ```
  [Your comments here]
  ```

### Finding 4: {{GENE_NAME_4}} - {{FINDING_TYPE_4}}
- **Finding ID:** {{FINDING_ID_4}}
- **Confidence Score:** FDR = {{FDR_VALUE_4}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_4}}
- **Validation Status:**
  - [ ] CONFIRMED
  - [ ] UNCERTAIN
  - [ ] INCORRECT
- **Comments:**
  ```
  [Your comments here]
  ```

### Finding 5: {{GENE_NAME_5}} - {{FINDING_TYPE_5}}
- **Finding ID:** {{FINDING_ID_5}}
- **Confidence Score:** FDR = {{FDR_VALUE_5}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_5}}
- **Validation Status:**
  - [ ] CONFIRMED
  - [ ] UNCERTAIN
  - [ ] INCORRECT
- **Comments:**
  ```
  [Your comments here]
  ```

### Finding 6: {{GENE_NAME_6}} - {{FINDING_TYPE_6}}
- **Finding ID:** {{FINDING_ID_6}}
- **Confidence Score:** FDR = {{FDR_VALUE_6}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_6}}
- **Validation Status:**
  - [ ] CONFIRMED
  - [ ] UNCERTAIN
  - [ ] INCORRECT
- **Comments:**
  ```
  [Your comments here]
  ```

### Finding 7: {{GENE_NAME_7}} - {{FINDING_TYPE_7}}
- **Finding ID:** {{FINDING_ID_7}}
- **Confidence Score:** FDR = {{FDR_VALUE_7}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_7}}
- **Validation Status:**
  - [ ] CONFIRMED
  - [ ] UNCERTAIN
  - [ ] INCORRECT
- **Comments:**
  ```
  [Your comments here]
  ```

### Finding 8: {{GENE_NAME_8}} - {{FINDING_TYPE_8}}
- **Finding ID:** {{FINDING_ID_8}}
- **Confidence Score:** FDR = {{FDR_VALUE_8}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_8}}
- **Validation Status:**
  - [ ] CONFIRMED
  - [ ] UNCERTAIN
  - [ ] INCORRECT
- **Comments:**
  ```
  [Your comments here]
  ```

### Finding 9: {{GENE_NAME_9}} - {{FINDING_TYPE_9}}
- **Finding ID:** {{FINDING_ID_9}}
- **Confidence Score:** FDR = {{FDR_VALUE_9}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_9}}
- **Validation Status:**
  - [ ] CONFIRMED
  - [ ] UNCERTAIN
  - [ ] INCORRECT
- **Comments:**
  ```
  [Your comments here]
  ```

### Finding 10: {{GENE_NAME_10}} - {{FINDING_TYPE_10}}
- **Finding ID:** {{FINDING_ID_10}}
- **Confidence Score:** FDR = {{FDR_VALUE_10}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_10}}
- **Validation Status:**
  - [ ] CONFIRMED
  - [ ] UNCERTAIN
  - [ ] INCORRECT
- **Comments:**
  ```
  [Your comments here]
  ```

---

## SECTION 3: CLINICAL GUIDELINE COMPLIANCE ⭐ REQUIRED

### NCCN Guidelines Alignment
Review treatment recommendations against current NCCN (National Comprehensive Cancer Network) guidelines for this cancer type.

- [ ] **ALIGNED** - All recommendations follow NCCN guidelines
- [ ] **PARTIAL** - Some recommendations deviate from guidelines (explain below)
- [ ] **NOT_ALIGNED** - Significant deviations from NCCN guidelines (explain below)

**Deviations from NCCN Guidelines (if any):**
```
[List specific guideline deviations and justification]

Example:
- Recommendation for off-label alpelisib use: Justified by PIK3CA E545K mutation and platinum resistance despite Category 2B NCCN rating
```

### Institutional Protocol Alignment
Review treatment recommendations against institutional formulary and clinical protocols.

- [ ] **ALIGNED** - All recommendations follow institutional protocols
- [ ] **PARTIAL** - Some recommendations deviate from protocols (explain below)
- [ ] **NOT_ALIGNED** - Significant deviations from institutional protocols (explain below)

**Deviations from Institutional Protocols (if any):**
```
[List specific protocol deviations and justification]

Example:
- Bevacizumab continuation: Not typically used in platinum-resistant setting per institutional protocol, but justified by absence of GI perforation risk factors
```

---

## SECTION 4: QUALITY FLAGS ASSESSMENT

Review automated quality flags raised during analysis. Assess whether each flag is acceptable or requires action before approval.

### Flag 1: {{FLAG_CHECK_NAME_1}}
- **Severity:** {{FLAG_SEVERITY_1}} (critical / warning / info)
- **Message:** {{FLAG_MESSAGE_1}}
- **Recommendation:** {{FLAG_RECOMMENDATION_1}}
- **Reviewer Assessment:**
  - [ ] **ACCEPTABLE** - Flag noted but not concerning, does not affect clinical validity
  - [ ] **REQUIRES_ACTION** - Flag must be addressed before approval
- **Comments:**
  ```
  [Explain your assessment of this quality flag]

  Example:
  Minimum region size of 35 spots is acceptable. While below ideal threshold of 50, the strong effect sizes (log2FC > 3) and very low FDR values (< 1e-10) provide sufficient confidence.
  ```

### Flag 2: {{FLAG_CHECK_NAME_2}}
- **Severity:** {{FLAG_SEVERITY_2}}
- **Message:** {{FLAG_MESSAGE_2}}
- **Recommendation:** {{FLAG_RECOMMENDATION_2}}
- **Reviewer Assessment:**
  - [ ] ACCEPTABLE
  - [ ] REQUIRES_ACTION
- **Comments:**
  ```
  [Your assessment here]
  ```

### Flag 3: {{FLAG_CHECK_NAME_3}}
- **Severity:** {{FLAG_SEVERITY_3}}
- **Message:** {{FLAG_MESSAGE_3}}
- **Recommendation:** {{FLAG_RECOMMENDATION_3}}
- **Reviewer Assessment:**
  - [ ] ACCEPTABLE
  - [ ] REQUIRES_ACTION
- **Comments:**
  ```
  [Your assessment here]
  ```

**If no quality flags were raised:**
☑️ No quality flags - all automated checks passed

---

## SECTION 5: TREATMENT RECOMMENDATIONS REVIEW

Review each treatment recommendation for appropriateness given patient context, molecular profile, and clinical guidelines.

### Recommendation 1: {{THERAPY_NAME_1}}
- **Recommendation ID:** {{REC_ID_1}}
- **Molecular Target:** {{MOLECULAR_TARGET_1}}
- **Expected Efficacy:** {{EXPECTED_EFFICACY_1}}
- **Agreement:**
  - [ ] **AGREE** - Recommendation is appropriate for this patient
  - [ ] **DISAGREE** - Recommendation is not appropriate (explain below)
- **Comments:**
  ```
  [Your assessment of this recommendation]

  Example:
  AGREE. PI3K inhibitor + PARP inhibitor combination is appropriate given PIK3CA E545K mutation and BRCA1 germline status. Patient meets eligibility criteria for clinical trial NCT12345678.
  ```

### Recommendation 2: {{THERAPY_NAME_2}}
- **Recommendation ID:** {{REC_ID_2}}
- **Molecular Target:** {{MOLECULAR_TARGET_2}}
- **Expected Efficacy:** {{EXPECTED_EFFICACY_2}}
- **Agreement:**
  - [ ] AGREE
  - [ ] DISAGREE
- **Comments:**
  ```
  [Your assessment here]
  ```

### Recommendation 3: {{THERAPY_NAME_3}}
- **Recommendation ID:** {{REC_ID_3}}
- **Molecular Target:** {{MOLECULAR_TARGET_3}}
- **Expected Efficacy:** {{EXPECTED_EFFICACY_3}}
- **Agreement:**
  - [ ] AGREE
  - [ ] DISAGREE
- **Comments:**
  ```
  [Your assessment here]
  ```

### Recommendation 4: {{THERAPY_NAME_4}}
- **Recommendation ID:** {{REC_ID_4}}
- **Molecular Target:** {{MOLECULAR_TARGET_4}}
- **Expected Efficacy:** {{EXPECTED_EFFICACY_4}}
- **Agreement:**
  - [ ] AGREE
  - [ ] DISAGREE
- **Comments:**
  ```
  [Your assessment here]
  ```

---

## SECTION 6: ATTESTATION & SIGNATURE ⭐ REQUIRED

By completing this review, you attest to the following:

- [ ] **I have reviewed all findings in this report** (10 molecular findings validated)
- [ ] **I have assessed clinical guideline compliance** (NCCN + institutional protocols)
- [ ] **I have evaluated quality flags** (all automated QC flags assessed)
- [ ] **My decision reflects my clinical judgment** based on expertise and patient context
- [ ] **I understand this review is part of the medical record** (10-year retention, HIPAA-compliant)

**Reviewer Information:**

**Full Name:** {{REVIEWER_NAME}}
**Credentials:** {{CREDENTIALS}} (e.g., MD, Gynecologic Oncology)
**Email:** {{EMAIL}}
**Role:** {{ROLE}} (oncologist / pathologist / radiologist / bioinformatician)
**Date:** {{DATE}} (YYYY-MM-DD)
**Time:** {{TIME}} (HH:MM timezone)

**Digital Signature:**
*[Auto-generated SHA-256 hash by citl_submit_review.py - leave blank]*

---

## SECTION 7: REVISION INSTRUCTIONS (Required if status = REVISE or REJECT)

**Complete this section ONLY if you selected REVISE or REJECT in Section 1.**

### Issues to Address

List specific issues that must be corrected before resubmission:

1. **Issue 1:**
   ```
   [Specific issue description]

   Example:
   Exclude tumor_necrotic region (35 spots) due to insufficient sample size. FDR values for DEGs in this region are unreliable.
   ```

2. **Issue 2:**
   ```
   [Specific issue description]
   ```

3. **Issue 3:**
   ```
   [Specific issue description]
   ```

### Re-analysis Parameters

Specify parameter adjustments for re-running the analysis:

**Regions to exclude:**
```
[List region IDs to exclude from analysis]

Example:
- tumor_necrotic (n=35 spots)
- stroma_artifact (n=28 spots)
```

**FDR threshold adjustment:**
```
[Specify new FDR threshold if needed]

Example:
- Change from FDR < 0.05 to FDR < 0.01 for more stringent criteria
```

**Minimum spots per region:**
```
[Specify minimum sample size threshold]

Example:
- Increase from 30 to 50 spots minimum per region
```

**Additional tests required:**
```
[List any additional analyses needed]

Example:
- Validate TP53 mutation with orthogonal sequencing method (Sanger)
- Request immunohistochemistry confirmation for Ki67 proliferation rate
```

### Resubmission Timeline

**Expected resubmission date:** {{RESUBMISSION_DATE}} (YYYY-MM-DD)

**Estimated time for corrections:** {{ESTIMATED_TIME}} (hours/days/weeks)

---

## NOTES FOR REVIEWERS

**This form structure:**
- Sections 1, 2, 3, 6: **REQUIRED** for all reviews
- Section 4: Required if quality flags exist, otherwise skip
- Section 5: Required if treatment recommendations exist
- Section 7: **REQUIRED** if status is REVISE or REJECT

**Converting to JSON:**
1. Complete all required sections in this markdown template
2. Use provided conversion script: `python scripts/convert_review_to_json.py`
3. OR manually structure as JSON following `shared/schemas/citl_review_schema.json`
4. Submit using: `python tools/reports/citl_submit_review.py --patient-id {{PATIENT_ID}} --review-file <path>`

**Typical completion time by status:**
- **APPROVE:** 20-25 minutes (straightforward validation)
- **REVISE:** 25-35 minutes (requires detailed revision instructions)
- **REJECT:** 30-40 minutes (requires detailed error documentation and escalation plan)

**Questions or issues?**
- Consult the [CitL Workflow Guide](CITL_WORKFLOW_GUIDE.md) for detailed instructions
- Contact bioinformatics team for technical questions
- Contact clinic administrator for workflow questions

---

**Document Version:** 1.0
**Last Updated:** 2026-01-13
**Part of:** Precision Medicine MCP - Clinician-in-the-Loop Validation Workflow
