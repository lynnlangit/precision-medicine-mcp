# Bias Audit Checklist

**Purpose:** Practical checklist for conducting algorithmic bias audits in precision medicine workflows

**Target Audience:** Clinicians, Researchers, Quality Assurance Teams

**Time Required:** 3-5 days per audit

**Last Updated:** 2026-01-12

---

## Quick Start

This checklist guides you through three phases of bias auditing:
1. **Pre-Analysis** (Day 1): Prepare data and documentation
2. **During Analysis** (Days 2-3): Run audit script and analyze results
3. **Post-Analysis** (Days 4-5): Document findings and implement mitigations

**Prerequisites:**
- Access to de-identified patient data (≥100 patients per ancestry group)
- Python environment with audit script installed
- Reference dataset documentation
- Clinical domain expertise

---

## Pre-Analysis Checklist

### Data Preparation

- [ ] **Collect de-identified patient data**
  - Genomic data (VCF files, gene expression matrices)
  - Clinical data (FHIR resources, de-identified)
  - Spatial transcriptomics (expression matrices, coordinates)
  - Demographics (self-reported ancestry, NOT race-as-biology)

- [ ] **Verify sample size**
  - Minimum 100 patients per ancestry group for statistical power
  - If <100 for any group, note as limitation in report

- [ ] **Check de-identification compliance**
  - All 18 HIPAA identifiers removed
  - No dates more specific than year
  - Geographic subdivisions larger than state
  - Use HIPAA-compliant de-identification tool if needed

---

### Reference Dataset Documentation

- [ ] **Inventory all reference datasets used**
  - Genomics: ClinVar, COSMIC, gnomAD, GTEx, etc.
  - Spatial: SingleR, Human Cell Atlas, custom references
  - Clinical: Local EHR norms, MIMIC-IV, etc.

- [ ] **Document ancestry distribution for each dataset**
  - Create table: Dataset | European % | African % | Latino % | Asian % | Other %
  - Flag datasets with <10% representation for any major group
  - Example:
    ```
    Dataset    | European | African | Latino | Asian | Other | Flag
    -----------|----------|---------|--------|-------|-------|------
    ClinVar    | 70%      | 10%     | 10%    | 5%    | 5%    | ⚠️
    gnomAD     | 43%      | 21%     | 14%    | 9%    | 8%    | ✅
    GTEx       | 85%      | 10%     | 5%     | 0%    | 0%    | ⚠️⚠️
    All of Us  | 40%      | 25%     | 20%    | 8%    | 7%    | ✅
    ```

- [ ] **Check reference dataset versions**
  - gnomAD: v4.1 (latest as of Jan 2026)
  - ClinVar: Monthly updates (use current month)
  - GTEx: v8 (check if v9 released)
  - All of Us: Quarterly releases

- [ ] **Identify datasets requiring updates**
  - If reference is >1 year old, check for updates
  - If newer diverse dataset available (e.g., All of Us), consider switching

---

### Ethical Review

- [ ] **Confirm appropriate use of ancestry data**
  - ✅ Appropriate: Ancestry for genomic variant interpretation
  - ❌ Inappropriate: Race as biological proxy (e.g., kidney function)
  - ✅ Appropriate: Socioeconomic ZIP code for resource access
  - ❌ Inappropriate: Insurance status influencing clinical recommendations

- [ ] **Verify IRB approval (if applicable)**
  - Bias audit on existing de-identified data: May not require new IRB
  - Check institutional policy
  - Document IRB determination

- [ ] **Confirm no PHI in analysis**
  - Audit uses aggregate statistics only
  - No individual patient re-identification
  - Stratification by demographic groups is permitted under HIPAA

---

### Workflow Documentation

- [ ] **Document current workflow**
  - Which tools are used (mcp-spatialtools, mcp-multiomics, etc.)
  - Which reference datasets feed each tool
  - Which features/variables are inputs to algorithms
  - Where ancestry/demographic data is used

- [ ] **Identify potential bias entry points**
  - Data collection: Are certain populations undersampled?
  - Feature selection: Are demographic proxies used (ZIP code, insurance)?
  - Model training: Was training data diverse?
  - Clinical interpretation: Are clinicians aware of limitations?

- [ ] **Review past bias incidents (if any)**
  - Have clinicians reported suspected bias?
  - Were certain ancestries disproportionately affected?
  - What mitigations were implemented?

---

### Tool Setup

- [ ] **Install Python dependencies**
  ```bash
  cd /path/to/spatial-mcp
  pip install -r requirements.txt
  ```

- [ ] **Verify audit script is runnable**
  ```bash
  python infrastructure/audit/audit_bias.py --help
  ```

- [ ] **Test on sample data**
  ```bash
  python infrastructure/audit/audit_bias.py \
    --workflow patientone \
    --genomics-data data/genomics/sample.vcf \
    --output /tmp/test_audit.html
  ```

- [ ] **Prepare output directory**
  ```bash
  mkdir -p reports/bias_audits/$(date +%Y)
  ```

---

### Stakeholder Communication

- [ ] **Notify clinical team**
  - Inform clinicians that bias audit is being conducted
  - Request input on suspected biases or edge cases
  - Schedule review meeting for audit results

- [ ] **Notify compliance team**
  - Hospital quality assurance
  - Privacy officer (confirm de-identification)
  - Legal (if needed for regulatory submission)

- [ ] **Set expectations**
  - Audit will take 3-5 days
  - Results may require workflow changes
  - Mitigations will be implemented before production deployment

---

## During Analysis Checklist

### Data Representation Analysis

- [ ] **Calculate ancestry distribution of reference datasets**
  - Use documented percentages from Pre-Analysis phase
  - Compare to US population baseline:
    * European: 60%
    * Hispanic/Latino: 18%
    * African/Black: 13%
    * Asian: 6%
    * Native American: 1%
    * Multiracial: 2%

- [ ] **Flag underrepresented groups**
  - <10% representation = HIGH RISK
  - <5% representation = CRITICAL RISK
  - 0% representation = UNACCEPTABLE

- [ ] **Document findings**
  - Create "Reference Dataset Diversity Report" table
  - For each dataset, note risk level
  - Prioritize datasets for replacement (critical risk first)

---

### Algorithm Fairness Testing

- [ ] **Stratify test data by ancestry**
  - European/White
  - African/Black/African American
  - Hispanic/Latino
  - Asian (specify: East Asian, South Asian, etc.)
  - Native American
  - Multiracial/Other

- [ ] **Calculate performance metrics per group**
  - Accuracy, sensitivity, specificity
  - Positive predictive value (PPV), negative predictive value (NPV)
  - For continuous outcomes: RMSE, MAE
  - For rankings: AUC-ROC

- [ ] **Test for statistical significance**
  - Use chi-square test for categorical outcomes
  - Use t-test for continuous outcomes
  - Flag disparities >10% as requiring review
  - Flag disparities >20% as critical

- [ ] **Calculate fairness metrics**
  - **Demographic parity:** P(ŷ=1|A=a) equal across groups?
  - **Equalized odds:** Equal TPR and FPR across groups?
  - **Calibration:** Predicted probabilities match observed rates?
  - **Clinical utility parity:** Equal benefit (QALYs, lives saved)?

- [ ] **Document fairness violations**
  - Which metric is violated?
  - Which groups are affected?
  - What is the magnitude of disparity?
  - Example:
    ```
    Metric: Equalized Odds (Sensitivity)
    European: 92%
    African: 80%
    Disparity: 12 percentage points ⚠️
    Conclusion: Requires mitigation
    ```

---

### Output Stratification Analysis

- [ ] **Stratify treatment recommendations by ancestry**
  - Count aggressive vs. conservative treatments by group
  - Run chi-square test for independence
  - Flag if p < 0.05 (statistically significant disparity)

- [ ] **Analyze confidence scores by group**
  - Calculate mean confidence per ancestry
  - Flag if >10% difference between groups
  - Lower confidence for underrepresented groups is EXPECTED (indicates uncertainty)
  - Example:
    ```
    Group      | Mean Confidence | Interpretation
    -----------|-----------------|----------------
    European   | 0.92            | High confidence
    African    | 0.78            | Medium confidence ⚠️ (14% lower due to limited reference data)
    Asian      | 0.70            | Medium confidence ⚠️⚠️ (22% lower)
    ```

- [ ] **Check for demographic proxy influence**
  - Identify proxy features: ZIP code, insurance type, language preference
  - Check feature importance scores
  - Run model with/without proxies
  - If performance changes >5%, proxies are influential (BAD)
  - Example:
    ```
    Feature          | Importance | Flag
    -----------------|------------|------
    BRCA1 status     | 0.45       | ✅ (clinical)
    CA-125 level     | 0.32       | ✅ (clinical)
    ZIP code         | 0.08       | ⚠️⚠️ (proxy - should be 0)
    Insurance type   | 0.05       | ⚠️ (proxy - should be 0)
    ```

- [ ] **Stratify by socioeconomic indicators (if appropriate)**
  - Medicaid vs. commercial insurance (should NOT affect clinical recommendations)
  - Urban vs. rural (should only affect resource access, NOT treatment)
  - Education level (NOT collected in most clinical workflows)

---

### Run Audit Script

- [ ] **Execute bias audit script**
  ```bash
  python infrastructure/audit/audit_bias.py \
    --workflow patientone \
    --genomics-data data/genomics/patient_variants.vcf \
    --clinical-data data/fhir/patients_deidentified.json \
    --spatial-data data/spatial/spatial_expression.csv \
    --demographics data/demographics/ancestry_labels.csv \
    --output reports/bias_audit_$(date +%Y-%m-%d).html \
    --min-representation 0.10 \
    --max-disparity 0.10
  ```

- [ ] **Verify script completes without errors**
  - Check for Python exceptions
  - Verify all input files loaded correctly
  - Confirm output files generated

- [ ] **Review automated findings**
  - Open HTML report in browser
  - Read executive summary (Pass/Fail)
  - Review flagged issues
  - Check recommendations

---

### Manual Validation

- [ ] **Spot-check algorithmic results**
  - Manually review 10-20 cases per ancestry group
  - Do predictions make clinical sense?
  - Are confidence scores appropriate?
  - Are disclaimers present for underrepresented groups?

- [ ] **Validate stratified metrics**
  - Re-calculate fairness metrics manually for 1 group
  - Verify automated script is correct
  - Cross-check with clinical domain knowledge

- [ ] **Review edge cases**
  - Patients with rare variants
  - Patients from underrepresented ancestries
  - Patients with conflicting data (e.g., discordant ClinVar classifications)

- [ ] **Check for data quality issues**
  - Are there missing values?
  - Are ancestry labels accurate (self-reported preferred)?
  - Are there sample size imbalances?

---

### Generate Visualizations

- [ ] **Create fairness metrics bar chart**
  - X-axis: Ancestry groups
  - Y-axis: Sensitivity (or other metric)
  - Highlight disparities >10%

- [ ] **Create reference dataset diversity heatmap**
  - Rows: Datasets (ClinVar, gnomAD, GTEx, etc.)
  - Columns: Ancestries (European, African, Latino, Asian, Other)
  - Color: Green (≥20%), Yellow (10-20%), Red (<10%)

- [ ] **Create confidence score distribution plot**
  - Separate distributions per ancestry group
  - Overlay to show differences

- [ ] **Include visualizations in HTML report**

---

## Post-Analysis Checklist

### Document Findings

- [ ] **Write executive summary**
  - Overall assessment: PASS / CONDITIONAL PASS / FAIL
  - Risk level: LOW / MEDIUM / HIGH / CRITICAL
  - Number of findings: X warnings, Y critical issues
  - Recommended action: Deploy / Deploy with mitigations / Do not deploy

- [ ] **Categorize findings by severity**
  - **Critical:** >20% disparity, 0% representation, proxy features used
  - **High:** 10-20% disparity, <5% representation
  - **Medium:** 5-10% disparity, <10% representation
  - **Low:** <5% disparity, disclaimers missing

- [ ] **Provide specific examples**
  - "African patients have 12% lower sensitivity (80% vs. 92% for European)"
  - "GTEx has 0% Asian representation → flagged in all reports"
  - "ZIP code has 8% feature importance → removed from model"

- [ ] **Map findings to mitigation strategies**
  - For each finding, specify mitigation (see Mitigation section below)

---

### Risk Assessment

- [ ] **Evaluate overall bias risk**
  - Genomics risk: LOW / MEDIUM / HIGH
  - Gene expression risk: LOW / MEDIUM / HIGH
  - Clinical risk: LOW / MEDIUM / HIGH
  - Spatial transcriptomics risk: LOW / MEDIUM / HIGH
  - Multiomics risk: LOW / MEDIUM / HIGH

- [ ] **Determine deployment readiness**
  - If all risks LOW → Deploy
  - If any risk MEDIUM → Deploy with mitigations + quarterly monitoring
  - If any risk HIGH → Implement mitigations first, re-audit
  - If any risk CRITICAL → Do not deploy until resolved

- [ ] **Estimate mitigation effort**
  - Days to implement mitigations
  - Resources required (personnel, compute, data)
  - Timeline for re-audit

- [ ] **Document limitations**
  - Known biases that cannot be fully mitigated
  - Reference datasets that cannot be replaced (no diverse alternative)
  - Acknowledge in reports and patient communications

---

### Implement Mitigations

**For Underrepresented Reference Data:**

- [ ] **Add diverse reference datasets**
  - Replace ClinVar → ClinVar + gnomAD + All of Us
  - Replace GTEx → GTEx + TOPMed + Human Cell Atlas
  - Document in workflow

- [ ] **Add ancestry caveats to reports**
  - Example: "GTEx has 0% Asian representation. Results may not accurately reflect normal variation in Asian ancestries."
  - Display prominently in clinician reports

- [ ] **Adjust thresholds for underrepresented groups**
  - Increase fold-change cutoff (log2FC > 2.0 instead of 1.5)
  - Widen confidence intervals
  - Require manual review

---

**For Algorithmic Fairness Violations:**

- [ ] **Retrain with fairness constraints**
  - Use fairlearn library for fairness-aware training
  - Constrain demographic parity or equalized odds
  - Validate improved fairness metrics

- [ ] **Post-hoc calibration by group**
  - Train separate calibration models per ancestry
  - Validate calibration curves

- [ ] **Adjust decision thresholds by group**
  - If sensitivity is lower for African patients, lower threshold for that group
  - Document rationale

---

**For Demographic Proxy Features:**

- [ ] **Remove proxy features**
  - Drop ZIP code, insurance type, language preference from clinical models
  - Retain only for resource allocation (not clinical decisions)

- [ ] **Retrain model without proxies**
  - Verify performance does not degrade significantly (<5% accuracy loss acceptable)
  - If performance degrades >5%, investigate why (may indicate confounding)

- [ ] **Add audit log for proxy usage**
  - If proxies must be used (e.g., resource allocation), log every use
  - Enable review of inappropriate proxy usage

---

**For Confidence Score Disparities:**

- [ ] **Implement ancestry-aware confidence scoring**
  - Reduce confidence by X% if patient ancestry has <10% representation
  - Flag for manual review if confidence <0.70

- [ ] **Add uncertainty quantification**
  - Display confidence intervals, not just point estimates
  - Example: "Risk: 65% (95% CI: 55%-75%)"

---

### Clinical Validation

- [ ] **Schedule review meeting with clinical team**
  - Present audit findings
  - Discuss clinical significance of disparities
  - Get input on mitigations

- [ ] **Obtain sign-off from experts**
  - Medical geneticist (genomics findings)
  - Oncologist (treatment recommendations)
  - Bioinformatician (technical validation)
  - Bioethicist (fairness assessment)

- [ ] **Document expert recommendations**
  - Do experts agree with audit findings?
  - Are mitigations sufficient?
  - Are there additional concerns?

- [ ] **Decide on deployment**
  - Go / No-Go decision
  - If Go: Implement mitigations first or concurrent monitoring?
  - If No-Go: Timeline for re-audit

---

### Update Documentation

- [ ] **Update HIPAA compliance docs**
  - Add reference to bias audit in §8 (new section)
  - Cross-link to ETHICS_AND_BIAS.md

- [ ] **Update Operations Manual**
  - Add bias audit to monthly compliance checklist
  - Include in incident response procedures

- [ ] **Update Admin Guide**
  - Add "Bias Audit" section with instructions
  - Document how to run audit script

- [ ] **Update workflow documentation**
  - Note which diverse reference datasets are used
  - Document ancestry caveats
  - Update tool descriptions with bias mitigations

---

### Archive Audit Report

- [ ] **Store audit report**
  ```bash
  mkdir -p reports/bias_audits/$(date +%Y)
  mv reports/bias_audit_$(date +%Y-%m-%d).html \
     reports/bias_audits/$(date +%Y)/$(date +%Y-%m-%d)_initial_audit/
  ```

- [ ] **Generate compliance summary**
  - One-page PDF for compliance officers
  - Key findings, risk level, mitigations

- [ ] **Set retention policy**
  - Keep for 10 years (HIPAA alignment)
  - Include in annual compliance reviews

- [ ] **Update audit log**
  ```json
  {
    "audit_date": "2026-01-12",
    "auditor": "Jane Doe, PhD",
    "workflow": "PatientOne Complete Analysis",
    "risk_level": "MEDIUM",
    "findings": 2,
    "mitigations_implemented": 3,
    "next_audit": "2026-04-12"
  }
  ```

---

### Schedule Next Audit

- [ ] **Determine next audit date**
  - Quarterly audits: 3 months from now
  - Triggered audits: After workflow changes
  - Annual comprehensive audit: 12 months from now

- [ ] **Set calendar reminders**
  - 2 weeks before: Prepare data
  - 1 week before: Review checklist
  - Audit date: Run script and analyze

- [ ] **Assign responsibilities**
  - Who will run the audit?
  - Who will review results?
  - Who will implement mitigations?

- [ ] **Update audit schedule**
  ```markdown
  ## Bias Audit Schedule
  - 2026-01-12: Initial Audit (COMPLETE)
  - 2026-04-12: Quarterly Audit Q2 (SCHEDULED)
  - 2026-07-12: Quarterly Audit Q3 (SCHEDULED)
  - 2026-10-12: Quarterly Audit Q4 (SCHEDULED)
  - 2027-01-12: Annual Comprehensive Audit (SCHEDULED)
  ```

---

### Communicate Results

- [ ] **Present to clinical team**
  - Summarize findings (non-technical language)
  - Explain mitigations
  - Answer questions

- [ ] **Notify compliance officers**
  - Share audit report
  - Confirm compliance with ethical AI standards

- [ ] **Update stakeholders**
  - Hospital administration
  - Quality assurance
  - Privacy office

- [ ] **External communication (if required)**
  - FDA (if part of pre-market submission)
  - IRB (if research use)
  - Publications (if publishing methodology)

---

## Common Pitfalls to Avoid

### Data Issues

- ❌ **Small sample sizes:** <100 patients per group → Insufficient statistical power
  - ✅ Solution: Pool data across institutions or extend data collection period

- ❌ **Mislabeled ancestry:** Administrative race vs. self-reported ancestry
  - ✅ Solution: Use self-reported ancestry from FHIR Patient.extension (US Core)

- ❌ **Survivorship bias:** Only analyzing patients who completed treatment
  - ✅ Solution: Include all patients, stratify by completion status

---

### Analysis Issues

- ❌ **Cherry-picking metrics:** Only reporting fairness metrics that look good
  - ✅ Solution: Report all standard metrics (demographic parity, equalized odds, calibration)

- ❌ **Ignoring small disparities:** "5% difference isn't clinically meaningful"
  - ✅ Solution: Aggregate small disparities across population can have large impact

- ❌ **Post-hoc rationalization:** Finding excuses for bias instead of mitigating
  - ✅ Solution: Implement mitigations first, then re-assess

---

### Mitigation Issues

- ❌ **Band-aid solutions:** Adding disclaimers without addressing root cause
  - ✅ Solution: Replace biased reference datasets, retrain models

- ❌ **Over-correction:** Enforcing demographic parity when base rates differ
  - ✅ Solution: Use appropriate fairness metric (calibration for risk scores)

- ❌ **Ignoring trade-offs:** Fairness improvements may reduce overall accuracy
  - ✅ Solution: Document trade-offs, get clinical approval

---

## Quick Reference: Thresholds

| Metric | Threshold | Action |
|--------|-----------|--------|
| **Representation** | <10% | HIGH RISK - Find diverse alternative |
| **Representation** | <5% | CRITICAL RISK - Do not use |
| **Fairness Disparity** | >10% | Implement mitigation |
| **Fairness Disparity** | >20% | Critical - Do not deploy |
| **Proxy Importance** | >5% | Remove feature, retrain |
| **Confidence Difference** | >15% | Add ancestry-specific disclaimers |
| **Sample Size** | <100 per group | Extend data collection |

---

## Checklist Summary

**Pre-Analysis (Day 1):**
- ☑ Data collected and de-identified
- ☑ Reference datasets documented
- ☑ Ethical review confirmed
- ☑ Tools set up and tested
- ☑ Stakeholders notified

**During Analysis (Days 2-3):**
- ☑ Data representation analyzed
- ☑ Fairness metrics calculated
- ☑ Output stratification completed
- ☑ Audit script executed
- ☑ Manual validation performed

**Post-Analysis (Days 4-5):**
- ☑ Findings documented
- ☑ Risk assessment completed
- ☑ Mitigations implemented
- ☑ Clinical validation obtained
- ☑ Reports archived and communicated

---

**For detailed technical guidance, see:**
- [ETHICS_AND_BIAS.md](ETHICS_AND_BIAS.md) - Comprehensive framework
- [PATIENTONE_BIAS_AUDIT.md](PATIENTONE_BIAS_AUDIT.md) - Concrete demonstration

**Questions?** File an issue on GitHub or consult your institutional bioethics committee.

---

**Document Status:** ✅ Complete (Week 1 Deliverable)
**Last Updated:** 2026-01-12
**Version:** 1.0
