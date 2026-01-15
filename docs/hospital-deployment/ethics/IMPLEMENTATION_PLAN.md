# Ethics & Bias Framework - Implementation Plan (Revised)

**Priority:** 1 (Colleague Feedback)
**Target Audiences:** Clinicians, Researchers (US Hospitals)
**Estimated Timeline:** 3-4 weeks
**Status:** Planning Phase
**Scope:** US Healthcare Context

---

## Executive Summary

**Goal:** Add dedicated "Ethics & Bias" section to repository demonstrating how to audit precision medicine workflows for algorithmic bias across diverse populations.

**Why:** Trust is the primary barrier to clinical AI adoption. This framework aligns PatientOne demo with US healthcare AI standards for ethical medical AI (WHO, FDA).

**Current Gap:** Repository has strong HIPAA compliance and data privacy coverage, but lacks systematic approach to detecting and mitigating algorithmic bias in precision medicine.

**Focus:** US hospital deployment with manual, batch-based bias audits using recommended diverse reference datasets.

---

## Scope

### What We're Adding

1. **Ethics & Bias Framework** - Comprehensive guide to ethical AI in precision medicine
2. **Bias Audit Methodology** - Step-by-step process for detecting bias
3. **PatientOne Bias Audit** - Concrete demonstration on existing workflow
4. **Bias Detection Tools** - Python utilities for ongoing monitoring
5. **Audit Checklist** - Practical checklist for clinicians/researchers

### What We're Auditing For

**Genomics Bias (Ancestral Populations):**
- BRCA1/BRCA2 variant interpretation: Are databases Euro-centric?
- Gene expression reference ranges: Do they account for ancestral variation?
- Pathway enrichment databases: Are they representative?
- Drug-gene interaction data: Tested in diverse genetic backgrounds?

**Clinical Bias (Socioeconomic Factors):**
- Insurance status: Does it influence treatment recommendations?
- Geographic location: Used as proxy for socioeconomic status?
- Race/ethnicity coding: Appropriate use (ancestry for genomics) vs. inappropriate (proxy for biology)?
- Language barriers: Are non-English speakers considered?

**Spatial Transcriptomics Bias:**
- Reference cell type signatures: From diverse populations?
- Deconvolution algorithms: Work equally well across tissue types?

**Multiomics Bias:**
- PDX models: Representative of patient populations?
- Drug sensitivity data: Diverse genetic backgrounds?

---

## Deliverables

### Phase 1: Core Documentation (Week 1)

#### 1. `docs/hospital-deployment/ethics/ETHICS_AND_BIAS.md` (~800 lines)

**Purpose:** Comprehensive ethical AI framework for precision medicine

**Sections:**
1. **Introduction** (100 lines)
   - Why ethics & bias matter in precision medicine
   - Trust as barrier to adoption
   - Legal/regulatory landscape
   - Relationship to HIPAA compliance

2. **US Healthcare AI Standards & Frameworks** (150 lines)
   - FDA AI/ML-Based Software as Medical Device (SaMD) guidance (2021)
   - WHO Ethics and Governance of AI for Health (2021) - Global best practices
   - NIH All of Us Research Program diversity requirements
   - 21st Century Cures Act - Advancing health equity through technology
   - AMA Code of Medical Ethics Opinion 2.3.2 - Physicians' use of AI
   - HIPAA Privacy Rule and algorithmic fairness considerations

3. **Types of Bias in Precision Medicine** (200 lines)
   - **Data bias**: Representation bias, selection bias, measurement bias
   - **Algorithmic bias**: Model training bias, feature selection bias
   - **Interpretation bias**: Clinical decision bias, confirmation bias
   - Examples specific to:
     * Genomics (variant calling, ancestry inference)
     * Spatial transcriptomics (cell type deconvolution)
     * Multiomics (pathway enrichment, drug target prediction)
     * Clinical data (treatment recommendations)

4. **Bias Audit Methodology** (200 lines)
   - **Step 1:** Data representation analysis
     * Demographic stratification
     * Ancestry distribution in genomic datasets
     * Geographic/socioeconomic coverage
   - **Step 2:** Algorithm fairness testing
     * Fairness metrics (demographic parity, equalized odds, calibration)
     * Performance stratification by subgroup
   - **Step 3:** Output stratification analysis
     * Treatment recommendations by demographic
     * Confidence intervals by subgroup
   - **Step 4:** Clinical validation
     * Expert review across populations
     * Real-world performance monitoring

5. **Fairness Metrics** (100 lines)
   - Demographic parity
   - Equalized odds (equal FPR/TPR across groups)
   - Calibration (predicted probabilities match observed rates)
   - Clinical utility parity (equal benefit across groups)
   - When to use which metric

6. **Transparency & Explainability** (100 lines)
   - Model explainability (SHAP values, feature importance)
   - Decision provenance (audit logs)
   - Uncertainty quantification
   - Disclaimers and limitations

7. **Mitigation Strategies** (100 lines)
   - Data augmentation for underrepresented groups
   - Fairness-aware training (fairness constraints)
   - Post-hoc calibration
   - Human-in-the-loop validation
   - Use diverse reference datasets (see §8)

8. **Recommended Diverse Reference Datasets** (150 lines)

   **Genomics & Variant Interpretation:**
   - **gnomAD (Genome Aggregation Database)**:
     * 76,156 genomes, 125,748 exomes
     * Diverse ancestry: African/African American (21%), Latino/Admixed American (14%), East Asian (9%), South Asian (8%), Other (5%), European (non-Finnish) (43%)
     * Use for: Population allele frequencies, pathogenicity assessment
     * URL: https://gnomad.broadinstitute.org/

   - **All of Us Research Program**:
     * 245,000+ participants (as of 2024)
     * 80% from underrepresented groups in biomedical research
     * Use for: Reference ranges, variant interpretation, phenotype associations
     * URL: https://www.researchallofus.org/

   - **ClinGen (Clinical Genome Resource)**:
     * Diverse variant classifications
     * Expert-curated pathogenicity assertions
     * Use for: Variant interpretation, gene-disease associations
     * URL: https://www.clinicalgenome.org/

   - **TOPMed (Trans-Omics for Precision Medicine)**:
     * 180,000+ whole genomes
     * Diverse US population
     * Use for: Rare variant discovery, ancestry-specific analysis
     * URL: https://www.nhlbiwgs.org/

   - **1000 Genomes Project**:
     * 2,504 genomes from 26 populations
     * Comprehensive ancestry representation
     * Use for: Structural variant analysis, ancestry inference
     * URL: https://www.internationalgenome.org/

   **Gene Expression & Transcriptomics:**
   - **GTEx (Genotype-Tissue Expression)**:
     * 85% European ancestry (limitation - acknowledge in reports)
     * 948 donors, 54 tissue types
     * Use for: Expression reference ranges (with ancestry caveat)
     * URL: https://gtexportal.org/

   - **Human Cell Atlas**:
     * Single-cell references from diverse populations
     * Use for: Cell type identification, deconvolution
     * URL: https://www.humancellatlas.org/

   **Clinical Data:**
   - **MIMIC-IV (Medical Information Mart for Intensive Care)**:
     * De-identified ICU data
     * Diverse patient demographics
     * Use for: Clinical phenotype validation
     * URL: https://mimic.mit.edu/

   **Best Practices:**
   - Use multiple reference datasets to cross-validate
   - Document ancestry distribution of each reference
   - Flag results when patient ancestry not well-represented
   - Update references annually as new diverse datasets emerge

9. **Batch Bias Auditing Approach** (100 lines)

   **Manual Execution Model:**
   - Run bias audits quarterly or after major workflow changes
   - Generate HTML reports for stakeholder review
   - Store audit history for compliance documentation
   - Manual review and approval before workflow updates

   **Audit Schedule:**
   - **Initial Audit:** Before production deployment
   - **Quarterly Audits:** Every 3 months during production
   - **Triggered Audits:** After workflow changes, new data sources, or bias incidents

   **Workflow:**
   ```bash
   # Step 1: Collect workflow data (genomics, clinical, spatial)
   # Step 2: Run bias audit script
   python infrastructure/audit/audit_bias.py --workflow patientone --output reports/

   # Step 3: Review HTML report
   open reports/bias_audit_2026-01-12.html

   # Step 4: Document findings and actions
   # Step 5: Implement mitigations if needed
   # Step 6: Archive report for compliance
   ```

   **Report Retention:**
   - Keep all audit reports for 10 years (HIPAA compliance alignment)
   - Include in annual compliance reviews
   - Make available for external audits

---

#### 2. `docs/hospital-deployment/ethics/BIAS_AUDIT_CHECKLIST.md` (~300 lines)

**Purpose:** Practical checklist for conducting bias audits

**Sections:**

1. **Pre-Analysis Checklist** (100 lines)
   - [ ] Dataset demographics documented
   - [ ] Ancestral population representation >10% for each major group
   - [ ] Socioeconomic factors recorded (if applicable)
   - [ ] Reference databases reviewed for diversity
   - [ ] Model training data stratified by subgroup
   - [ ] Known limitations documented
   - [ ] IRB approval for bias analysis (if applicable)

2. **During Analysis Checklist** (100 lines)
   - [ ] Stratify results by ancestry/ethnicity
   - [ ] Stratify results by socioeconomic indicators
   - [ ] Check for performance disparities >10% between groups
   - [ ] Validate reference ranges across populations
   - [ ] Review feature importance for demographic proxies
   - [ ] Test edge cases in underrepresented groups
   - [ ] Document all stratified analyses

3. **Post-Analysis Checklist** (100 lines)
   - [ ] Bias metrics calculated and documented
   - [ ] Disparities >10% flagged for review
   - [ ] Mitigation strategies implemented (if bias found)
   - [ ] Results validated by domain experts
   - [ ] Audit report generated with findings
   - [ ] Recommendations incorporated into workflow
   - [ ] Continuous monitoring plan established

---

#### 3. `docs/hospital-deployment/ethics/PATIENTONE_BIAS_AUDIT.md` (~500 lines)

**Purpose:** Concrete demonstration of bias audit on PatientOne workflow

**Sections:**

1. **Overview** (50 lines)
   - PatientOne profile: 63-year-old woman, Stage IV HGSOC, platinum-resistant
   - Workflow scope: Genomics, spatial transcriptomics, multiomics, clinical data
   - Audit objectives: Identify potential biases, document findings, propose mitigations

2. **Genomics Bias Analysis** (150 lines)

   **2.1 BRCA1/BRCA2 Variant Interpretation**
   - **Current Data Source:** ClinVar, COSMIC databases
   - **Bias Check:** Are pathogenic variants primarily from European studies?
   - **PatientOne Finding:** BRCA1 variant (c.5266dupC) well-studied in European populations, limited data in others
   - **Potential Impact:** Pathogenicity classification may be less certain in non-European ancestries
   - **Recommended Diverse References:**
     * **gnomAD:** Check population-specific allele frequencies (21% African/African American, 14% Latino, 9% East Asian)
     * **ClinGen:** Expert-curated variant classifications with ancestry context
     * **All of Us:** Cross-reference with 80% underrepresented groups data
   - **Mitigation:**
     * Flag variants with <5 studies in patient's ancestry
     * Recommend genetic counseling for ancestry-specific interpretation
     * Always check gnomAD for population-specific frequencies

   **2.2 Gene Expression Reference Ranges**
   - **Current Data Source:** GTEx (Genotype-Tissue Expression project)
   - **Bias Check:** GTEx is 85% European ancestry, 10% African, 5% other
   - **PatientOne Finding:** Differential expression analysis uses GTEx normal tissue baseline
   - **Potential Impact:** Reference ranges may not reflect true variation in underrepresented ancestries
   - **Recommended Diverse References:**
     * **GTEx:** Use with explicit ancestry caveat in reports
     * **TOPMed:** 180,000+ genomes from diverse US population for validation
     * **Human Cell Atlas:** Single-cell references from diverse populations
   - **Mitigation:**
     * Document GTEx ancestry distribution in all reports (85% European)
     * Cross-validate with TOPMed data when patient ancestry differs
     * Apply larger thresholds (log2FC >2 instead of >1.5) for conservatism

   **2.3 Pathway Enrichment Databases**
   - **Current Data Source:** KEGG, Reactome, GO databases
   - **Bias Check:** Are pathway definitions universal or population-specific?
   - **PatientOne Finding:** Pathways based on aggregate human data, not population-stratified
   - **Potential Impact:** Pathway relevance may vary by ancestry
   - **Recommended Diverse References:**
     * **KEGG + Reactome + GO:** Cross-validate across all three databases
     * **All of Us Workbench:** Check pathway associations in diverse cohorts
   - **Mitigation:**
     * Cross-validate with multiple databases (KEGG + Reactome + GO)
     * Flag pathways with >30% of genes showing ancestry-specific expression
     * Reference All of Us data for pathway validation

3. **Clinical Bias Analysis** (150 lines)

   **3.1 Insurance Status & Treatment Recommendations**
   - **Data Source:** FHIR CoverageResource
   - **Bias Check:** Does insurance type affect treatment recommendations?
   - **PatientOne Finding:** Treatment recommendations based ONLY on molecular/clinical data (BRCA status, platinum resistance, tumor markers)
   - **Verification:** Removed `Coverage` resource from analysis prompts
   - **Result:** ✅ No insurance bias detected

   **3.2 Geographic & Socioeconomic Factors**
   - **Data Source:** FHIR Patient address
   - **Bias Check:** Is zip code used as proxy for treatment eligibility?
   - **PatientOne Finding:** Address used ONLY for healthcare provider coordination
   - **Verification:** Reviewed mcp-epic tool calls - no zip code filtering
   - **Result:** ✅ No geographic bias detected

   **3.3 Race/Ethnicity Coding**
   - **Data Source:** FHIR Patient.extension (US Core Race/Ethnicity)
   - **Bias Check:** Is race used appropriately?
   - **PatientOne Finding:**
     * Race/ethnicity recorded as "White/European ancestry"
     * Used for genomic variant interpretation (appropriate)
     * NOT used for treatment eligibility (appropriate)
   - **Best Practice:** Use ancestry for genomics, NOT race as biology
   - **Result:** ✅ Appropriate use confirmed

4. **Spatial Transcriptomics Bias Analysis** (100 lines)

   **4.1 Cell Type Reference Signatures**
   - **Data Source:** SingleR reference datasets
   - **Bias Check:** Are reference cell types from diverse tissues?
   - **PatientOne Finding:** Using generic immune cell references (pan-tissue)
   - **Potential Impact:** May miss ovarian cancer-specific cell states
   - **Mitigation:**
     * Use ovarian-specific reference when available (GSE146026)
     * Validate deconvolution with H&E pathology review

   **4.2 Spatial Autocorrelation Algorithms**
   - **Data Source:** Moran's I implementation
   - **Bias Check:** Does algorithm perform equally across tissue types?
   - **PatientOne Finding:** Moran's I is tissue-agnostic statistical method
   - **Result:** ✅ No tissue-specific bias

5. **Multiomics Bias Analysis** (50 lines)

   **5.1 PDX Model Representativeness**
   - **Data Source:** Xenograft PDX models
   - **Bias Check:** Are PDX models representative?
   - **PatientOne Finding:** PDX models are from diverse ovarian cancer patients
   - **Limitation:** PDX models may not capture immune microenvironment
   - **Mitigation:** Combine with patient's own spatial transcriptomics

6. **Summary & Recommendations** (50 lines)

   **Biases Detected:**
   1. BRCA variant databases: Euro-centric (MEDIUM risk)
   2. GTEx reference ranges: 85% European (MEDIUM risk)
   3. Cell type references: Generic, not cancer-specific (LOW risk)

   **Biases NOT Detected:**
   1. Insurance status: Not used in recommendations ✅
   2. Geographic location: Not used in treatment decisions ✅
   3. Race/ethnicity: Appropriately used for genomics only ✅

   **Recommendations:**
   1. Add ancestry diversity warnings to genomic reports
   2. Flag variants with limited non-European data
   3. Use ancestry-matched reference data when available
   4. Document reference database ancestry distributions
   5. Implement continuous monitoring for bias drift

---

### Phase 2: Bias Detection Tools (Week 2)

#### 4. `shared/utils/bias_detection.py` (~400 lines)

**Purpose:** Reusable Python utilities for bias detection

**Functions:**

```python
def check_dataset_representation(
    df: pd.DataFrame,
    demographic_col: str,
    min_representation: float = 0.10
) -> Dict[str, Any]:
    """
    Check if each demographic group has minimum representation.

    Args:
        df: DataFrame with data
        demographic_col: Column with demographic labels
        min_representation: Minimum fraction required (default 10%)

    Returns:
        {
            "underrepresented_groups": [...],
            "representation": {"group1": 0.45, "group2": 0.08, ...},
            "meets_threshold": True/False
        }
    """

def calculate_fairness_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    groups: np.ndarray
) -> Dict[str, Dict[str, float]]:
    """
    Calculate fairness metrics stratified by group.

    Returns demographic parity, equalized odds, calibration.
    """

def flag_demographic_proxy_features(
    feature_importance: Dict[str, float],
    proxy_features: List[str]
) -> List[str]:
    """
    Flag if demographic proxy features (zip code, language, etc.)
    have high feature importance.
    """

def audit_reference_database_diversity(
    database_metadata: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Audit ancestry diversity in reference databases
    (e.g., ClinVar, GTEx, gnomAD).
    """

def stratify_results_by_group(
    results: pd.DataFrame,
    group_col: str,
    metrics: List[str]
) -> pd.DataFrame:
    """
    Stratify analysis results by demographic group.
    Calculate performance metrics per group.
    """
```

**Tests:** `tests/test_bias_detection.py` (~200 lines)

---

#### 5. `infrastructure/audit/audit_bias.py` (~300 lines)

**Purpose:** Standalone script to run bias audit on workflows

**Usage:**
```bash
# Audit PatientOne workflow
python infrastructure/audit/audit_bias.py \
  --workflow patientone \
  --genomics-data data/genomics/patient001.vcf \
  --clinical-data data/fhir/patient001.json \
  --output reports/bias_audit_patient001.html

# Audit with custom thresholds
python infrastructure/audit/audit_bias.py \
  --workflow patientone \
  --min-representation 0.15 \
  --max-disparity 0.10 \
  --output reports/bias_audit_custom.html
```

**Output:** HTML report with:
- Demographic representation summary
- Fairness metrics by group
- Flagged disparities
- Mitigation recommendations
- Visualizations (bar charts, heatmaps)

---

### Phase 3: Integration & Cross-References (Week 3)

#### 6. Update Existing Documentation

**6.1 Update `docs/hospital-deployment/../compliance/hipaa.md`**
- Add new section: "§8. Ethical AI & Bias Mitigation"
- Cross-reference to `docs/hospital-deployment/ethics/ETHICS_AND_BIAS.md`
- Position ethics as complementary to privacy compliance

**6.2 Update `docs/hospital-deployment/OPERATIONS_MANUAL.md`**
- Add bias audit to monthly compliance checklist
- Include bias monitoring in incident response procedures

**6.3 Update `docs/hospital-deployment/ADMIN_GUIDE.md`**
- Add "Bias Audit Dashboard" to monitoring section
- Document how to run `audit_bias.py` script

**6.4 Update `docs/EXECUTIVE_SUMMARY.md`**
- Add "Ethical AI & Bias Mitigation" to key features
- Highlight alignment with WHO/FDA/EU AI Act standards
- Emphasize trust-building for clinical adoption

**6.5 Update `test-docs/patient-one-scenario/README.md`**
- Add reference to PATIENTONE_BIAS_AUDIT.md
- Document expected findings and mitigations

**6.6 Update `architecture/README.md`**
- Add "Ethics & Bias Framework" to architecture overview
- Link to detailed documentation

---

### Phase 4: Validation & Reporting (Week 4)

#### 7. Run PatientOne Bias Audit

1. Execute `audit_bias.py` on PatientOne workflow
2. Generate comprehensive bias audit report
3. Document all findings (positive and negative)
4. Validate mitigations are effective

#### 8. Create Example Reports

- `reports/patientone_bias_audit_2026-01-12.html` - Full HTML report
- `reports/patientone_bias_summary.csv` - CSV with metrics
- `reports/patientone_bias_findings.md` - Markdown summary for docs

---

## Implementation Order

### Week 1: Core Documentation
**Day 1-2:** Write ETHICS_AND_BIAS.md (comprehensive framework)
**Day 3:** Write BIAS_AUDIT_CHECKLIST.md (practical checklist)
**Day 4-5:** Write PATIENTONE_BIAS_AUDIT.md (concrete demonstration)

### Week 2: Tools & Automation
**Day 1-2:** Implement `bias_detection.py` utilities
**Day 3:** Write tests for bias detection
**Day 4-5:** Create `audit_bias.py` script and test

### Week 3: Integration
**Day 1-2:** Update hospital deployment docs (HIPAA, OPERATIONS_MANUAL, ADMIN_GUIDE)
**Day 3:** Update EXECUTIVE_SUMMARY and architecture docs
**Day 4:** Update PatientOne README and test docs
**Day 5:** Review all cross-references for consistency

### Week 4: Validation & Documentation
**Day 1-2:** Run manual bias audit on PatientOne, document findings
**Day 3:** Generate example HTML reports for stakeholder review
**Day 4:** Create summary documentation and audit archive
**Day 5:** Internal review and polish (no external validation at this stage)

---

## Success Criteria

1. ✅ Comprehensive ethics framework documented (ETHICS_AND_BIAS.md) aligned with US healthcare standards
2. ✅ Practical audit checklist for US hospital clinicians/researchers
3. ✅ PatientOne workflow audited with concrete findings using diverse reference datasets
4. ✅ Python utilities for manual batch bias detection
5. ✅ Standalone audit script with HTML report generation (manual execution)
6. ✅ Integration with existing HIPAA compliance documentation
7. ✅ Alignment with FDA and WHO AI standards (US healthcare focus)
8. ✅ Clear mitigations for identified biases with recommended diverse datasets
9. ✅ Transparency about limitations and data source ancestry distributions
10. ✅ Quarterly batch auditing framework established
11. ✅ Recommended diverse reference datasets documented (gnomAD, All of Us, ClinGen, TOPMed, 1000 Genomes)

---

## Adjustments from Original Plan

**✅ Implemented:**
1. **US Healthcare Focus:** Removed EU AI Act, kept FDA/WHO standards
2. **Manual Execution:** Batch auditing with quarterly schedule, no CI/CD automation
3. **Diverse Reference Datasets:** Specific recommendations for gnomAD, All of Us, ClinGen, TOPMed, 1000 Genomes, Human Cell Atlas
4. **No External Validation:** Internal review only at this stage
5. **Batch Processing:** HTML reports for manual review, 10-year retention

**Scope Adjustments:**
- Target: US hospitals only
- Standards: FDA SaMD, WHO, NIH All of Us, AMA Ethics
- Execution: Manual quarterly audits
- Tools: Standalone Python script (no dashboards)
- Validation: Internal only (no ethics board review yet)

---

## Remaining Questions for Review

1. **Audit Frequency:** Is quarterly sufficient, or should we audit monthly initially?
2. **Report Format:** Is HTML sufficient, or do we need PDF export for compliance?
3. **Stakeholder Training:** Should we add training materials for running audits?
4. **Integration Priority:** Which existing doc should we update first (HIPAA compliance or Operations Manual)?

---

## References

**US Healthcare AI Standards:**
1. FDA AI/ML-Based SaMD Action Plan (2021) - https://www.fda.gov/medical-devices/software-medical-device-samd/artificial-intelligence-and-machine-learning-aiml-enabled-medical-devices
2. AMA Code of Medical Ethics Opinion 2.3.2 - Physicians' use of AI - https://www.ama-assn.org/delivering-care/ethics/augmented-intelligence-ai-medicine
3. 21st Century Cures Act - https://www.fda.gov/regulatory-information/selected-amendments-fdc-act/21st-century-cures-act

**Global Best Practices:**
4. WHO Ethics and Governance of AI for Health (2021) - https://www.who.int/publications/i/item/9789240029200

**Diverse Reference Datasets:**
5. NIH All of Us Research Program - https://allofus.nih.gov/
6. gnomAD (Genome Aggregation Database) - https://gnomad.broadinstitute.org/
7. ClinGen (Clinical Genome Resource) - https://www.clinicalgenome.org/
8. TOPMed (Trans-Omics for Precision Medicine) - https://www.nhlbiwgs.org/
9. 1000 Genomes Project - https://www.internationalgenome.org/
10. GTEx (Genotype-Tissue Expression) - https://gtexportal.org/
11. Human Cell Atlas - https://www.humancellatlas.org/

**Academic Literature:**
12. Mehrabi et al., "A Survey on Bias and Fairness in Machine Learning" (2021)
13. Obermeyer et al., "Dissecting racial bias in an algorithm used to manage the health of populations" Science 2019
14. Popejoy & Fullerton, "Genomics is failing on diversity" Nature 2016
15. Sirugo et al., "The Missing Diversity in Human Genetic Studies" Cell 2019

---

**Next Step:** Review this revised plan. Once approved, I'll begin implementation with Week 1 documentation (ETHICS_AND_BIAS.md).

**Plan Status:** ✅ Revised per user feedback - Ready for approval
