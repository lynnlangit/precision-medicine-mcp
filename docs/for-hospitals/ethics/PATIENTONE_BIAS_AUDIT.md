# PatientOne Bias Audit Demonstration

**Purpose:** Concrete demonstration of bias audit methodology on PatientOne precision medicine workflow

**Patient:** PAT001-OVC-2025 (63-year-old woman, Stage IV HGSOC, platinum-resistant)

**Audit Date:** 2026-01-12

**Auditor:** Precision Medicine MCP Team

**Status:** DEMONSTRATION (Educational Example)

---

## Executive Summary

**Overall Assessment:** ✅ CONDITIONAL PASS

**Risk Level:** MEDIUM

**Key Findings:**
- 3 biases detected (2 MEDIUM risk, 1 LOW risk)
- 0 critical biases
- 5 biases checked and not detected (insurance, geographic, race coding, spatial algorithms, PDX models)

**Recommended Action:** Deploy with implemented mitigations + quarterly monitoring

**Mitigations Required:**
1. Add ancestry diversity warnings to all genomic reports
2. Flag BRCA variants with <5 studies in patient ancestry
3. Document GTEx 85% European composition in expression reports
4. Cross-validate with TOPMed/All of Us for non-European patients
5. Use ovarian-specific cell type references (GSE146026)

---

## Table of Contents

1. [Workflow Overview](#workflow-overview)
2. [Audit Objectives](#audit-objectives)
3. [Genomics Bias Analysis](#genomics-bias-analysis)
4. [Clinical Bias Analysis](#clinical-bias-analysis)
5. [Spatial Transcriptomics Bias Analysis](#spatial-transcriptomics-bias-analysis)
6. [Multiomics Bias Analysis](#multiomics-bias-analysis)
7. [Summary & Recommendations](#summary--recommendations)

---

## Workflow Overview

### Patient Profile

**Demographics:**
- **Age:** 63 years
- **Sex:** Female
- **Ancestry:** Self-reported European (White)
- **Diagnosis:** Stage IV High-Grade Serous Ovarian Carcinoma (HGSOC)
- **Clinical Status:** Platinum-resistant, post-chemotherapy

**Clinical History:**
- Diagnosed 2023 (age 61)
- Treatment: Carboplatin + Paclitaxel + Bevacizumab (2023-2024)
- Progression: Platinum-resistant disease after 6 cycles
- Current: Exploring precision medicine options

---

### Workflow Scope

**Data Modalities Analyzed:**
1. **Genomics:** BRCA1/BRCA2 sequencing, germline variants
2. **Spatial Transcriptomics:** Tumor gene expression, spatial patterns
3. **Multiomics:** RNA-Protein-Phospho integration from PDX models
4. **Clinical:** FHIR data (conditions, medications, observations)

**MCP Servers Used:**
- mcp-epic (clinical EHR data)
- mcp-fgbio (genomic reference data)
- mcp-spatialtools (spatial analysis)
- mcp-multiomics (multi-omics integration)
- mcp-openimagedata (H&E and MxIF histology)
- mcp-tcga (cohort comparison)

**Total Tool Calls:** 27 across 6 servers

**Analysis Duration:** 127 seconds

**Cost:** $1.45 (Claude Sonnet 4.5)

---

## Audit Objectives

### Primary Objectives

1. **Identify potential biases** across all data modalities
2. **Quantify bias severity** (LOW, MEDIUM, HIGH, CRITICAL)
3. **Document findings** with specific examples
4. **Propose mitigations** using diverse reference datasets

### Audit Questions

**Genomics:**
- Are BRCA variant databases Euro-centric?
- Do gene expression references represent diverse populations?
- Are pathway databases universal or population-specific?

**Clinical:**
- Does insurance status influence treatment recommendations?
- Is geographic location used inappropriately?
- Is race/ethnicity used correctly (ancestry for genomics, NOT biology)?

**Spatial:**
- Are cell type references from diverse populations?
- Do spatial algorithms work equally across tissue types?

**Multiomics:**
- Are PDX models representative of patient populations?
- Do drug sensitivity predictions generalize?

---

## Genomics Bias Analysis

### 1. BRCA1/BRCA2 Variant Interpretation

#### Current Data Sources

**Primary References:**
- **ClinVar:** Variant pathogenicity classifications
- **COSMIC:** Cancer-specific mutations
- **gnomAD:** Population allele frequencies

**Ancestry Distribution:**
```
Database  | European | African | Latino | Asian | Other | Source
----------|----------|---------|--------|-------|-------|--------
ClinVar   | ~70%     | ~10%    | ~10%   | ~5%   | ~5%   | Estimated
COSMIC    | ~75%     | ~8%     | ~8%    | ~6%   | ~3%   | Estimated
gnomAD    | 43%      | 21%     | 14%    | 9%    | 8%    | Published
```

---

#### PatientOne Finding

**Variant:** BRCA1 c.5266dupC (p.Gln1756Profs*74)

**ClinVar Classification:** Pathogenic

**Evidence:**
- 50+ studies in European populations
- Founder mutation in Ashkenazi Jewish population (1.2% allele frequency)
- Well-characterized: frameshift leading to premature stop codon
- Strong association with breast/ovarian cancer risk

**gnomAD Frequencies:**
```
Population              | Allele Frequency | Interpretation
------------------------|------------------|----------------
European (non-Finnish)  | 0.0001          | Rare
Ashkenazi Jewish        | 0.012           | Common (founder)
African/African American| 0.00003         | Very rare
Latino/Admixed American | 0.00008         | Very rare
East Asian              | 0.00001         | Extremely rare
South Asian             | 0.00002         | Extremely rare
```

---

#### Bias Assessment

**For PatientOne (European Ancestry):**
- ✅ **APPROPRIATE:** Variant is well-studied in patient's ancestry
- ✅ **CONFIDENCE:** HIGH (50+ European studies)
- ✅ **CLASSIFICATION:** Pathogenic (supported by gnomAD data)

**For Non-European Patients:**
- ⚠️ **LIMITED DATA:** <5 studies in African, Latino, Asian ancestries
- ⚠️ **LOWER CONFIDENCE:** Pathogenicity extrapolated from European data
- ⚠️ **POTENTIAL IMPACT:** May be classified as VUS (variant of uncertain significance) in other ancestries

---

#### Bias Severity

**Risk Level:** 🟡 MEDIUM

**Justification:**
- ClinVar is 70% European → Adequate for European patients
- gnomAD has 43% European, 21% African, 14% Latino → Good diversity
- However, variant-specific studies are Euro-centric (50+ European vs. <5 others)

---

#### Recommended Mitigations

**✅ Implemented in PatientOne Report:**

1. **Cross-reference with gnomAD:**
   - Always check population-specific allele frequencies
   - Flag if frequency >1% in any population (likely benign)

2. **Document study count by ancestry:**
   ```
   Variant: BRCA1 c.5266dupC
   Studies by Ancestry:
     - European: 50+ studies (HIGH confidence)
     - African: 2 studies (LOW confidence)
     - Latino: 1 study (LOW confidence)
     - Asian: 1 study (LOW confidence)
   ```

3. **Add ancestry-specific disclaimer:**
   ```
   ═══════════════════════════════════════════════════════════════
   ⚠️ ANCESTRY CONSIDERATION

   This variant has been extensively studied in European populations
   (50+ studies) but has limited data in non-European ancestries
   (<5 studies each).

   For patients of non-European ancestry, pathogenicity classification
   may have greater uncertainty. Genetic counseling is recommended.

   Recommended Validation:
   - Check gnomAD for population-specific frequencies
   - Reference All of Us (80% underrepresented groups)
   - Consider ClinGen expert panel review
   ═══════════════════════════════════════════════════════════════
   ```

**📊 Recommended Diverse References:**
- **gnomAD v4.1:** Use for population-specific frequencies (21% African, 14% Latino)
- **ClinGen:** Expert-curated variant assessments with ancestry context
- **All of Us:** Cross-validate with 80% underrepresented groups (when available)

---

### 2. Gene Expression Reference Ranges

#### Current Data Source

**Primary Reference:**
- **GTEx (Genotype-Tissue Expression) v8**
  - 948 deceased donors
  - 54 tissue types (including ovary)
  - **Ancestry Distribution:** 85% European, 10% African American, 5% Other, 0% Asian

---

#### PatientOne Finding

**Analysis:** Differential gene expression in ovarian tumor vs. GTEx normal ovary

**Top Upregulated Genes:**
```
Gene   | Patient FC | GTEx Normal | Interpretation
-------|------------|-------------|----------------
MUC16  | +8.2       | Low         | CA-125, tumor marker
TP53   | +6.5       | Medium      | Tumor suppressor (mutated)
BRCA1  | -3.8       | High        | DNA repair (germline mutation)
```

**Reference Range Issue:**
- GTEx normal ovary expression is **85% European**
- For PatientOne (European), this is APPROPRIATE
- For Asian patients, GTEx has **0% representation**

---

#### Bias Assessment

**For PatientOne (European Ancestry):**
- ✅ **APPROPRIATE:** Patient ancestry matches GTEx majority
- ✅ **CONFIDENCE:** HIGH for differential expression calls

**For Non-European Patients:**
- ⚠️ **GTEx BIAS:** 85% European, 0% Asian
- ⚠️ **POTENTIAL IMPACT:** Normal expression ranges may not generalize
- ⚠️ **EXAMPLE:** If baseline MUC16 expression varies by ancestry, fold-change interpretation may be incorrect

---

#### Bias Severity

**Risk Level:** 🟡 MEDIUM

**Justification:**
- GTEx is the gold-standard tissue expression atlas
- 85% European composition is known limitation
- No comparably comprehensive diverse alternative currently available
- Risk mitigated by using multiple references and larger fold-change thresholds

---

#### Recommended Mitigations

**✅ Implemented in PatientOne Report:**

1. **Document GTEx ancestry distribution:**
   ```
   ═══════════════════════════════════════════════════════════════
   ⚠️ REFERENCE DATA LIMITATION

   Gene expression reference ranges are from GTEx v8:
     - European (non-Finnish): 85%
     - African American: 10%
     - Other: 5%
     - Asian: 0%

   Patient Ancestry: European → APPROPRIATE match

   For Asian patients: Results may not accurately reflect normal
   variation. Validation with TOPMed or All of Us recommended.
   ═══════════════════════════════════════════════════════════════
   ```

2. **Apply larger fold-change thresholds for non-European patients:**
   - Standard: log2FC > 1.5 (2.8-fold change)
   - Non-European: log2FC > 2.0 (4-fold change) for conservatism

3. **Cross-validate with diverse references when patient ancestry differs:**
   - **TOPMed:** 180,000+ genomes, diverse US population
   - **Human Cell Atlas:** Single-cell references, global diversity
   - **All of Us:** Validate expression ranges in underrepresented groups

**📊 Recommended Diverse References:**
- **GTEx v8:** Use with explicit 85% European caveat
- **TOPMed RNA-seq:** Validate expression ranges in diverse populations
- **Human Cell Atlas:** Ovarian single-cell data for cell-type-specific expression
- **All of Us Genomic Data:** Cross-validate expression patterns (when available)

---

### 3. Pathway Enrichment Databases

#### Current Data Sources

**Primary References:**
- **KEGG (Kyoto Encyclopedia of Genes and Genomes)**
- **Reactome**
- **GO (Gene Ontology)**

**Ancestry Consideration:**
- Pathway definitions based on **aggregate human data**
- Not stratified by ancestry
- Assumption: Pathways are universal across populations

---

#### PatientOne Finding

**Enriched Pathways (from differential expression):**
```
Pathway                          | p-value  | Genes | Database
---------------------------------|----------|-------|----------
DNA repair (BRCA pathway)        | 1.2e-08  | 12    | KEGG
p53 signaling                    | 3.4e-07  | 15    | Reactome
Apoptosis                        | 8.9e-06  | 18    | GO
Cell cycle checkpoint            | 2.1e-05  | 14    | KEGG
```

---

#### Bias Assessment

**Current Understanding:**
- ✅ Core biological pathways (DNA repair, apoptosis, cell cycle) are **conserved across ancestries**
- ⚠️ Pathway **activity levels** may vary by ancestry due to:
  - Genetic modifiers (ancestry-specific variants affecting pathway genes)
  - Gene expression differences (baseline expression varies by population)
  - Environmental factors (diet, environment interacting with genetics)

**For PatientOne:**
- ✅ **APPROPRIATE:** DNA repair and p53 pathways are universal
- ✅ **CONFIDENCE:** HIGH for pathway enrichment findings

**Potential Concern:**
- ⚠️ **DRUG RESPONSE PATHWAYS:** Pharmacogenomic pathways may show ancestry-specific variation
- ⚠️ **EXAMPLE:** CYP450 metabolism varies significantly by ancestry → Drug dosing implications

---

#### Bias Severity

**Risk Level:** 🟢 LOW

**Justification:**
- Core cancer pathways (DNA repair, p53, apoptosis) are biologically conserved
- Pathway enrichment is statistical test, not absolute truth
- Cross-validation across multiple databases (KEGG + Reactome + GO) mitigates single-source bias

---

#### Recommended Mitigations

**✅ Implemented in PatientOne Report:**

1. **Cross-validate across multiple pathway databases:**
   - Use KEGG **AND** Reactome **AND** GO
   - Flag pathways found in only 1 database (lower confidence)
   - Prioritize pathways found in all 3 (higher confidence)

2. **Flag pathways with known ancestry-specific variation:**
   - Pharmacogenomic pathways (CYP450, drug metabolism)
   - Immune response pathways (HLA variation)
   - Nutrient metabolism (lactose tolerance, alcohol metabolism)

3. **Validate pathway associations in All of Us (when available):**
   - Check if pathway-phenotype associations hold in diverse cohorts
   - Example: "BRCA pathway enrichment → ovarian cancer risk" validated across ancestries?

**📊 Recommended Diverse References:**
- **KEGG + Reactome + GO:** Cross-validate across all three
- **All of Us Workbench:** Validate pathway associations in diverse cohorts (when available)
- **PharmGKB:** Ancestry-specific drug-gene interaction data

---

## Clinical Bias Analysis

### 1. Insurance Status & Treatment Recommendations

#### Audit Objective

**Question:** Does insurance status influence treatment recommendations?

**Ethical Requirement:** Insurance status should NEVER affect clinical recommendations (only billing and access to resources).

---

#### Data Source

**FHIR Resources:**
- `Patient` - Demographics
- `Condition` - Diagnosis (Stage IV HGSOC)
- `Observation` - Lab values (CA-125, BRCA status)
- `MedicationStatement` - Treatment history
- `Coverage` - Insurance information

---

#### PatientOne Finding

**Treatment Recommendation:** Based on:
1. **BRCA1 mutation:** Indicates PARP inhibitor sensitivity
2. **Platinum-resistant:** Rules out platinum re-challenge
3. **CA-125 elevated:** Confirms active disease

**Insurance Status:** NOT used in recommendation algorithm

**Verification:**
```python
# Reviewed mcp-epic tool calls
tool_calls = [
    "get_patient_demographics",      # ✅ Name, age, ancestry
    "get_patient_conditions",         # ✅ Stage IV HGSOC diagnosis
    "get_patient_observations",       # ✅ CA-125, BRCA status
    "get_patient_medications"         # ✅ Treatment history
]

# Coverage resource NOT queried
assert "get_patient_coverage" not in tool_calls  # ✅ PASS
```

**Treatment Recommendation:**
```
Recommended Next-Line Therapy:
1. PARP inhibitor (Olaparib or Niraparib) - BRCA1 mutation indicates sensitivity
2. Clinical trial enrollment - Platinum-resistant HGSOC
3. Consider immunotherapy - PD-L1 status pending

Insurance Status: NOT CONSIDERED in recommendation
```

---

#### Bias Assessment

**Result:** ✅ **NO BIAS DETECTED**

**Evidence:**
- Insurance/coverage data NOT queried during analysis
- Treatment recommendations based solely on:
  - Molecular markers (BRCA1, CA-125)
  - Clinical status (platinum-resistant, Stage IV)
  - Evidence-based guidelines (NCCN ovarian cancer guidelines)

---

#### Bias Severity

**Risk Level:** 🟢 NONE

**Conclusion:** Insurance bias is **NOT present** in PatientOne workflow.

---

### 2. Geographic & Socioeconomic Factors

#### Audit Objective

**Question:** Is geographic location (ZIP code) used as proxy for socioeconomic status to influence treatment eligibility?

**Ethical Requirement:** ZIP code should ONLY be used for:
- ✅ Resource access (finding nearest clinic, pharmacy)
- ✅ Logistics (scheduling, transportation assistance)
- ❌ NOT for treatment eligibility or clinical recommendations

---

#### Data Source

**FHIR Patient.address:**
```json
{
  "address": [{
    "use": "home",
    "line": ["456 Oak Street"],
    "city": "Boston",
    "state": "MA",
    "postalCode": "02108",
    "country": "US"
  }]
}
```

---

#### PatientOne Finding

**Address Usage:**
- ✅ Used for healthcare provider coordination (finding local oncologist)
- ✅ Used for pharmacy selection (nearest specialty pharmacy for PARP inhibitors)
- ❌ **NOT** used in treatment recommendation algorithm

**Verification:**
```python
# Check if ZIP code is in feature importance
feature_importance = {
    "BRCA1_status": 0.45,      # ✅ Clinical
    "CA125_level": 0.32,       # ✅ Clinical
    "platinum_resistant": 0.18, # ✅ Clinical
    "stage": 0.05,             # ✅ Clinical
    "postal_code": 0.0         # ✅ NOT USED
}

assert feature_importance["postal_code"] == 0.0  # ✅ PASS
```

---

#### Bias Assessment

**Result:** ✅ **NO BIAS DETECTED**

**Evidence:**
- ZIP code NOT in algorithm feature set
- No correlation between geographic location and treatment recommendations
- Address used appropriately (resource coordination only)

---

#### Bias Severity

**Risk Level:** 🟢 NONE

**Conclusion:** Geographic/socioeconomic bias is **NOT present** in PatientOne workflow.

---

### 3. Race/Ethnicity Coding

#### Audit Objective

**Question:** Is race/ethnicity used appropriately?

**Ethical Requirements:**
- ✅ **Appropriate:** Ancestry for genomic variant interpretation (genetic context)
- ❌ **Inappropriate:** Race as biological proxy (e.g., assuming kidney function differs by race)

---

#### Data Source

**FHIR Patient.extension (US Core Race/Ethnicity):**
```json
{
  "extension": [{
    "url": "http://hl7.org/fhir/us/core/StructureDefinition/us-core-race",
    "extension": [{
      "url": "ombCategory",
      "valueCoding": {
        "system": "urn:oid:2.16.840.1.113883.6.238",
        "code": "2106-3",
        "display": "White"
      }
    }, {
      "url": "text",
      "valueString": "European ancestry"
    }]
  }]
}
```

---

#### PatientOne Finding

**Ancestry Usage:**
1. **Genomic Variant Interpretation:**
   - ✅ APPROPRIATE: "European ancestry" used to select appropriate gnomAD population frequencies
   - ✅ APPROPRIATE: Checked if BRCA1 variant has sufficient studies in European populations (50+)

2. **Treatment Recommendation:**
   - ✅ APPROPRIATE: Ancestry NOT used to determine treatment eligibility
   - ✅ APPROPRIATE: No race-based treatment algorithms (e.g., no "use Drug X for White patients, Drug Y for Black patients")

3. **Clinical Calculations:**
   - ✅ APPROPRIATE: No race-based eGFR calculation (race coefficient removed)
   - ✅ APPROPRIATE: No race-based dosing adjustments (unless pharmacogenomic indication)

---

#### Bias Assessment

**Result:** ✅ **APPROPRIATE USE**

**Evidence:**
- Ancestry used ONLY for genomic context (matching patient to reference populations)
- Race NOT used as biological proxy
- Consistent with best practices: "Ancestry is genetic, race is social construct"

---

#### Best Practice Validation

**✅ Ancestry for Genomics (CORRECT):**
```
Patient: European ancestry
Variant: BRCA1 c.5266dupC
Reference: gnomAD European population (43% of database)
Interpretation: Pathogenic (consistent with European studies)
```

**❌ Race as Biology (WOULD BE INCORRECT):**
```
# EXAMPLE OF INAPPROPRIATE USE (NOT in PatientOne)
Patient: White race
eGFR Calculation: Multiply creatinine clearance by 1.0 (White) vs. 1.2 (Black)
Conclusion: Race-based medical algorithm (INAPPROPRIATE)
```

**PatientOne does NOT use race-based clinical algorithms.** ✅

---

#### Bias Severity

**Risk Level:** 🟢 NONE

**Conclusion:** Race/ethnicity coding is **APPROPRIATE** in PatientOne workflow.

---

## Spatial Transcriptomics Bias Analysis

### 1. Cell Type Reference Signatures

#### Current Data Source

**Primary Reference:**
- **SingleR default references:** Generic immune cell signatures
  - Source: Bulk RNA-seq from multiple tissue types
  - Ancestry: Not explicitly documented (assumed diverse but unknown)
  - Tissue specificity: Pan-tissue (not ovarian-specific)

---

#### PatientOne Finding

**Cell Type Deconvolution Results:**
```
Cell Type                | Proportion | Confidence
-------------------------|------------|------------
Tumor cells (epithelial) | 68%        | HIGH
T cells (CD8+)           | 12%        | MEDIUM
Fibroblasts (CAFs)       | 15%        | HIGH
Macrophages              | 3%         | MEDIUM
B cells                  | 2%         | LOW
```

**Concern:**
- Generic immune references may NOT capture ovarian cancer-specific cell states
- Example: Tumor-associated macrophages (TAMs) in ovary may differ from TAMs in lung
- Ovarian cancer has unique tumor microenvironment (ascites, peritoneal seeding)

---

#### Bias Assessment

**For PatientOne (European Ancestry):**
- ⚠️ **MODERATE CONCERN:** Using generic (not ovarian-specific) references
- ⚠️ **POTENTIAL IMPACT:** May misclassify ovarian-specific cell states
- ✅ **MITIGATION:** Cross-validate with H&E pathology review

**For Non-European Patients:**
- ⚠️ **UNKNOWN ANCESTRY:** SingleR reference ancestry not documented
- ⚠️ **POTENTIAL IMPACT:** If references are Euro-centric, deconvolution accuracy may vary

---

#### Bias Severity

**Risk Level:** 🟡 LOW-MEDIUM

**Justification:**
- Cell types (epithelial, immune, fibroblast) are biologically conserved
- However, cell **states** (e.g., exhausted T cells, M1 vs. M2 macrophages) may be tissue- and ancestry-specific
- Risk partially mitigated by pathology review

---

#### Recommended Mitigations

**✅ Implemented in PatientOne Workflow:**

1. **Use ovarian-specific reference when available:**
   - **GSE146026:** Single-cell RNA-seq of human ovarian tumors
   - Contains: Tumor cells, T cells, CAFs, macrophages, B cells, dendritic cells
   - Ovarian-specific cell states (e.g., ascites-associated macrophages)

2. **Cross-validate with H&E pathology review:**
   - Pathologist estimates: "70-80% tumor cells, moderate immune infiltration"
   - Deconvolution estimates: "68% tumor cells, 12% T cells"
   - ✅ **CONCORDANCE:** High agreement between deconvolution and pathology

3. **Document reference ancestry (when known):**
   - For Human Cell Atlas references, document donor demographics
   - Flag if patient ancestry is underrepresented in reference

**📊 Recommended Diverse References:**
- **GSE146026:** Ovarian cancer single-cell atlas (tissue-specific)
- **Human Cell Atlas:** Ovarian tissue data (global diversity)
- **Custom references:** Validate with local patient population when possible

---

### 2. Spatial Autocorrelation Algorithms

#### Audit Objective

**Question:** Do spatial algorithms (Moran's I, Geary's C) work equally well across tissue types and ancestries?

---

#### PatientOne Finding

**Spatial Autocorrelation Analysis:**
```
Gene    | Moran's I | p-value | Interpretation
--------|-----------|---------|----------------
MUC16   | 0.72      | <0.001  | Strong spatial clustering
TP53    | 0.58      | <0.001  | Moderate clustering
BRCA1   | 0.35      | 0.003   | Weak clustering
```

**Algorithm:** Moran's I (standard spatial statistics)

---

#### Bias Assessment

**Result:** ✅ **NO BIAS DETECTED**

**Justification:**
- Moran's I is a **mathematical algorithm** (not trained on data)
- Tissue-agnostic: Works on any spatial coordinate system
- Ancestry-agnostic: No population-specific parameters

**Formula (Ancestry-Independent):**
```
Moran's I = (N/W) * Σ(w_ij * (x_i - x̄) * (x_j - x̄)) / Σ(x_i - x̄)²

Where:
  N = number of spatial locations
  W = sum of spatial weights
  w_ij = spatial weight between locations i and j
  x_i = gene expression at location i
  x̄ = mean expression across all locations
```

---

#### Bias Severity

**Risk Level:** 🟢 NONE

**Conclusion:** Spatial autocorrelation algorithms are **NOT biased** by ancestry or tissue type.

---

## Multiomics Bias Analysis

### 1. PDX Model Representativeness

#### Current Data Source

**PDX (Patient-Derived Xenograft) Models:**
- Source: Ovarian cancer patients (diverse cohort)
- Modalities: RNA (gene expression), Protein (TMT proteomics), Phospho (phosphoproteomics)
- Use: Drug sensitivity prediction, upstream regulator analysis

---

#### PatientOne Finding

**Multiomics Integration:**
- Identified: 42 RNA-Protein-Phospho associations
- Top pathway: PI3K/AKT signaling (dysregulated)
- Predicted drug targets: PI3K inhibitors, AKT inhibitors

**PDX Cohort Diversity:**
- Ovarian cancer PDX models from multiple institutions
- Represents diverse patient populations (assumed, but not explicitly validated)

---

#### Bias Assessment

**Concern:**
- ⚠️ **PDX LIMITATION:** Xenografts lack human immune microenvironment
- ⚠️ **ANCESTRY UNKNOWN:** PDX model donor demographics not always documented
- ⚠️ **POTENTIAL IMPACT:** Drug response may differ in patients with intact immune system

**For PatientOne:**
- ✅ **MITIGATION:** Combined PDX predictions with patient's own spatial transcriptomics
  - PDX: Tumor-intrinsic biology (RNA, protein, phospho)
  - Patient spatial: Tumor microenvironment (immune infiltration)

---

#### Bias Severity

**Risk Level:** 🟡 LOW-MEDIUM

**Justification:**
- PDX models are valuable despite limitations
- Combining PDX data with patient-specific spatial data mitigates bias
- Drug predictions should be validated in clinical trials (diverse cohorts)

---

#### Recommended Mitigations

**✅ Implemented in PatientOne Workflow:**

1. **Combine PDX with patient-specific data:**
   - PDX: Tumor-intrinsic drug sensitivity
   - Patient spatial: Microenvironment context (immune infiltration, CAFs)
   - Integration: More comprehensive than PDX alone

2. **Validate drug predictions with clinical trial data:**
   - Check if PI3K inhibitor trials included diverse populations
   - Example: KEYNOTE-100 (pembrolizumab in ovarian cancer) → 23% non-White patients

3. **Document PDX limitations:**
   ```
   ═══════════════════════════════════════════════════════════════
   ⚠️ PDX MODEL LIMITATION

   Drug sensitivity predictions are based on Patient-Derived Xenograft
   (PDX) models, which have the following limitations:

   1. Lack of human immune microenvironment
   2. PDX donor demographics not always documented
   3. Drug response may differ in patients with intact immunity

   Mitigation: Combined PDX predictions with patient's spatial
   transcriptomics data to capture tumor microenvironment context.

   Recommendation: Validate predictions in clinical trials with
   diverse patient populations.
   ═══════════════════════════════════════════════════════════════
   ```

**📊 Recommended Diverse References:**
- **PDX cohorts:** Ensure diverse donor representation when available
- **Clinical trial data:** Cross-validate with diverse trial populations
- **Patient-specific spatial data:** Complement PDX with patient's own microenvironment

---

## Summary & Recommendations

### Biases Detected

#### 1. BRCA Variant Databases: Euro-Centric (MEDIUM Risk)

**Finding:**
- ClinVar/COSMIC are ~70% European
- BRCA1 c.5266dupC has 50+ European studies, <5 in other ancestries

**Impact:**
- ✅ PatientOne (European) → Well-studied, HIGH confidence
- ⚠️ Non-European patients → Lower confidence, may be VUS

**Mitigation:**
- Flag variants with <5 studies in patient ancestry
- Cross-reference gnomAD (43% European, 21% African, 14% Latino)
- Recommend genetic counseling for VUS in underrepresented ancestries

**Status:** ✅ MITIGATED (disclaimers added, gnomAD cross-reference implemented)

---

#### 2. GTEx Reference Ranges: 85% European (MEDIUM Risk)

**Finding:**
- GTEx is 85% European, 10% African, 0% Asian
- Gene expression reference ranges may not generalize

**Impact:**
- ✅ PatientOne (European) → Appropriate match
- ⚠️ Asian patients → 0% representation, HIGH uncertainty

**Mitigation:**
- Document GTEx 85% European composition in all reports
- Apply larger fold-change thresholds for non-European patients (log2FC > 2.0)
- Cross-validate with TOPMed/Human Cell Atlas for non-European patients

**Status:** ✅ MITIGATED (disclaimers added, validation strategy documented)

---

#### 3. Cell Type References: Generic, Not Cancer-Specific (LOW Risk)

**Finding:**
- Using generic immune cell references (not ovarian-specific)
- May miss ovarian cancer-specific cell states

**Impact:**
- ⚠️ Moderate concern for all patients (ancestry-independent)
- Potential misclassification of tissue-specific cell states

**Mitigation:**
- Use ovarian-specific references (GSE146026)
- Cross-validate with H&E pathology review

**Status:** ✅ MITIGATED (ovarian-specific references recommended, pathology concordance confirmed)

---

### Biases NOT Detected

#### 1. Insurance Status ✅

**Finding:** Insurance/coverage data NOT used in treatment recommendations

**Evidence:** Coverage resource NOT queried, feature importance = 0%

**Status:** ✅ NO BIAS

---

#### 2. Geographic Location ✅

**Finding:** ZIP code NOT used for treatment eligibility

**Evidence:** ZIP code feature importance = 0%, used only for resource coordination

**Status:** ✅ NO BIAS

---

#### 3. Race/Ethnicity Coding ✅

**Finding:** Ancestry used appropriately (genomics only, NOT as biological proxy)

**Evidence:** No race-based clinical algorithms (eGFR, dosing, treatment)

**Status:** ✅ APPROPRIATE USE

---

#### 4. Spatial Algorithms ✅

**Finding:** Moran's I is mathematical algorithm, ancestry- and tissue-agnostic

**Evidence:** Algorithm has no population-specific parameters

**Status:** ✅ NO BIAS

---

#### 5. PDX Models ✅ (with mitigation)

**Finding:** PDX models lack immune microenvironment, ancestry unknown

**Mitigation:** Combined with patient-specific spatial data

**Status:** ✅ MITIGATED

---

### Overall Risk Assessment

**Genomics:** 🟡 MEDIUM
- BRCA databases Euro-centric
- Mitigated by gnomAD cross-reference and disclaimers

**Gene Expression:** 🟡 MEDIUM
- GTEx 85% European, 0% Asian
- Mitigated by ancestry documentation and validation strategy

**Clinical:** 🟢 LOW
- No insurance or geographic bias
- Appropriate use of ancestry

**Spatial:** 🟢 LOW
- Generic cell references (tissue, not ancestry issue)
- Mitigated by ovarian-specific references

**Multiomics:** 🟢 LOW-MEDIUM
- PDX limitations acknowledged
- Mitigated by combining with patient spatial data

**OVERALL:** 🟡 MEDIUM (acceptable with implemented mitigations)

---

### Deployment Decision

**✅ RECOMMENDED: Deploy with Mitigations + Quarterly Monitoring**

**Conditions:**
1. All mitigations implemented (ancestry disclaimers, diverse references)
2. Quarterly bias audits to monitor for drift
3. Genetic counseling available for VUS cases
4. Continuous monitoring of fairness metrics

---

### Next Steps

**Immediate (Before Production):**
- [x] Add ancestry diversity warnings to all genomic reports
- [x] Document GTEx 85% European composition
- [x] Implement gnomAD cross-reference
- [x] Use ovarian-specific cell type references

**Quarterly (Every 3 Months):**
- [ ] Re-run fairness metrics on new patient data
- [ ] Check for updated diverse reference datasets
- [ ] Validate no bias drift over time

**Annual (Every 12 Months):**
- [ ] Comprehensive bias audit (all 4 steps)
- [ ] Update reference datasets (gnomAD, GTEx, All of Us)
- [ ] Review clinical validation with expert panel

---

**For methodology details, see:**
- [ETHICS_AND_BIAS.md](ETHICS_AND_BIAS.md) - Comprehensive framework
- [BIAS_AUDIT_CHECKLIST.md](BIAS_AUDIT_CHECKLIST.md) - Practical checklist
- [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md) - Week-by-week plan

**Questions?** File an issue on GitHub or consult institutional bioethics committee.

---

**Document Status:** ✅ Complete (Week 1 Deliverable)
**Last Updated:** 2026-01-12
**Version:** 1.0
