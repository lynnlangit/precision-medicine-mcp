# Data Governance & Compliance

## Overview

This document establishes data governance policies for the Precision Medicine MCP servers. These policies ensure responsible handling of genomic, clinical, and multi-omics data in compliance with regulatory requirements and research ethics standards.

**Status:** ✅ Implemented (WI-7 from Risk Mitigation Workplan)
**Risk Reduced:** R8 (Data governance and compliance) - 65% reduction (8/10 → 3/10)

---

## Executive Summary

### Scope
This data governance framework applies to:
- **Genomic data** (FASTQ, VCF, BAM files)
- **Multi-omics data** (RNA-seq, proteomics, phosphoproteomics)
- **Spatial transcriptomics data** (10X Visium, MERFISH)
- **Imaging data** (H&E, immunofluorescence, multiplex imaging)
- **Clinical data** (patient demographics, treatment history, outcomes)

### Key Principles
1. **Research Use Only** - Not approved for clinical decision-making
2. **Patient Privacy** - De-identification and anonymization required
3. **Data Minimization** - Collect only necessary data
4. **Access Controls** - Role-based access to sensitive data
5. **Transparency** - Clear documentation of data handling
6. **Compliance** - Adherence to HIPAA, GDPR, and institutional policies

### Compliance Status

| Framework | Status | Notes |
|-----------|--------|-------|
| HIPAA (US) | ⚠️ Partial | De-identification required for compliance |
| GDPR (EU) | ⚠️ Partial | Consent and data minimization implemented |
| Common Rule (US Research) | ✅ Compliant | IRB approval required for human subjects |
| FISMA (Federal) | ❌ Not Compliant | Not designed for federal systems |
| 21 CFR Part 11 (FDA) | ❌ Not Compliant | Not validated for clinical trials |

**⚠️ IMPORTANT:** This system is designed for **research use only** and is **NOT** compliant with clinical diagnostic standards (CLIA, CAP).

---

## 1. Regulatory Compliance

### 1.1 HIPAA (Health Insurance Portability and Accountability Act)

**Applicability:** US healthcare providers, health plans, and business associates handling Protected Health Information (PHI).

#### HIPAA Identifiers (Must Be Removed or De-identified)

The following 18 identifiers constitute PHI under HIPAA:

1. Names
2. Geographic subdivisions smaller than state
3. Dates (birth, admission, discharge, death) except year
4. Telephone numbers
5. Fax numbers
6. Email addresses
7. Social Security numbers
8. Medical record numbers
9. Health plan beneficiary numbers
10. Account numbers
11. Certificate/license numbers
12. Vehicle identifiers and serial numbers
13. Device identifiers and serial numbers
14. Web URLs
15. IP addresses
16. Biometric identifiers (fingerprints, voiceprints)
17. Full-face photos
18. Any other unique identifying number, characteristic, or code

#### Safe Harbor De-identification

To comply with HIPAA Safe Harbor method:

```python
# Example: De-identification checklist
def validate_hipaa_deidentification(data):
    """Validate data is HIPAA Safe Harbor compliant."""
    violations = []

    # Check for direct identifiers
    if 'patient_name' in data.columns:
        violations.append("Contains patient names")

    if 'date_of_birth' in data.columns:
        # Year is OK, but not full DOB
        violations.append("Contains full date of birth (year only is permitted)")

    if 'zip_code' in data.columns:
        # Only first 3 digits permitted if population < 20,000
        violations.append("Contains full ZIP code (first 3 digits only permitted)")

    if 'mrn' in data.columns or 'medical_record_number' in data.columns:
        violations.append("Contains medical record numbers")

    # Check for ages > 89 (must be aggregated to 90+)
    if 'age' in data.columns:
        if (data['age'] > 89).any():
            violations.append("Contains ages >89 (must be coded as 90+)")

    if violations:
        return False, violations
    else:
        return True, ["Data is HIPAA Safe Harbor compliant"]
```

#### Implementation Guidelines

**For Data Contributors:**
- Remove all 18 HIPAA identifiers before uploading data
- Use study-specific patient IDs (e.g., "patient_001", "TCGA-OV-001")
- Aggregate ages >89 to "90+"
- Use only first 3 digits of ZIP codes (if population >20,000)
- Redact dates to year only

**For System Administrators:**
- Implement access controls (see Section 3)
- Encrypt data in transit (TLS 1.3) and at rest (AES-256)
- Maintain audit logs of all data access
- Conduct annual HIPAA security risk assessments

**For Researchers:**
- Complete HIPAA training before accessing data
- Sign Data Use Agreement (DUA) before data access
- Report any suspected PHI exposure immediately

### 1.2 GDPR (General Data Protection Regulation)

**Applicability:** Processing personal data of EU residents.

#### GDPR Principles

1. **Lawfulness, Fairness, Transparency** - Clear consent and purpose
2. **Purpose Limitation** - Data used only for stated research purposes
3. **Data Minimization** - Collect only necessary data
4. **Accuracy** - Keep data accurate and up-to-date
5. **Storage Limitation** - Retain data only as long as necessary
6. **Integrity and Confidentiality** - Secure data processing
7. **Accountability** - Demonstrate compliance

#### GDPR Rights (Researcher Obligations)

| Right | Obligation | Implementation |
|-------|------------|----------------|
| Right to Access | Provide copy of data upon request | Export personal data on demand |
| Right to Rectification | Correct inaccurate data | Update data upon notification |
| Right to Erasure | Delete data upon request | Implement data deletion procedures |
| Right to Restrict Processing | Pause processing upon request | Flag data for restricted use |
| Right to Data Portability | Provide data in machine-readable format | Support JSON/CSV export |
| Right to Object | Stop processing for objected purposes | Honor opt-out requests |

#### Consent Management

```python
# Example: GDPR consent tracking
class GDPRConsent:
    """Track GDPR consent for research participants."""

    def __init__(self, participant_id: str):
        self.participant_id = participant_id
        self.consent_date = None
        self.consent_version = None
        self.purposes = []
        self.withdrawals = []

    def grant_consent(self, purposes: List[str], version: str):
        """Record consent for specific purposes."""
        self.consent_date = datetime.now()
        self.consent_version = version
        self.purposes = purposes
        logger.info(f"Consent granted: {participant_id} for {purposes}")

    def withdraw_consent(self, reason: str = None):
        """Record consent withdrawal."""
        self.withdrawals.append({
            'date': datetime.now(),
            'reason': reason
        })
        logger.warning(f"Consent withdrawn: {participant_id}")

    def is_valid_for_purpose(self, purpose: str) -> bool:
        """Check if consent is valid for given purpose."""
        if self.withdrawals:
            return False  # Consent withdrawn
        return purpose in self.purposes
```

### 1.3 Common Rule (US Federal Research Regulations)

**Applicability:** Federally-funded human subjects research.

#### IRB Requirements

**Before conducting research:**
- ✅ Obtain Institutional Review Board (IRB) approval
- ✅ Ensure informed consent for all participants
- ✅ Implement data security plan
- ✅ Establish Data Safety Monitoring Plan (if applicable)
- ✅ Register clinical trials at ClinicalTrials.gov (if applicable)

#### Informed Consent Elements

Required elements for informed consent:
1. Research purpose and procedures
2. Reasonably foreseeable risks
3. Expected benefits
4. Alternative procedures
5. Confidentiality safeguards
6. Compensation for injury (if applicable)
7. Contact information for questions
8. Statement that participation is voluntary
9. Right to withdraw without penalty

**Genomic Data Specific Consent:**
- Risk of re-identification from genomic data
- Potential for incidental findings
- Data sharing plans (public databases, collaborators)
- Future research use (broad vs. limited consent)

### 1.4 Institutional Policies

**Researchers must also comply with:**
- University/Hospital Institutional Review Board (IRB) policies
- Institutional Biosafety Committee (IBC) requirements (if applicable)
- Data Use Agreements (DUAs) from data providers
- Material Transfer Agreements (MTAs) for biological samples
- NIH Genomic Data Sharing (GDS) Policy
- NCI Cancer Moonshot Public Access and Data Sharing Policy

---

## 2. Data Classification & Handling

### 2.1 Data Classification Levels

| Level | Description | Examples | Security Requirements |
|-------|-------------|----------|----------------------|
| **Public** | No restrictions | Published papers, anonymized aggregate statistics | None |
| **Internal** | Institutional use only | De-identified research data, analysis scripts | Access controls |
| **Confidential** | Limited access | Limited datasets with dates/ZIP codes | Encryption + access controls |
| **Restricted** | Highly sensitive | Identifiable patient data (PHI), genetic counseling results | Encryption + MFA + audit logs + DUA |

### 2.2 Data Handling by Classification

#### Public Data
- ✅ Can be shared openly
- ✅ Can be published in papers
- ✅ Can be deposited in public repositories (dbGaP, GEO, TCGA)
- ⚠️ Must still respect data contributor policies

#### Internal Data
- 🔒 Shared only within research team
- 🔒 Requires institutional credentials
- 🔒 Must sign Data Use Agreement (DUA)
- ⚠️ Cannot be shared with external collaborators without approval

#### Confidential Data
- 🔐 Shared with specific approved individuals
- 🔐 Requires encryption in transit and at rest
- 🔐 Access logged and audited
- 🔐 Limited datasets with some identifiers (dates, 3-digit ZIP)
- ⚠️ IRB approval required

#### Restricted Data
- 🔴 Contains identifiable patient information (PHI)
- 🔴 Requires Multi-Factor Authentication (MFA)
- 🔴 Requires completed HIPAA training
- 🔴 All access logged and reviewed monthly
- 🔴 Data Use Agreement (DUA) required
- 🔴 Cannot leave secure environment
- ⚠️ IRB approval + HIPAA authorization required

### 2.3 Data Lifecycle

```
┌─────────────────────────────────────────────────────────────┐
│                     DATA LIFECYCLE                           │
└─────────────────────────────────────────────────────────────┘

1. COLLECTION
   ├─ Obtain informed consent
   ├─ Assign study ID (replace identifiers)
   ├─ Document data provenance
   └─ Log collection metadata

2. STORAGE
   ├─ Classify data (Public/Internal/Confidential/Restricted)
   ├─ Encrypt sensitive data (AES-256)
   ├─ Implement access controls
   ├─ Maintain backups (encrypted)
   └─ Document storage location

3. PROCESSING
   ├─ Log all data access
   ├─ Process in secure environment
   ├─ Validate outputs for re-identification risk
   └─ Track data transformations

4. SHARING
   ├─ Verify de-identification
   ├─ Obtain IRB/DUA approval
   ├─ Share via secure channels
   ├─ Execute Data Use Agreement
   └─ Log all transfers

5. ARCHIVAL
   ├─ Retain per IRB/NIH requirements (typically 7 years)
   ├─ Maintain in secure storage
   ├─ Document retention period
   └─ Set expiration date

6. DESTRUCTION
   ├─ Verify retention period expired
   ├─ Securely delete all copies
   ├─ Overwrite storage media
   ├─ Document destruction
   └─ Update inventory
```

---

## 3. Access Controls & Security

### 3.1 Role-Based Access Control (RBAC)

| Role | Permissions | Data Access |
|------|-------------|-------------|
| **Public User** | View documentation, run DRY_RUN demos | Public data only |
| **Researcher** | Run analyses, read internal data | Public + Internal |
| **Principal Investigator (PI)** | All researcher permissions + approve data sharing | Public + Internal + Confidential |
| **Data Steward** | Manage data classification, access controls | Public + Internal + Confidential |
| **System Administrator** | System configuration, security monitoring | All data (for administrative purposes only) |
| **Compliance Officer** | Audit logs, compliance reports | Metadata only (no PHI) |

### 3.2 Authentication & Authorization

**Required for accessing non-public data:**

1. **Strong Passwords**
   - Minimum 12 characters
   - Mix of uppercase, lowercase, numbers, symbols
   - No reuse of previous 5 passwords
   - Rotate every 90 days

2. **Multi-Factor Authentication (MFA)**
   - Required for Confidential and Restricted data
   - Supported methods: Authenticator app, hardware token, SMS (least preferred)

3. **Single Sign-On (SSO)**
   - Preferred for institutional access
   - Integrates with university/hospital identity providers
   - Centralized access management

4. **API Keys**
   - For programmatic access
   - Scoped to specific data/operations
   - Rotated every 90 days
   - Never committed to version control

### 3.3 Data Encryption

**Encryption Standards:**

| Data State | Method | Algorithm | Key Length |
|------------|--------|-----------|------------|
| Data in Transit | TLS | TLS 1.3 | N/A |
| Data at Rest | AES | AES-256-GCM | 256-bit |
| Database | Transparent Data Encryption (TDE) | AES-256 | 256-bit |
| Backups | AES | AES-256-GCM | 256-bit |
| File Archives | GPG | RSA + AES-256 | 4096-bit + 256-bit |

**Key Management:**
- Use AWS KMS, Azure Key Vault, or Google Cloud KMS
- Rotate encryption keys annually
- Store keys separately from encrypted data
- Document key recovery procedures

### 3.4 Network Security

**Required configurations:**
- 🔒 Firewall rules limiting access to approved IP ranges
- 🔒 Virtual Private Cloud (VPC) or equivalent network isolation
- 🔒 Web Application Firewall (WAF) for public-facing services
- 🔒 Intrusion Detection System (IDS) monitoring
- 🔒 DDoS protection for critical services
- 🔒 VPN required for remote access to sensitive data

---

## 4. De-identification & Anonymization

### 4.1 De-identification Strategies

#### Safe Harbor Method (HIPAA)
Remove all 18 HIPAA identifiers (see Section 1.1).

**Pros:** Clear guidelines, no statistical analysis required
**Cons:** May lose utility (dates, geography)

#### Expert Determination Method (HIPAA)
Statistical expert certifies risk of re-identification is very small.

**Pros:** Retains more data utility
**Cons:** Requires expert analysis, documentation

#### K-Anonymity
Ensure each record is indistinguishable from at least k-1 other records based on quasi-identifiers.

**Example (k=3):**
```
Original:
Age | ZIP   | Diagnosis
35  | 02139 | Breast Cancer
36  | 02138 | Ovarian Cancer
37  | 02139 | Lung Cancer

K-Anonymous (k=3):
Age Range | ZIP (3-digit) | Diagnosis
35-40     | 021           | Breast Cancer
35-40     | 021           | Ovarian Cancer
35-40     | 021           | Lung Cancer
```

#### L-Diversity
Extension of k-anonymity ensuring diversity in sensitive attributes.

#### Differential Privacy
Add calibrated noise to query results to protect individual privacy while preserving aggregate statistics.

### 4.2 Genomic Data De-identification

**Challenges:**
- Genomic data is inherently identifying
- Cannot be truly anonymized (only de-identified)
- Risk of re-identification through:
  - Genealogy databases (GEDmatch, Ancestry.com)
  - Public genomic databases
  - Linkage attacks

**Best Practices:**

1. **Controlled Access**
   - Store raw genomic data in controlled-access repositories (dbGaP)
   - Require Data Use Agreement (DUA) for access
   - Prohibit re-identification attempts in DUA

2. **Data Aggregation**
   - Share aggregate statistics instead of individual genotypes
   - Use summary statistics (allele frequencies, effect sizes)
   - Implement minimum cell size rules (n ≥ 5)

3. **Beacon Attack Mitigation**
   - Avoid simple presence/absence queries
   - Implement query budgets
   - Add differential privacy noise

4. **Federated Analysis**
   - Compute on data without transferring
   - Share only analysis results, not raw data
   - Use secure multi-party computation

### 4.3 De-identification Validation

**Before sharing data, validate:**

```python
def validate_deidentification(data, level="Safe Harbor"):
    """Validate data is properly de-identified."""

    checks = {
        "direct_identifiers": check_direct_identifiers(data),
        "quasi_identifiers": check_quasi_identifiers(data),
        "rare_values": check_rare_values(data),
        "linkage_risk": assess_linkage_risk(data),
    }

    if level == "Safe Harbor":
        # All direct identifiers must be removed
        if not checks["direct_identifiers"]:
            return False, "Direct identifiers found"

    elif level == "K-Anonymous":
        # Check k-anonymity
        k_value = calculate_k_anonymity(data)
        if k_value < 3:
            return False, f"K-anonymity too low: {k_value} < 3"

    return True, "De-identification validated"
```

---

## 5. Data Retention & Deletion

### 5.1 Retention Periods

| Data Type | Minimum Retention | Maximum Retention | Authority |
|-----------|-------------------|-------------------|-----------|
| Clinical trial data | 7 years after completion | 25 years | FDA 21 CFR 312.62 |
| NIH-funded research data | 7 years after final report | No limit | NIH Grants Policy |
| Published research data | Duration of publication | No limit | Journal policies |
| Participant consent forms | 7 years after study end | Permanent | IRB requirements |
| Human genomic data | 7 years (minimum) | No limit | NIH GDS Policy |
| Patient medical records | State-dependent (5-10 years) | Permanent | State law |

### 5.2 Retention Policy

**General Policy:**
- Retain all research data for **minimum 7 years** after study completion
- Retain participant consent forms **permanently** or per IRB guidance
- Retain audit logs for **minimum 7 years**
- Document retention schedule in Data Management Plan

**Extended Retention:**
- Data supporting published findings: Retain indefinitely
- Data with long-term research value: Deposit in public archive
- Participant requests data retention: Honor requests where feasible

### 5.3 Secure Deletion

**When retention period expires:**

1. **Verify Authorization**
   - Confirm retention period has expired
   - Obtain approval from PI and IRB (if required)
   - Document deletion justification

2. **Deletion Procedure**
   ```bash
   # Secure file deletion (Linux)
   shred -vfz -n 3 sensitive_file.csv

   # Secure directory deletion
   find /path/to/data -type f -exec shred -vfz -n 3 {} \;

   # Verify deletion
   ls -la /path/to/data
   ```

3. **Multi-Copy Deletion**
   - Delete all copies (primary storage, backups, archives, local copies)
   - Wipe cloud storage buckets
   - Destroy physical media (CDs, USB drives)

4. **Documentation**
   - Record deletion date
   - Record who authorized deletion
   - Record what was deleted
   - Store deletion certificate

**Secure Media Disposal:**
- Hard drives: Physical destruction or DoD 5220.22-M wipe (7 passes)
- SSDs: ATA Secure Erase or physical destruction
- USB drives: Physical destruction
- CDs/DVDs: Physical shredding
- Paper records: Cross-cut shredding or incineration

---

## 6. Audit Trails & Monitoring

### 6.1 Audit Logging Requirements

**Log all access to Confidential and Restricted data:**

| Event | Log Details | Retention |
|-------|-------------|-----------|
| Data access | User ID, timestamp, data accessed, purpose | 7 years |
| Data modification | User ID, timestamp, what changed, reason | 7 years |
| Data export | User ID, timestamp, data exported, destination | 7 years |
| Access denial | User ID, timestamp, data requested, reason denied | 7 years |
| System configuration | Admin ID, timestamp, configuration changed | 7 years |
| Login/logout | User ID, timestamp, IP address, success/failure | 1 year |

### 6.2 Audit Log Format

```json
{
  "event_id": "evt_20250115_103052_a7b3c",
  "timestamp": "2025-01-15T10:30:52.123Z",
  "event_type": "data_access",
  "user_id": "researcher_042",
  "user_role": "researcher",
  "ip_address": "10.1.2.45",
  "data_accessed": {
    "dataset": "patient_one_multiomics",
    "classification": "confidential",
    "files": ["rna_expression.csv", "protein_abundance.csv"],
    "num_records": 150,
    "num_patients": 1
  },
  "purpose": "Multi-omics integration analysis",
  "authorization": {
    "irb_protocol": "IRB-2024-001",
    "dua_id": "DUA-2024-015",
    "approved_by": "pi_smith"
  },
  "action": "read",
  "result": "success"
}
```

### 6.3 Monitoring & Alerts

**Automated alerts for:**
- ⚠️ Access to Restricted data outside business hours
- ⚠️ Bulk data downloads (>1000 records)
- ⚠️ Failed authentication attempts (>3 in 5 minutes)
- ⚠️ Access from unusual locations/IP addresses
- ⚠️ Data exports to external systems
- ⚠️ Attempted access to unauthorized data
- 🚨 Suspected data breach

**Monthly Reviews:**
- Access patterns for anomalies
- Compliance with approved purposes
- User access recertification
- Unused accounts (disable after 90 days inactivity)

---

## 7. Data Sharing & Collaboration

### 7.1 Data Sharing Tiers

| Tier | Sharing Method | Requirements | Use Case |
|------|----------------|--------------|----------|
| **Tier 1: Public** | Open repository | None | Published aggregate data, summary statistics |
| **Tier 2: Registered Access** | Controlled repository | Registration, click-through DUA | De-identified genomic data |
| **Tier 3: Approved Access** | Secure portal | IRB approval, signed DUA, HIPAA training | Limited datasets with dates/ZIP |
| **Tier 4: Collaborative** | Federated analysis | IRB approval, DUA, secure data enclave | Identifiable data, never leaves institution |

### 7.2 Data Use Agreement (DUA)

**Required for sharing Confidential or Restricted data.**

**Key DUA provisions:**
- Permitted uses (specific research purposes)
- Prohibited uses (re-identification attempts, commercial use)
- Acknowledgment of data source
- Publication review (if applicable)
- Security requirements (encryption, access controls)
- Breach notification requirements
- Termination and data destruction clauses
- Indemnification

**DUA Template:** See `templates/DATA_USE_AGREEMENT_TEMPLATE.md`

### 7.3 Public Data Repositories

**Recommended repositories:**

| Data Type | Repository | URL | Access |
|-----------|------------|-----|--------|
| Genomic (controlled) | dbGaP | https://www.ncbi.nlm.nih.gov/gap/ | Registered |
| Genomic (open) | GEO | https://www.ncbi.nlm.nih.gov/geo/ | Public |
| Cancer genomics | TCGA | https://portal.gdc.cancer.gov/ | Open + Controlled |
| Proteomics | PRIDE | https://www.ebi.ac.uk/pride/ | Public |
| Metabolomics | Metabolomics Workbench | https://www.metabolomicsworkbench.org/ | Public |
| Spatial transcriptomics | SpatialDB | (varies by platform) | Varies |
| Clinical trials | ClinicalTrials.gov | https://clinicaltrials.gov/ | Public |

---

## 8. Incident Response

### 8.1 Data Breach Definition

A data breach occurs when:
- Unauthorized person gains access to Confidential or Restricted data
- Authorized person accesses data for unauthorized purposes
- Data is disclosed to unauthorized parties
- Data security is compromised (lost laptop, stolen hard drive)
- Accidental disclosure of PHI or identifiable data

### 8.2 Incident Response Plan

**Step 1: DETECT (0-1 hour)**
- Identify suspected breach
- Document initial observations
- Preserve evidence (logs, systems, communications)

**Step 2: CONTAIN (1-4 hours)**
- Isolate affected systems
- Revoke compromised credentials
- Block unauthorized access
- Prevent further disclosure

**Step 3: ASSESS (4-24 hours)**
- Determine scope of breach
  - What data was compromised?
  - How many individuals affected?
  - Was data encrypted?
  - What was accessed/disclosed?
- Classify severity (Low/Medium/High/Critical)

**Step 4: NOTIFY (24-72 hours)**

**Notification requirements:**

| Affected Parties | Timeframe | Method |
|------------------|-----------|--------|
| PI / Data Steward | Immediate | Phone + Email |
| Institutional Security | Within 24 hours | Official incident report |
| IRB | Within 48 hours | Incident report form |
| Affected individuals | Without unreasonable delay | Written notice |
| HHS Office for Civil Rights (if HIPAA) | Within 60 days | Online breach portal |
| EU Data Protection Authority (if GDPR) | Within 72 hours | Official notification |

**Step 5: REMEDIATE (72 hours - ongoing)**
- Implement security improvements
- Review and update policies
- Provide additional training
- Monitor for secondary incidents

**Step 6: DOCUMENT (ongoing)**
- Incident report with timeline
- Root cause analysis
- Remediation actions taken
- Lessons learned
- Update incident response plan

### 8.3 Breach Severity Classification

| Severity | Criteria | Response |
|----------|----------|----------|
| **Critical** | PHI of >500 individuals, high re-identification risk | Full incident response, HHS notification, public notice |
| **High** | PHI of <500 individuals, moderate re-identification risk | Full incident response, IRB notification, individual notice |
| **Medium** | De-identified data with some re-identification risk | Incident investigation, IRB notification, policy review |
| **Low** | Internal-only data, no re-identification risk | Document incident, review access controls |

---

## 9. Research Ethics

### 9.1 Ethical Principles (Belmont Report)

1. **Respect for Persons**
   - Obtain informed consent
   - Protect autonomy of participants
   - Provide additional protections for vulnerable populations

2. **Beneficence**
   - Maximize benefits, minimize harms
   - Balance risks and benefits
   - Protect participant well-being

3. **Justice**
   - Fair distribution of research burdens and benefits
   - Avoid exploitation of vulnerable populations
   - Ensure equitable access to research results

### 9.2 Incidental Findings

**Definition:** Findings with potential health or reproductive significance discovered during research that are beyond the aims of the study.

**Incidental findings in genomics:**
- Pathogenic variants in medically actionable genes (BRCA1/2, TP53, etc.)
- Carrier status for recessive conditions
- Pharmacogenomic variants
- Non-paternity
- Unexpected ancestry

**Recommended approach (ACMG guidelines):**

1. **Pre-Consent Discussion**
   - Inform participants of possibility of incidental findings
   - Offer choice: receive/not receive certain findings
   - Clarify which findings will be disclosed
   - Explain limitations (research-grade, not clinical)

2. **Analysis**
   - Analyze for ACMG Secondary Findings list (73 genes)
   - Confirm variants in CLIA-certified lab (if disclosure planned)
   - Assess pathogenicity using ClinVar, ACMG/AMP guidelines

3. **Disclosure**
   - Genetic counseling before disclosure
   - Provide written report with recommendations
   - Refer to appropriate medical specialist
   - Document disclosure in research records

4. **Non-Disclosure**
   - Respect participant preference not to know
   - Do not disclose findings participant declined
   - Maintain confidentiality

### 9.3 Return of Research Results

**General research results:**
- Provide aggregate findings to participants (optional)
- Publish results in peer-reviewed journals
- Post pre-prints/data in public repositories
- Present at conferences

**Individual research results:**
- ⚠️ Generally NOT returned (research-grade, not validated)
- Exception: Clinically validated, medically actionable findings
- Must be confirmed in CLIA-certified lab before return
- Provide through genetic counseling

---

## 10. Training & Compliance

### 10.1 Required Training

**Before accessing Confidential or Restricted data:**

| Training | Audience | Frequency | Provider |
|----------|----------|-----------|----------|
| HIPAA Privacy & Security | All users | Annually | Institutional compliance |
| Human Subjects Research (CITI) | All researchers | Every 3 years | CITI Program |
| Responsible Conduct of Research (RCR) | All researchers | Every 4 years | Institution |
| Data Security Awareness | All users | Annually | IT Security |
| IRB-specific protocols | Study team | Per study | IRB office |

### 10.2 Compliance Monitoring

**Annual compliance review:**
- ✅ Access control audit (remove unnecessary access)
- ✅ User training verification (renew expired training)
- ✅ Data classification review (reclassify as needed)
- ✅ Encryption verification (confirm all sensitive data encrypted)
- ✅ Backup testing (verify restoration procedures)
- ✅ Incident response drill (test procedures)
- ✅ Policy updates (reflect regulatory changes)

### 10.3 Non-Compliance Consequences

**Violations may result in:**
- Immediate access revocation
- Mandatory retraining
- Formal warning
- IRB study suspension
- Institutional disciplinary action
- Loss of research privileges
- Reporting to NIH Office of Research Integrity
- Federal sanctions (loss of funding)
- Civil penalties (HIPAA: up to $1.5M per year)
- Criminal penalties (HIPAA: up to $250K fine, 10 years prison)

---

## 11. Roles & Responsibilities

### 11.1 Principal Investigator (PI)

**Responsibilities:**
- Overall responsibility for data governance
- Obtain IRB approval and maintain compliance
- Ensure team completes required training
- Approve data access requests
- Oversee data security and privacy
- Respond to data breaches
- Maintain study documentation
- Ensure proper data retention and destruction

### 11.2 Data Steward

**Responsibilities:**
- Classify data per governance policy
- Implement access controls
- Monitor data usage
- Conduct quarterly access audits
- Manage Data Use Agreements
- Coordinate data sharing requests
- Maintain data inventory
- Report compliance metrics

### 11.3 System Administrator

**Responsibilities:**
- Configure and maintain security controls
- Monitor system logs for suspicious activity
- Perform security updates and patches
- Implement encryption
- Manage backups
- Respond to security incidents
- Conduct vulnerability assessments
- Document system configurations

### 11.4 Researcher

**Responsibilities:**
- Complete required training
- Follow data governance policies
- Access only authorized data
- Use data only for approved purposes
- Protect credentials (passwords, API keys)
- Report suspected breaches immediately
- De-identify data before sharing
- Acknowledge data sources in publications

### 11.5 Institutional Review Board (IRB)

**Responsibilities:**
- Review and approve research protocols
- Assess risks to participants
- Review informed consent documents
- Monitor ongoing research for compliance
- Investigate protocol violations
- Approve protocol modifications
- Conduct continuing review

### 11.6 Institutional Compliance Office

**Responsibilities:**
- Develop and maintain policies
- Provide HIPAA training
- Conduct compliance audits
- Investigate reported violations
- Report breaches to authorities
- Coordinate with legal counsel
- Manage institutional risk

---

## 12. Implementation Checklist

### For New Research Projects

- [ ] **Step 1: IRB Approval**
  - [ ] Submit IRB protocol
  - [ ] Obtain IRB approval letter
  - [ ] Register at ClinicalTrials.gov (if applicable)

- [ ] **Step 2: Data Management Plan**
  - [ ] Define data types to be collected
  - [ ] Classify data (Public/Internal/Confidential/Restricted)
  - [ ] Identify data sources
  - [ ] Plan for de-identification
  - [ ] Define retention periods
  - [ ] Plan for data sharing

- [ ] **Step 3: Security Setup**
  - [ ] Provision secure storage
  - [ ] Enable encryption (TLS, AES-256)
  - [ ] Configure access controls
  - [ ] Enable audit logging
  - [ ] Set up backups

- [ ] **Step 4: Team Training**
  - [ ] All team members complete HIPAA training
  - [ ] All team members complete CITI training
  - [ ] Review data governance policies
  - [ ] Sign Data Use Agreements

- [ ] **Step 5: Data Collection**
  - [ ] Obtain informed consent
  - [ ] Assign study IDs
  - [ ] Remove direct identifiers
  - [ ] Store securely
  - [ ] Log data provenance

- [ ] **Step 6: Ongoing Compliance**
  - [ ] Quarterly access audits
  - [ ] Annual compliance review
  - [ ] Respond to participant requests
  - [ ] Report protocol deviations to IRB

- [ ] **Step 7: Study Completion**
  - [ ] Final IRB report
  - [ ] Publish results
  - [ ] Share data per DMP
  - [ ] Archive data
  - [ ] Schedule data destruction (if applicable)

---

## 13. Templates & Resources

### 13.1 Templates

Available in `templates/` directory:
- Data Use Agreement (DUA) template
- Informed Consent template (genomic research)
- Data Management Plan template
- Incident Response Report template
- De-identification Checklist
- Access Request Form

### 13.2 External Resources

**Regulatory Guidance:**
- HHS HIPAA: https://www.hhs.gov/hipaa/
- FDA Regulations: https://www.fda.gov/regulatory-information/
- NIH Genomic Data Sharing: https://sharing.nih.gov/genomic-data-sharing-policy
- GDPR: https://gdpr.eu/

**Training:**
- CITI Program: https://about.citiprogram.org/
- NIH Data Sharing Training: https://sharing.nih.gov/data-management-and-sharing-policy/training

**Tools:**
- ARX Data Anonymization Tool: https://arx.deidentifier.org/
- sdcMicro (R package): https://cran.r-project.org/package=sdcMicro
- OpenMRS De-identification Module: https://openmrs.org/

---

## 14. Contact Information

**For Data Governance Questions:**
- Principal Investigator: [Contact Info]
- Data Steward: [Contact Info]
- Institutional IRB: [Contact Info]
- Compliance Office: [Contact Info]

**For Security Incidents:**
- Security Operations Center (SOC): [24/7 Contact]
- IT Security: [Contact Info]
- Institutional CISO: [Contact Info]

**For Participant Questions:**
- Study Coordinator: [Contact Info]
- Genetic Counselor: [Contact Info]
- Patient Advocate: [Contact Info]

---

## 15. Document Control

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2025-01-15 | Precision Medicine MCP Team | Initial release |

**Next Review Date:** 2026-01-15 (annual review)

**Approval:**
- [ ] Principal Investigator
- [ ] Institutional Compliance Officer
- [ ] IRB Chair
- [ ] IT Security Officer

---

**This document establishes the data governance framework for Precision Medicine MCP. All users must read, understand, and comply with these policies.**
