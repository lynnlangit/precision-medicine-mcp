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

| Framework | Status | Documentation |
|-----------|--------|---------------|
| HIPAA (US) | ✅ Compliant | **[See HIPAA Compliance Guide](hipaa.md)** |
| GDPR (EU) | ⚠️ Partial | See Section 1.1 below |
| Common Rule (US Research) | ✅ Compliant | See Section 1.2 below |
| FISMA (Federal) | ❌ Not Compliant | Not designed for federal systems |
| 21 CFR Part 11 (FDA) | ❌ Not Compliant | Not validated for clinical trials |

**⚠️ IMPORTANT:** This system is designed for **research use only** and is **NOT** compliant with clinical diagnostic standards (CLIA, CAP).

---

## 1. Regulatory Compliance

### 1.1 GDPR (General Data Protection Regulation)

**Applicability:** Processing personal data of EU residents.

**📚 For HIPAA compliance, see: [HIPAA Compliance Guide](hipaa.md)**

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

    def withdraw_consent(self, reason: str = None):
        """Record consent withdrawal."""
        self.withdrawals.append({
            'date': datetime.now(),
            'reason': reason
        })

    def is_valid_for_purpose(self, purpose: str) -> bool:
        """Check if consent is valid for given purpose."""
        if self.withdrawals:
            return False
        return purpose in self.purposes
```

### 1.2 Common Rule (US Federal Research Regulations)

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

### 1.3 Institutional Policies

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

**📚 For detailed HIPAA de-identification requirements, see: [HIPAA Compliance Guide](hipaa.md#de-identification-validation)**

Removes all 18 HIPAA identifiers.

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
- Risk of re-identification through genealogy databases, linkage attacks

**Best Practices:**

1. **Controlled Access** - Store in controlled-access repositories (dbGaP), require DUA
2. **Data Aggregation** - Share aggregate statistics instead of individual genotypes
3. **Beacon Attack Mitigation** - Avoid simple presence/absence queries, add noise
4. **Federated Analysis** - Compute on data without transferring, share only results

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

### 5.2 Secure Deletion

**When retention period expires:**

1. **Verify Authorization** - Confirm period expired, obtain PI/IRB approval
2. **Deletion Procedure**
   ```bash
   # Secure file deletion
   shred -vfz -n 3 sensitive_file.csv
   ```
3. **Multi-Copy Deletion** - Delete all copies (primary, backups, archives, local)
4. **Documentation** - Record deletion date, authorization, what was deleted

**Secure Media Disposal:**
- Hard drives: Physical destruction or DoD 5220.22-M wipe (7 passes)
- SSDs: ATA Secure Erase or physical destruction
- USB drives: Physical destruction
- CDs/DVDs: Physical shredding
- Paper records: Cross-cut shredding or incineration

---

## 6. Audit Trails & Monitoring

**📚 For detailed HIPAA audit requirements, see: [HIPAA Compliance Guide](hipaa.md#audit-controls)**

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

### 6.2 Monitoring & Alerts

**Automated alerts for:**
- ⚠️ Access to Restricted data outside business hours
- ⚠️ Bulk data downloads (>1000 records)
- ⚠️ Failed authentication attempts (>3 in 5 minutes)
- ⚠️ Access from unusual locations/IP addresses
- 🚨 Suspected data breach

---

## 7. Data Sharing & Collaboration

### 7.1 Data Sharing Tiers

| Tier | Sharing Method | Requirements | Use Case |
|------|----------------|--------------|----------|
| **Tier 1: Public** | Open repository | None | Published aggregate data |
| **Tier 2: Registered Access** | Controlled repository | Registration, click-through DUA | De-identified genomic data |
| **Tier 3: Approved Access** | Secure portal | IRB approval, signed DUA, HIPAA training | Limited datasets with dates/ZIP |
| **Tier 4: Collaborative** | Federated analysis | IRB approval, DUA, secure enclave | Identifiable data, never leaves institution |

### 7.2 Data Use Agreement (DUA)

**Required for sharing Confidential or Restricted data.**

**Key DUA provisions:**
- Permitted uses (specific research purposes)
- Prohibited uses (re-identification attempts)
- Security requirements
- Breach notification requirements
- Data destruction clauses

### 7.3 Public Data Repositories

| Data Type | Repository | URL | Access |
|-----------|------------|-----|--------|
| Genomic (controlled) | dbGaP | https://www.ncbi.nlm.nih.gov/gap/ | Registered |
| Genomic (open) | GEO | https://www.ncbi.nlm.nih.gov/geo/ | Public |
| Cancer genomics | TCGA | https://portal.gdc.cancer.gov/ | Open + Controlled |
| Proteomics | PRIDE | https://www.ebi.ac.uk/pride/ | Public |

---

## 8. Incident Response

**📚 For detailed HIPAA incident response procedures, see: [HIPAA Compliance Guide](hipaa.md#incident-response)**

### 8.1 Data Breach Definition

A data breach occurs when:
- Unauthorized person gains access to Confidential or Restricted data
- Data is disclosed to unauthorized parties
- Accidental disclosure of PHI or identifiable data

### 8.2 Incident Response Plan

**Step 1: DETECT (0-1 hour)** - Identify and document
**Step 2: CONTAIN (1-4 hours)** - Isolate affected systems
**Step 3: ASSESS (4-24 hours)** - Determine scope
**Step 4: NOTIFY (24-72 hours)** - Notify affected parties
**Step 5: REMEDIATE (72 hours+)** - Implement improvements
**Step 6: DOCUMENT (ongoing)** - Record lessons learned

---

## 9. Research Ethics

### 9.1 Ethical Principles (Belmont Report)

1. **Respect for Persons** - Informed consent, protect autonomy
2. **Beneficence** - Maximize benefits, minimize harms
3. **Justice** - Fair distribution of burdens and benefits

### 9.2 Incidental Findings

**Definition:** Findings with potential health significance discovered beyond study aims.

**Recommended approach:**
1. Pre-consent discussion
2. Analysis (ACMG Secondary Findings list)
3. Disclosure (genetic counseling)
4. Non-disclosure (respect preference not to know)

---

## 10. Training & Compliance

### 10.1 Required Training

**Before accessing Confidential or Restricted data:**

| Training | Frequency | Provider |
|----------|-----------|----------|
| HIPAA Privacy & Security | Annually | Institutional compliance |
| Human Subjects Research (CITI) | Every 3 years | CITI Program |
| Responsible Conduct of Research | Every 4 years | Institution |
| Data Security Awareness | Annually | IT Security |

---

## 11. Roles & Responsibilities

### Principal Investigator (PI)
- Overall responsibility for data governance
- Obtain IRB approval and maintain compliance
- Approve data access requests

### Data Steward
- Classify data per governance policy
- Implement access controls
- Conduct quarterly audits

### Researcher
- Complete required training
- Access only authorized data
- Report suspected breaches immediately

---

## Related Documentation

- **[HIPAA Compliance](hipaa.md)** - Detailed HIPAA compliance requirements
- **[Risk Assessment](../compliance/risk-assessment.md)** - Risk mitigation strategies
- **[Disclaimers](../compliance/disclaimers.md)** - Legal disclaimers

---

**Last Updated:** 2026-01-13

**This document establishes the data governance framework for Precision Medicine MCP. All users must read, understand, and comply with these policies.**
