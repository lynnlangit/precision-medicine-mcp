# For Educators & Students

**This section is for:** Professors, instructors, students, and anyone teaching or learning precision medicine, bioinformatics, or computational biology.

---

## What You Can Accomplish Here

- ✅ **Teach precision medicine** with real-world case studies (100% synthetic data)
- ✅ **Demonstrate bioinformatics workflows** (differential expression, pathway enrichment, spatial analysis)
- ✅ **Practice AI-assisted research** (natural language → bioinformatics pipelines)
- ✅ **Learn cost-effectively** (DRY_RUN mode ~$0.32 per complete analysis)
- ✅ **Publish educational materials** (no patient privacy concerns)

---

## Why This Is Perfect for Teaching

### 100% Synthetic Data
- **No patient privacy concerns** - Safe for classroom use, public demos, publications
- **No IRB approval needed** - Freely shareable with students
- **Realistic but synthetic** - Clinically plausible scenarios without real PHI
- **PatientOne dataset** - Stage IV ovarian cancer case with 5 modalities

### Low Cost for Education
- **DRY_RUN mode:** ~$0.32 per complete analysis (synthetic data, no API costs)
- **Classroom pricing:** 20 students × $0.32 = **$6.40 per class session**
- **Semester course:** 15 sessions × 20 students = **$96 total**
- **Compare to:** Commercial platforms ($500-2,000 per student)

### Comprehensive Coverage
All major bioinformatics domains in one platform:
- Clinical data (FHIR standards)
- Genomics (VCF variant calling)
- Multi-omics (RNA, protein, phospho integration)
- Spatial transcriptomics (Visium analysis)
- Imaging (H&E, multiplex IF)
- AI orchestration (natural language queries)
- Cloud deployment (GCP, Docker, serverless)

### Well-Documented
- Step-by-step tutorials with expected outputs
- Example prompts for common analyses
- Architecture diagrams showing data flow
- Statistical methods fully explained
- Reproducible workflows

---

## Quick Start for Instructors

### 1. Set Up Class Environment (30-60 minutes)

**Option A: Cloud-Based (Recommended for large classes)**
```bash
# Deploy to GCP with DRY_RUN mode (synthetic data)
./scripts/deployment/deploy_to_gcp.sh --development

# Students access via web interface (Streamlit UI)
# No local installation needed
```

**Option B: Local Installation (Best for small classes/workshops)**
```bash
# Students install Claude Desktop locally
# Follow installation guide: docs/getting-started/installation.md
# Time: 10-15 minutes per student
```

### 2. Assign PatientOne Case Study (25-35 minutes)

**Case:** PAT001-OVC-2025 - 58-year-old female, Stage IV HGSOC, platinum-resistant

**Learning objectives:**
- Interpret clinical FHIR data
- Analyze genomic variants (TP53, BRCA1)
- Integrate multi-omics datasets
- Perform spatial transcriptomics analysis
- Generate treatment recommendations

**Student deliverable:**
- Written report with:
  - Clinical summary
  - Genomic findings
  - Pathway analysis results
  - Treatment recommendations with rationale
  - Visualizations (plots, heatmaps)

**See:** [PatientOne Guide](../test-docs/patient-one-scenario/README.md)

### 3. Assess Student Work

**Rubric categories:**
- Clinical interpretation (20%) - Correct understanding of FHIR data
- Genomic analysis (20%) - Variant identification and interpretation
- Pathway analysis (20%) - Correct statistical methods and interpretation
- Integration (20%) - Synthesis across modalities
- Presentation (20%) - Clear visualizations and reporting

**Example outputs:** [PatientOne Results](../test-docs/patient-one-scenario/architecture/outputs.md)

---

## Quick Start for Students

### 1. Install & Setup (10-15 minutes)

Follow: [Installation Guide](../getting-started/installation.md)

**What you'll install:**
- Claude Desktop (local) OR access Streamlit UI (cloud)
- MCP servers (configured by instructor)

### 2. Try Your First Analysis (5 minutes)

**Prompt:**
```
What clinical data is available for PatientOne (PAT001-OVC-2025)?
```

**Expected output:**
- Demographics: 58-year-old female
- Diagnosis: Stage IV High-Grade Serous Ovarian Carcinoma (HGSOC)
- Treatment history: Platinum-based chemotherapy (carboplatin + paclitaxel)
- Biomarkers: Elevated CA-125
- Genetic: BRCA1 germline variant

### 3. Complete PatientOne Workflow (25-35 minutes)

Follow: [PatientOne Guide](../test-docs/patient-one-scenario/README.md)

**What you'll learn:**
- Load and interpret clinical data
- Analyze genomic variants
- Integrate multi-omics datasets
- Perform spatial pathway enrichment
- Generate treatment recommendations

**Cost:** ~$0.32 (using DRY_RUN mode with synthetic data)

---

## Educational Topics Covered

### Bioinformatics Fundamentals
- **Data formats:** FHIR, VCF, CSV matrices, Visium spatial format
- **Quality control:** VCF validation, expression matrix QC
- **Normalization:** Batch correction, log transformation
- **Visualization:** Heatmaps, volcano plots, spatial maps

### Statistical Methods
- **Differential expression:** Mann-Whitney U test
- **Multiple testing:** Benjamini-Hochberg FDR correction
- **Pathway enrichment:** Fisher's exact test
- **Spatial statistics:** Moran's I spatial autocorrelation
- **Meta-analysis:** Stouffer's Z-score method
- **Effect sizes:** Log2 fold change, odds ratios

### Precision Medicine Workflows
- **Clinical → Genomic integration:** FHIR + VCF analysis
- **Multi-modal analysis:** RNA + Protein + Spatial
- **Treatment matching:** Pathway → Drug recommendations
- **Biomarker discovery:** Differential expression → Validation
- **Tumor microenvironment:** Spatial transcriptomics + Imaging

### AI & Orchestration
- **Natural language queries:** Prompting best practices
- **Tool composition:** Chaining multiple bioinformatics tools
- **Error handling:** Interpreting and fixing errors
- **Workflow design:** Building reproducible pipelines

### Cloud & DevOps
- **Serverless deployment:** GCP Cloud Run
- **Containerization:** Docker basics
- **API design:** RESTful APIs, MCP protocol
- **Cost optimization:** DRY_RUN mode, caching strategies

---

## Classroom Activities

### Activity 1: PatientOne Case Study (1-2 class sessions)

**Format:** Individual or group (2-3 students)

**Time:** 25-35 minutes analysis + 30-60 minutes report writing

**Assignment:**
1. Analyze PatientOne using platform
2. Write clinical report with treatment recommendations
3. Present findings to class (5-minute presentation)

**Learning outcomes:**
- Interpret multi-modal cancer data
- Apply statistical methods correctly
- Communicate findings to non-technical audience

### Activity 2: Compare Real vs. Synthetic Data (1 class session)

**Format:** Class discussion

**Time:** 50-75 minutes

**Assignment:**
1. Run PatientOne in DRY_RUN mode (synthetic data)
2. Review real data examples (instructor provides)
3. Discuss differences and limitations

**Learning outcomes:**
- Understand data synthesis methods
- Recognize limitations of synthetic data
- Appreciate importance of real data validation

### Activity 3: Design Custom Workflow (1-2 class sessions)

**Format:** Group project (3-4 students)

**Time:** 2-3 hours

**Assignment:**
1. Choose research question (e.g., "Identify biomarkers for platinum resistance")
2. Design workflow using available tools
3. Test with PatientOne data
4. Present workflow design and results

**Learning outcomes:**
- Plan bioinformatics analyses
- Chain tools into workflows
- Document reproducible methods

### Activity 4: Extend Platform (Advanced, 2-3 class sessions)

**Format:** Individual project

**Time:** 4-8 hours

**Assignment:**
1. Identify missing functionality (e.g., "Add metabolomics analysis")
2. Implement new server using boilerplate template
3. Write tests (≥35% coverage)
4. Document new tools

**Learning outcomes:**
- Software engineering best practices
- API design and implementation
- Testing and documentation

**See:** [Add New Modality Server Guide](../for-developers/ADD_NEW_MODALITY_SERVER.md)

---

## Course Integration Ideas

### Undergraduate Courses

**Bioinformatics 101:**
- Module 1: Introduction to genomic data formats (VCF, FASTA)
- Module 2: Differential expression analysis
- Module 3: Pathway enrichment
- **Use PatientOne:** As running example throughout course

**Precision Medicine:**
- Week 1-2: Clinical data standards (FHIR)
- Week 3-4: Genomic variant interpretation
- Week 5-6: Multi-omics integration
- Week 7-8: Spatial transcriptomics
- **Final project:** Complete PatientOne analysis with report

**Data Science for Biology:**
- Topic: Multi-modal data integration
- Dataset: PatientOne (5 modalities)
- Methods: PCA, clustering, meta-analysis
- **Lab:** Integrate PatientOne RNA + Protein + Spatial data

### Graduate Courses

**Computational Biology:**
- Advanced topic: AI orchestration in bioinformatics
- Assignment: Design and implement custom workflow
- **Capstone:** Extend platform with new server (metabolomics, single-cell, etc.)

**Translational Research:**
- Case study: Ovarian cancer precision medicine
- Data: PatientOne multi-modal dataset
- Outcome: Treatment recommendation report
- **Presentation:** Molecular tumor board simulation

**Cloud Computing for Bioinformatics:**
- Architecture: MCP protocol, serverless computing
- Deployment: GCP Cloud Run, Docker
- Cost optimization: DRY_RUN mode, caching
- **Project:** Deploy and scale bioinformatics pipeline

---

## Assessment Ideas

### Formative Assessments

**1. Quick Checks (5-10 minutes each)**
- "What is a VCF file? What information does it contain?"
- "Explain the difference between p-value and FDR-corrected p-value"
- "Why do we use batch correction in spatial transcriptomics?"

**2. Prompt Engineering Exercises (15-20 minutes)**
- "Write a prompt to identify upregulated genes in tumor vs. normal"
- "Critique this prompt: [bad example]. How would you improve it?"

**3. Error Debugging (20-30 minutes)**
- Provide error message, ask students to diagnose and fix
- Example: "FileNotFoundError: /data/patient-001/spatial/data.csv"

### Summative Assessments

**1. PatientOne Analysis Report (Individual, 2-3 hours)**

**Grading rubric (100 points):**
- Clinical summary (15 pts)
- Genomic analysis (20 pts)
- Multi-omics integration (20 pts)
- Spatial analysis (20 pts)
- Treatment recommendations (15 pts)
- Visualizations (10 pts)

**2. Literature Review + Replication (Group, 4-6 hours)**

**Assignment:**
- Choose paper using spatial transcriptomics
- Replicate key analysis using platform
- Compare results (synthetic vs. published data)
- Write methods section for reproducibility

**3. Custom Workflow Design (Group, 8-12 hours)**

**Assignment:**
- Define research question
- Design 5-step workflow
- Implement and test with PatientOne
- Document workflow (README, example prompts)
- Present to class (10-minute presentation)

---

## Example Syllabi

### Undergraduate: Introduction to Precision Medicine (8 weeks)

| Week | Topic | Activity | PatientOne Module |
|------|-------|----------|-------------------|
| 1 | Clinical data standards | Explore FHIR resources | Clinical summary |
| 2 | Genomic variants | VCF analysis | Identify TP53, BRCA1 |
| 3 | RNA-seq basics | Differential expression | Upregulated genes |
| 4 | Pathway analysis | Fisher's exact test | Enriched pathways |
| 5 | Multi-omics integration | Stouffer meta-analysis | RNA + Protein |
| 6 | Spatial transcriptomics | Moran's I, cell deconvolution | Spatial patterns |
| 7 | Imaging & AI | Cell segmentation | H&E analysis |
| 8 | Final project | Complete PatientOne report | All modules |

**Learning outcomes:**
- Interpret multi-modal cancer data
- Apply statistical methods correctly
- Generate treatment recommendations
- Communicate findings effectively

### Graduate: Computational Methods in Cancer Biology (12 weeks)

| Week | Topic | Reading | Assignment |
|------|-------|---------|------------|
| 1-2 | AI orchestration in bioinformatics | MCP protocol spec | Setup environment |
| 3-4 | Statistical methods | Mann-Whitney, FDR, Fisher's | PatientOne differential expression |
| 5-6 | Multi-omics integration | HAllA paper | Integrate PatientOne modalities |
| 7-8 | Spatial transcriptomics | Visium tutorial | PatientOne spatial analysis |
| 9-10 | Cloud deployment | GCP Cloud Run docs | Deploy custom server |
| 11-12 | Final project | Literature review | Extend platform with new feature |

**Learning outcomes:**
- Design and implement bioinformatics workflows
- Apply advanced statistical methods
- Deploy scalable cloud infrastructure
- Contribute to open-source bioinformatics

---

## Teaching Resources

### Instructor Materials
- **[PatientOne Guide](../test-docs/patient-one-scenario/README.md)** - Complete walkthrough
- **[Example Outputs](../test-docs/patient-one-scenario/architecture/outputs.md)** - Expected results
- **[Statistical Methods](../for-researchers/ANALYSIS_WORKFLOWS.md)** - Detailed method explanations
- **[Cost Analysis](../../operations/cost-and-budget.md)** - Budgeting for classroom use

### Student Resources
- **[Installation Guide](../getting-started/installation.md)** - Setup instructions
- **[Quick Start Demo](../demos/NINETY_SECOND_PITCH.md)** - 90-second overview
- **[Architecture Overview](../for-developers/ARCHITECTURE.md)** - System design
- **[Troubleshooting Guide](../../hospital-deployment/RUNBOOKS/)** - Common issues

### Additional Materials
- **Lecture slides:** [To be added - community contributions welcome]
- **Video tutorials:** [To be added - community contributions welcome]
- **Problem sets:** [To be added - community contributions welcome]

---

## Instructor Support

### Getting Help
- **Technical questions:** GitHub Discussions
- **Teaching strategies:** [instructor-mailing-list placeholder]
- **Bug reports:** GitHub Issues
- **Feature requests:** GitHub Issues with "education" tag

### Contributing Educational Materials
We welcome contributions!
- Lecture slides
- Problem sets
- Video tutorials
- Assessment rubrics
- Course syllabi

**How to contribute:** See [CONTRIBUTING.md](../for-developers/CONTRIBUTING.md)

### Pilot Program
Interested in piloting this platform in your course?

**We offer:**
- Free GCP credits for classroom use
- Priority support during semester
- Featured on platform website
- Co-authorship on educational papers

**Contact:** [email placeholder]

---

## Frequently Asked Questions

### "Is this approved for classroom use?"
**A:** Yes! PatientOne data is 100% synthetic with no patient identifiers. No IRB approval needed for educational use.

### "How much does it cost for a class of 20 students?"
**A:** DRY_RUN mode (synthetic data): ~$6.40 per class session (20 students × $0.32)
Production mode (real data): Requires institutional GCP account and varies by usage.

### "Can students work on their own computers?"
**A:** Yes with Claude Desktop (local installation). OR use cloud-based Streamlit UI (no local install needed).

### "What if students find bugs or have questions?"
**A:** Students can:
1. Check documentation first
2. Ask in GitHub Discussions
3. Email instructor
4. Open GitHub Issue for bugs

### "Can I modify the PatientOne case study?"
**A:** Yes! Synthetic data is fully customizable. See [Data Generation Guide](../../data/patient-data/README.md) for how to create custom cases.

### "Is this suitable for high school students?"
**A:** Yes for advanced high school (AP Biology, AP Computer Science). Requires basic biology knowledge and comfort with technology. Instructor guidance recommended.

### "What if my institution blocks Claude Desktop?"
**A:** Use cloud-based Streamlit UI instead. Deployed on GCP, accessible via browser.

---

## Next Steps for Instructors

1. **Try the demo** → [PatientOne Guide](../test-docs/patient-one-scenario/README.md) (25-35 min)
2. **Review example outputs** → [Results](../test-docs/patient-one-scenario/architecture/outputs.md)
3. **Estimate costs** → [Cost Analysis](../../operations/cost-and-budget.md)
4. **Plan integration** → Choose activities and timeline for your course
5. **Request pilot** → Contact us for classroom support

## Next Steps for Students

1. **Set up environment** → [Installation Guide](../getting-started/installation.md) (10-15 min)
2. **Try first analysis** → Simple clinical query (5 min)
3. **Complete PatientOne** → Full workflow (25-35 min)
4. **Explore advanced topics** → Custom workflows, server development

---

**Related Resources:**
- 🔬 [Researcher Guide](../for-researchers/README.md) - For bioinformatics details
- 💻 [Developer Guide](../for-developers/README.md) - For extending the platform
- 💰 [Funding Information](../for-funders/README.md) - For grant applications
- 🏠 [Back to Main Documentation](../README.md)

---

**Last Updated:** 2026-01-14
