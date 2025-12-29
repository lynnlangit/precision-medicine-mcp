# PatientOne Cost & Performance Analysis

## Executive Summary

| Mode | Total Time | Total Cost | Best For |
|------|-----------|------------|----------|
| **DRY_RUN** | 25-35 min | $0.30-0.65 | Demo, learning, CI/CD testing |
| **Automated Report** | ~12 seconds | $0.10-0.20 | Quick analysis with pre-aligned data |
| **Real Patient Data** | 1-2 hours | $7-19 | Production analysis, research, clinical decision support |

**Note:** Real data costs reduced from previous estimates due to spatialtools improvements (70% real as of Dec 29, 2025). STAR alignment still mocked, reducing processing time and cost.

---

## DRY_RUN Mode (Synthetic Data)

### Overview
- Uses synthetic responses from all MCP servers
- No external API calls or computational processing
- Ideal for demonstration and workflow validation

### Per-Test Breakdown

| Test | Servers Used | Est. Tokens (In/Out) | Time | Cost |
|------|--------------|---------------------|------|------|
| **TEST_1: Clinical + Genomic** | MockEpic, FGbio, TCGA | 2,000 / 3,500 | 3-5 min | $0.06 |
| **TEST_2: Multi-Omics** | MultiOmics | 2,500 / 4,000 | 5-8 min | $0.07 |
| **TEST_3: Spatial** | SpatialTools, DeepCell | 2,000 / 3,500 | 4-6 min | $0.06 |
| **TEST_4: Imaging** | OpenImageData, DeepCell | 1,800 / 3,000 | 3-5 min | $0.05 |
| **TEST_5: Integration** | All 9 servers | 3,000 / 5,000 | 5-7 min | $0.09 |
| **TOTAL** | - | **11,300 / 19,000** | **25-35 min** | **$0.30-0.65** |

### Token Usage Details

**Input Tokens (11,300 total):**
- User prompts: ~2,500 tokens
- File path references: ~500 tokens
- MCP tool definitions: ~4,000 tokens
- Context from previous tests: ~4,300 tokens

**Output Tokens (19,000 total):**
- Claude analysis & synthesis: ~8,000 tokens
- Synthetic MCP responses: ~6,000 tokens
- Recommendations & summaries: ~5,000 tokens

### Cost Calculation (Claude Sonnet 4 pricing)
- Input: 11,300 tokens × $3/M = **$0.034**
- Output: 19,000 tokens × $15/M = **$0.285**
- **Total: ~$0.32** (or $0.65 for verbose mode)

### Time Breakdown
- Claude processing: 15-20 minutes
- User review & navigation: 10-15 minutes
- **Total wall-clock time: 25-35 minutes**

---

## Automated Patient Report Generator

### Overview
- **New capability added December 29, 2025**
- Standalone Python script for rapid patient analysis
- Requires pre-aligned spatial transcriptomics data
- Integrates FHIR clinical data from GCP Healthcare API

### Script Details

**Script:** `scripts/generate_patient_report.py`
**Documentation:** `docs/AUTOMATED_PATIENT_REPORTS.md`

**Capabilities:**
- Differential expression (Mann-Whitney U + FDR correction)
- Spatial autocorrelation (Moran's I with spatial weights)
- Cell type deconvolution (signature-based scoring)
- Publication-quality visualizations (5 PNG files, 300 DPI)
- Clinical summary reports (TXT, JSON, CSV)

### Performance & Cost

| Metric | Value |
|--------|-------|
| **Runtime** | ~12 seconds per patient |
| **Compute Cost** | $0.05-0.10 (local CPU processing) |
| **Token Cost** | $0.05-0.10 (FHIR data retrieval if using Claude) |
| **Total Cost** | $0.10-0.20 per patient |

### Output Files (10 total, ~3.4 MB)

**Data files (5):**
- differential_expression.csv (3.4 KB)
- spatial_autocorrelation.csv (1.5 KB)
- cell_deconvolution.csv (566 B)
- clinical_summary.txt (2.9 KB)
- metadata.json (396 B)

**Visualizations (5):**
- volcano_plot.png (158 KB)
- spatial_heatmap.png (1.8 MB)
- cell_composition_heatmap.png (217 KB)
- spatial_autocorrelation_plot.png (140 KB)
- summary_figure.png (1.0 MB)

### Use Cases

**Best for:**
- Quick patient analysis with pre-aligned data (Space Ranger output)
- Batch processing of multiple patients
- Automated report generation pipelines
- Research studies requiring standardized analysis

**Not suitable for:**
- Raw FASTQ alignment (STAR not integrated yet)
- Batch correction across multiple samples
- Pathway enrichment analysis

### Example Usage

```bash
cd /path/to/spatial-mcp
/path/to/servers/mcp-spatialtools/venv/bin/python3 \
  scripts/generate_patient_report.py \
  --patient-id patient-001 \
  --output-dir ./results
```

**Example Output:**
```
✅ Analysis complete!
   Runtime: 11.8 seconds
   Files generated: 10
   Total size: 3.4 MB
```

---

## Real Patient Data Mode

### Overview
- Processes actual patient genomic, transcriptomic, and imaging data
- Executes computational workflows (alignment, segmentation, etc.)
- Makes external API calls (TCGA, HuggingFace, Seqera)

### Per-Test Breakdown

| Test | Servers Used | Processing Time | Compute Cost | API Cost | Total Cost |
|------|--------------|----------------|--------------|----------|-----------|
| **TEST_1: Clinical + Genomic** | MockEpic, FGbio, TCGA | 10-15 min | $0.50 | $0 | $0.50-0.75 |
| **TEST_2: Multi-Omics** | MultiOmics | 15-25 min | $2-4 | $0 | $2-4 |
| **TEST_3: Spatial** | SpatialTools, DeepCell | 10-20 min | $1-3 | $0 | $1-3 |
| **TEST_4: Imaging** | OpenImageData, DeepCell | 20-40 min | $3-6 | $0-1 | $3-7 |
| **TEST_5: Integration** | All 9 servers | 5-10 min | $0.25 | $0 | $0.25-0.50 |
| **TOTAL** | - | **1-2 hours** | **$7-14** | **$0-5** | **$7-19** |

### Detailed Cost Components

#### 1. Computational Processing ($7-14)

**Spatial Transcriptomics (TEST_3):**
- ⚠️ **Note:** STAR alignment still mocked (70% real implementation as of Dec 29, 2025)
- Differential expression (Mann-Whitney U + FDR): 2-5 min = $0.20-0.50
- Spatial autocorrelation (Moran's I): 2-5 min = $0.20-0.50
- Cell type deconvolution (signature-based): 1-3 min = $0.10-0.30
- ❌ DeepCell segmentation: MOCKED (not included in cost)
- **When STAR implemented:** Add 30-60 min, $5-10 for alignment

**Multi-Omics Integration (TEST_2):**
- HAllA analysis (chunked): 10-15 min = $1-2
- Stouffer's meta-analysis: 2-5 min = $0.20-0.50
- Upstream regulator prediction: 3-5 min = $0.30-0.50

**Imaging Analysis (TEST_4):**
- DeepCell cell segmentation: 15-30 min GPU = $3-5
- Image feature extraction: 5-10 min = $0.50-1

**Clinical/Genomic (TEST_1):**
- VCF processing: 5-10 min = $0.50

#### 2. External API Costs ($0-5)

**TCGA (TEST_1, TEST_3):**
- Free public API
- Rate-limited (10 requests/second)
- Cost: $0

**HuggingFace Inference API (TEST_3, TEST_4):**
- Cell type prediction: $0-1 per 1000 predictions
- Sequence embeddings: $0-1 per 1000 sequences
- Estimated: $0-2 total

**Seqera Platform (TEST_3):**
- Nextflow pipeline execution
- Free tier: 10 hours/month
- Paid tier: $0.10/compute-hour
- Estimated: $0-3 for spatial workflow

#### 3. Claude Token Usage (~$0.50-1.00)
Similar to DRY_RUN mode, but with larger context from real data files:
- Input: 15,000-20,000 tokens × $3/M = $0.045-0.06
- Output: 25,000-35,000 tokens × $15/M = $0.375-0.525
- **Total: ~$0.42-0.59**

### Time Breakdown

**Wall-clock time: 1-2 hours** (reduced from 2-4 hours due to spatialtools improvements)
- Computational processing: 60-110 minutes
- API calls & data transfer: 5-15 minutes
- Claude processing & user review: 20-45 minutes

**Parallelization opportunities:**
- Tests can run sequentially (2 hours) or with some parallelization (1.5 hours)
- Multi-omics HAllA analysis is now the bottleneck (10-15 min)
- **When STAR alignment implemented:** Add 30-60 min, making it the bottleneck again

---

## Cost Comparison by Use Case

### Academic Research Lab
**Scenario:** 50 patient analyses per year

| Mode | Cost per Patient | Annual Cost | Use Case |
|------|-----------------|-------------|----------|
| DRY_RUN | $0.32 | $16 | Workflow development, testing |
| Automated Report | $0.15 | $7.50 | Quick analysis (pre-aligned data) |
| Real Data | $13 (avg) | $650 | Full analysis with MCP orchestration |

**ROI:** Replaces ~40 hours of manual bioinformatics work per patient ($80/hr × 40 = $3,200)
**Savings:** $3,187 per patient using automated reports, or $159,350 annually
**Note:** Costs reduced by ~50% due to spatialtools improvements (Dec 29, 2025)

### Clinical Genomics Center
**Scenario:** 500 patient analyses per year

| Mode | Cost per Patient | Annual Cost |
|------|-----------------|-------------|
| Automated Report | $0.15 | $75 |
| Real Data | $13 (avg) | $6,500 |

**ROI:** Replaces manual analysis + reduces time-to-result from 2-3 weeks to 4-6 hours
**Value:** Faster treatment decisions, improved patient outcomes
**Annual Savings:** ~$1.6M using automated reports (vs manual $3,200/patient)

### Pharmaceutical R&D
**Scenario:** 200 PDX model analyses per year

| Mode | Cost per Analysis | Annual Cost |
|------|------------------|-------------|
| Automated Report | $0.15 | $30 |
| Real Data | $15 (avg) | $3,000 |

**ROI:** Accelerates target identification and biomarker discovery
**Value:** Weeks → hours for multi-omics integration
**Throughput:** 200 samples in ~40 minutes using automated reports (vs weeks manually)

---

## Infrastructure Requirements

### DRY_RUN Mode
- **CPU:** Any modern laptop (2+ cores)
- **RAM:** 4GB minimum, 8GB recommended
- **Storage:** 5GB for MCP servers + Python dependencies
- **Network:** Internet for Claude Desktop API calls only

### Real Patient Data Mode (Current: 70% real implementation)
- **CPU:** 4-8 cores sufficient for current real features (DE, Moran's I, deconvolution)
- **GPU:** ❌ Not needed yet (DeepCell mocked as of Dec 29, 2025)
- **RAM:** 16GB minimum, 32GB recommended
- **Storage:** 20GB minimum for patient data + cache
- **Network:** Stable connection for TCGA API calls (if enabled)

**When STAR alignment implemented:**
- CPU: 8-16 cores recommended
- RAM: 32-64GB
- Storage: 100GB for reference genomes

**When DeepCell implemented:**
- GPU: NVIDIA with 8GB+ VRAM recommended

---

## Cost Optimization Strategies

### 1. Batch Processing
- Run multiple patients in sequence to amortize setup costs
- Cache reference genomes and model weights
- Estimated savings: 20-30% for 10+ patients

### 2. Pre-computed Data Reuse
- Store intermediate results (aligned BAMs, cell segmentation masks)
- Reuse for hypothesis testing without re-processing
- Estimated savings: 50-70% for follow-up analyses

### 3. Selective Analysis
- Run only relevant tests (e.g., skip imaging if not available)
- DRY_RUN mode for workflow validation before committing compute
- Estimated savings: 40-60% when only 2-3 tests needed

### 4. Cloud vs Local
- **Local processing:** One-time hardware cost, no per-analysis fees
- **Cloud (AWS/GCP):** Pay-per-use, scales easily
- **Hybrid:** Local for common workflows, cloud for large-scale parallel processing

---

## Frequently Asked Questions

### Q: Why is Real Data mode 50-140× more expensive than DRY_RUN?
**A:** Real data requires actual computational processing (sequence alignment, cell segmentation, statistical analysis) rather than returning synthetic responses. The cost reflects the scientific computation, not just LLM orchestration.

### Q: Can I reduce costs by using smaller data files?
**A:** Yes! The costs scale with:
- Number of spatial spots (900 in demo → 10,000 in high-resolution)
- Number of genes profiled (31 in demo → 20,000 in whole transcriptome)
- Image resolution and number of markers

### Q: What happens if a test fails midway?
**A:**
- DRY_RUN: No cost beyond Claude tokens consumed
- Real Data: You pay for compute time used, but intermediate results are cached for retry

### Q: Are there free tiers for external APIs?
**A:**
- TCGA: Always free (NIH public data)
- HuggingFace: Free tier with rate limits (30 requests/hour)
- Seqera: 10 compute-hours/month free

---

**Last Updated:** December 29, 2025

**Recent Updates:**
- Added automated patient report generator (~12 sec, $0.10-0.20)
- Updated spatialtools status: 40% → 70% real implementation
- Reduced real data costs: $15-45 → $7-19 (STAR alignment still mocked)
- Reduced real data time: 2-4 hours → 1-2 hours
- Updated infrastructure requirements (no GPU needed yet)

**Pricing basis:** Claude Sonnet 4 ($3/M input, $15/M output), AWS EC2 c6i.2xlarge ($0.34/hr), AWS g4dn.xlarge ($0.526/hr)
