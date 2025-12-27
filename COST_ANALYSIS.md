# PatientOne Cost & Performance Analysis

## Executive Summary

| Mode | Total Time | Total Cost | Best For |
|------|-----------|------------|----------|
| **DRY_RUN** | 25-35 min | $0.30-0.65 | Demo, learning, CI/CD testing |
| **Real Patient Data** | 2-4 hours | $15-45 | Production analysis, research, clinical decision support |

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
| **TEST_3: Spatial** | SpatialTools, DeepCell | 45-90 min | $8-15 | $0-2 | $8-17 |
| **TEST_4: Imaging** | OpenImageData, DeepCell | 20-40 min | $3-6 | $0-1 | $3-7 |
| **TEST_5: Integration** | All 9 servers | 5-10 min | $0.25 | $0 | $0.25-0.50 |
| **TOTAL** | - | **2-4 hours** | **$14-30** | **$0-5** | **$15-45** |

### Detailed Cost Components

#### 1. Computational Processing ($14-30)

**Spatial Transcriptomics (TEST_3):**
- STAR alignment (900 spots): 30-60 min on 8-core CPU = $5-10
- UMI counting: 5-10 min = $0.50-1
- Spatial autocorrelation: 5-10 min = $0.50-1
- DeepCell segmentation: 10-20 min GPU = $2-3

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

**Wall-clock time: 2-4 hours**
- Computational processing: 95-180 minutes
- API calls & data transfer: 5-15 minutes
- Claude processing & user review: 20-45 minutes

**Parallelization opportunities:**
- Tests can run sequentially (4 hours) or with some parallelization (2.5 hours)
- STAR alignment is the bottleneck (30-60 min)

---

## Cost Comparison by Use Case

### Academic Research Lab
**Scenario:** 50 patient analyses per year

| Mode | Cost per Patient | Annual Cost | Use Case |
|------|-----------------|-------------|----------|
| DRY_RUN | $0.32 | $16 | Workflow development, testing |
| Real Data | $25 (avg) | $1,250 | Actual patient analysis |

**ROI:** Replaces ~40 hours of manual bioinformatics work per patient ($80/hr × 40 = $3,200)
**Savings:** $2,950 per patient or $147,500 annually

### Clinical Genomics Center
**Scenario:** 500 patient analyses per year

| Mode | Cost per Patient | Annual Cost |
|------|-----------------|-------------|
| Real Data | $25 (avg) | $12,500 |

**ROI:** Replaces manual analysis + reduces time-to-result from 2-3 weeks to 1 day
**Value:** Faster treatment decisions, improved patient outcomes

### Pharmaceutical R&D
**Scenario:** 200 PDX model analyses per year

| Mode | Cost per Analysis | Annual Cost |
|------|------------------|-------------|
| Real Data | $30 (avg, more spatial data) | $6,000 |

**ROI:** Accelerates target identification and biomarker discovery
**Value:** Weeks → days for multi-omics integration

---

## Infrastructure Requirements

### DRY_RUN Mode
- **CPU:** Any modern laptop (2+ cores)
- **RAM:** 4GB minimum, 8GB recommended
- **Storage:** 5GB for MCP servers + Python dependencies
- **Network:** Internet for Claude Desktop API calls only

### Real Patient Data Mode
- **CPU:** 8-16 core server recommended for STAR alignment
- **GPU:** Optional but recommended for DeepCell (NVIDIA with 8GB+ VRAM)
- **RAM:** 32GB minimum, 64GB recommended
- **Storage:** 100GB minimum for reference genomes + patient data
- **Network:** Stable connection for TCGA/HuggingFace API calls

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

**Last Updated:** December 27, 2025
**Pricing basis:** Claude Sonnet 4 ($3/M input, $15/M output), AWS EC2 c6i.2xlarge ($0.34/hr), AWS g4dn.xlarge ($0.526/hr)
