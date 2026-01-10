# Spatial Transcriptomics Architecture

**Status:** Production - 2 servers deployed (mcp-spatialtools, mcp-openimagedata)
**Last Updated:** January 10, 2026

---

## System Overview

```mermaid
graph TB
    subgraph Input["📁 Input Data"]
        CSV1[coordinates.csv<br/>900 spots x,y]
        CSV2[expression.csv<br/>31 genes]
        CSV3[annotations.csv<br/>6 regions]
    end

    subgraph MCP["🔧 mcp-spatialtools<br/>(14 tools)"]
        Load[Load CSV Data]

        subgraph Analysis["Analysis Tools"]
            DE[Differential<br/>Expression]
            SA[Spatial<br/>Autocorrelation]
            PE[Pathway<br/>Enrichment]
        end

        subgraph Viz["Visualization Tools"]
            VH[Spatial<br/>Heatmap]
            VG[Gene×Region<br/>Heatmap]
            VR[Region<br/>Chart]
            VM[Moran's I<br/>Plot]
        end

        Bridge[Bridge Tool]
    end

    subgraph Output["📊 Outputs"]
        Stats[Statistical Results<br/>p-values, fold changes]
        Paths[Enriched Pathways<br/>PI3K/AKT/mTOR]
        Imgs[PNG Visualizations<br/>4 plots]
    end

    Multi[mcp-multiomics<br/>Integration]

    CSV1 --> Load
    CSV2 --> Load
    CSV3 --> Load

    Load --> DE
    Load --> SA
    Load --> PE
    Load --> VH
    Load --> VG
    Load --> VR
    Load --> VM

    DE --> Stats
    SA --> Stats
    PE --> Paths
    VH --> Imgs
    VG --> Imgs
    VR --> Imgs
    VM --> Imgs

    Load --> Bridge
    Bridge --> Multi

    classDef inputStyle fill:#e1f5ff,stroke:#0288d1,stroke-width:2px
    classDef mcpStyle fill:#fff3e0,stroke:#f57c00,stroke-width:2px
    classDef outputStyle fill:#f1f8e9,stroke:#689f38,stroke-width:2px
    classDef analysisStyle fill:#fce4ec,stroke:#c2185b,stroke-width:2px
    classDef vizStyle fill:#f3e5f5,stroke:#7b1fa2,stroke-width:2px

    class Input inputStyle
    class MCP mcpStyle
    class Output outputStyle
    class Analysis analysisStyle
    class Viz vizStyle
    class Multi mcpStyle
```

---

## Quick Navigation

### Core Documentation
- **[OVERVIEW.md](OVERVIEW.md)** - System architecture and design principles
- **[SERVERS.md](SERVERS.md)** - All 9 MCP servers (4 deployed, 5 future)
- **[DEPLOYMENT.md](DEPLOYMENT.md)** - GCP Cloud Run deployment

### Workflows
- **[CSV_WORKFLOW.md](CSV_WORKFLOW.md)** - ⭐ **Current implementation** - Pre-processed tabular data (PatientOne tests)
- **[FASTQ_WORKFLOW.md](FASTQ_WORKFLOW.md)** - Raw sequencing alignment (implemented, not tested)

### Reference
- **[TOOL_REFERENCE.md](TOOL_REFERENCE.md)** - All 14 spatialtools tools (10 analysis + 4 visualization)
- **[GLOSSARY.md](GLOSSARY.md)** - Terms and definitions

---

## What This Is

Spatial transcriptomics component of the Precision Medicine MCP system. Processes gene expression data with spatial coordinates for cancer analysis.

**Current implementation:** CSV/tabular workflow (900 spots, 31 genes, 6 regions)
**Server:** mcp-spatialtools (14 tools: 10 analysis + 4 visualization)

---

## Quick Start

**For users:** Read [CSV_WORKFLOW.md](CSV_WORKFLOW.md) → Run [PatientOne TEST_3_SPATIAL.txt](../../tests/manual_testing/PatientOne-OvarianCancer/implementation/TEST_3_SPATIAL.txt)

**For developers:** See [OVERVIEW.md](OVERVIEW.md) for architecture

**For deployers:** See [DEPLOYMENT.md](DEPLOYMENT.md) for GCP procedures

---

**See also:** [Main Architecture](../README.md) | [PatientOne Workflow](../../tests/manual_testing/PatientOne-OvarianCancer/README.md)
