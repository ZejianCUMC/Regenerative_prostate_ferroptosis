# Prostate Cancer Ferroptosis – scRNA-seq & Bulk RNA-seq Analysis

This repository contains all R scripts used for bioinformatics analysis in our study on ferroptosis and luminal cell plasticity in prostate cancer under androgen deprivation therapy (ADT).

## Overview

The analysis pipeline covers:

1. **Single-cell RNA-seq** – CellHash demultiplexing, Seurat clustering, doublet removal, and luminal cell extraction for five in-house mouse samples (ML001, ML002, MW3, MW4, MW8)
2. **Public dataset reanalysis** – Reproduction of Karthaus et al. 2020 (*Science*) results for mouse castration time-course and human prostate atlas
3. **Luminal signature scoring** – LumA/LumP activity scores (bin-matched null distribution) + integration with a machine-learning classifier
4. **Luminal signature visualization** – Figure panels (Figs 1e, 2, 3e, 4, 5i) comparing cohorts across genetic perturbations and drug treatments
5. **GSEA** – Cell-death pathway enrichment in luminal cells, mCRPC ENZA responders, and CSPC vs CRPC tumor cells
6. **Stromal/mesenchymal analysis** – Ptn/Igf1 expression dynamics across castration time points, AR stratification
7. **Bulk RNA-seq** – DESeq2 time-course analysis, Pre vs Post ADT (GSE48403), GCT/CLS export for GSEA

## Repository Structure

```
prostate_ferroptosis/
├── R/
│   ├── 00_setup.R                      # Libraries, color palettes, marker gene sets
│   ├── 01_data_loading_hto_demux.R     # 10X loading + HTODemux for ML001/002/MW3/4/8
│   ├── 02_clustering_doublet_removal.R # Clustering, DoubletFinder, luminal extraction
│   ├── 03_karthaus_science_reanalysis.R# Mouse + Human public dataset re-analysis
│   ├── 04_luminal_signature_scoring.R  # LumA/LumP scores + classifier integration
│   ├── 05_luminal_signature_visualization.R # Figure panels (scatter plots)
│   ├── 06_gsea_cell_death.R            # GSEA: luminal, ENZA, CSPC vs CRPC
│   ├── 07_stromal_mesenchymal_analysis.R # Ptn/Igf1 dynamics in stromal cells
│   └── 08_bulk_rnaseq.R                # DESeq2 + GCT/CLS export
├── data/                               # Input data (not tracked by git)
│   ├── Luis_marker.xlsx                # Reference marker gene set
│   ├── Pathway/                        # GMT pathway files
│   └── GSE48403_ADT_RNA/               # Human pre/post ADT RNA-seq
├── output/                             # Generated figures and result tables
└── README.md
```

## Dependencies

All scripts are written in **R** using the **ir** (R) kernel. Key packages:

| Package | Version | Purpose |
|---|---|---|
| Seurat | ≥5.0 | Single-cell analysis |
| harmony | - | Batch integration |
| DoubletFinder | - | Doublet detection |
| clusterProfiler | - | GSEA |
| enrichplot | - | GSEA visualization |
| DESeq2 | - | Bulk RNA-seq DE |
| biomaRt | - | Mouse↔Human ortholog mapping |
| ggplot2 / patchwork / ggpubr | - | Visualization |
| data.table | - | Fast I/O |

Install with:
```r
# Bioconductor packages
BiocManager::install(c("DESeq2", "clusterProfiler", "enrichplot", "org.Hs.eg.db",
                       "org.Mm.eg.db", "MAST", "biomaRt"))

# CRAN packages
install.packages(c("Seurat", "harmony", "ggpubr", "patchwork",
                   "ggnewscale", "data.table", "readxl", "tibble"))

# GitHub packages
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
```

## Usage

Run scripts in order from within the project root:

```r
source("R/00_setup.R")
source("R/01_data_loading_hto_demux.R")
source("R/02_clustering_doublet_removal.R")
# ... and so on
```

Steps 01–02 require access to raw 10X Genomics output directories.  
Steps 03–08 can be run independently given the intermediate `.RDS` objects.

## Data Availability

- **Public mouse data**: [GSE146811](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146811) (Karthaus et al. 2020 *Science*)
- **Human prostate atlas**: Available from the Karthaus et al. 2020 *Science* supplementary portal
- **Pre/post ADT bulk RNA-seq**: [GSE48403](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48403)
- **mCRPC CSPC/CRPC tumor cells**: [GSE264573](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264573)

In-house scRNA-seq data (ML001, ML002, MW3, MW4, MW6, MW7, MW8) will be deposited to GEO upon publication.

## Contact

For questions, please contact the corresponding author or open a GitHub Issue.
