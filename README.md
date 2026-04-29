# Prostate Ferroptosis – scripts

This repository contains all scripts used for bioinformatics analysis in our study on dual pathway regulation of castration response and ferroptosis in the prostate epithelium


## Overview

The analysis pipeline covers:

1. **Single-cell RNA-seq** – CellHash demultiplexing, Seurat clustering, doublet removal, and luminal cell extraction for five in-house mouse samples (ML001, ML002, MW3, MW4, MW8)
2. **Public dataset reanalysis** – Reproduction of Karthaus et al. 2020 (*Science*) results for mouse castration time-course and several human prostate cancer datasets
3. **Luminal signature scoring** – LumA/LumP activity scores (bin-matched null distribution) + integration with a machine-learning classifier
4. **Luminal signature visualization** – Figure panels (Figs 1e, 2, 3e, 4, 5i) comparing cohorts across genetic perturbations and drug treatments
5. **GSEA** – Cell-death pathway enrichment in luminal cells, mCRPC ENZA responders, and CSPC vs CRPC tumor cells
6. **Stromal/mesenchymal analysis** – Ptn/Igf1 expression dynamics across castration time points, AR stratification
7. **Bulk RNA-seq** – DESeq2 time-course analysis, Pre vs Post ADT (GSE48403), GCT/CLS export for GSEA
8. **VIPER** - pyVIPER workflow to infer protein activity

## Repository Structure

```
prostate_ferroptosis/
├── scripts/
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
│   ├── Pathway/                        # GMT pathway files
│   └── GSE48403_ADT_RNA/               # Human pre/post ADT RNA-seq
├── output/                             # Generated figures and result tables, also deposited in FigShare
└── README.md
```

## R Dependencies

| Package | Version |
|---|---|
| Seurat | 5.1.0 |
| harmony | 1.2.1 | 
| DoubletFinder | - | 
| clusterProfiler | 4.14.3 |
| enrichplot | 1.26.2 |
| DESeq2 | 1.46.0 | 
| biomaRt | 2.62.0 | 
| ggplot2 | 4.0.1 | 
| data.table | 1.16.4 | 

```

## Usage

Run scripts in order from within the project root:

```r
source("R/00_setup.R")
source("R/01_data_loading_hto_demux.R")
source("R/02_clustering_doublet_removal.R")
# ... and so on
```

## Public Data

- **Public mouse data**: [GSE146811](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146811) (Karthaus et al. 2020 *Science*)
- **Pre/post ADT bulk RNA-seq**: [GSE48403](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48403)
- **mCRPC CSPC/CRPC tumor cells**: [GSE264573](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264573)

In-house mouse prostate scRNA-seq data (GSE295132 and GSE295388) will be released upon publication.

## Contact

For questions, please contact with Weiping Li (wl2761@cumc.columbia.edu), Zejian Wang (zw2994@cumc.columbia.edu) and Michael Shen (ms3625@cumc.columbia.edu) if you have any questions rated with this record
