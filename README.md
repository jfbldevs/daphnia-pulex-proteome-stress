# Integrative Bioinformatics Analysis of *Daphnia pulex* Proteome Under Anthropogenic Stress

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Language: R](https://img.shields.io/badge/Language-R%20≥4.3-blue.svg)](https://www.r-project.org/)

> **Paper:** Integrative bioinformatics analysis reveals disrupted metabolic pathways and key hub proteins in the *Daphnia pulex* proteome under anthropogenic stress
>
> **Authors:** Juan Alejandro Norambuena, Lisandra Herrera-Belén, Patricia Poblete-Grant, Jorge F. Beltrán*

## Overview

This repository contains the reproducible analysis pipeline for the systems biology study of *Daphnia pulex* proteomic reorganization in response to anthropogenic stress in Northern Patagonian lakes (Chile). The analysis compares proteomic profiles from Lake Icalma (oligotrophic reference) and Lake Llanquihue (anthropized).

The pipeline integrates:
- **Differential Expression Analysis** (t-test + BH correction)
- **Gene Ontology Enrichment** (Over-Representation Analysis)
- **Gene Set Enrichment Analysis** (GSEA)
- **Metabolic Pathway Classification** (keyword-based from UniProt annotations)
- **Protein-Protein Interaction Network Analysis** (STRING database + igraph)

## Project Structure

```
Paper_Systems_Biology/
├── data/
│   ├── proteinGroups.csv                 # Raw MaxQuant output
│   └── proteinGroups_clean.csv           # Cleaned data for analysis
├── scripts/
│   ├── 00_Setup.R                        # Package installation & configuration
│   ├── 01_Differential_Expression.R      # Data prep, DEA, volcano/MA/heatmap
│   ├── 02_GO_Enrichment_ORA.R            # GO annotation retrieval & ORA
│   ├── 03_GSEA_Analysis.R                # Gene Set Enrichment Analysis
│   ├── 04_Pathway_Analysis.R             # Metabolic pathway classification
│   ├── 05_PPI_Network.R                  # PPI network construction & analysis
│   ├── Run_All_Pipeline.R                # Master script to run all steps
│   ├── Aux_GSEA_Combined_Figure.R        # Combined GSEA figure
│   ├── Aux_PPI_Figures.R                 # Focused PPI network figures
│   └── Aux_Audit_Verification.R          # End-to-end data verification
├── results/                              # Generated tables and data files
├── figures/                              # Generated publication figures
├── .gitignore
├── LICENSE
└── README.md
```

## Requirements

### R Version
- R ≥ 4.3.0

### R Packages

**CRAN:**
tidyverse, readxl, writexl, ggplot2, pheatmap, RColorBrewer, ggrepel,
gridExtra, scales, igraph, ggraph, tidygraph, httr, jsonlite,
VennDiagram, UpSetR, cowplot, viridis, here

**Bioconductor:**
clusterProfiler, org.Dm.eg.db, GO.db, DOSE, enrichplot, AnnotationDbi,
pathview, STRINGdb, biomaRt

> Packages are automatically installed by `00_Setup.R` on first run.

## How to Run

### Complete Pipeline (recommended)

```r
# From the project root directory:
source("scripts/Run_All_Pipeline.R")
```

### Step-by-Step

```r
# 1. Set working directory to scripts/
setwd("path/to/Paper_Systems_Biology/scripts")

# 2. Run each script in order:
source("00_Setup.R")
source("01_Differential_Expression.R")
source("02_GO_Enrichment_ORA.R")
source("03_GSEA_Analysis.R")
source("04_Pathway_Analysis.R")
source("05_PPI_Network.R")

# 3. Optional: Generate auxiliary figures
source("Aux_GSEA_Combined_Figure.R")
source("Aux_PPI_Figures.R")

# 4. Optional: Verify data integrity
source("Aux_Audit_Verification.R")
```

## Expected Outputs

### Tables (in `results/`)

| File | Description |
|------|-------------|
| `Table1_Differential_Expression.csv` | Full differential expression results |
| `Table2_GSEA_Enrichment.csv` | GSEA summary across GO categories |
| `Table3_Pathway_Analysis.csv` | Metabolic pathway classification summary |
| `Table4_Hub_Proteins.csv` | Hub protein identification and functions |
| `TableS1_DEPs_Full_List.csv` | Complete DEPs list with annotations |
| `TableS2–S4_GSEA_*.csv` | GSEA results by GO category (BP, MF, CC) |
| `TableS5_Pathway_Proteins.csv` | Per-protein pathway assignments |
| `TableS6_PPI_Centrality.csv` | Network centrality metrics for all nodes |
| `TableS7_PPI_Interactions.csv` | PPI edge list with confidence scores |

### Figures (in `figures/`)

| File | Description |
|------|-------------|
| `Fig1A_Volcano_Plot` | Differential expression volcano plot |
| `Fig1B_MA_Plot` | MA plot of expression differences |
| `Fig1C_Heatmap_DEPs` | Hierarchical clustering heatmap |
| `Fig2_GSEA_Enrichment` | Combined GSEA results |
| `Fig3A_Metabolic_Pathways` | Metabolic pathway bar chart |
| `Fig3B_Pathway_Heatmap` | Pathway expression heatmap |
| `Fig4A_PPI_Network_Metabolic` | Full PPI network with metabolic functions |
| `Fig4B_PPI_Energy_Metabolism` | Focused energy metabolism subnetwork |

### Reproducibility

Session information is automatically captured at `results/session_info.txt` when running the full pipeline.

## Data

Proteomic data is from Norambuena et al. (2026). The raw MaxQuant output (`proteinGroups.csv`) contains LFQ intensities for 3 biological replicates per lake.

## Notes

- **Internet required:** Scripts 02 and 05 query the UniProt and STRING APIs. Cached results are reused if available.
- **Runtime:** The full pipeline takes approximately 10–20 minutes depending on API response times.

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

## Citation

If you use this code, please cite:

> Norambuena, J.A., Herrera-Belén, L., Poblete-Grant, P., & Beltrán, J.F. (2026). Integrative bioinformatics analysis reveals disrupted metabolic pathways and key hub proteins in the *Daphnia pulex* proteome under anthropogenic stress.
