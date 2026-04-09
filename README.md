# Cellular Signatures of Melanocortin Pathway Genes Across the Locus Coeruleus

This repository contains the analysis code accompanying the manuscript:

> **Cellular Signatures of Melanocortin Pathway Genes Across the Locus Coeruleus**  
> Basak et al., *Acta Neuropathologica Communications* (accepted)

---

## Overview

We characterize melanocortin pathway gene expression in the human locus coeruleus (LC) using single-nucleus RNA sequencing, Visium spatial transcriptomics, and RNAscope HiPlex. The code in this repository reproduces the snRNA-seq and Visium analyses reported in the paper. RNAscope image analysis was performed in HALO (Indica Labs) and is not included here.

---

## Repository Structure

```
.
├── README.md
├── data/
│   └── README.md                 # Data access instructions and sample metadata
├── resources/
│   ├── gene_lists/
│   │   └── mito_ribo_genes.rds   # Mitochondrial/ribosomal gene list for filtering
│   └── IUPHAR/
│       ├── README.md             # Source and version information
│       ├── GPCRTargets.csv       # IUPHAR GPCR gene list
│       └── GtP_to_HGNC_mapping.csv
├── R/
│   ├── apply_qc_thresholds.R     # Batch-aware QC filtering
│   ├── qc_plots_by_sample.R      # QC visualization
│   └── DotPlot_allow_dups.R      # DotPlot with duplicate gene support
├── 01_basak_snRNAseq.Rmd         # Basak snRNA-seq: loading, QC, integration, annotation
├── 02_meta_integration.Rmd       # Three-dataset integration, neuronal subtypes, DE
├── 03_basak_visium.Rmd           # Basak Visium: loading, QC, section processing
├── 04_visium_integrated.Rmd      # Basak + Weber Visium integration and annotation
├── figures_main.Rmd              # Code for main figures 1–2
└── figures_supplementary.Rmd     # Code for supplementary figures 1–3, 7–8
```

All analysis scripts read from `data/` and `resources/` and write intermediate
objects to `output/` (generated automatically at runtime).

---

## Execution Order

Scripts are designed to be run in numbered order:

1. `01_basak_snRNAseq.Rmd`
2. `02_meta_integration.Rmd`
3. `03_basak_visium.Rmd`
4. `04_visium_integrated.Rmd`
5. `figures_main.Rmd`
6. `figures_supplementary.Rmd`

---

## Data Availability

Raw data (Cell Ranger / Space Ranger output matrices) and processed Seurat objects
are deposited in NCBI GEO. See `data/README.md` for accession numbers,
sample metadata, and download instructions.

Published datasets from Weber et al. (2024) and Siletti et al. (2023) are available
from the original publications. Expected download paths are documented in
`data/README.md`.

---

## Software Requirements

| Package        | Version  |
|----------------|----------|
| R              | 4.3.2    |
| Seurat         | 5.0.1    |
| SeuratWrappers | ≥ 0.3.2  |
| SeuratObject   | 5.0.1    |
| dplyr          | 1.1.4    |
| ggplot2        | 3.5.1    |
| edgeR          | ≥ 3.44   |
| MAST           | ≥ 1.28   |
| ggrepel        | ≥ 0.9    |
| gridExtra      | ≥ 2.3    |
| patchwork      | ≥ 1.2    |
| tibble         | ≥ 3.2    |
| tidyr          | ≥ 1.3    |
| RColorBrewer   | ≥ 1.1    |
