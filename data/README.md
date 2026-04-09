# Data Availability

## GEO Accessions

Raw data and processed Seurat objects are deposited in NCBI GEO:

| Dataset | Accession |
|---------|-----------|
| Basak snRNA-seq | [GSE327442](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327442) |
| Basak Visium | [GSE327441](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327441) |

The deposit contains:
- Raw count matrices (Cell Ranger / Space Ranger output) for all Basak samples
- Four processed Seurat objects:
  - `basak_snRNA_seurat.rds` — 2,071 nuclei, all cell types, QC-passed and annotated (used in `01_basak_snRNAseq.Rmd`)
  - `meta_snRNA_allcells_seurat.rds` — 46,417 nuclei, all cell types, three-dataset integration (used in `02_meta_integration.Rmd`)
  - `meta_snRNA_neurons_seurat.rds` — 33,375 neuronal nuclei, 68 subtypes annotated (used in `02_meta_integration.Rmd`)
  - `meta_visium_seurat.rds` — 32,875 Visium spots, Basak + Weber integrated and annotated (used in `04_visium_integrated.Rmd`)

---

## Basak Dataset — Sample Metadata

### snRNA-seq (Supplementary Table 2)

Raw data folder: `data/snRNA/`

| GEO Label / Folder | Donor ID | Age | Sex | Anatomical Level | Cell Count (post-QC) | DBH+ Neurons |
|---|---|---|---|---|---|---|
| `LC001A_filtered_feature_bc_matrix` | Donor_1 | 74 | F | Rostral | 455 | 41 |
| `LC001B_filtered_feature_bc_matrix` | Donor_1 | 74 | F | Rostral | 485 | 41 |
| `LC002A_filtered_feature_bc_matrix` | Donor_1 | 74 | F | Caudal  | 555 | 140 |
| `LC002B_filtered_feature_bc_matrix` | Donor_1 | 74 | F | Caudal  | 576 | 149 |

Tissue sourced from the New York Brain Bank, Columbia University.

### Visium Spatial Transcriptomics (Supplementary Table 3)

Raw data folder: `data/Visium/`

| GEO Label / Folder | Donor ID | Age | Sex | Anatomical Level | Visium Spot Count | NE Neuron Spots |
|---|---|---|---|---|---|---|
| `Upper_Rostral` | Donor_2 | 76 | F | Rostral, section 1 of 2 | 3,085 | 201 |
| `Lower_Rostral` | Donor_2 | 76 | F | Rostral, section 2 of 2 | 2,846 | 205 |
| `Upper_Caudal`  | Donor_2 | 76 | F | Caudal, section 1 of 2  | 3,328 | 185 |
| `Lower_Caudal`  | Donor_2 | 76 | F | Caudal, section 2 of 2  | 3,236 | 203 |

---

## Directory Structure

After downloading all data, the `data/` folder should have the following structure:

```
data/
├── README.md
├── snRNA/
│   ├── LC001A_filtered_feature_bc_matrix/
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── LC001B_filtered_feature_bc_matrix/
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── LC002A_filtered_feature_bc_matrix/
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   └── LC002B_filtered_feature_bc_matrix/
│       ├── barcodes.tsv.gz
│       ├── features.tsv.gz
│       └── matrix.mtx.gz
├── Visium/
│   ├── Upper_Rostral/
│   │   ├── filtered_feature_bc_matrix.h5
│   │   └── spatial/
│   │       ├── aligned_fiducials.jpg
│   │       ├── scalefactors_json.json
│   │       ├── tissue_hires_image.png
│   │       ├── tissue_lowres_image.png
│   │       └── tissue_positions.csv
│   ├── Lower_Rostral/   (same structure)
│   ├── Upper_Caudal/    (same structure)
│   └── Lower_Caudal/    (same structure)
└── external/
    ├── Weber/
    │   ├── weber_snRNA_seurat.rds      # snRNA-seq Seurat object (20,191 nuclei)
    │   └── weber_visium_seurat.rds     # Visium Seurat object (20,380 spots)
    └── Siletti/
        └── Pons_PnRF_PB_seurat.rds    # 24,155 nuclei
```

---

## Published Datasets

The following published datasets are integrated in this analysis.
Download them from the original sources and place them at the paths shown above.

### Weber et al. 2024 (snRNA-seq + Visium)

**Citation:** Weber LM, Divecha HR et al. *eLife* (2024). DOI: https://doi.org/10.7554/eLife.84628

**Data access:** Install the `WeberDivechaLCdata` R/Bioconductor package (available in Bioconductor ≥ 3.16):

```r
BiocManager::install("WeberDivechaLCdata")
```

GitHub and installation instructions: https://github.com/lmweber/WeberDivechaLCdata

#### Downloading and formatting the Weber snRNA-seq data

```r
library(SingleCellExperiment)
library(WeberDivechaLCdata)
library(Seurat)

# Load SingleCellExperiment object
sce <- WeberDivechaLCdata_singleNucleus()

# Convert gene IDs from Ensembl to gene symbols
gene_symbols <- rowData(sce)$gene_name
ens_ids      <- rownames(sce)
gene_symbols_clean <- ifelse(
  is.na(gene_symbols) | gene_symbols == "",
  ens_ids,
  gene_symbols
)
rownames(sce) <- make.unique(gene_symbols_clean)
rowData(sce)$ensembl_id <- ens_ids

# Convert to Seurat
weber_snRNA_seurat <- as.Seurat(sce, counts = "counts", data = "logcounts", project = "weber_snRNAseq")
weber_snRNA_seurat[["RNA"]] <- weber_snRNA_seurat[["originalexp"]]
DefaultAssay(weber_snRNA_seurat) <- "RNA"
weber_snRNA_seurat[["originalexp"]] <- NULL

# Save
saveRDS(weber_snRNA_seurat, file = "data/external/Weber/weber_snRNA_seurat.rds")
```

The resulting object contains 20,191 nuclei across 33,352 features with `counts` and `data` layers,
and PCA/UMAP reductions pre-computed by Weber et al.

#### Downloading and formatting the Weber Visium data

```r
library(SpatialExperiment)
library(WeberDivechaLCdata)
library(Seurat)

# Load SpatialExperiment object
spe <- WeberDivechaLCdata_Visium()

# Convert gene IDs from Ensembl to gene symbols
new_rownames <- rowData(spe)$gene_name
new_rownames[is.na(new_rownames)] <- rownames(spe)[is.na(new_rownames)]
rownames(spe) <- make.unique(new_rownames)

# Extract spatial coordinates
coords_df <- as.data.frame(spatialCoords(spe))
rownames(coords_df) <- colnames(spe)

# Convert to Seurat and add coordinates
weber_visium_seurat <- as.Seurat(spe, counts = "counts", data = "logcounts")
weber_visium_seurat$x_coord <- coords_df$pxl_col_in_fullres
weber_visium_seurat$y_coord <- coords_df$pxl_row_in_fullres
weber_visium_seurat <- RenameAssays(weber_visium_seurat, assay.name = "originalexp", new.assay.name = "RNA")

# Integrate across samples (RPCA)
weber_visium_seurat[["RNA"]] <- split(weber_visium_seurat[["RNA"]], f = weber_visium_seurat$sample_part_id)
weber_visium_seurat <- NormalizeData(weber_visium_seurat) |>
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA()
weber_visium_seurat <- IntegrateLayers(weber_visium_seurat, method = RPCAIntegration, k.weight = 50)
weber_visium_seurat <- RunUMAP(weber_visium_seurat, dims = 1:30, reduction = "integrated.dr")
weber_visium_seurat <- FindNeighbors(weber_visium_seurat, reduction = "integrated.dr", dims = 1:30, k.param = 20)
weber_visium_seurat <- FindClusters(weber_visium_seurat, resolution = seq(0.1, 1.6, by = 0.2))
weber_visium_seurat[["RNA"]] <- JoinLayers(weber_visium_seurat[["RNA"]])

# Save
saveRDS(weber_visium_seurat, file = "data/external/Weber/weber_visium_seurat.rds")
```

The resulting object contains 20,380 spots across 23,728 features.

---

### Siletti et al. 2023 (snRNA-seq)

**Citation:** Siletti K et al. *Science* (2023). DOI: https://doi.org/10.1126/science.add7046

**Source:** CZ CELLxGENE Discover:  
https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443

Only the parabrachial/pontine reticular formation section was used (Dissection: Pons → PnRF → PB;
24,155 nuclei; 10x 3' v3; *Homo sapiens*; tissue type: normal). This section contains the
DBH-expressing NE neuron cluster used in the meta-analysis.

#### Downloading and formatting the Siletti data

Download `Pons_PnRF_PB.h5ad` from CELLxGENE and convert to Seurat as follows:

```r
library(anndata)
library(Matrix)
library(Seurat)

# Load h5ad file
PnRF_PB_data <- read_h5ad("Pons_PnRF_PB.h5ad")

# Extract gene symbols
var_df       <- as.data.frame(PnRF_PB_data$var)
ensembl_ids  <- rownames(var_df)
gene_symbols <- as.character(var_df[["feature_name"]])
gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- ensembl_ids
gene_symbols <- make.unique(gene_symbols)

# Build counts matrix (transpose to genes x cells for Seurat)
counts <- Matrix::t(PnRF_PB_data$X)
rownames(counts) <- gene_symbols
colnames(counts) <- PnRF_PB_data$obs_names

# Metadata
meta <- as.data.frame(PnRF_PB_data$obs)
rownames(meta) <- PnRF_PB_data$obs_names

# Create Seurat object and normalize
PnRF_PB_seurat <- CreateSeuratObject(counts = counts, meta.data = meta)
PnRF_PB_seurat <- NormalizeData(PnRF_PB_seurat, verbose = FALSE)

# Transfer UMAP from h5ad
umap_mat <- as.matrix(PnRF_PB_data$obsm$X_UMAP)
rownames(umap_mat) <- colnames(PnRF_PB_seurat)
PnRF_PB_seurat[["umap"]] <- CreateDimReducObject(
  embeddings = umap_mat,
  key        = "UMAP_",
  assay      = "RNA"
)

# Save
saveRDS(PnRF_PB_seurat, file = "data/external/Siletti/Pons_PnRF_PB_seurat.rds")
```

The resulting object contains 24,155 nuclei across 58,232 features.
