# Single-cell RNA-seq analysis of GSE125527 (Ulcerative Colitis)

Seurat-based analysis of the public scRNA-seq dataset
[GSE125527](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125527):
immune cells from peripheral blood (PBMC) and rectal mucosa of ulcerative
colitis patients and healthy controls.

The workflow covers quality control, normalization, dimensionality reduction,
clustering, marker detection, cell-type annotation, and a comparison of
differential expression (diseased vs. healthy) at the single-cell level versus
the donor-level pseudobulk level.

## Dataset

| Item | Value |
|------|-------|
| Cells (before QC) | 68,627 |
| Donors | 16 (7 UC patients `U1`–`U7`, 9 healthy controls `C1`–`C9`) |
| Tissues | PBMC (39,255 cells), rectum (29,372 cells) |
| Metadata columns | `cell_id`, `patient_assignment`, `tissue_assignment`, `disease_assignment`, `celltype` |

## Repository layout

```
.
├── Snakefile                 # Snakemake workflow that runs the R script
├── config/config.yaml        # List of expected output files
├── envs/snakemake_env.yaml   # Conda environment used by Snakemake
├── scripts/script.R          # Full analysis
├── data/                     # GEO input files (tracked with Git LFS)
├── results/                  # Figures and tables produced by the analysis
├── renv.lock                 # R package versions (renv)
├── Dockerfile                # RStudio Server image (rocker/rstudio 4.4.3)
└── docker-compose.yml
```

## Input data

Download the following files from the GEO series page and place them in `data/`:

- `GSE125527_UMI_cell_table_sparse.csv` – UMI counts in triplet format (gene index, cell index, count). Not tracked in the repository because of its size.
- `GSE125527_gene_id_rownames.csv` – gene names (matrix rows)
- `GSE125527_cell_id_colnames.csv` – cell barcodes (matrix columns)
- `GSE125527_cell_metadata.csv` – per-cell metadata

The three smaller files are included in the repository through Git LFS; run
`git lfs pull` after cloning.

## Analysis steps

1. **Matrix construction** – build a sparse gene × cell matrix and a Seurat object with the cell metadata.
2. **QC** – keep cells with 200 < nFeature_RNA < 4000, nCount_RNA < 20,000 and mitochondrial content < 10%.
3. **Normalization and feature selection** – LogNormalize (scale factor 10,000); 2,000 variable genes (vst).
4. **Scaling and PCA** – scale all genes, PCA on the variable genes, elbow plot.
5. **Clustering and UMAP** – kNN graph on PCs 1–10, Louvain clustering at resolution 0.2 (8 clusters), UMAP on PCs 1–10.
6. **Markers** – `FindAllMarkers` (positive markers only), top-10 heatmap, canonical marker violin and feature plots.
7. **Annotation** – manual assignment of cluster labels.
8. **Differential expression, diseased vs. healthy, in the cluster labelled "CD8 T"**
   - single-cell level: Wilcoxon test with `FindMarkers`
   - pseudobulk level: counts aggregated per donor × condition × cell type with `AggregateExpression`, tested with DESeq2
   - genes are split into those with p < 0.05 in both analyses, only in the single-cell analysis, and only in the pseudobulk analysis.

## Outputs (`results/`)

| File | Content |
|------|---------|
| `qc.jpg`, `scatter_plot.jpg` | QC metrics before filtering |
| `variable_features.jpg` | Highly variable genes |
| `PCA_plot_dim.jpg`, `DIM_PCA.jpg`, `DIM_PCA_500*.jpg`, `elbowplot.jpg` | PCA diagnostics |
| `UMAP_plot1.jpg` | UMAP coloured by cluster |
| `seurat.markers.csv`, `heatmap_top10.jpg` | Cluster markers |
| `violin1.jpg`, `violin2.jpg`, `feature1.jpg` | Canonical marker expression |
| `dim_plot_new.jpg`, `final_clustering.jpg` | Annotated UMAP |
| `common_genes.csv` | Genes with p < 0.05 in both single-cell and pseudobulk DE (932 genes) |
| `sc_only_genes.csv` | Genes with p < 0.05 only in the single-cell DE (5,121 genes) |
| `bulk_only_genes.csv` | Genes with p < 0.05 only in the pseudobulk DE (57 genes) |
| `CD8_Violin_by_*.jpg`, `features_in_sc_only_by*.jpg` | Expression of selected DE genes by condition and by donor |

## Running the analysis

### Snakemake + Conda

```bash
snakemake --use-conda --cores 4
```

The workflow writes its outputs to `output/`.

### Docker (RStudio Server)

```bash
docker compose up --build
```

Then open `http://localhost:8787` (user `rstudio`) and restore packages with
`renv::restore()`. The image targets `linux/arm64`; change the `--platform`
line in the `Dockerfile` for other architectures.

### Plain R

```r
source("scripts/script.R")
```

The script must be run from the repository root.

## Requirements

- R 4.4.3
- Seurat 5.3.0, SeuratObject, Matrix, dplyr, ggplot2
- DESeq2 (Bioconductor) for the pseudobulk test
- SeuratData (loaded by the script)
- ~16 GB RAM recommended (`ScaleData` is run on all genes for ~68k cells)

## Reference

Boland BS, He Z, Tsai MS, et al. Heterogeneity and clonal relationships of
adaptive immune cells in ulcerative colitis revealed by single-cell analyses.
*Science Immunology* 2020;5(50):eabb4432.

## License

See [LICENSE](LICENSE).
