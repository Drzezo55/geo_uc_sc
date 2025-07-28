# Activate renv environment
#source("renv/activate.R")

# Ensure remotes is installed and loaded (Conda will handle installation)
if (!requireNamespace("remotes", quietly = TRUE)) {
  # This part should ideally not be reached if r-remotes is in environment.yaml
  # but as a failsafe, you could put install.packages("remotes") here.
  # For now, let's assume Conda installs it.
  stop("The 'remotes' package is required but not found. Please ensure it's in your Conda environment.")
}

# Install SeuratData from GitHub
if (!requireNamespace("SeuratData", quietly = TRUE)) {
  remotes::install_github("satijalab/seurat-data")
}
#install.packages("Matrix")
# Load required libraries
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
#install.packages("remotes")  # if not installed
#remotes::install_github("satijalab/seurat-data")
library(SeuratData)
# Set seed for reproducibility
set.seed(42)
# Load UMI sparse matrix (triplet format)
umi_table <- read.csv("data/GSE125527_UMI_cell_table_sparse.csv")
# Load row and column names
gene_names <- read.csv("data/GSE125527_gene_id_rownames.csv", header = FALSE)[[1]]
cell_ids <- read.csv("data/GSE125527_cell_id_colnames.csv", header = FALSE)[[1]]
# Check: column = cell, row = gene
head(umi_table)
# Rename columns for clarity
colnames(umi_table) <- c("gene", "cell", "count")
# Construct sparse matrix
sparse_counts <- sparseMatrix(
  i = umi_table$gene,
  j = umi_table$cell,
  x = umi_table$count
)
# Assign row and column names
rownames(sparse_counts) <- gene_names
colnames(sparse_counts) <- cell_ids
# Check matrix shape
dim(sparse_counts)
# metdata 
# Load metadata
metadata <- read.csv("data/GSE125527_cell_metadata.csv")
rownames(metadata) <- metadata$cell_id
# Ensure metadata matches cell order
metadata <- metadata[colnames(sparse_counts), ]
# Create Seurat object
seurat <- CreateSeuratObject(counts = sparse_counts, meta.data = metadata, project = "GSE125527")
head(seurat@meta.data)
# qc
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
qc <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "output/qc.jpg", height = 7, width = 12, plot = qc, quality = 50)

#scatter
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
scatter_plot <- plot1 + plot2
ggsave(filename = "output/scatter_plot.jpg", height = 7, width = 12, plot = scatter_plot, quality = 50)
seurat <- subset(seurat,
                 subset = nFeature_RNA > 200 &
                   nFeature_RNA < 4000 &
                   nCount_RNA < 20000 &
                   percent.mt < 10)
#norm
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
variable_features <- plot1 + plot2
ggsave(filename = "output/variable_features.jpg", height = 7, width = 12, plot = variable_features, quality = 50)

# scale and dim
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
# Examine and visualize PCA results a few different ways
print(seurat[["pca"]], dims = 1:5, nfeatures = 5)
PCA_plot_dim <- VizDimLoadings(seurat, dims = 1:2, reduction = "pca")
ggsave(filename = "output/PCA_plot_dim.jpg", height = 7, width = 12, plot = PCA_plot_dim, quality = 50)

DIM_PCA <- DimPlot(seurat, reduction = "pca") + NoLegend()
ggsave(filename = "output/DIM_PCA.jpg", height = 7, width = 12, plot = DIM_PCA, quality = 50)

DIM_PCA_500 <- DimHeatmap(seurat, dims = 1, cells = 500, balanced = TRUE)
ggsave(filename = "output/DIM_PCA_500.jpg", height = 7, width = 12, plot = DIM_PCA_500, quality = 50)

DIM_PCA_500_15 <- DimHeatmap(seurat, dims = 1:15, cells = 500, balanced = TRUE)
ggsave(filename = "output/DIM_PCA_500_15.jpg", height = 7, width = 12, plot = DIM_PCA_500_15, quality = 50)

elbowplot <- ElbowPlot(seurat)
ggsave(filename = "output/elbowplot.jpg", height = 7, width = 12, plot = elbowplot, quality = 50)
# clustering 
seurat <- FindNeighbors(seurat, dims = 1:10)
seurat <- FindClusters(seurat, resolution = 0.2)
# non linear 
seurat <- RunUMAP(seurat, dims = 1:10)
UMAP_plot1 <- DimPlot(seurat, reduction = "umap", label = TRUE)
ggsave(filename = "output/UMAP_plot1.jpg", height = 7, width = 12, plot = UMAP_plot1, quality = 50)

# marker 
# find all markers of cluster 4
cluster4.markers <- FindMarkers(seurat, ident.1 = 4)
head(cluster4.markers, n = 5)
# find all markers distinguishing cluster 3 from clusters 0 and 1,2,7
cluster3.markers <- FindMarkers(seurat, ident.1 = 3, ident.2 = c(0, 1,2,6))
levels(seurat)
head(cluster3.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE)
seurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
write.csv(seurat.markers, "output/seurat.markers.csv", row.names = TRUE)

violin1 <- VlnPlot(seurat, features = c("IL7R", "CCR7", "S100A4", "CD8A"))
ggsave(filename = "output/violin1.jpg", height = 7, width = 12, plot = violin1, quality = 50)

# plot raw counts as well
violin2 <- VlnPlot(seurat, features = c("MS4A1", "NKG7", "GNLY"), slot = "counts", log = TRUE)
ggsave(filename = "output/violin2.jpg", height = 7, width = 12, plot = violin2, quality = 50)

feature1 <- FeaturePlot(seurat, features = c("IL7R", "CCR7", "S100A4", "CD8A","MS4A1", "NKG7", "GNLY"))
ggsave(filename = "output/feature1.jpg", height = 7, width = 12, plot = feature1, quality = 50)

seurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
heatmap_top10 <- DoHeatmap(seurat, features = top10$gene) + NoLegend()
ggsave(filename = "output/heatmap_top10.jpg", height = 7, width = 12, plot = heatmap_top10, quality = 50)

# assign
new.cluster.ids <- c("CD8 T", "Memory CD4 T", "Naive CD4 T", "B", "NK", "CD14+ Mono", "Memory CD4 T", "Memory CD4 T")
names(new.cluster.ids) <- levels(seurat)
seurat <- RenameIdents(seurat, new.cluster.ids)
seurat$seurat_annotations <- Idents(seurat)
dim_plot_new <- DimPlot(seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(filename = "output/dim_plot_new.jpg", height = 7, width = 12, plot = dim_plot_new, quality = 50)

library(ggplot2)
final_clustering <- DimPlot(seurat, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "output/final_clustering.jpg", height = 7, width = 12, plot = final_clustering, quality = 50)


seurat$seurat_annotations <- Idents(seurat)
#now we will compare the cells across conditions
seurat$new_idents <- paste(seurat$seurat_annotations, seurat$disease_assignment, sep = "_")
Idents(seurat) <- "new_idents"
# Set the future global size limit to 8 GB
options(future.globals.maxSize = 8000 * 1024^2)  # 8 GB , u need it due to the internal function of paralellization, otherwise U will get error as it will exceed the default of the memory size
tcell <- FindMarkers(seurat, ident.1 = "CD8 T_diseased", ident.2 = "CD8 T_healthy", verbose = FALSE)
head(tcell, n = 10)
# aggregate 
# pseudobulk the counts based on donor-condition-celltype
pseudo_seurat <- AggregateExpression(seurat, assays = "RNA", return.seurat = T, group.by = c("patient_assignment","disease_assignment",
                                                                                             "seurat_annotations"))
tail(Cells(pseudo_seurat))
pseudo_seurat$new_idents <- paste(pseudo_seurat$seurat_annotations, pseudo_seurat$disease_assignment, sep = "_")
# each 'cell' is a donor-condition-celltype pseudobulk profile
Idents(pseudo_seurat) <- "new_idents"
bulk.tcr <- FindMarkers(object = pseudo_seurat, 
                            ident.1 = "CD8 T_diseased", 
                            ident.2 = "CD8 T_healthy",
                            test.use = "DESeq2")
head(bulk.tcr, n = 15)
# compare the DE P-values between the single-cell level and the pseudobulk level results
names(bulk.tcr) <- paste0(names(bulk.tcr), ".bulk")
bulk.tcr$gene <- rownames(bulk.tcr)

names(tcell) <- paste0(names(tcell), ".sc")
tcell$gene <- rownames(tcell)

merge_dat <- merge(tcell, bulk.tcr, by = "gene")
merge_dat <- merge_dat[order(merge_dat$p_val.bulk), ]

# Number of genes that are marginally significant in both; marginally significant only in bulk; and marginally significant only in single-cell
common <- merge_dat$gene[which(merge_dat$p_val.bulk < 0.05 & 
                                 merge_dat$p_val.sc < 0.05)]
only_sc <- merge_dat$gene[which(merge_dat$p_val.bulk > 0.05 & 
                                  merge_dat$p_val.sc < 0.05)]
only_bulk <- merge_dat$gene[which(merge_dat$p_val.bulk < 0.05 & 
                                    merge_dat$p_val.sc > 0.05)]
print(paste0('# Common: ',length(common)))
print(paste0('# Only in bulk: ',length(only_bulk)))
print(paste0('# Only in single-cell: ',length(only_sc)))
# export the data frames 
common_df <- merge_dat[merge_dat$gene%in%common,]
common_df <- common_df[order(common_df$p_val.sc, decreasing = TRUE),]
write.csv(common_df, "output/common_genes.csv", row.names = FALSE)
only_sc_df <- merge_dat[merge_dat$gene%in%only_sc,]
only_sc_df <- only_sc_df[order(only_sc_df$p_val.sc, decreasing = TRUE),]
write.csv(only_sc, "output/sc_only_genes.csv", row.names = FALSE)
only_bulk_df <- merge_dat[merge_dat$gene%in%only_bulk,]
only_bulk_df <- only_bulk_df[order(only_bulk_df$p_val.sc, decreasing = TRUE),]
write.csv(only_bulk, "output/bulk_only_genes.csv", row.names = FALSE)





# create a new column to annotate sample-condition-celltype in the single-cell dataset
seurat$id.disease <- paste0(seurat$patient_assignment, "-", seurat$disease_assignment)
# generate violin plot 
seurat$new_idents <- paste(seurat$seurat_annotations, seurat$disease_assignment, sep = "_")
Idents(seurat) <- "new_idents"
print(merge_dat[merge_dat$gene%in%common[1:2],c('gene','p_val.sc','p_val.bulk')])
CD8_Violin_by_disease <- VlnPlot(seurat, features = common[1:10], idents = c("CD8 T_diseased", "CD8 T_healthy"), group.by = "disease_assignment") 
ggsave(filename = "output/CD8_Violin_by_disease.jpg", height = 7, width = 12, plot = CD8_Violin_by_disease, quality = 50)

CD8_Violin_by_patient <- VlnPlot(seurat, features = common[1:2], idents = c("CD8 T_diseased", "CD8 T_healthy"), group.by = "patient_assignment", ncol = 1) # to ensure robustness
ggsave(filename = "output/CD8_Violin_by_patient.jpg", height = 7, width = 12, plot = CD8_Violin_by_patient, quality = 50)

print(merge_dat[merge_dat$gene%in%c('TMCO3','MFAP1'),c('gene','p_val.sc','p_val.bulk')]) # DE only under sc

features_in_sc_only_bydisease <- VlnPlot(seurat, features <- c('PPIB','POMC'), idents = c("CD8 T_diseased", "CD8 T_healthy"), group.by = "disease_assignment") 
ggsave(filename = "output/features_in_sc_only_bydisease.jpg", height = 7, width = 12, plot = features_in_sc_only_bydisease, quality = 50)


features_in_sc_only_bypatient <- VlnPlot(seurat, features <- c('PPIB','POMC'), idents = c("CD8 T_diseased", "CD8 T_healthy"), group.by = "patient_assignment", ncol = 1) # to ensure robustness
ggsave(filename = "output/features_in_sc_only_bypatient.jpg", height = 7, width = 12, plot = features_in_sc_only_bypatient, quality = 50)










