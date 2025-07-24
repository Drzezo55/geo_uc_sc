BiocManager::install(c("GEOquery", "DESeq2", "pheatmap"), ask = FALSE)

# Load libraries
library(GEOquery)
library(DESeq2)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(ggfortify)

df <- read.delim("data/GSE243625_RNAseq_colon.txt", check.names = FALSE)
df
row.names(df) <- df$`Gene symbol`
# the above will return error due to 1-mar 2-mar as they are are read by r as two similar values
rownames(df) <- make.unique(df$`Gene symbol`) # this will append .1 , .2
df$`Gene symbol` <- NULL

# create a metadata 
samples <- colnames(df)
condition <- ifelse(grepl("^N", samples), "Control",
                    ifelse(grepl("^I", samples), "Inactive_UC",
                           ifelse(grepl("^A", samples), "Active_UC", NA)))

# Filter out any columns not assigned
keep <- !is.na(condition)
df <- df[, keep]
condition <- condition[keep]

# Create metadata
coldata <- data.frame(row.names = colnames(df),
                      condition = factor(condition, levels = c("Control", "Inactive_UC", "Active_UC")))
# in the above code levels should be ordered the same as the order of ifelse loops in the condition vector
all(colnames(df) == rownames(coldata))  # should return TRUE to check if all the same 

# Transpose the data (samples as rows)
log_expr <- log2(df + 1)
pca <- prcomp(t(log_expr), scale. = FALSE)

# Plot PCA
library(ggplot2)
autoplot(pca, data = coldata, colour = "condition")

autoplot(pca, data = coldata, colour = "condition") +
  ggtitle("PCA of Gene Expression") +
  theme_minimal()

# to compare the active versus control
# Subset samples
keep <- coldata$condition %in% c("Active_UC", "Control")
df_sub <- df[, keep]
coldata_sub <- coldata[keep, , drop = FALSE]

# create the desing matrix
library(edgeR)
library(limma)

# Create a DGEList
dge <- DGEList(counts = df_sub)
# Use 'none' normalization because TPM/FPKM are already normalized
dge <- calcNormFactors(dge, method = "none")
# Design matrix
design <- model.matrix(~0 + coldata_sub$condition)
colnames(design) <- levels(coldata_sub$condition)
# Apply voom transformation
v <- voom(dge, design, plot = TRUE)
# Fit the linear model
fit <- lmFit(v, design)
# Define contrast: Active UC vs Control
contrast_matrix <- makeContrasts(Active_vs_Control = Active_UC - Control, levels = design)
# Fit contrast
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
# Get results
results <- topTable(fit2, coef = "Active_vs_Control", number = Inf)
head(results)
# Filter for significant DEGs
deg_filtered <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
library(EnhancedVolcano)
# volcano plot
EnhancedVolcano(results,
                lab = rownames(results),
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Active UC vs Control')

# heatmap
top_degs <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
top_genes <- rownames(head(top_degs[order(top_degs$adj.P.Val), ], 50))
# Extract expression of top DEGs
heatmap_data <- df_sub[top_genes, ]
# Scale (Z-score across samples)
heatmap_scaled <- t(scale(t(heatmap_data)))
# Annotation for sample groups
annotation_col <- data.frame(
  Condition = coldata_sub$condition
)
rownames(annotation_col) <- rownames(coldata_sub)
#heatmap plotting
library(pheatmap)

pheatmap(heatmap_scaled,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 10,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         scale = "none",  # already scaled above
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Top 50 DEGs: Active UC vs Control")











