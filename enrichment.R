BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"))
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(tidyverse)
# Example: top DEGs from limma
deg_genes <- rownames(results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ])
# Map gene symbols to Entrez IDs
gene_df <- bitr(deg_genes, fromType = "SYMBOL",
                toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_ids <- gene_df$ENTREZID
# check for the unmapped gene symbols and their percent
unmapped_genes <- setdiff(deg_genes, gene_df$SYMBOL)
# Get the top 50 gene names
top50_genes <- rownames(deg_filtered)[1:50]
# Check which unmapped genes are in the top 50
unmapped_in_top50 <- intersect(unmapped_genes, top50_genes)
unmapped_in_top50
# proceed to GO 
ego <- enrichGO(gene = entrez_ids,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",       # BP: biological process
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)

# View top terms
head(ego)
dotplot(ego, showCategory = 20) + ggtitle("GO Enrichment (BP)")
emapplot(ego)

library(enrichplot)

ego_sim <- pairwise_termsim(ego)  # compute term similarity
emapplot(ego_sim, showCategory = 30)  # plot top 30 terms

library(dplyr)
top_terms <- ego@result %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::select(Description, p.adjust) %>%
  head(20)
print(top_terms)
dotplot(ego, showCategory = 30) + ggtitle("Top GO Enrichment Terms")
barplot(ego, showCategory = 20, title = "Top GO Terms")
gene_fc_vector <- deg_filtered$logFC
names(gene_fc_vector) <- deg_filtered$ENTREZID
# Create mapping vector
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Hs.eg.db)
cnetplot(ego, categorySize = "pvalue", foldChange = gene_fc_vector)
cnetplot(ego, categorySize = "pvalue")

# proceed to KEGG
ekegg <- enrichKEGG(
  gene = entrez_ids,
  organism = 'hsa',       # 'hsa' = Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
dotplot(ekegg, showCategory = 20) + ggtitle("KEGG Pathway Enrichment")
library(enrichplot)
kegg_sim <- pairwise_termsim(ekegg)
cnetplot(ekegg, categorySize="pvalue", foldChange=gene_fc_vector)
emapplot(kegg_sim, showCategory=20)
library(org.Hs.eg.db)
library(AnnotationDbi)

gene_map <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = names(gene_fc_vector),
                                  columns = "SYMBOL",
                                  keytype = "ENTREZID")

# Replace names with symbols (optional)
names(gene_fc_vector) <- gene_map$SYMBOL[match(names(gene_fc_vector), gene_map$ENTREZID)]
library(org.Hs.eg.db)
library(AnnotationDbi)
cnetplot(ekegg, categorySize="pvalue", foldChange=gene_fc_vector)

# prepare for gsea 

# deg_filtered is your DE results data frame with 'logFC' and 'symbol'
library(clusterProfiler)
library(org.Hs.eg.db)
deg_filtered$SYMBOL <- row.names(deg_filtered)
row.names(deg_filtered) <- NULL
mapped_genes <- bitr(deg_filtered$SYMBOL, fromType="SYMBOL",
                     toType="ENTREZID", OrgDb=org.Hs.eg.db)
# Merge back to get ENTREZID and logFC
deg_filtered_mapped <- merge(deg_filtered, mapped_genes, by.x="SYMBOL", by.y="SYMBOL")
gene_list <- deg_filtered_mapped$logFC
names(gene_list) <- deg_filtered_mapped$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

library(clusterProfiler)

gsea_res <- gseKEGG(geneList = gene_list,
                    organism = "hsa",
                    nPerm = 1000,            # number of permutations
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    verbose = FALSE)

# Example: plot the first pathway in results
gseaplot2(gsea_res, geneSetID = gsea_res$ID[1], title = gsea_res$Description[1])
dotplot(gsea_res, showCategory=20)

emapplot(pairwise_termsim(gsea_res), showCategory=20)

ridgeplot(gsea_res, showCategory = 20)
heatplot(gsea_res, showCategory = 10, foldChange = gene_list)
gseaplot2(gsea_res, 
          geneSetID = c(gsea_res$ID[1:5]), 
          title = "Top 5 KEGG Pathways")









