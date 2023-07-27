setwd("/srv/GT/analysis/adhideb/p3360-Ferdinand/snSeurat_multiSample/short/")

library(Seurat)
library(tidyverse)

# Load annotated rds object
#snIntegrated <- readRDS("snMulti_short_annotated.rds")
snIntegrated=snMulti_CHD_timePoint_annotated
DimPlot(snIntegrated, reduction = "umap", split.by = "age", label = TRUE)
# Create new variable combining condition and cell type info
snIntegrated@meta.data$cellTypeCond <-  paste(snIntegrated$age, snIntegrated$cellType, sep = "_")

# Update cell idents
Idents(snIntegrated) <- "cellTypeCond"

# Define cell types for which DEGs are going to computed
cellType <- c("APCs", "FIPs")

# Intialize empty list to dtore DEGs
diffGenes <- list()

# For each defined cell type, compute DEGs
for(eachCluster in cellType)
{
  # Update the conditions for DE analysis
  markersEach <- try(FindMarkers(snIntegrated, ident.1=paste0("31 wk_", eachCluster),
                                 ident.2=paste0("18 wk_", eachCluster), test.use = "wilcox"))
  
  # Skip cell types with few cells
  if(class(markersEach) != "try-error"){
    diffGenes[[eachCluster]] <- as_tibble(markersEach, rownames="gene")
  }
}

# Combine cell type specific DEGs
diffGenes <- bind_rows(diffGenes, .id="cellType")

# Calculate pct difference
# pct.1: cell fraction expressing DE gene in cond 1
# pct.2: cell fraction expressing DE gene in cond 1
diff_pct = abs(diffGenes$pct.1-diffGenes$pct.2)
diffGenes$diff_pct <- diff_pct
diffGenes <- diffGenes[order(diffGenes$diff_pct, decreasing = TRUE),]
rownames(diffGenes) <- NULL

# Only keep the significant DEGs
diffGenes <- subset(diffGenes, p_val_adj < 0.05)

write.table(diffGenes, "LHFD_vs_LCHD_DEGs.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.csv(diffGenes, "31wk_vs_18wk.csv", sep = ",")

VlnPlot(snAPCs_long, c("nCount_RNA", "nFeature_RNA"))







