
# This R script has been used to generate results used in 
# Figures 2,5 and associated extended data

library(SingleR)
library(xlsx)
library(AUCell)
library(Seurat)
library(pheatmap)
library(cowplot)
library(ggplot2)
library(harmony)
library(RColorBrewer)
library(clustree)
library(dplyr)
library(ezRun)

# Following code chunk has been illustrated for weight loss mouse snRNA samples
# Same code chunk has been used for YoYo mouse snRNA samples

# Load each mouse snRNA sample
c <- readRDS("../snC.rds")
cc <- readRDS("../snCC.rds")
h <- readRDS("../snH.rds")
hc <- readRDS("../snHC.rds")
ccc <- readRDS("../snCCC.rds")
hh <- readRDS("../snHH.rds")
hhc <- readRDS("../snHCC.rds")

# Create list of all Seurat objects
snSampleList <- list(c, cc, ccc,
                     h, hc, hh, hhc)

# Perform SCT on each object regressing out the mito and ribo genes
for (i in 1:length(snSampleList)) {
  snSampleList[[i]] <- SCTransform(snSampleList[[i]], vars.to.regress = c("percent.mt", "percent.ribo"), verbose = FALSE)
}

# Integrate SCT normalized objects
features <- SelectIntegrationFeatures(object.list = snSampleList)

snSampleList <- PrepSCTIntegration(object.list = snSampleList, anchor.features = features)

int.anchors <- FindIntegrationAnchors(object.list = snSampleList, normalization.method = "SCT",
                                      anchor.features = features)
snIntegrated <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT")

# Check assay
names(snIntegrated@assays)

snIntegrated@active.assay

# Run standard workflow for visualization and clustering
snIntegrated <- ScaleData(snIntegrated, verbose = FALSE)
snIntegrated <- RunPCA(snIntegrated, verbose = FALSE)
ElbowPlot(snIntegrated, ndims = 50)

snIntegrated <- RunUMAP(snIntegrated, reduction = "pca", dims = 1:30)
snIntegrated <- RunTSNE(snIntegrated, reduction = "pca", dims = 1:30)
snIntegrated <- FindNeighbors(snIntegrated, reduction = "pca", dims = 1:30)

# Select a range of resolutions
resolution.range <- seq(from = 0.1, to = 1, by = 0.1)

# Find clusters using a range of resolutions
snIntegrated <- Seurat::FindClusters(object = snIntegrated, resolution = resolution.range)

# Clustering overview
clustree(snIntegrated)

snIntegrated <- FindClusters(snIntegrated, resolution = 0.4)

# Visualize cell clusters of integrated samples
DimPlot(snIntegrated, reduction = "umap", label = TRUE)
DimPlot(snIntegrated, reduction = "tsne", label = TRUE)

DimPlot(snIntegrated, reduction = "pca", group.by = "condition")
DimPlot(snIntegrated, reduction = "tsne", group.by = "condition")
DimPlot(snIntegrated, reduction = "umap", group.by = "condition")

DimPlot(snIntegrated, reduction = "umap", split.by = "condition", ncol = 3)
DimPlot(snIntegrated, reduction = "tsne", split.by = "condition", ncol = 3)

# Visualize QC metrics across cell clusters
VlnPlot(snIntegrated, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), pt.size = 0, ncol=2)

# Find cell cluster specific marker genes
DefaultAssay(snIntegrated) <- "RNA"
markers <- FindAllMarkers(snIntegrated, 
                          only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, return.thresh = 0.01)

saveRDS(snIntegrated, "../epiAT_weightLoss.rds")
