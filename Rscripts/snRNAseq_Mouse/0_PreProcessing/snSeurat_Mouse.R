
# This R script has been used to generate results used in 
# Figures 3,6 and associated extended data

library(Seurat)
library(dplyr)
library(tibble)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(SingleR)
library(RColorBrewer)
library(scDblFinder)

# Following code chunk has been illustrated for single mouse snRNA sample
# Same code chunk has been used for remaining all mouse snRNA samples

## QC and pre-processing

# Load filtered count matrix from 10X Cell ranger
rawCounts <- Read10X(data.dir = "../mmepiAT_HC/filtered_feature_bc_matrix")

# Initialize Seurat object with raw filtered data
snData <- CreateSeuratObject(counts = rawCounts, project = "HC", min.cells = 3, min.features = 200)

# Detect doublets
sce <- as.SingleCellExperiment(snData)
sce <- scDblFinder(sce)

# Transform sce object into Seurat
xx <- as.Seurat(sce, counts = "counts", data="logcounts")
zz <- xx@meta.data
zz <- zz[,c("ident" ,"scDblFinder.class","scDblFinder.score","scDblFinder.weighted",
            "scDblFinder.cxds_score")]

yy <- snData@meta.data

meta <- merge(yy, zz, by="row.names")
meta <- column_to_rownames(meta, "Row.names")

snData@meta.data <- meta

table(snData$scDblFinder.class)

# Remove doublets
snData <- subset(snData, scDblFinder.class != "doublet")

snData

# QC on MT and Ribo genes
snData[["percent.mt"]] <- PercentageFeatureSet(snData, pattern = "^mt-")
snData[["percent.ribo"]] <- PercentageFeatureSet(snData, pattern = "^Rp[sl]")

# Visualize QC metrics
VlnPlot(snData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2)

# Visualize feature-feature relationships
FeatureScatter(snData, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filter cells with > 5% MT counts
# Filter cells with unique feature counts < 300 or > 3000
# Filter cells with unique UMI counts > 12,000
snData_sub <- subset(snData,subset = nFeature_RNA > 300 & nFeature_RNA < 3000 & nCount_RNA < 12000 & percent.mt < 5)

VlnPlot(snData_sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2)
FeatureScatter(snData_sub, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

snData_sub

# Log normalize data
snData_sub <- NormalizeData(snData_sub, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify variable features
snData_sub<-FindVariableFeatures(snData_sub, selection.method = "vst", nfeatures = 2000)

# Identify 20 most variable genes
top20 <- head(VariableFeatures(snData_sub), 20)

# Plot variable features
plot1<-VariableFeaturePlot(snData_sub)
plot2<-LabelPoints(plot = plot1, points = top20, repel = F, xnudge = 0, ynudge = 0)
plot2

# Scale data
all.genes<-rownames(snData_sub)
snData_sub<-ScaleData(snData_sub, features = all.genes)

# Perform linear dimensional reduction
snData_sub <- RunPCA(snData_sub, features = VariableFeatures(object = snData_sub))

# Visualize PCA results
DimPlot(snData_sub, reduction = "pca")

ElbowPlot(snData_sub, ndims=50)

# Run non-linear dimensional reduction
snData_sub <- RunUMAP(snData_sub, reduction = "pca", dims = 1:30)
snData_sub <- RunTSNE(snData_sub, reduction = "pca", dims = 1:30)
snData_sub <- FindNeighbors(snData_sub, reduction = "pca", dims = 1:30)

# Find cell clusters
snData_sub <- FindClusters(snData_sub, resolution = 0.5)

# Visualize cell clusters
p1 <- DimPlot(snData_sub, reduction = "umap", label = TRUE)
p2 <- DimPlot(snData_sub, reduction = "tsne", label = TRUE)
p1 + p2

# Visualize QC metrics across cell clusters
VlnPlot(snData_sub, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)

# Add sample metadata
snData_sub$condition <- "HC"

################################################################################

## Cell cluster annotation

# Load Emont et al dataset for epiAT
mmepiAT <- readRDS("../mouse_epiAT.rds") 

# Select common variable features
features <- intersect(mmepiAT$integrated@var.features,snData_sub$RNA@var.features)

# Reference mapping via transfer anchors
epiAT.anchors <- FindTransferAnchors(reference = mmepiAT, query = snData_sub,
                                     dims = 1:30, reference.reduction = "pca", features = features)

# Label transfer
predictions <- TransferData(anchorset = epiAT.anchors, refdata = mmepiAT$ct3,
                            dims = 1:30)
snData_sub <- AddMetaData(snData_sub, metadata = predictions)

# Query mapping on to UMAP embeddings of reference data set
mmepiAT <- RunUMAP(mmepiAT, dims = 1:30, reduction = "pca", return.model = TRUE)
snData_sub <- MapQuery(anchorset = epiAT.anchors, reference = mmepiAT, query = snData_sub,
                       refdata = list(celltype = "ct3"), reference.reduction = "pca", reduction.model = "umap")

snData_sub <- AddMetaData(object = snData_sub, metadata = MappingScore(epiAT.anchors,ndim=30), col.name = "mapping.score")

# Visualize predicted cell type annotations
p1 <- DimPlot(snData_sub, reduction = "umap",group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE)  + ggtitle("Query transfered labels") + NoLegend()

p2 <- DimPlot(snData_sub, reduction = "tsne",group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE)  + ggtitle("Query transfered labels") + NoLegend()

p1 + p2

# Visualize mapping and cell type prediction scores
FeaturePlot(snData_sub,features = "mapping.score")
FeaturePlot(snData_sub,features = "predicted.celltype.score")

saveRDS(snData_sub, "../epiAT_HC.rds")

