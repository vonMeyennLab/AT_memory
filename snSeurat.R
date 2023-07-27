setwd("/srv/GT/analysis/adhideb/p3360-Ferdinand/")

library(Seurat)
library(dplyr)
library(tibble)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(SingleR)
library(RColorBrewer)
library(scDblFinder)

# Load the dataset
rawCounts <- Read10X(data.dir = "/srv/GT/analysis/adhideb/p3360-Ferdinand/snSeurat/S019_4_mmepiAT_HC/filtered_feature_bc_matrix")

# Initialize Seurat object with raw filtered data
snData <- CreateSeuratObject(counts = rawCounts, project = "HC", min.cells = 3, min.features = 200)

# Doublet detection
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

snData <- subset(snData, scDblFinder.class != "doublet")

snData

# QC on MT and Ribo genes
snData[["percent.mt"]] <- PercentageFeatureSet(snData, pattern = "^mt-")
snData[["percent.ribo"]] <- PercentageFeatureSet(snData, pattern = "^Rp[sl]")

# Visualize QC metrics
VlnPlot(snData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2)

# Visualize feature-feature relationships
FeatureScatter(snData, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filter cells with >6% MT counts
# Filter cells with unique feature counts <300 or >3000
snData_sub <- subset(snData,subset = nFeature_RNA > 300 & nFeature_RNA < 3000 & nCount_RNA < 40000 & percent.mt < 5)

VlnPlot(snData_sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2)
FeatureScatter(snData_sub, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Normalize the data
snData_sub <- NormalizeData(snData_sub, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify variable features
snData_sub<-FindVariableFeatures(snData_sub, selection.method = "vst", nfeatures = 2000)

# Identify 20 most variable genes
top20<-head(VariableFeatures(snData_sub), 20)

# Plot variable features
plot1<-VariableFeaturePlot(snData_sub)
plot2<-LabelPoints(plot = plot1, points = top20, repel = F, xnudge = 0, ynudge = 0)
plot2

# Scale the data
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
snData_sub <- FindClusters(snData_sub, resolution = 0.5)

# Visualization of cell clusters
p1 <- DimPlot(snData_sub, reduction = "umap", label = TRUE)
p2 <- DimPlot(snData_sub, reduction = "tsne", label = TRUE)
p1 + p2

VlnPlot(snData_sub, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
FeaturePlot(snData_sub, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
FeaturePlot(snData_sub, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), reduction = "tsne")

##################################################################################################################################

snData_sub$group <- "HC"
snData_sub$group2 <- "short"
snData_sub$diet <- "HFD_CHD"
snData_sub$age <- "20 wk"
snData_sub$timePoint <- "t2"

##################################################################################################################################

mmRef <- (MouseRNAseqData <- MouseRNAseqData())

singler.results <- SingleR(
  test = GetAssayData(snData_sub), ref = mmRef,
  labels = mmRef$label.fine)

snData_sub[["SingleR.labels"]] <- singler.results$labels

##################################################################################################################################

mmepiAT <- readRDS("mouse_epiAT.rds")

# Metadata info
summary(mmepiAT@meta.data)
table(mmepiAT$name) # Mm_EPI / Mm_ING / Mm_POV
table(mmepiAT$sample) # Mm_EPI / Mm_ING / Mm_POV
table(mmepiAT$depot) # ING / PG
table(mmepiAT$diet) # Chow / HC
table(mmepiAT$animal) # Chow / HC
table(mmepiAT$feeding) # Chow / HC
table(mmepiAT$sex) # female / male
table(mmepiAT$cell_type) # Detailed sub populations
table(mmepiAT$cell_type2) # Broad populations
table(mmepiAT$ct3) # Broad populations with adipocyte sub populations

DimPlot(mmepiAT)

##################################################################################################################################

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

p1 <- DimPlot(mmepiAT, reduction = "umap", group.by = "ct3", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(snData_sub, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

FeaturePlot(snData_sub,features = "mapping.score")
FeaturePlot(snData_sub,features = "predicted.celltype.score")

p1 <- DimPlot(snData_sub, reduction = "umap",group.by = "predicted.celltype", label = TRUE,
        label.size = 3, repel = TRUE)  + ggtitle("Query transfered labels") + NoLegend()

p2 <- DimPlot(snData_sub, reduction = "tsne",group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE)  + ggtitle("Query transfered labels") + NoLegend()

p1 + p2

saveRDS(snData_sub, "snSeurat/S019_4_mmepiAT_HC/snHC.rds")
