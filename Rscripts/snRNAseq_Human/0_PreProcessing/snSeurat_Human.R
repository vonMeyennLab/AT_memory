
# This R script has been used to generate results used in 
# Figure 1 and associated extended data

library(Seurat)
library(dplyr)
library(tibble)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(SingleR)
library(RColorBrewer)
library(scDblFinder)

# Following code chunk has been illustrated for single human snRNA sample
# Same code chunk has been used for remaining all human snRNA samples

## SoupX correction of ambient RNA

# Load Cell Ranger output
scData <- load10X('../omAT_Lean_MTSS')

## Visualize Adipocyte specific genes for sanity check (tSNE embeddings from CellRanger)
adipoGenes <- c("ADIPOQ", "PLIN1", "LIPE", "GPAM")

plotMarkerMap(scData, adipoGenes[1]) + ggtitle(adipoGenes[1]) + plotMarkerMap(scData, adipoGenes[2]) + ggtitle(adipoGenes[2])
plotMarkerMap(scData, adipoGenes[3]) + ggtitle(adipoGenes[3]) + plotMarkerMap(scData, adipoGenes[4]) + ggtitle(adipoGenes[4])

plotMarkerMap(scData, adipoGenes) + ggtitle("Adipocyte markers")

# Calculate contamination fraction
scData <- autoEstCont(scData) 

# List most highly expressed genes
head(scData$soupProfile[order(scData$soupProfile$est,decreasing=TRUE),],n=20)

# Visuailze observed to expected expression distribution of marker genes
plotMarkerDistribution(scData)

# Remove background contamination from count matrix
out <- adjustCounts(scData)

# Visualize gene expression changes after correction
# Default: fraction of counts identified as contamination
plotChangeMap(scData,out, adipoGenes[1]) + ggtitle(adipoGenes[1]) + plotChangeMap(scData,out, adipoGenes[2]) + ggtitle(adipoGenes[2])
plotChangeMap(scData,out, adipoGenes[3]) + ggtitle(adipoGenes[3]) + plotChangeMap(scData,out, adipoGenes[4]) + ggtitle(adipoGenes[4])

# Binary: cells expressed or not
# Counts: corrected and uncorrected counts
plotChangeMap(scData,out, adipoGenes) + ggtitle("Adipocyte markers")
plotChangeMap(scData,out, adipoGenes, dataType = "counts") + ggtitle("Adipocyte markers")
plotChangeMap(scData,out, adipoGenes, dataType = "binary") + ggtitle("Adipocyte markers")

################################################################################

## QC and pre-processing

# Initialize Seurat object with SoupX corrected count data
snData <- CreateSeuratObject(counts = out, project = "omAT_Lean_MTSS", min.cells = 3, min.features = 200)

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
snData[["percent.mt"]] <- PercentageFeatureSet(snData, pattern = "^MT-")
snData[["percent.ribo"]] <- PercentageFeatureSet(snData, pattern = "^RP[SL]")

# Visualize QC metrics
VlnPlot(snData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2)

# Visualize feature-feature relationships
FeatureScatter(snData, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filter cells with > 5% MT counts
# Filter cells with unique feature counts < 300 or > 4000 # 6000 for scAT NEFA samples
# Filter cells with unique UMI counts > 15,000 # 25,000 for scAT NEFA samples
snData_sub <- subset(snData,subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 5)

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
snData_sub$tissue <- "omAT"
snData_sub$status <- "lean"
snData_sub$cohort <- "MTSS"

################################################################################

## Cell cluster annotation

# Load Emont et al dataset for omAT
# For scAT samples, scAt refrence has been used from Emont et al.
hsAT <- readRDS("../human_visAT_all.rds") 

# Select only Caucasian individuals
hsAT <- subset(hsAT, ethnicity == "Caucasian")

# Select common variable features
features <- intersect(hsAT$integrated@var.features,snData_sub$RNA@var.features)

# Reference mapping via transfer anchors
epiAT.anchors <- FindTransferAnchors(reference = hsAT, query = snData_sub,
                                     dims = 1:30, reference.reduction = "pca", features = features)

# Label transfer
predictions <- TransferData(anchorset = epiAT.anchors, refdata = hsAT$cell_type2,
                            dims = 1:30)
snData_sub <- AddMetaData(snData_sub, metadata = predictions)

# Query mapping on to UMAP embeddings of reference data set
hsAT <- RunUMAP(hsAT, dims = 1:30, reduction = "pca", return.model = TRUE)
snData_sub <- MapQuery(anchorset = epiAT.anchors, reference = hsAT, query = snData_sub,
                       refdata = list(celltype = "cell_type2"), reference.reduction = "pca", reduction.model = "umap")

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

saveRDS(snData_sub, "../omAT_Lean_MTSS.rds")

