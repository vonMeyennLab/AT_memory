setwd("/srv/GT/analysis/adhideb/p3360-Ferdinand/")

library(SingleR)
library(xlsx)
library(AUCell)
library(Seurat)
library(pheatmap)
library(cowplot)
library(ggplot2)
library(RColorBrewer)

lchd <- readRDS("/srv/GT/analysis/adhideb/p3360-Ferdinand/snSeurat/S075_6_mmepiAT_LCHD/snLCHD.rds")
lcc <- readRDS("/srv/GT/analysis/adhideb/p3360-Ferdinand/snSeurat/S019_1_mmepiAT_LCC/snLCC.rds")
lcch <- readRDS("/srv/GT/analysis/adhideb/p3360-Ferdinand/snSeurat/S075_2_mmepiAT_LCCH/snLCCH.rds")

#######################################################################################################################################

# Create list of all Seurat objects
snSampleList <- list(lchd, lcc, lcch)

for (i in 1:length(snSampleList)) {
  snSampleList[[i]] <- FindVariableFeatures(snSampleList[[i]], selection.method = "vst", 
                                            nfeatures = 2000, verbose = FALSE)
}

hvgs_per_dataset <- lapply(snSampleList, function(x) {
  x@assays$RNA@var.features
})
venn::venn(hvgs_per_dataset, opacity = 0.4, zcolor = (scales::hue_pal())(3), cexsn = 1, 
           cexil = 1, lwd = 1, col = "white", frame = F, borders = NA)

# Select variable features that are common across data sets for integration
features <- SelectIntegrationFeatures(object.list = snSampleList)

# Identify integration anchors
sampleAnchors <- FindIntegrationAnchors(object.list = snSampleList, anchor.features = features, dims = 1:30,
                                      scale = FALSE, reduction = "cca")

# Integrate data sets
scIntegrated <- IntegrateData(anchorset = sampleAnchors, dims = 1:30, new.assay.name = "CCA")

# Check assay
names(scIntegrated@assays)
scIntegrated@active.assay
scIntegrated$condition <- paste0(scIntegrated$timePoint, "_", scIntegrated$group)

# Run the standard workflow for visualization and clustering
scIntegrated <- ScaleData(scIntegrated, verbose = FALSE)
scIntegrated <- RunPCA(scIntegrated, npcs = 30, verbose = FALSE)
scIntegrated <- RunUMAP(scIntegrated, reduction = "pca", dims = 1:30)
scIntegrated <- RunTSNE(scIntegrated, reduction = "pca", dims = 1:30)
scIntegrated <- FindNeighbors(scIntegrated, reduction = "pca", dims = 1:30)
scIntegrated <- FindClusters(scIntegrated, resolution = 0.4)

DimPlot(scIntegrated, reduction = "umap", label = TRUE)
DimPlot(scIntegrated, reduction = "tsne", label = TRUE)

DimPlot(scIntegrated, reduction = "tsne", group.by = "condition")
DimPlot(scIntegrated, reduction = "umap", group.by = "condition")

DimPlot(scIntegrated, reduction = "umap", split.by = "condition")
DimPlot(scIntegrated, reduction = "tsne", split.by = "condition")

# Update assay
DefaultAssay(scIntegrated) <- "RNA"

DimPlot(scIntegrated, reduction = "tsne", group.by = "SingleR.labels", split.by = "condition")
DimPlot(scIntegrated, reduction = "umap", group.by = "SingleR.labels", split.by = "condition")

DimPlot(scIntegrated, reduction = "umap", group.by = "predicted.celltype", split.by = "condition", label = TRUE, label.size = 3, repel = TRUE)
DimPlot(scIntegrated, reduction = "tsne", group.by = "predicted.celltype", split.by = "condition", label = TRUE, label.size = 3, repel = TRUE)

#######################################################################################################################################

mmRef <- (MouseRNAseqData <- MouseRNAseqData())

singler.results <- SingleR(
  test = GetAssayData(scIntegrated), ref = mmRef,
  labels = mmRef$label.fine)

scIntegrated[["SingleR.labels"]] <- singler.results$labels

#######################################################################################################################################

mmepiAT <- readRDS("mouse_epiAT.rds")

DimPlot(mmepiAT)

# Select common variable features
features <- intersect(mmepiAT$integrated@var.features,scIntegrated$CCA@var.features)

# Reference mapping via transfer anchors
epiAT.anchors <- FindTransferAnchors(reference = mmepiAT, query = scIntegrated,
                                     dims = 1:30, reference.reduction = "pca", features = features)

# Label transfer
predictions <- TransferData(anchorset = epiAT.anchors, refdata = mmepiAT$ct3,
                            dims = 1:30)
scIntegrated <- AddMetaData(scIntegrated, metadata = predictions)

# Query mapping on to UMAP embeddings of reference data set
mmepiAT <- RunUMAP(mmepiAT, dims = 1:30, reduction = "pca", return.model = TRUE)
scIntegrated <- MapQuery(anchorset = epiAT.anchors, reference = mmepiAT, query = scIntegrated,
                       refdata = list(celltype = "ct3"), reference.reduction = "pca", reduction.model = "umap")

scIntegrated <- AddMetaData(object = scIntegrated, metadata = MappingScore(epiAT.anchors,ndim=30), col.name = "mapping.score")

p1 <- DimPlot(mmepiAT, reduction = "umap", group.by = "ct3", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(scIntegrated, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

FeaturePlot(scIntegrated,features = "mapping.score")
FeaturePlot(scIntegrated,features = "predicted.celltype.score")

p1 <- DimPlot(scIntegrated, reduction = "umap",group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE)  + ggtitle("Query transfered labels") + NoLegend()

p2 <- DimPlot(scIntegrated, reduction = "tsne",group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE)  + ggtitle("Query transfered labels") + NoLegend()

p1 + p2

DimPlot(scIntegrated, reduction = "umap",group.by = "predicted.celltype", label = TRUE,
        label.size = 3, repel = TRUE, split.by = "condition")

saveRDS(scIntegrated, "snLCHD_LCC_LCCH.rds")
#scIntegrated <- readRDS("snHFD_HC_YoYo.rds")
