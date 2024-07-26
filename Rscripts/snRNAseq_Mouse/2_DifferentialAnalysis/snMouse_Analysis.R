#### This code for snRNAseq plots (after integration) of mouse epiAT data used in Figures 3 and 6 and associated Extended Data. It also contains the "explanation analysis" from Figure 6.


###load packages
library(cowplot)
library(ezRun)
library(ggplot2)
library(kableExtra)
library(Matrix)
library(plotly)
library(RColorBrewer)
library(scater)
library(Seurat)
library(SingleR)
library(tidyverse)
library(SCpubr)
library(ggdist)
library(biomaRt)
library(dplyr)
library(pheatmap)
library(Rsubread)
library(scales)
library(tibble)
library(tidyverse)
library(ggfortify)
library(ggplot2)
library(openxlsx)
library(readxl)
library(janitor)
library(ggpubr)

###load object
### this object is already integrated and QCed
scIntegrated <- readRDS("../snMulti.rds")
DefaultAssay(scIntegrated) <- "RNA"

pvalue_allMarkers <- 0.01
pvalue_all2allMarkers <- 0.01

###assign new names to clusters depending on plot

specific.cellTypes.2 <- c("APCs", "NPVMs", "LAMs", "Adipocytes", "MesoCs", "FIPS", "PVMs", "EpiCs","EndoCs", "Monocytes", "Tcells", "Bcells", "LECs", "P-LAMs", "SMCs", "Adipocytes_2","DCs","MastCs")

clubbed.cell.types <- c("APCs", "Macrophages", "Macrophages", "Adipocytes", "MesoCs", "FIPS", "Macrophages", "EpiCs","EndoCs", "Monocytes", "Tcells", "Bcells", "LECs", "Macrophages", "SMCs", "Adipocytes_2","DCs","MastCs")


Idents(scIntegrated) <- scIntegrated$specific_cellTypes
scIntegrated$specific_cellTypes2 <- Idents(scIntegrated)
names(specific.cellTypes.2) <- levels(scIntegrated$specific_cellTypes)
Idents(scIntegrated) <- scIntegrated$specific_cellTypes2
scIntegrated <-RenameIdents(scIntegrated, specific.cellTypes.2)
scIntegrated$specific_cellTypes2 <- Idents(scIntegrated)

Idents(scIntegrated) <- scIntegrated$specific_cellTypes

scIntegrated$major_celltypes <- Idents(scIntegrated)
names(clubbed.cell.types) <- levels(scIntegrated$specific_cellTypes)
Idents(scIntegrated) <- scIntegrated$major_celltypes
scIntegrated <-RenameIdents(scIntegrated, clubbed.cell.types)
scIntegrated$major_celltypes <- Idents(scIntegrated)

### define color palette
cellTypes2 <- as.character(sort(unique(scIntegrated$specific_cellTypes2)))
xx <- as.character(brewer.pal(12, "Set3"))
yy <- as.character(brewer.pal(8, "Set2"))
zz <- as.character(brewer.pal(6, "Set2"))
colPalette2 <- c(xx, yy, zz)
colPalette2 <- colPalette2[c(1:22)]
names(colPalette2) <- cellTypes2

cellTypes <- as.character(sort(unique(scIntegrated$major_celltypes)))
xx <- as.character(brewer.pal(12, "Set3"))
yy <- as.character(brewer.pal(8, "Set2"))
zz <- as.character(brewer.pal(6, "Set2"))
colPalette <- c(xx, yy, zz)
colPalette <- colPalette[c(1,3:18)]
names(colPalette) <- cellTypes

### plot UMAP

UMAP_broad <- SCpubr::do_DimPlot(scIntegrated,
                   legend.position = "bottom",
                   #order = rev(cellTypes),
                   colors.use = colPalette,
                   label = TRUE,
                   label.color = "black",
                   font.size = 28)



### Proportions barplot
Proportions <- SCpubr::do_BarPlot(scIntegrated, 
                                  group.by = "major_celltypes",
                                  split.by = "condition",
                                  legend.position = "right",
                                  colors.use = colPalette,
                                  position = "fill",
                                  ylab = "Cell type frequency",
                                  legend.title = "Cell types",
                                  flip = FALSE, font.size = 18)

## save proportions
Prop_specific <-prop.table(table(Idents(scIntegrated),scIntegrated$orig.ident), margin = 2)
Prop_specific<- Prop_specific[, c(2, 1, 5, 4, 3, 7,6)]
colnames(Prop_specific) <- c("C", "CC","CCC", "H", "HC", "HH", "HHC")
write.csv(Prop_broad, ".../Proportions.broad.csv")

### DotPlot for markers
Single_Markers <- c("Pdgfra", "Adgre1","Plin1", "Gpm6a", "Cd55","Ank3","Mecom","Ccr2", "Themis", "Pax5",  "Prox1" ,"Myh11", "Clec9a",  "Cpa3")

SCpubr::do_DotPlot(sample = scIntegrated, 
                   features = Single_Markers, font.size = 16, cluster.idents = F,
                   viridis_color_map = "A", use_viridis = T, 
                   flip = T)

#### Macrophages
### subset main object 
scMacro <- subset(scIntegrated, idents = "Macrophages")
Idents(scMacro) <- scMacro$specific_cellTypes2

### plot UMAP
SCpubr::do_DimPlot(scMacro,
                         legend.position = "bottom",
                         #order = rev(cellTypes),
                         colors.use = colPalette2,
                         label = TRUE,
                         label.color = "black",
                         font.size = 30)

### Proportions; percentage of total macrophage content was added manually to the plot
SCpubr::do_BarPlot(scMacro, 
                   group.by = "specific_cellTypes2",
                   split.by = "condition",
                   legend.position = "bottom",
                   #plot.title = "Relative cell type proportions",
                   colors.use = colPalette2,
                   position = "fill",
                   ylab = "Macrophage subclass frequency of total macrophages",
                   legend.title = "Subclasses",
                   flip = FALSE, font.size = 24)

### DotPlot for Markers
Mac_Markers <- c("Ms4a7","Cd74", "Trem2", "Plin2","Cd163" , "Lyve1","Top2a", "Kif15")

SCpubr::do_DotPlot(sample = scMacro, 
                   features = Mac_Markers, font.size = 16, cluster.idents = F,
                   viridis_color_map = "A", use_viridis = T, 
                   flip = T)


### Differential analysis per cell cluster (every comparison is done separately)

# Create new variable combining condition and cell type info
scIntegrated@meta.data$cellTypeCond <-  paste(scIntegrated$group, scIntegrated$specific_cellTypes, sep = "_")

# Update cell idents
Idents(scIntegrated) <- "cellTypeCond"

# Define cell types for which DEGs are going to computed
cellType <- c(levels(scIntegrated$clubbed.cellTypes))

# Intialize empty list to dtore DEGs
diffGenes <- list()

# For each defined cell type, compute DEGs
for(eachCluster in cellType)
{
  # Update the conditions for DE analysis
  markersEach <- try(FindMarkers(scIntegrated, ident.1=paste0("H_", eachCluster),
                                 ident.2=paste0("C_", eachCluster), test.use = "wilcox"))
  
  # Skip cell types with few cells
  if(class(markersEach) != "try-error"){
    diffGenes[[eachCluster]] <- as_tibble(markersEach, rownames="gene")
  }
}

# Combine cell type specific DEGs
diffGenes <- bind_rows(diffGenes, .id="specific_cellTypes")

# Calculate pct difference
# pct.1: cell fraction expressing DE gene in cond 1
# pct.2: cell fraction expressing DE gene in cond 1
diff_pct = abs(diffGenes$pct.1-diffGenes$pct.2)
diffGenes$diff_pct <- diff_pct
diffGenes <- diffGenes[order(diffGenes$diff_pct, decreasing = TRUE),]
rownames(diffGenes) <- NULL

# Only keep the significant DEGs and only abs(logFC) > 0.5 and p.adj < 0.05 (can be adjusted if needed)
diffGenes <- subset(diffGenes, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
diffGenes<- diffGenes %>%
  mutate(regulation = case_when(
    diffGenes$avg_log2FC >0 ~ "Up",
    diffGenes$avg_log2FC<0 ~ "Down"
  ))

### save files
### this code counts the DEGs for up and down per cell type and generates a csv summary
write_csv((tabyl(diffGenes, specific_cellTypes2, regulation)), "/HvC_DEGs_summary.csv")
### this code lists the DEGs, logFC, p.adj etc per cell type
write_csv(diffGenes, ".../HvC_DEGs.csv")


#### Transcriptional Memory by cell type
## read Seurat DE output per cell type (example)
H_vs_C_sigDEGs <- read_csv(".../H_vs_C_sigDEGs.csv")
HC_vs_CC_sigDEGs <- read_csv("...HC_vs_CC_sigDEGs.csv")
## factorize cell type
H_vs_C_sigDEGs$specific_cellTypes <- as.factor(H_vs_C_sigDEGs$specific_cellTypes)
HC_vs_CC_sigDEGs$specific_cellTypes <- as.factor(HC_vs_CC_sigDEGs$specific_cellTypes)

## split in up and down
H_vs_C_sigDEGs_up <- subset(H_vs_C_sigDEGs, regulation == "Up")
HC_vs_CC_sigDEGs_up <- subset(HC_vs_CC_sigDEGs, regulation == "Up")
H_vs_C_sigDEGs_down <- subset(H_vs_C_sigDEGs, regulation == "Down")
HC_vs_CC_sigDEGs_down <- subset(HC_vs_CC_sigDEGs, regulation == "Down")

cellTypes <- levels(H_vs_C_sigDEGs$specific_cellTypes)
## initialize lists
recov.down <- list()
recov.up <- list()
non.recov.down <- list()
non.recov.up <- list()

### get overlaps in numbers
for(i in cellTypes)
{
  
  non.recov.up[[i]] <- length(intersect(H_vs_C_sigDEGs_up[H_vs_C_sigDEGs_up$specific_cellTypes ==i,]$gene, HC_vs_CC_sigDEGs_up[HC_vs_CC_sigDEGs_up$specific_cellTypes ==i,]$gene))
  
}

for(i in cellTypes)
{
  
  non.recov.down[[i]] <- length(intersect(H_vs_C_sigDEGs_down[H_vs_C_sigDEGs_down$specific_cellTypes ==i,]$gene, HC_vs_CC_sigDEGs_down[HC_vs_CC_sigDEGs_down$specific_cellTypes ==i,]$gene))
  
}

for(i in cellTypes)
{
  
  recov.down[[i]] <- length(setdiff(H_vs_C_sigDEGs_down[H_vs_C_sigDEGs_down$specific_cellTypes ==i,]$gene,HC_vs_CC_sigDEGs_down[HC_vs_CC_sigDEGs_down$specific_cellTypes ==i,]$gene))
  
}

for(i in cellTypes)
{
  
  recov.up[[i]] <- length(setdiff(H_vs_C_sigDEGs_up[H_vs_C_sigDEGs_up$specific_cellTypes ==i,]$gene, HC_vs_CC_sigDEGs_up[HC_vs_CC_sigDEGs_up$specific_cellTypes ==i,]$gene))
  
}

### This code counts the number of DEGs that are retained or not per cell type
summary_counts <- as.data.frame(rbind(recov.up, recov.down, non.recov.down, non.recov.up))
fwrite(summary_counts,".../Counts_Short_Memory.csv")

### get gene names 

for(i in cellTypes)
{
   non.recov.up[[i]] <- print(intersect(H_vs_C_sigDEGs_up[H_vs_C_sigDEGs_up$specific_cellTypes ==i,]$gene, 
                                        HC_vs_CC_sigDEGs_up[HC_vs_CC_sigDEGs_up$specific_cellTypes ==i,]$gene))
}

for(i in cellTypes)
{
  non.recov.down[[i]] <- print(intersect(H_vs_C_sigDEGs_down[H_vs_C_sigDEGs_down$specific_cellTypes ==i,]$gene, 
                                         HC_vs_CC_sigDEGs_down[HC_vs_CC_sigDEGs_down$specific_cellTypes ==i,]$gene))
}

for(i in cellTypes)
{
  recov.down[[i]] <- print(setdiff(H_vs_C_sigDEGs_down[H_vs_C_sigDEGs_down$specific_cellTypes ==i,]$gene,
                                   HC_vs_CC_sigDEGs_down[HC_vs_CC_sigDEGs_down$specific_cellTypes ==i,]$gene))
}

for(i in cellTypes)
{
  recov.up[[i]] <- print(setdiff(H_vs_C_sigDEGs_up[H_vs_C_sigDEGs_up$specific_cellTypes ==i,]$gene,
                                 HC_vs_CC_sigDEGs_up[HC_vs_CC_sigDEGs_up$specific_cellTypes ==i,]$gene))
}

summary_genes <- as.data.frame(rbind(recov.up, recov.down, non.recov.down, non.recov.up))
fwrite(summary_genees,".../Genes_Short_Memory.csv")

## calculate percentages in excel and load sheet
Recovery_Short <- read_csv(".../Counts_Short_Memory.csv")
# Upregulated memory; pie chart done per cell type; exchange cell type with another for other pie charts

ggpie(Recovery_Short[5:6,], "Adipocytes", label = paste0(round(100*(Recovery_long[5:6,]$Adipocytes)),"%"),
        lab.pos = "in", lab.font = "black",
        fill = "state", color = "black",
        palette = c("#FF9300", "#069153")) + ggtitle("Adipocytes", subtitle = paste0("n = ", Recovery_long[9,]$Adipocytes))

# Downregulated memory, pie chart done per cell type

ggpie(Recovery_short[7:8,], "Adipocytes", label = paste0(round(100*(Recovery_short[7:8,]$Adipocytes)),"%"),
      lab.pos = "in", lab.font = "black",
      fill = "state", color = "black",
      palette = c("#FF9300", "#069153")) + ggtitle("Adipocytes", subtitle = paste0("n = ", Recovery_short[10,]$Adipocytes))



### Violinplots of adipocytes

## make gene lists
a <- c("Ctsc", "Maob", "Adcy8","Dnah3")
b<- c("Myo1b", "Abca8a", "Gsta3", "Cyp2e1")
c <- c("Gpam", "Acacb", "Cdk8", "Gulp1")
d <- c("Cd74", "Ctss", "Runx2", "Lyz2")
e<- c("Ctsd", "Apobec1", "Mmp12", "Gpnmb")
f <- c("Tyrobp", "Ltc4s", "Mgp", "Tcim")

## plot
pa <- VlnPlot(scIntegrated, features=c(a, e), split.by="condition", 
              stack=T, flip=T, idents = "Adipocytes", 
              cols = c("#7A81FF", "#8E67FF", "#8D2AC5", "#FF9300", "#069153", "#FF6A02", "#009192"))


pb <- VlnPlot(scIntegrated, features=c(b, d), split.by="condition", 
              stack=T, flip=T, idents = "Adipocytes", 
              cols = c("#7A81FF", "#8E67FF", "#8D2AC5", "#FF9300", "#069153", "#FF6A02", "#009192"))


pc <- VlnPlot(scIntegrated, features=c(c, f), split.by="condition", 
              stack=T, flip=T, idents = "Adipocytes", 
              cols = c("#7A81FF", "#8E67FF", "#8D2AC5", "#FF9300", "#069153", "#FF6A02", "#009192"))


pa|pb|pc


#### GSEA

Recovery_Short_Genes <- read_csv(".../Genes_Short_Memory.csv")

## load EnrichR

# Set necessary enrichR global options. This is copied from EnrichR code to avoid having to load the package.
suppressMessages({
  options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")
  options(enrichR.live = TRUE)
  options(modEnrichR.use = TRUE)
  options(enrichR.sites.base.address = "https://maayanlab.cloud/")
  options(enrichR.sites = c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr"))
  
  # Set the search to mammalian genes.
  enrichR::setEnrichrSite(site = "Enrichr")
  
  websiteLive <- TRUE
  dbs <- enrichR::listEnrichrDbs()
  # Get all the possible databases to query.
  dbs <- sort(dbs$libraryName)
})
# Choose the datasets to query against. Pick one or run all
dbs_use <- c("Reactome_2022",
             "GO_Cellular_Component_2021", 
              "WikiPathways_2019_Mouse" )


## Retrieve the enriched terms and plot them for the retained terms; done per cell type and per comparison

# Adipocytes
enriched_terms <- enrichR::enrichr(unlist(Recovery_Short_Genes$Adipocytes[4]), dbs_use)
## save all hits with pVal etc in excel 
write.xlsx(enriched_terms, ".../Adipo_short_non_recov_up_terms.xlsx")
## plot circular plots
SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms, plot.title = NULL, plot.subtitle = "Upregulated non restored in adipocytes that remain from H to HC",  text_labels_size = 5, colors.use = c("red", "#FFFACD"), nchar_wrap = 20, nterms = 8)

#############################################################################################################################################

#### This code was used to generate the plots from Fig. 6 and associated Extended Data

###load object
scIntegrated <- readRDS("../snMulti_YoYo.rds")
DefaultAssay(scIntegrated) <- "RNA"

pvalue_allMarkers <- 0.01
pvalue_all2allMarkers <- 0.01

### rename clusters

Idents(scIntegrated_yoyo) <- scIntegrated_yoyo$seurat_clusters

specific.cellTypes.2 <- c("APCs", "LAMs", "Adipocytes", "NPVMs","PVMs", "MesoCs", "EpiCs", "FIPS","EndoCs", "Monocytes", "Tcells", "EndoCs", "LECs", "P-LAMs", "Bcells", "SMCs","DCs", "SMCs", "EpiCs", "MastCs")
scIntegrated_yoyo$specific_cellTypes2 <- Idents(scIntegrated_yoyo)
names(specific.cellTypes.2) <- levels(scIntegrated_yoyo$seurat_clusters)
Idents(scIntegrated_yoyo) <- scIntegrated_yoyo$specific_cellTypes2
scIntegrated_yoyo <-RenameIdents(scIntegrated_yoyo, specific.cellTypes.2)
scIntegrated_yoyo$specific_cellTypes2 <- Idents(scIntegrated_yoyo)

### set color palette to match other colors
cellTypes2 <- as.character(sort(unique(scIntegrated_yoyo$specific_cellTypes2)))
xx <- as.character(brewer.pal(12, "Set3"))
yy <- as.character(brewer.pal(8, "Set2"))
zz <- as.character(brewer.pal(6, "Set2"))
colPalette2 <- c("#8DD3C7", "#80B1D3", "#FDB462", "#FB8072", "#D9D9D9","#B3DE69", "#BC80BD", "#FCCDE5", "#CCEBC5", "#FFED6F", "#66C2A5", "#8DA0CB", "#E78AC3","#FC8D62","#A6D854", "#E5C494", "#B3B3B3")
names(colPalette2) <- cellTypes2

### UMAP

SCpubr::do_DimPlot(scIntegrated_yoyo,
                   legend.position = "bottom",
                   #order = rev(cellTypes),
                   colors.use = colPalette2,
                   label = TRUE,
                   label.color = "black",
                   font.size = 24, group.by = "specific_cellTypes2", na.value = "grey99", split.by = "condition2")

### Proportions (use colPalette from above)

SCpubr::do_BarPlot(scIntegrated_yoyo, 
                   group.by = "cellTypes",
                   split.by = "condition2",
                   legend.position = "right",
                   #plot.title = "Relative cell type proportions",
                   colors.use = colPalette,
                   position = "fill",
                   ylab = "Cell type frequency",
                   
                   legend.title = "Cell types",
                   flip = FALSE, font.size = 18)

### Macrophages

Idents(scIntegrated_yoyo) <- scIntegrated_yoyo$specific_cellTypes2
sc_yoyo_macro <- subset(scIntegrated_yoyo, idents = c("LAMs", "PVMs", "NPVMs", "P-LAMs"))


SCpubr::do_BarPlot(sc_yoyo_macro, 
                   group.by = "specific_cellTypes2",
                   split.by = "condition2",
                   legend.position = "right",
                   #plot.title = "Relative cell type proportions",
                   colors.use = colPalette2,
                   position = "fill",
                   ylab = "Cell type frequency",
                   
                   legend.title = "Cell types",
                   flip = FALSE, font.size = 18)

### differentially expressed genes
# Create new variable combining condition and cell type info
scIntegrated_yoyo@meta.data$cellTypeCond <-  paste(scIntegrated_yoyo$condition2, scIntegrated_yoyo$specific_cellTypes2, sep = "_")

# Update cell idents
Idents(scIntegrated_yoyo) <- "cellTypeCond"

# Define cell types for which DEGs are going to computed
cellType <- levels(scIntegrated_yoyo$specific_cellTypes2)

# Intialize empty list to store DEGs
diffGenes <- list()

# For each defined cell type, compute DEGs
for(eachCluster in cellType)
{
  # Update the conditions for DE analysis
  markersEach <- try(FindMarkers(scIntegrated_yoyo, ident.1=paste0("HCH_", eachCluster),
                                 ident.2=paste0("CCH_", eachCluster), test.use = "wilcox"))
  
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
diffGenes <- subset(diffGenes, p_val_adj < 0.05, abs(avg_log2FC) >= 0.5)
diffGenes<- diffGenes %>%
  mutate(regulation = case_when(
    diffGenes$avg_log2FC >0 ~ "Up",
    diffGenes$avg_log2FC<0 ~ "Down"
  ))

## save files 
write.csv(diffGenes, ".../HCH_vs_CCH_DEGs.csv", sep = ",")
write_csv((tabyl(diffGenes, cellType, regulation)), ".../HCH_vs_CCH_DEGs_summary.csv")

##################################################################################################
#### explanation pie charts in Figure 6

### transcriptional regulation at HC time point

### load DE hits from post weight loss and after rebound

HC_vs_CC_sigDEGs <- read_csv(".../HC_vs_CC_sigDEGs.csv")
HC_Adipo <- subset(HC_vs_CC_sigDEGs, HC_vs_CC_sigDEGs$specific_cellTypes == "Adipocytes")
HCH_vs_CCH_DEGs <- read_csv(".../HCH_vs_CCH_DEGs.csv")
Yoyo_Adipo <- subset(HCH_vs_CCH_DEGs, HCH_vs_CCH_DEGs$cellType == "Adipocytes") 

## Up
HC_up_Yoyo_up <- list(HC_upregulated =  subset(HC_Adipo, HC_Adipo$regulation == "Up")$gene, Upregulated_Yoyo = subset(Yoyo_Adipo, Yoyo_Adipo$regulation == "Up")$gene)
ggvenn(HC_up_Yoyo_up) 

## Down
HC_down_Yoyo_down <- list(HC_upregulated =  subset(HC_Adipo, HC_Adipo$regulation == "Down")$gene, Upregulated_Yoyo = subset(Yoyo_Adipo, Yoyo_Adipo$regulation == "Down")$gene)
ggvenn(HC_down_Yoyo_down) 

###write into a txt file

### transcriptional memory
###load transcriptional memory DEGs
Recovery_short_DEGs <- read_csv(".../Genes_Short_Memory.csv")
HCH_vs_CCH_DEGs <- read_csv(".../HCH_vs_CCH_DEGs.csv")
Yoyo_Adipo <- subset(HCH_vs_CCH_DEGs, HCH_vs_CCH_DEGs$cellType == "Adipocytes")

## Up
non.recov.up.up <- list(Memory_upregulated = unlist(str_split(Recovery_short_DEGs$Adipocytes[4], ", ")), Upregulated_Yoyo=subset(Yoyo_Adipo, Yoyo_Adipo$regulation == "Up")$gene)
ggvenn(non.recov.up.up) 

## Down
non.recov.down.down <- list(Memory_downregulated = unlist(str_split(Recovery_short_DEGs$Adipocytes[3], ", ")), downregulated_Yoyo=subset(Yoyo_Adipo, Yoyo_Adipo$regulation == "Down")$gene)
ggvenn(non.recov.down.down) 

###write into a text file

#### epigenetics at HC (= memory)

##### Enhancers

### H3K4me1
### load Peaks and annotate to genes 
H3K4me1_sigDP_HC_vs_CC <- read_csv(".../H3K4me1/HC_vs_CC/H3K4me1_sigDP_HC_vs_CC.csv")
H3K4me1 <- rbind(H3K4me1_sigDP_HC_vs_CC)
H3K4me1 <- makeGRangesFromDataFrame(H3K4me1,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=FALSE,
                                    seqinfo=NULL,
                                    seqnames.field=c("chromosome", "chrom",
                                                     "chr", "chromosome_name"),
                                    start.field="start",
                                    end.field=c("end", "stop"),
                                    strand.field="strand",
                                    starts.in.df.are.0based=FALSE)
H3K4me1 <- annotatePeak(H3K4me1, tssRegion=c(-2000, 2000),  TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")
H3K4me1 <- as.data.frame(H3K4me1)

### ATAC
### load Peaks and annotate to genes
ATAC_sigDP_HC_vs_CC <- read_csv(".../ATAC/HC_vs_CC/ATAC_sigDP_HC_vs_CC.csv")
ATAC <- rbind(ATAC_sigDP_HC_vs_CC)
ATAC <- makeGRangesFromDataFrame(ATAC,
                                 keep.extra.columns=TRUE,
                                 ignore.strand=FALSE,
                                 seqinfo=NULL,
                                 seqnames.field=c("chromosome", "chrom",
                                                  "chr", "chromosome_name"),
                                 start.field="start",
                                 end.field=c("end", "stop"),
                                 strand.field="strand",
                                 starts.in.df.are.0based=FALSE)
ATAC <- annotatePeak(ATAC, tssRegion=c(-2000, 2000),  TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")
ATAC <- as.data.frame(ATAC)

### H3K27ac
### load Peaks and annotate to genes 
H3K27ac_sigDP_HC_vs_CC <- read_csv(".../H3K27ac/HC_vs_CC/H3K27ac_sigDP_HC_vs_CC.csv")
H3K27ac <- rbind(H3K27ac_sigDP_HC_vs_CC)
H3K27ac <- makeGRangesFromDataFrame(H3K27ac,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=FALSE,
                                    seqinfo=NULL,
                                    seqnames.field=c("chromosome", "chrom",
                                                     "chr", "chromosome_name"),
                                    start.field="start",
                                    end.field=c("end", "stop"),
                                    strand.field="strand",
                                    starts.in.df.are.0based=FALSE)
H3K27ac <- annotatePeak(H3K27ac, tssRegion=c(-2000, 2000),  TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")
H3K27ac <- as.data.frame(H3K27ac)


##### Promoters (are arleady annotated and collapsed at gene level)
### load files 
### H3K4me3
H3K4me3_sigDEGPr_HC_vs_CC <- read_csv(".../H3K4me3_sigDEGPr_HC_vs_CC.csv")
H3K4me3 <- rbind( H3K4me3_sigDEGPr_HC_vs_CC)
H4K3me3 <- as.data.frame(H3K4me3)

### H3K27me3
H3K27me3_sigDEGPr_HC_vs_CC <- read_csv(".../H3K27me3_sigDEGPr_HC_vs_CC.csv")
H3K27me3 <- rbind(H3K27me3_sigDEGPr_HC_vs_CC)
H3K27me3 <- as.data.frame(H3K27me3)

### H3K27ac
H3K27ac_sigDEGPr_HC_vs_CC <- read_csv(".../H3K27ac_sigDEGPr_HC_vs_CC.csv")
H3K27ac_pr <- rbind(H3K27ac_sigDEGPr_HC_vs_CC)
H3K27ac_pr <- as.data.frame(H3K27ac_pr)



### Calculate overlaps

short_up_final_Yoyo <- list(Upregulated_Yoyo = subset(Yoyo_Adipo, Yoyo_Adipo$regulation == "Up")$gene, 
                            H3K4me1_up =intersect(subset(Yoyo_Adipo, Yoyo_Adipo$regulation == "Up")$gene,subset(H3K4me1, H3K4me1$regulation == "up")$SYMBOL), 
                            
                            ATAC_up = intersect(subset(Yoyo_Adipo, Yoyo_Adipo$regulation == "Up")$gene, subset(ATAC, ATAC$regulation == "up")$SYMBOL), 
                            H3K27ac_up =intersect(subset(Yoyo_Adipo, Yoyo_Adipo$regulation == "Up")$gene, subset(H3K27ac, H3K27ac$regulation == "up")$SYMBOL), 
                            H3K4me3 = intersect(subset(Yoyo_Adipo, Yoyo_Adipo$regulation == "Up")$gene, subset(H3K4me3, H3K4me3$regulation == "up")$external_gene_name),  
                            H3K27me3 = intersect(subset(Yoyo_Adipo, Yoyo_Adipo$regulation == "Up")$gene, subset(H3K27me3, H3K27me3$regulation == "down")$external_gene_name),
                            H3K27ac_pr = intersect(subset(Yoyo_Adipo, Yoyo_Adipo$regulation == "Up")$gene, subset(H3K27ac_pr, H3K27ac_pr$regulation == "up")$external_gene_name))

expla_up <- unique(unlist(short_up_final_Yoyo[2:7]))
write.csv(expla_up, ".../Yoyo_Up_genes.explained.csv")

non_expla_up <- print(setdiff(unlist(short_up_final_Yoyo[1]), expla_up))

write.csv(non_expla_up, ".../Yoyo_Up_genes.non.explained.csv")

venn::venn(short_up_final_Yoyo, zcolor = "style")+ title("Can epigenetics explain upregulation in yoyo?")


short_down_final_Yoyo <- list(downregulated_Yoyo = subset(Yoyo_Adipo, Yoyo_Adipo$regulation == "Down")$gene, 
                              H3K4me1_down =intersect(subset(Yoyo_Adipo, Yoyo_Adipo$regulation == "Down")$gene, subset(H3K4me1, H3K4me1$regulation == "down")$SYMBOL), 
                              
                              ATAC_down = intersect(subset(Yoyo_Adipo, Yoyo_Adipo$regulation == "Down")$gene, subset(ATAC, ATAC$regulation == "down")$SYMBOL), 
                              H3K27ac_down =intersect(subset(Yoyo_Adipo, Yoyo_Adipo$regulation == "Down")$gene, subset(H3K27ac, H3K27ac$regulation == "down")$SYMBOL), 
                              H3K4me3 = intersect(subset(Yoyo_Adipo, Yoyo_Adipo$regulation == "Down")$gene, subset(H3K4me3, H3K4me3$regulation == "down")$external_gene_name),  
                              H3K27me3 = intersect(subset(Yoyo_Adipo, Yoyo_Adipo$regulation == "Down")$gene, subset(H3K27me3, H3K27me3$regulation == "up")$external_gene_name),
                              H3K27ac_pr = intersect(subset(Yoyo_Adipo, Yoyo_Adipo$regulation == "Down")$gene, subset(H3K27ac_pr, H3K27ac_pr$regulation == "down")$external_gene_name))

expla_down <- unique(unlist(short_down_final_Yoyo[2:7]))
write.csv(expla_down, ".../Yoyo_down_genes.explained.csv")

non_expla_down <- print(setdiff(unlist(short_down_final_Yoyo[1]), expla_down))

write.csv(non_expla_down, ".../Yoyo_down_genes.non.explained.csv")

venn::venn(short_down_final_Yoyo, zcolor = "style") + title("Can epigenetics explain downregulation in yoyo?")

#### write into a txt file and combine with txt files from other analysis


##### plot pie charts

### read files 

Yoyo_Percent_Overlaps <- read_excel(".../Yoyo_Percent_Overlaps.xlsx", 
                                    col_types = c("numeric", "text", "text", 
                                                  "text", "numeric", "numeric"))

ggpie(subset(Yoyo_Percent_Overlaps, Yoyo_Percent_Overlaps$Category == "Up" & Yoyo_Percent_Overlaps$Modality == "Epigenetics_HC"), "Percent", label = paste0(round(100*(subset(Yoyo_Percent_Overlaps, Yoyo_Percent_Overlaps$Category == "Up" & Yoyo_Percent_Overlaps$Modality == "Epigenetics_HC")$Percent)),"%"),lab.pos = "in", lab.font = "black",
            fill = "Explained", color = "black",
            palette = c("#0047AB", "#D2042D")) 



