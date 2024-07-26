#### This code can be used to generate results and plots found in Figure 1 and associated extended data

### load libraries
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
library(readxl)
library(biomaRt)
library(enrichR)
require(openxlsx)
library(janitor)
library(ggpubr)

### load integrated object for pooled omAT/scAT
### add metadata to assign donors 

scIntegrated <- readRDS("snHumanAT.rds")
DefaultAssay(scIntegrated) <- "RNA"

pvalue_allMarkers <- 0.01
pvalue_all2allMarkers <- 0.01


### add donor information if needed to have consistent names
scIntegrated$Donor <- scIntegrated$best_singlet
Idents(scIntegrated) <- scIntegrated$Donor
Donor <- c("Lean_D4","Lean_D1","Lean_D3", "Lean_D2","Lean_D0","Res_T0_D4","Res_T0_D0","Res_T0_D2","Res_T0_D1", "Res_T0_D3", "Res_T1_D0","Res_T1_D1","Res_T1_D4","Res_T1_D2","Res_T1_D3")
names(Donor) <- levels(Idents(scIntegrated))
Donor
scIntegrated <-RenameIdents(scIntegrated, Donor)
scIntegrated$Donor <- Idents(scIntegrated)
scIntegrated$Donor <- factor(scIntegrated$Donor, levels = c("Lean_D0","Lean_D1", "Lean_D2",
                                                            "Lean_D3","Lean_D4",
                                                            "Res_T0_D0", "Res_T0_D1", 
                                                            "Res_T0_D2",
                                                            "Res_T0_D3","Res_T0_D4",
                                                            "Res_T1_D0",
                                                            "Res_T1_D1",
                                                            "Res_T1_D2",
                                                            "Res_T1_D3","Res_T1_D4"))

### To do per person DE one needs to "merge' the control; this is used for Extended Data Figure 1 heatmaps
Idents(scIntegrated) <- scIntegrated$best_singlet
scIntegrated$Pseudo <- scIntegrated$best_singlet
Idents(scIntegrated) <- scIntegrated$Pseudo
Pseudo <- c("Lean", "Lean", "Lean", "Lean", "Lean","Res_T0_D4","Res_T0_D0","Res_T0_D2","Res_T0_D1", "Res_T0_D3", "Res_T1_D0","Res_T1_D1","Res_T1_D4","Res_T1_D2","Res_T1_D3")
names(Pseudo) <- levels(Idents(scIntegrated))
Pseudo
scIntegrated <-RenameIdents(scIntegrated, Pseudo)
scIntegrated$Pseudo <- Idents(scIntegrated)
scIntegrated$Pseudo <- factor(scIntegrated$Pseudo, levels = c("Lean",
                                                              "Res_T0_D0", "Res_T0_D1", 
                                                              "Res_T0_D2",
                                                              "Res_T0_D3","Res_T0_D4",
                                                              "Res_T1_D0",
                                                              "Res_T1_D1",
                                                              "Res_T1_D2",
                                                              "Res_T1_D3",
                                                              "Res_T1_D4"))

### order cell types
Idents(scIntegrated) <- scIntegrated$cellType
cellTypes_2 <- c("Adipo", "APCs1","APCs2","FAPs","Macro","DCs","Tcells","Bcells","EndoVCs","EndoSCs","EndoACs","LECs","MesoCs","SMCs","Peri","MastCs","NeurCs1","NeurCs2")
names(cellTypes_2)<-levels(Idents(scIntegrated))
scIntegrated <- RenameIdents(scIntegrated, cellTypes_2)
scIntegrated$cellType <- Idents(scIntegrated)
scIntegrated$cellType <- factor(scIntegrated$cellType, levels = c("Adipo", "APCs1","APCs2","FAPs", "Macro","DCs","Tcells","Bcells","EndoVCs","EndoSCs","EndoACs","LECs","MesoCs","SMCs","Peri","MastCs","NeurCs1", "NeurCs2"))


### Plots
## load colour palette
colPalette <- 
  c("Adipo" =  "#8DD3C7",
    "APCs1" = "#FFFFB3",
    "APCs2" = "#BEBADA",
    "FAPs" = "#FB8072",
    "DCs"   = "#FDB462",
    "Macro" = "#80B1D3",
    # "Macro2" = "#B3B3B3",
    # "Macro3" = "#FCFD55",
    "Tcells" = "#B3DE69",
    "EndoVCs" = "#D9D9D9",
    "EndoSCs" = "#BC80BD",
    "EndoACs" = "#CCEBC5",
    "Peri" =  "#8DA0CB",
    "MastCs" = "#E78AC3",
    "LECs" = "#FFED6F",
    "SMCs"   =  "#FC8D62",
    "Bcells" = "#FCCDE5",
    "NeurCs1" = "#A6D854",
    "NeurCs2" =  "#FFD92F",
    # "ASDCs" = "#E5C494",
    "MesoCs" = "#66C2A5")

## UMAP
Idents(scIntegrated) <- scIntegrated$cellType
SCpubr::do_DimPlot(scIntegrated,
                   legend.position = "bottom",
                   group.by = "cellType",
                   #order = rev(cellTypes),
                   colors.use = colPalette,
                   label = TRUE,
                   label.color = "black",
                   font.size = 24)

## stacked barplot by condition
SCpubr::do_BarPlot(scIntegrated, 
                   group.by = "cellType",
                   split.by = "condition",
                   legend.position = "right",
                   #plot.title = "Relative cell type proportions",
                   colors.use = colPalette,
                   position = "fill",
                   ylab = "Cell type frequency",
                   legend.title = "Cell types",
                   flip = FALSE, font.size = 13)

## stacked barplot by Donor
SCpubr::do_BarPlot(scIntegrated, 
                   group.by = "cellType",
                   split.by = "Donor",
                   legend.position = "right",
                   #plot.title = "Relative cell type proportions",
                   colors.use = colPalette,
                   position = "fill",
                   ylab = "Cell type frequency",
                   legend.title = "Cell types",
                   flip = FALSE, font.size = 13)

### write proportions/cell numbers into table
pt <- table(scIntegrated$cellType, scIntegrated$Donor)

###DotPlot for markers

Single_Markers <- c("GPAM","DCN", "CRISPLD2","MFAP5","PDZRN4", "MRC1","TLR2","THEMIS", "BANK1", "MECOM", "CADM2", "TPO", "NKAIN2","MMRN1","PTPRQ", "MYH11","STEAP4",  "CPA3",  "NRXN1")

SCpubr::do_DotPlot(sample = scIntegrated, 
                   features = Single_Markers, font.size = 16,
                   #cluster.idents = F,
                   viridis.palette = "A", use_viridis = T, 
                   flip = T, dot.scale = 7)

### QC plots by cellType

VlnPlot(scIntegrated, features= c("nFeature_RNA", "nCount_RNA"), 
        stack=F, flip=T, pt.size = 0,
        cols = colPalette)

## QC plots by condition
Idents(scIntegrated) <- scIntegrated$condition
VlnPlot(scIntegrated, features= c("nFeature_RNA", "nCount_RNA"), 
        stack=F, flip=T, pt.size = 0,
        cols = c("#8A307F","#F47A60","#316879"))

### QC plots by donor
Idents(scIntegrated) <- scIntegrated$Donor
VlnPlot(scIntegrated, features= c("nFeature_RNA", "nCount_RNA"), 
        stack=F, flip=T, pt.size = 0,
        cols = c("#8A307F","#8A307F","#8A307F","#8A307F","#8A307F", "#F47A60","#F47A60","#F47A60","#F47A60","#F47A60","#316879","#316879","#316879","#316879","#316879"))


Idents(scIntegrated) <- scIntegrated$cellType

#### DE analysis
### this analysis can be done per donor if one chooses Pseudo instead of condition when createing cellTypeCond

# Create new variable combining condition and cell type info
scIntegrated@meta.data$cellTypeCond <-  paste(scIntegrated$condition, scIntegrated$cellType, sep = "_")

# Update cell idents
Idents(scIntegrated) <- "cellTypeCond"

# Define cell types for which DEGs are going to computed
cellType <- c(levels(scIntegrated$cellType))

# Intialize empty list to dtore DEGs
diffGenes <- list()
levels(scIntegrated$condition) #check exact names of conditions

# For each defined cell type, compute DEGs
for(eachCluster in cellType)
{
  # Update the conditions for DE analysis
  markersEach <- try(FindMarkers(scIntegrated, ident.1=paste0("res_t0_", eachCluster),
                                 ident.2=paste0("lean_", eachCluster), test.use = "wilcox"))
  
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

# Only keep the significant DEGs and only abs(logFC) > 1
diffGenes <- subset(diffGenes, p_val_adj < 0.01 & abs(avg_log2FC) > 0.5)
diffGenes<- diffGenes %>%
  mutate(regulation = case_when(
    diffGenes$avg_log2FC >0 ~ "Up",
    diffGenes$avg_log2FC<0 ~ "Down"
  ))

## write diffgenes per cell type and the number of diff genes per cell type
write_csv((tabyl(diffGenes, cellType, regulation)),".../T0_vs_Lean_DEGs_count.csv")
write_csv(diffGenes, ".../T0_vs_Lean_DEGs.csv")


##### Transcriptional Retention analysis

##memory (for individuals or for bulk per cell type)
# load files
t0_L_R_vs_t0_L_C_DEGs <- read_csv(".../T0_vs_Lean.csv")
t1_L_R_vs_t0_L_C_DEGs <- read_csv(".../T1vs_Lean.csv")


t0_L_R_vs_t0_L_C_DEGs$cellType <- as.factor(t0_L_R_vs_t0_L_C_DEGs$cellType)
t1_L_R_vs_t0_L_C_DEGs$cellType <- as.factor(t1_L_R_vs_t0_L_C_DEGs$cellType)

H_vs_C_sigDEGs_up <- subset(t0_L_R_vs_t0_L_C_DEGs, regulation == "Up")
HC_vs_CC_sigDEGs_up <- subset(t1_L_R_vs_t0_L_C_DEGs, regulation == "Up")
H_vs_C_sigDEGs_down <- subset(t0_L_R_vs_t0_L_C_DEGs, regulation == "Down")
HC_vs_CC_sigDEGs_down <- subset(t1_L_R_vs_t0_L_C_DEGs, regulation == "Down")

cellTypes <- levels(t0_L_R_vs_t0_L_C_DEGs$cellType)

### Counts
recov.down <- list()
recov.up <- list()
non.recov.down <- list()
non.recov.up <- list()


for(i in cellTypes)
{
  non.recov.up[[i]] <- length(intersect(H_vs_C_sigDEGs_up[H_vs_C_sigDEGs_up$cellType ==i,]$gene, 
                                        HC_vs_CC_sigDEGs_up[HC_vs_CC_sigDEGs_up$cellType ==i,]$gene))
}

for(i in cellTypes)
{
  non.recov.down[[i]] <- length(intersect(H_vs_C_sigDEGs_down[H_vs_C_sigDEGs_down$cellType ==i,]$gene, 
                                          HC_vs_CC_sigDEGs_down[HC_vs_CC_sigDEGs_down$cellType ==i,]$gene))
}

for(i in cellTypes)
{
  recov.down[[i]] <- length(setdiff(H_vs_C_sigDEGs_down[H_vs_C_sigDEGs_down$cellType ==i,]$gene, 
                                    HC_vs_CC_sigDEGs_down[HC_vs_CC_sigDEGs_down$cellType ==i,]$gene))
}

for(i in cellTypes)
{
  recov.up[[i]] <- length(setdiff(H_vs_C_sigDEGs_up[H_vs_C_sigDEGs_up$cellType ==i,]$gene, 
                                  HC_vs_CC_sigDEGs_up[HC_vs_CC_sigDEGs_up$cellType ==i,]$gene))
}


summary__counts <- as.data.frame(rbind(recov.up, recov.down, non.recov.down, non.recov.up))
fwrite(summary__counts, ".../Memory_counts.csv")
### these files are used to compile the summary file used for barplots 

### DEGs
recov.down <- list()
recov.up <- list()
non.recov.down <- list()
non.recov.up <- list()

for(i in cellTypes)
{
  non.recov.up[[i]] <- c(intersect(H_vs_C_sigDEGs_up[H_vs_C_sigDEGs_up$cellType ==i,]$gene, 
                                   HC_vs_CC_sigDEGs_up[HC_vs_CC_sigDEGs_up$cellType ==i,]$gene))
}

for(i in cellTypes)
{
  non.recov.down[[i]] <- c(intersect(H_vs_C_sigDEGs_down[H_vs_C_sigDEGs_down$cellType ==i,]$gene, 
                                     HC_vs_CC_sigDEGs_down[HC_vs_CC_sigDEGs_down$cellType ==i,]$gene))
}

for(i in cellTypes)
{
  recov.down[[i]] <- c(setdiff(H_vs_C_sigDEGs_down[H_vs_C_sigDEGs_down$cellType ==i,]$gene, 
                               HC_vs_CC_sigDEGs_down[HC_vs_CC_sigDEGs_down$cellType ==i,]$gene))
}

for(i in cellTypes)
{
  recov.up[[i]] <- c(setdiff(H_vs_C_sigDEGs_up[H_vs_C_sigDEGs_up$cellType ==i,]$gene, 
                             HC_vs_CC_sigDEGs_up[HC_vs_CC_sigDEGs_up$cellType ==i,]$gene))
}

summary_DEGs <- as.data.frame(rbind(recov.up, recov.down, non.recov.down, non.recov.up))
fwrite(summary_DEGs, ".../Memory_DEGs.csv")

### GSEA terms for hits 
Memory_Up_human <-summary_DEGs[4,]
Memory_Down_human <- summary_DEGs[3,]

suppressMessages({
  options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")
  options(enrichR.live = TRUE)
  options(modEnrichR.use = TRUE)
  options(enrichR.sites.base.address = "https://maayanlab.cloud/")
  options(enrichR.sites = c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr"))
  
  # Set the search to Human genes.
  enrichR::setEnrichrSite(site = "Enrichr")
  
  websiteLive <- TRUE
  dbs <- enrichR::listEnrichrDbs()
  # Get all the possible databases to query.
  dbs <- sort(dbs$libraryName)
})

dbs_use <- c("WikiPathways_2019_Human")

enriched_terms <- enrichR::enrichr(unlist(Memory_Up_human$Adipo), dbs_use)
write.xlsx(enriched_terms,file = ".../Adipo_resp_nonrecov_up_terms.xlsx")

SCpubr::do_TermEnrichmentPlot(enriched_terms =enriched_terms, plot.title = NULL, plot.subtitle = "Upregulated non restored in adipocytes that remain after WL",  text_labels_size = 4, 
                              colors.use = c("red", "#FFFACD"), nchar_wrap = 15, nterms = 8)

enriched_terms <- enrichR::enrichr(unlist(Memory_Down_human$Adipo), dbs_use)
write.xlsx(enriched_terms,file = ".../Adipo_resp_nonrecov_down_terms.xlsx")

SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms, plot.title = NULL, plot.subtitle = "Downregulated non restored in adipocytes that remain after WL",  text_labels_size = 4, 
                              colors.use = c("red", "#FFFACD"), nchar_wrap = 15, nterms =8)





  
  