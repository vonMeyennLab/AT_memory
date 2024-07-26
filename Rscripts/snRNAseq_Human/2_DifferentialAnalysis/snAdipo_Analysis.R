#### This code can be used to generate results and plots in Figure 1 and associated figures for integrated adipocytes

#### load libraries
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
#library(xlsx)
library(enrichR)
require(openxlsx)
library(janitor)
library(ggpubr)

### load object of integrated adipocytes from omAT or scAT

scIntegrated <- readRDS(".../snAdipo.rds")
DefaultAssay(scIntegrated) <- "RNA"

pvalue_allMarkers <- 0.01
pvalue_all2allMarkers <- 0.01

Idents(scIntegrated) <- scIntegrated$cellType
Idents(scIntegrated) <- scIntegrated$cellType

### UMAP
SCpubr::do_DimPlot(scIntegrated,
                         legend.position = "bottom",
                         group.by = "cohort",
                         #split.by = "cohort",
                         #order = rev(cellTypes),
                         colors.use = c(L1 = "#c6d7eb", L2= "#d9a5b3"),
                         label = F,
                         label.color = "black",
                         font.size = 24,
                         na.value = "white")

#### perform differential analysis
# Define cell types for which DEGs are going to computed
cellType <- c(levels(Idents(scIntegrated)))

# Create new variable combining condition and cell type info
scIntegrated@meta.data$cellTypeCond <-  paste(scIntegrated$condition, scIntegrated$cellType, sep = "_")

# Update cell idents
Idents(scIntegrated) <- scIntegrated$cellTypeCond

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

# Only keep the significant DEGs and only abs(logFC) > 0.5 and at least 10% of cells need to express it in either or both. Adjust p-value according to what is needed 
diffGenes <- subset(diffGenes, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5 & (diffGenes$pct.1 > 0.1 | diffGenes$pct.2 > 0.1))
diffGenes<- diffGenes %>%
  mutate(regulation = case_when(
    diffGenes$avg_log2FC >0 ~ "Up",
    diffGenes$avg_log2FC<0 ~ "Down"
  ))


### remove lncRNAs or other nonannotated ORFs because we noticed that in some cohorts we obtain more reads for such kind of genes than in others. Likely this is a technical effect. 
new_d <- diffGenes[grep("-AS", diffGenes$gene, invert = TRUE), ]
new_d<- new_d[grep("LINC", new_d$gene, invert = TRUE), ]
new_d<- new_d[grep("MIR", new_d$gene, invert = TRUE), ]
new_d<- new_d[grep("NEAT", new_d$gene, invert = TRUE), ]
new_d<- new_d[grep("XIST", new_d$gene, invert = TRUE), ]
new_d<- new_d[grep("MALAT", new_d$gene, invert = TRUE), ]
new_d<- new_d[grep("orf", new_d$gene, invert = TRUE), ]
new_d<- new_d[grep("MT-", new_d$gene, invert = TRUE), ]
pattern <- "^[A-Z]{2}\\d{6}\\.\\d$"
new_d<- new_d[grep(pattern, new_d$gene, invert = TRUE), ]

### write table with outputs: one with counts and one with DEGs
write_csv((tabyl(new_d, cellType, regulation)), ".../res_t0_vs_lean_noLNC.count.csv")
write_csv(new_d, ".../res_t0_vs_lean_noLNC.csv")

#### perform transcriptional retention analysis

###memory (for individuals or for bulk per cell type)
# load files
t0_L_R_vs_t0_L_C_DEGs <- read_csv(".../res_t0_vs_lean_noLNC.csv")
t1_L_R_vs_t0_L_C_DEGs <- read_csv(".../res_t1_vs_lean_noLNC.csv")

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


summary_counts <- as.data.frame(rbind(recov.up, recov.down, non.recov.down, non.recov.up))
fwrite(summary_responder_counts, ".../Memory_counts.csv")
### calculate percentages in excel or similar software

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

summary__DEGs <- as.data.frame(rbind(recov.up, recov.down, non.recov.down, non.recov.up))
fwrite(summary_DEGs, ".../Memory_DEGs.csv")

#### plot in pie chart

Memory_counts_sum <- read_csv(".../Memory_counts.csv")

p1 <- ggpie(Memory_counts_sum[7:8,], "Adipo", label = paste0(round(100*(Memory_counts_sum[7:8,]$Adipo)),"%"),
            lab.pos = "in", lab.font = "black",
            fill = "State", color = "black",
            palette = c("#F47A60", "#316879")) + ggtitle("omAT adipocytes", subtitle = paste0("n = ", Memory_counts_sum[10,]$Adipo))


p2 <- ggpie(Memory_counts_sum[5:6,], "Adipo", label = paste0(round(100*(Memory_counts_sum[5:6,]$Adipo)),"%"),
            lab.pos = "in", lab.font = "black",
            fill = "State", color = "black",
            palette = c("#F47A60", "#316879")) + ggtitle("omAT adipocytes", subtitle = paste0("n = ", Memory_counts_sum[9,]$Adipo))



p <- p1 | p2

p

#### ViolinPlot of memory DEGs

a <- c("IGF1", "LPIN1", "IDH1", "PDE3A")

VlnPlot(scIntegrated, features= c(a), split.by="condition", 
        stack=T, flip=T, idents = "Adipo", 
        cols =c("#8a307f", "#F47A60", "#316879"))


