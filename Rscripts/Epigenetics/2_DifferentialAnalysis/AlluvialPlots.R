##### This code can be used to generate alluvial plots for promoters or enhancers as found in Figures 4 and 5

#### Load libraries
library(cowplot)
library(ezRun)
library(ggplot2)
library(kableExtra)
library(Matrix)
library(plotly)
library(tidyverse)
library(ggalluvial)
library(biomaRt)
library(dplyr)
library(edgeR)
library(ggfortify)
library(Matrix)
library(matrixStats)
library(pheatmap)
library(RColorBrewer)
library(Rsubread)
library(scales)
library(tibble)
library(WGCNA)
library(ggrepel)
library(openxlsx)
library(readxl)
library(magrittr)
library(purrr)

#### Promoters

### H3K4me3 (for HC/H as an example)

## read DE hits (these contain gene names, up/down regulation, logFC etc.)
HC_vs_CC_DEGs <- read_csv(".../H3K4me3_sigDEGPr_HC_vs_CC.csv")
H_vs_C_DEGs <-  read_csv(".../H3K4me3_sigDEGPr_H_vs_C.csv")
#rename/add columns of up and down
HC_vs_CC_DEGs$score_HCvCC <- HC_vs_CC_DEGs$regulation
H_vs_C_DEGs$score_HvC <- H_vs_C_DEGs$regulation

## join the two datasets based on gene name (promoter counts were collapsed at gene name)

DE_overlap_s <- inner_join(H_vs_C_DEGs, HC_vs_CC_DEGs, by = "external_gene_name") 

## create dfs for promoters that are only significantly different in one comparison

H_vs_C_only <- subset(H_vs_C_DEGs, !(external_gene_name %in% HC_vs_CC_DEGs$external_gene_name))
HC_vs_CC_only <- subset(HC_vs_CC_DEGs, !(external_gene_name %in% H_vs_C_DEGs$external_gene_name))

### add a score 
H_vs_C_only <- H_vs_C_only %>%
  mutate(score_HCvCC = case_when(
    H_vs_C_only$logFC > 0.5 ~ "NS",
    H_vs_C_only$logFC < -0.5 ~ "NS"
  ))

HC_vs_CC_only <- HC_vs_CC_only %>%
  mutate(score_HvC = case_when(
    HC_vs_CC_only$logFC >0.5 ~ "NS",
    HC_vs_CC_only$logFC< -0.5 ~ "NS"
  ))

## overlap the dfs and create df that can be used for an alluvial plot

Overlap_s <- rbind(DE_overlap_s %>% dplyr::select(external_gene_name, score_HvC, score_HCvCC), 
                   H_vs_C_only %>% dplyr::select(external_gene_name, score_HvC, score_HCvCC), 
                   HC_vs_CC_only%>% dplyr::select(external_gene_name, score_HvC, score_HCvCC))
Overlap_s$combined_score <- paste0(Overlap_s$score_HvC, "_", Overlap_s$score_HCvCC)

### add a one for counting in the alluvial plot
Overlap_s <- Overlap_s %>%
  mutate(count = case_when(
    Overlap_s$score_HCvCC == "down" ~ 1,
    Overlap_s$score_HCvCC == "up"~ 1,
    Overlap_s$score_HCvCC == "NS" ~ 1
  ))

### save out put
write_csv(Overlap_s, ".../Overlap_short.csv")

### alluvial plot
ggplot(as.data.frame(Overlap_s),
       aes(y = count , axis1 = score_HvC, axis2 = score_HCvCC)) +
  geom_alluvium(aes(fill = combined_score)) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  geom_stratum(width = 1/12, fill = "grey50", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), label.size = 1)+
  scale_x_discrete(limits = c("HvC", "HCvCC_s"), expand = c(.05, .05))+
  scale_fill_brewer(type = "qual", palette = "Set1")+
  ggtitle("H3K4me3")+theme(plot.title = element_text(size = 30, face = "bold"))+ 
  scale_fill_manual(values=c("#0000FF","#6495ED" ,"honeydew3", "honeydew3","#FAA0A0", "#FF0000")) +
  theme(axis.text = element_text(size = 30), legend.text = element_text(size = 20), axis.title.y = element_text(size=20, face="bold"))


#### Enhancers

### H3K4me1 H/HC as example
## read DE hits; these contain genomic coordinates and annotated genes 
HC_vs_CC_DEGs_up <- read_csv(".../H3K4me1_HCvsCC_upPeaks.csv")
H_vs_C_DEGs_up <-  read_csv(".../H3K4me1_HvsC_upPeaks.csv")

HC_vs_CC_DEGs_down <- read_csv(".../H3K4me1_HCvsCC_downPeaks.csv")
H_vs_C_DEGs_down <-  read_csv(".../H3K4me1_HvsC_downPeaks.csv")

HC_vs_CC_DEGs <-rbind(HC_vs_CC_DEGs_up, HC_vs_CC_DEGs_down)
H_vs_C_DEGs <- rbind(H_vs_C_DEGs_up, H_vs_C_DEGs_down)

#rename/add columns of up and down
HC_vs_CC_DEGs$score_HCvCC <- HC_vs_CC_DEGs$regulation
H_vs_C_DEGs$score_HvC <- H_vs_C_DEGs$regulation

HC_vs_CC_DEGs$loci <- paste0(HC_vs_CC_DEGs$seqnames, ":", HC_vs_CC_DEGs$start, "-", HC_vs_CC_DEGs$end)
H_vs_C_DEGs$loci <- paste0(H_vs_C_DEGs$seqnames, ":", H_vs_C_DEGs$start, "-", H_vs_C_DEGs$end)                                                               

### join based on loci

DE_overlap_s <- inner_join(H_vs_C_DEGs, HC_vs_CC_DEGs, by = "loci") 
DE_overlap_s$SYMBOL <- DE_overlap_s$SYMBOL.y

H_vs_C_only <- subset(H_vs_C_DEGs, !(loci %in% HC_vs_CC_DEGs$loci))
HC_vs_CC_only <- subset(HC_vs_CC_DEGs, !(loci %in% H_vs_C_DEGs$loci))
### add a score 
H_vs_C_only <- H_vs_C_only %>%
  mutate(score_HCvCC = case_when(
    H_vs_C_only$logFC > 0.5 ~ "NS",
    H_vs_C_only$logFC < -0.5 ~ "NS"
  ))


HC_vs_CC_only <- HC_vs_CC_only %>%
  mutate(score_HvC = case_when(
    HC_vs_CC_only$logFC >0.5 ~ "NS",
    HC_vs_CC_only$logFC< -0.5 ~ "NS"
  ))

## overlap the dfs and create df that can be used for an alluvial plot
Overlap_s <- rbind(DE_overlap_s %>% dplyr::select(loci, score_HvC, score_HCvCC, SYMBOL),
                   H_vs_C_only %>% dplyr::select(loci, score_HvC, score_HCvCC, SYMBOL),
                   HC_vs_CC_only%>% dplyr::select(loci, score_HvC, score_HCvCC, SYMBOL))

Overlap_s$combined_score <- paste0(Overlap_s$score_HvC, "_", Overlap_s$score_HCvCC)

### add a 1 for counting in the alluvial plot
Overlap_s <- Overlap_s %>%
  mutate(count = case_when(
    Overlap_s$score_HCvCC == "down" ~ 1,
    Overlap_s$score_HCvCC == "up"~ 1,
    Overlap_s$score_HCvCC == "NS" ~ 1
  ))

### save outputs: one with information on both original data sets  and one for the alluvial plot
write_csv(DE_overlap_s,"...//Overlap_short.csv")
write_csv(Overlap_s, ".../Overlap_short_alluvial.csv")

### alluvial plot

ggplot(as.data.frame(Overlap_s),
       aes(y = count , axis1 = score_HvC, axis2 = score_HCvCC)) +
  geom_alluvium(aes(fill = combined_score)) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  geom_stratum(width = 1/12, fill = "grey50", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), label.size = 1)+
  scale_x_discrete(limits = c("HvC", "HCvCC_s"), expand = c(.05, .05))+
  scale_fill_brewer(type = "qual", palette = "Set1")+
  ggtitle("H3K4me1")+theme(plot.title = element_text(size = 30, face = "bold"))+ 
  scale_fill_manual(values=c("#0000FF","#6495ED" ,"honeydew3", "honeydew3","honeydew3","#FAA0A0", "#FF0000"))+
  theme(axis.text = element_text(size = 30), legend.text = element_text(size = 20), axis.title.y = element_text(size=20, face="bold")) 

#### when comparing HHC and HH we add one more colour: "#CF9FFF"