setwd("~/public/_Projects/AT_HFD_Memory/Epigenetic_memory/QC/CnT")

#.libPaths("/scratch/lib/R/4.1.0")

library(dplyr)
library(edgeR)
library(ggfortify)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(tibble)
library(tidyverse)

## Load sample data
metadata <- read.table("CnT_metadata.txt", header = TRUE, sep = "\t")
metadata <- subset(metadata, histone == "H3K4me1" & QC != "No" & timepoint == "AdipoERCre")

# Load raw peak counts
rawCounts <- data.frame(readRDS("H3K4me1_AdipoERCre_refC_filtPeakCounts.rds"))

samples <- colnames(rawCounts)
condition <- as.factor(metadata$condition)
experiment <- as.factor(metadata$experiment)

# Average expression for control reference
rawCounts[, "C_short_H3K4me1_AdipoERCre"] <- rowMeans(rawCounts[,1:3])

# Add sample information
group <- data.frame(condition, experiment, row.names = samples)
group["C_short_H3K4me1_AdipoERCre",] <- c("C", "short")

# Create DGEList object
y <- DGEList(counts=rawCounts, group=group$condition)

# Remove genes with low expression in > 50% samples
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Apply TMM normalization
y <- calcNormFactors(y)

# Extract normalized counts
normCounts <- edgeR::cpm(y, log=TRUE)

# Correlation coefficient
xx <- cor(normCounts, use="complete.obs", method = "spearman")
xx <- merge(xx, group, by="row.names")
xx <- xx[,c(1,24:26)]
xx$cond <- paste0(xx$condition, "_", xx$experiment)
xx[4,"cond"] <- "Control"

xx$cond <- factor(xx$cond, levels = c("Control", "C_short", "CC_short", "CC_long", "CCC_long",
                                      "H_short", "HC_short", "HH_long", "HHC_long"))

ggplot(xx) + 
  geom_bar(aes(cond, C_short_H3K4me1_AdipoERCre, fill = cond), 
           position = "dodge", stat = "summary", fun="mean") +
  #scale_fill_manual(values = my_col) + facet_wrap(. ~ variable, strip.position="bottom", ncol=4) +
  geom_point(aes(x=cond,y=C_short_H3K4me1_AdipoERCre,group=cond),
             position= position_dodge(0.9)) +
  labs(title = "H3K4me1",
       x="", y="Correlation coefficient") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face="bold"), 
        legend.title = element_text(size=10, face="italic"),
        #axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=10),
        strip.text = element_text(size=10, face="bold.italic", vjust = 2.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)


# Calculates mean, sd, se and IC
my_sum <- xx %>%
  group_by(cond) %>%
  summarise( 
    n=n(),
    mean=mean(C_short_H3K4me1_AdipoERCre),
    sd=sd(C_short_H3K4me1_AdipoERCre)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

ggplot(my_sum) + 
  geom_bar(aes(cond, mean, fill = cond), 
           position = "dodge", stat = "summary", fun="mean") +
  geom_errorbar( aes(x=cond, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.9) +
  #scale_fill_manual(values = my_col) +
  labs(title = "H3K4me1",
       x="", y="Correlation coefficient") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face="bold"), 
        legend.title = element_text(size=10, face="italic"),
        #axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=10),
        strip.text = element_text(size=10, face="bold.italic", vjust = 2.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)


write.csv(xx, "H3K4me1_refC_corrCoeff.csv", row.names = FALSE)
