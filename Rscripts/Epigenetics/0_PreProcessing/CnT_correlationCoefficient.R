
# This R script has been used to generate the following figures:
### Figure 4 and associated Extended Data

library(dplyr)
library(edgeR)
library(ggfortify)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(tibble)
library(tidyverse)

# Load sample metadata
# Metadata table has the following columns:
# ID | omics | replicate | diet | experiment | condition | timePoint | bamFilePath | peakFilePath
metadata <- read.table("../metadata.txt", header = TRUE, sep = "\t")

# Define samples and conditions
condition <- metadata$condition
samples <- metadata$ID

################################################################################

# Following code chunk has been illustrated for H3K4me1
# Same code chunk has been used for H3K4me3 and H3K27ac

# Load union peak based raw counts (for each histone mark)
# Peak regions are represented in (chr):(start)-(end) format
peakCounts <- readRDS("../filtPeakCounts.rds")

# Calculate avg peak counts for control reference (young chow diet mice)
# Create new column to store counts for hypothetically ideal control 
peakCounts[, "C_H3K4me1"] <- rowMeans(peakCounts[,1:3])

# Select relevant metadata info
group <- data.frame(condition, row.names = samples)

# Add control info to metadata
group["C_H3K4me1",] <- "Control"

# Create DGEList object
y <- DGEList(counts=rawCounts, group=group$condition)

# Remove peak regions with low expression in > 50% samples
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Apply TMM normalization
y <- calcNormFactors(y)

# Extract normalized peak counts
normCounts <- edgeR::cpm(y, log=TRUE)

# Calculate Pearson correlation coefficients
corCoeff <- cor(normCounts, use="complete.obs")

# Combine correlation coefficients with metadata
corCoeff <- merge(corCoeff, group, by="row.names")

# Define conditions as factors
corCoeff$condition <- factor(corCoeff$cond, levels = c("Control", "C", "CC_s", "CC_l", "CCC",
                                      "H", "HC", "HH", "HHC"))

# Calculate mean and sd for correlation coefficients
corCoeffDist <- corCoeff %>%
  group_by(condition) %>%
  summarise( 
    n=n(),
    mean=mean(C_H3K4me1),
    sd=sd(C_H3K4me1)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

# Visualize histone correlation coefficients of each condition against the hypothetical control
# Error bars indicate mean and sd
ggplot(corCoeffDist) + 
  geom_bar(aes(condition, mean, fill = condition), 
           position = "dodge", stat = "summary", fun="mean") +
  geom_errorbar( aes(x=condition, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.9) +
  labs(title = "H3K4me1",
       x="", y="Correlation coefficient") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face="bold"), 
        legend.title = element_text(size=10, face="italic"),
         axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=10),
        strip.text = element_text(size=10, face="bold.italic", vjust = 2.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)
