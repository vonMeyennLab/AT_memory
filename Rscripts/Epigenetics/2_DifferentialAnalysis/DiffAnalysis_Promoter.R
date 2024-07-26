##### This code can be used to generate data supporting the plots in Figures 4,5 and 6.

### Load libraries

library(biomaRt)
library(dplyr)
library(edgeR)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(tibble)
library(tidyverse)

#### Count QC and PCA

# Load sample metadata
# this metadata file contains ID, hPTM, replicate, condition, experiment and QC information that can be used to subset the metadata and then select samples
metadata <- read.table(".../CnT_metadata.txt", header = TRUE, sep = "\t")

# Select only info for the specific histone and for the correct GMO (time point of labelling)
metadata <- subset(metadata, histone == "H3K4me3" & QC != "No" & timepoint == "AdipoERCre")

# Load filtered peak counts. For peak/promoter quantification see other script. 
rawCounts <- readRDS(".../filtPeakPromCounts.rds")

rawCounts <- rawCounts[,metadata$ID]
rawCounts <- rawCounts[rowSums(rawCounts) > 0, ]

# Extract sample info
samples <- colnames(rawCounts)
histone <- as.factor(metadata$histone)
condition <- as.factor(metadata$condition)
experiment <- as.factor(metadata$experiment)

# Add sample information
group <- data.frame(condition, experiment, row.names = samples)


### use counts and perform normalization

# Create DGEList object
y <- edgeR::DGEList(counts=rawCounts, group=group$condition)

# Remove genes/promoters with low expression in > 50% samples
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Apply TMM normalization
y <- calcNormFactors(y)

# Extract normalized counts
normCounts <- edgeR::cpm(y, log=TRUE)


### Run PCA
pcDat  <- prcomp(t(normCounts), scale. = TRUE)

# PC1 vs PC2
autoplot(pcDat, data=group,
         colour="condition", shape="experiment", size=4, frame = TRUE) + theme_classic()

#### Differential analysis between two conditions (say H an C)

## Load sample data
metadata_sub <- subset(metadata, condition %in% c("C", "H"))

# Load raw peak counts
rawCounts_sub <- rawCounts[,metadata_sub$ID]
samples_sub <- colnames(rawCounts_sub)

# Add sample information
group_sub <- data.frame(condition=as.factor(metadata_sub$condition), row.names = samples_sub)

# Create DGEList object
y <- DGEList(counts=rawCounts_sub, group=group_sub$condition)

# Remove genes with low expression in > 50% samples
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Apply TMM normalization
y <- calcNormFactors(y)

# Extract normalized counts
normCounts_sub <- edgeR::cpm(y, log=TRUE)

# Construct design matrix
design <- model.matrix(~ 0 + group_sub$condition)
colnames(design) <- levels(group_sub$condition)

# Estimate dispersion
y <- estimateDisp(y, design, robust=TRUE)                  

# Estimate quasi-likelihood (QL) dispersions
fitGlm <- glmQLFit(y, design, robust = TRUE)

# Define contrast of interest
contr <- makeContrasts(H - C, levels=design)

# Apply QL F-test for differential expression
qlf <- glmQLFTest(fitGlm, contrast=contr)

# Extract HFD vs CHD DE results
res <- topTags(qlf, n = Inf, sort.by = "PValue")$table                         

# Annotate genes based on log2FC / adjusted p-value thresholds
res$regulation <- ""
res[which(res$logFC >= 0.5 & res$PValue <= 0.01),"regulation"] <- "up"
res[which(res$logFC <= -0.5 & res$PValue <= 0.01),"regulation"] <- "down"
res[which(abs(res$logFC) < 0.5 | res$PValue > 0.01),"regulation"] <- "non-significant"

# Hits in numbers
table(res$regulation)


# Select significant hits based on p. adj. and logFC
sigGenes <- subset(res, PValue < 0.01 & abs(logFC) > 0.5)

# Add gene names (promoters were already annotated and collapsed at gene level during quantification)
ensembl <- useMart('ensembl', dataset="mmusculus_gene_ensembl")

geneAnnotation <- getBM(attributes=c("entrezgene_id", "external_gene_name"), 
                        filters="entrezgene_id", values=rownames(sigGenes), mart=ensembl)

sigGenesPr <- merge(geneAnnotation, sigGenes, by.x="entrezgene_id", by.y="row.names")

# Save significant DEGs based on H3K18la promoter levels
write.csv(sigGenesPr, ".../H3K4me3_sigDEGPr_H_vs_C.csv")






