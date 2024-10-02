##### This code can be used to generate results and plots supporting Figure 4 and associated Extended Data

### load libraries

library(ChIPseeker)
library(ChIPpeakAnno)
library(GenomicRanges)
library(dplyr)
library(edgeR)
library(ggfortify)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(tibble)
library(tidyverse)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

### load counts and metadata

# Load sample metadata; this file contains information on hPTM, condition, QC and other experimental information, which can be used to select samples
metadata <- read.table(".../CnT_metadata.txt", header = TRUE, sep = "\t")

# Select only info for the specific histone; select specific samples for normalization and PCR clustering 
metadata <- subset(metadata, histone == "H3K4me1" & QC != "No" & timepoint == "AdipoERCre" & experiment == "short")

# Load raw peak counts; were quantified before
rawCounts <- readRDS(".../H3K4me1_filtPeakShortEnCounts.rds")

rawCounts <- rawCounts[,metadata$ID]

samples <- colnames(rawCounts)
histone <- as.factor(metadata$histone)
condition <- as.factor(metadata$condition)
experiment <- as.factor(metadata$experiment)
timePoint <- as.factor(metadata$timepoint)

# Add sample information
group <- data.frame(condition, row.names = samples)

# Create DGEList object
y <- DGEList(counts=rawCounts, group=group$condition)

# Remove genes with low expression in > 50% samples
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Apply TMM normalization
y <- calcNormFactors(y)

# Extract normalized counts
normCounts <- edgeR::cpm(y, log=TRUE)

### PCA plot

# Run PCA
pcDat  <- prcomp(t(normCounts), scale. = TRUE)

# PC1 vs PC2
autoplot(pcDat, data=group,
         colour="condition", size=6, frame = F) + theme_classic()

### if needed one can calculte the centres of mass and distances in the PCA (was done for revision round1 for all PCA plots)

pca_scores <- data.frame(pcDat$x, Condition = metadata$condition)

# Calculate centers of mass for each condition
centers_of_mass <- aggregate(. ~ Condition, data = pca_scores, mean)
## write as data.frame
write.csv(centers_of_mass, ".../centres_of_mass.csv")

# Calculate pairwise distances between centers of mass
distances <- dist(centers_of_mass[,-1]) # Exclude the 'Condition' column before calculating distances

# To view the distances in a more readable format
distance_matrix <- as.matrix(distances)
distance.df <- as.data.frame(distance_matrix)
rownames(distance.df) <- colnames(distance.df) <- centers_of_mass$Condition
write.csv(distance.df, ".../PC1PC2.csv")


### differential analysis between two conditions (H and C)
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

# Annotate genes based on log2FC / adjusted p-value thresholds; for enhancers we used FDR and not adj. p-value
res$regulation <- ""
res[which(res$logFC > 1 & res$FDR <= 0.05),"regulation"] <- "up"
res[which(res$logFC < -1 & res$FDR <= 0.05),"regulation"] <- "down"
res[which(abs(res$logFC) <= 1 | res$FDR > 0.05),"regulation"] <- "non-significant"

table(res$regulation)

##Select significant diff peaks based on p values
sig_peaks <- subset(res, FDR < 0.05 & abs(logFC) > 1)
sig_peaks <- rownames_to_column(sig_peaks, "loci")

sig_peaks <- separate(sig_peaks, "loci", c("chr", "pos"), sep=":")
sig_peaks <- separate(sig_peaks, "pos", c("start", "end"), sep="-")
sig_peaks$start <- as.integer(sig_peaks$start)
sig_peaks$end <- as.integer(sig_peaks$end)

# Save significant diff peaks data.frame as csv to check peaks in Seqmonk
write.csv(sig_peaks, ".../H3K4me1_sigDP_H_vs_C.csv")

# Update for peak fold enrichmnet
sig_peaks$chr <- gsub("chr", "", sig_peaks$chr)

# Split significant diff peaks into up and down
upPeaks <- subset(sig_peaks, regulation == "up")
downPeaks <- subset(sig_peaks, regulation == "down")

histDiffPeaks <- list("up" = upPeaks, 
                      "down" = downPeaks)

### Peak annotation

histDiffPeaks <- lapply(histDiffPeaks, function(x) {makeGRangesFromDataFrame(x,
                                                                             keep.extra.columns=TRUE,
                                                                             ignore.strand=FALSE,
                                                                             seqinfo=NULL,
                                                                             seqnames.field=c("chromosome", "chrom",
                                                                                              "chr", "chromosome_name"),
                                                                             start.field="start",
                                                                             end.field=c("end", "stop"),
                                                                             strand.field="strand",
                                                                             starts.in.df.are.0based=FALSE)})

histDiffPeaks <- lapply(histDiffPeaks, function(x) {renameSeqlevels(x, mapSeqlevels(seqlevels(x), "UCSC"))})

# Annotate the peaks
histDiffPeaks_anno <- lapply(histDiffPeaks, function(x) {annotatePeak(x, tssRegion=c(-2000, 2000), 
                                                                      TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")})

write.csv(data.frame(histDiffPeaks_anno$up), ".../H3K4me1_HvsC_upPeaks.csv")
write.csv(data.frame(histDiffPeaks_anno$down), ".../H3K4me1_HvsC_downPeaks.csv")




