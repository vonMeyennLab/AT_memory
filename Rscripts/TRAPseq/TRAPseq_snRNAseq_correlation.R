
# This R script has been used to generate the following figures:
### Figure S9C

library(Seurat)
library(ggplot2)
library(ggplot2)
library(ggpubr)
library(ezRun)
library(biomaRt)
library(corrplot)

# Load raw TRApseq counts
rawCounts_RNAseq <- read.table("../gene_counts_merged.txt", header = TRUE, sep = "\t")
rawCounts_RNAseq <- column_to_rownames(rawCounts_RNAseq, "Geneid")

# Rename samples
colnames(rawCounts_RNAseq) <- gsub("mmepiAT_", "", colnames(rawCounts_RNAseq))
colnames(rawCounts_RNAseq) <- gsub("TRAP.", "", colnames(rawCounts_RNAseq))
colnames(rawCounts_RNAseq) <- gsub("_GRCm38_hisat2", "", colnames(rawCounts_RNAseq))

# Remove genes with no count across all samples
rawCounts <- rawCounts_RNAseq[rowSums(rawCounts_RNAseq) > 0, ]

# Create condition columns for mean raw counts over replicates
rawCounts$H <- rowMeans(rawCounts[,grep("H", colnames(rawCounts_RNAseq))])
rawCounts$C <- rowMeans(rawCounts[,grep("C", colnames(rawCounts_RNAseq))])
rawCounts$HC <- rowMeans(rawCounts[,grep("HC", colnames(rawCounts_RNAseq))])
rawCounts$CC_s <- rowMeans(rawCounts[,grep("CC_s", colnames(rawCounts_RNAseq))])
rawCounts$HH <- rowMeans(rawCounts[,grep("HH", colnames(rawCounts_RNAseq))])
rawCounts$CC_l <- rowMeans(rawCounts[,grep("CC_l", colnames(rawCounts_RNAseq))])
rawCounts$HHC <- rowMeans(rawCounts[,grep("HHC", colnames(rawCounts_RNAseq))])
rawCounts$CCC <- rowMeans(rawCounts[,grep("CCC", colnames(rawCounts_RNAseq))])

TRAPseq_normCounts <- data.frame(edgeR::cpm(rawCounts, log = TRUE))
colnames(TRAPseq_normCounts) <- paste0("TRAPseq_", colnames(TRAPseq_normCounts))

TRAPseq_normCounts <- rownames_to_column(TRAPseq_normCounts, "geneID")

# Annotate counts with gene names
ensembl <- useMart('ensembl', dataset="mmusculus_gene_ensembl")

gene_mapping <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
                      filters = "ensembl_gene_id",
                      values = TRAPseq_normCounts$geneID, 
                      mart= ensembl)
colnames(gene_mapping) <- c("geneID", "geneName")

TRAPseq_normCounts <- merge(gene_mapping, TRAPseq_normCounts, by="geneID")

################################################################################

# Load annotated Seurat object
snIntegrated <- readRDS("../snMulti_short_annotated.rds")
DefaultAssay(snIntegrated) <- "RNA"

# Select sonly adipocytes
snIntegrated <- subset(snIntegrated, cellType == "Adipocytes")

# Calculate avg gene expression in adipocytes
snAd_meanGeneExp <- AverageExpression(snIntegrated, group.by = "condition")
snAd_meanGeneExp_df <- as.data.frame(snAd_meanGeneExp$RNA)
snAd_meanGeneExp_df <- snAd_meanGeneExp_df[,1:4]

# Remove genes with no count across all samples
avgRawCounts <- snAd_meanGeneExp_df[rowSums(snAd_meanGeneExp_df) > 0, ]

snRNA_avgNormCounts <- data.frame(edgeR::cpm(avgRawCounts, log = TRUE))

colnames(snRNA_avgNormCounts) <- gsub("t[1|2]_", "", colnames(snRNA_avgNormCounts))
colnames(snRNA_avgNormCounts) <- paste0("snRNA_", colnames(snRNA_avgNormCounts))

################################################################################

# Combine snRNAseq and TRAPseq data 
comb <- merge(TRAPseq_normCounts, snRNA_avgNormCounts, by.x="geneName", by.y="row.names")

# Make scatterplot to show Pearson correlations 
ggplot(comb,aes(y=snRNA_HFD,x=TRAPseq_HFD)) + geom_point(alpha = 0.3, colour="#FF9300") + 
  ggtitle('H')+xlab("Avg. norm. expr. (TRAPseq)") + 
  ylab("Avg. norm. expr. (snRNAseq)")+ theme_bw()

