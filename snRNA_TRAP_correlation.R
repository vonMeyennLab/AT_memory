setwd("/home/ghoshad/public/_Projects/AT_HFD_Memory/snRNAseq_TRAPSeq")

library(Seurat)
library(ggplot2)
library(ggplot2)
library(ggpubr)
library(ezRun)
library(biomaRt)
library(corrplot)

# Load raw RNAseq data
rawCounts_RNAseq <- read.table("~/public/_Sequencing/SEQ00086/aligned/counts/gene_counts_merged.txt", header = TRUE, sep = "\t")
rawCounts_RNAseq <- column_to_rownames(rawCounts_RNAseq, "Geneid")

# Rename samples
colnames(rawCounts_RNAseq) <- gsub("mmepiAT_", "", colnames(rawCounts_RNAseq))
colnames(rawCounts_RNAseq) <- gsub("TRAP.", "", colnames(rawCounts_RNAseq))
colnames(rawCounts_RNAseq) <- gsub("_GRCm38_hisat2", "", colnames(rawCounts_RNAseq))

# Remove genes with no count across all samples
rawCounts <- rawCounts_RNAseq[rowSums(rawCounts_RNAseq) > 0, ]

# Create condition columns for mean raw counts over replicates
rawCounts$HFD <- rowMeans(rawCounts[,grep("HFD", colnames(rawCounts_RNAseq))])
rawCounts$CHD <- rowMeans(rawCounts[,grep("CHD", colnames(rawCounts_RNAseq))])
rawCounts$HC <- rowMeans(rawCounts[,grep("HC", colnames(rawCounts_RNAseq))])
rawCounts$CC <- rowMeans(rawCounts[,grep("CC", colnames(rawCounts_RNAseq))])

TRAPseq_normCounts <- data.frame(edgeR::cpm(rawCounts, log = TRUE))
colnames(TRAPseq_normCounts) <- paste0("TRAPseq_", colnames(TRAPseq_normCounts))

TRAPseq_normCounts <- rownames_to_column(TRAPseq_normCounts, "geneID")

# Retrieve gene names using Ensembl gene id
ensembl <- useMart('ensembl', dataset="mmusculus_gene_ensembl")

gene_mapping <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
                      filters = "ensembl_gene_id",
                      values = TRAPseq_normCounts$geneID, 
                      mart= ensembl)
colnames(gene_mapping) <- c("geneID", "geneName")

TRAPseq_normCounts <- merge(gene_mapping, TRAPseq_normCounts, by="geneID")

###########################################################################################################################################

# Load annotated rds object
snIntegrated <- readRDS("/home/ghoshad/public/_Projects/AT_HFD_Memory/snSeurat_multiSample_old/snSeurat_multiSample/short/snMulti_short_annotated.rds")
DefaultAssay(snIntegrated) <- "RNA"

# Select same conditions as TRAPseq
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

###########################################################################################################################################

# Combine the single cell and bulk RNAseq data 
comb <- merge(TRAPseq_normCounts, snRNA_avgNormCounts, by.x="geneName", by.y="row.names")

# Make the scatterplot to show correlations along with their Pearson correlation coefficient
gh1 <- ggscatter(comb, x = "snRNA_HFD", y = "TRAPseq_HFD", 
          cor.coef = TRUE, cor.method = "spearman", color = "#2E8B57") +
  labs(title = "HFD",
       x="Avg. log normalized expression (snRNAseq)", y="Avg. log normalized expression (TRAPseq)") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face="bold"), 
        axis.text=element_text(size=10, face="bold"),
        strip.text.x=element_text(size=10, face="bold"),
        axis.title=element_text(size=11,face="bold"),
        legend.title = element_text(size=10, face="italic")) 

ggplot(comb,aes(y=snRNA_HFD,x=TRAPseq_HFD)) + geom_point(alpha = 0.3, colour="#FF9300") + 
  ggtitle('H')+xlab("Avg. norm. expr. (TRAPseq)") + 
  ylab("Avg. norm. expr. (snRNAseq)")+ theme_bw()

gh2 <- ggscatter(comb, x = "snRNA_HC", y = "TRAPseq_HC", 
          cor.coef = TRUE, cor.method = "spearman", color = "#2E8B57") +
  labs(title = "HC",
       x="Avg. log normalized expression (snRNAseq)", y="Avg. log normalized expression (TRAPseq)") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face="bold"), 
        axis.text=element_text(size=10, face="bold"),
        strip.text.x=element_text(size=10, face="bold"),
        axis.title=element_text(size=11,face="bold"),
        legend.title = element_text(size=10, face="italic")) 

ggplot(comb,aes(y=snRNA_HC,x=TRAPseq_HC)) + geom_point(alpha = 0.3, colour="#069153") + 
  ggtitle('HC')+xlab("Avg. norm. expr. (TRAPseq)") + 
  ylab("Avg. norm. expr. (snRNAseq)")+ theme_bw()

gh1 + gh2

gc1 <- ggscatter(comb, x = "snRNA_CHD", y = "TRAPseq_CHD", 
          cor.coef = TRUE, cor.method = "spearman", color = "#2E8B57") +
  labs(title = "CHD",
       x="Avg. log normalized expression (snRNAseq)", y="Avg. log normalized expression (TRAPseq)") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face="bold"), 
        axis.text=element_text(size=10, face="bold"),
        strip.text.x=element_text(size=10, face="bold"),
        axis.title=element_text(size=11,face="bold"),
        legend.title = element_text(size=10, face="italic")) 

ggplot(comb,aes(y=snRNA_CHD,x=TRAPseq_CHD)) + geom_point(alpha = 0.3, colour="#7A81FF") + 
  ggtitle('C')+xlab("Avg. norm. expr. (TRAPseq)") + 
  ylab("Avg. norm. expr. (snRNAseq)")+ theme_bw()

gc2 <- ggscatter(comb, x = "snRNA_CC", y = "TRAPseq_CC", 
          cor.coef = TRUE, cor.method = "spearman", color = "#2E8B57") +
  labs(title = "CC",
       x="Avg. log normalized expression (snRNAseq)", y="Avg. log normalized expression (TRAPseq)") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face="bold"), 
        axis.text=element_text(size=10, face="bold"),
        strip.text.x=element_text(size=10, face="bold"),
        axis.title=element_text(size=11,face="bold"),
        legend.title = element_text(size=10, face="italic")) 

ggplot(comb,aes(y=snRNA_CC,x=TRAPseq_CC)) + geom_point(alpha = 0.3, colour="#8E67FF") + 
  ggtitle('CC_s')+xlab("Avg. norm. expr. (TRAPseq)") + 
  ylab("Avg. norm. expr. (snRNAseq)")+ theme_bw()


gc1 + gc2

# Use of the mtcars data proposed by R
data <- cor(comb[,15:22], method="spearman")

# Build a Pannel of 100 colors with Rcolor Brewer
my_colors <- brewer.pal(9, "Spectral")
my_colors <- colorRampPalette(my_colors)(100)

# Order the correlation matrix
ord <- order(data[1, ])
data_ord <- data[ord, ord]
plotcorr(data_ord , col=my_colors[data_ord*50+50] , mar=c(1,1,1,1))


###########################################################################################################################################

# CHD
pheatmap(
  mat               = cor(comb[,c(3,7,11,19)], use="complete.obs", method = "spearman"),
  cutree_rows       = 2, 
  cutree_cols       = 2,
  silent            = F,
  color             = brewer.pal(n = 9, name = "YlOrRd"),
  fontsize_row      = 8, 
  fontsize_col      = 8,
  display_numbers   = FALSE, 
  cellwidth         = 40,
  cellheight        = 40,
  border_color      = "black")

# CC
pheatmap(
  mat               = cor(comb[,c(5,9,13,21)], use="complete.obs", method = "spearman"),
  cutree_rows       = 2, 
  cutree_cols       = 2,
  silent            = F,
  color             = brewer.pal(n = 9, name = "YlOrRd"),
  fontsize_row      = 8, 
  fontsize_col      = 8,
  display_numbers   = FALSE, 
  cellwidth         = 40,
  cellheight        = 40,
  border_color      = "black")

# HFD
pheatmap(
  mat               = cor(comb[,c(4,8,12,20)], use="complete.obs", method = "spearman"),
  cutree_rows       = 2, 
  cutree_cols       = 2,
  silent            = F,
  color             = brewer.pal(n = 9, name = "YlOrRd"),
  fontsize_row      = 8, 
  fontsize_col      = 8,
  display_numbers   = FALSE, 
  cellwidth         = 40,
  cellheight        = 40,
  border_color      = "black")

# CC
pheatmap(
  mat               = cor(comb[,c(6,10,14,22)], use="complete.obs", method = "spearman"),
  cutree_rows       = 2, 
  cutree_cols       = 2,
  silent            = F,
  color             = brewer.pal(n = 9, name = "YlOrRd"),
  fontsize_row      = 8, 
  fontsize_col      = 8,
  display_numbers   = FALSE, 
  cellwidth         = 40,
  cellheight        = 40,
  border_color      = "black")

pheatmap(
  mat               = cor(comb[,c(3:22)], use="complete.obs", method = "spearman"),
  cutree_rows       = 6, 
  cutree_cols       = 6,
  silent            = F,
  color             = brewer.pal(n = 9, name = "YlOrRd"),
  fontsize_row      = 8, 
  fontsize_col      = 8,
  display_numbers   = FALSE,
  border_color      = "black")
