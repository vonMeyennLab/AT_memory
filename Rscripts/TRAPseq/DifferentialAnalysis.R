#### This code can be used to analyse TRAPseq data (PCA plot and logFC extraction) used in Figure 4, 5 and associated data.

### load libraries

library(biomaRt)
library(dplyr)
library(edgeR)
library(ezRun)
library(ggfortify)
library(ggplot2)
library(Matrix)
library(matrixStats)
library(pheatmap)
library(RColorBrewer)
library(Rsubread)
library(scales)
library(tibble)
library(tidyverse)
library(WGCNA)
library(ggrepel)

### Load raw counts (are derived from feature counts from the RNAseq nextflow pipeline)

# Load raw RNAseq data and metadata to select samples
metadata <- read.table("TRAP_metadata.txt", header = TRUE, sep = "\t")
rawCounts <- read.table(".../ERCre_TRAP_rawCounts.txt", header = TRUE, sep = "\t")
#rawCounts <- column_to_rownames(rawCounts, "Geneid")

# Subset metadata
metadata <- subset(metadata, QC == "yes")

# Load filtered peak counts
rawCounts <- rawCounts[,metadata$ID]
rawCounts <- rawCounts[rowSums(rawCounts) > 0, ]

# Extract sample info
samples <- colnames(rawCounts)
condition <- as.factor(metadata$condition)
experiment <- as.factor(metadata$experiment)

# Add sample information
group <- data.frame(condition, experiment, row.names = samples)

### Normalize Counts

# Create DGEList object
y <- edgeR::DGEList(counts=rawCounts, group=group$condition)
# Remove genes with low expression in > 50% samples
# Remove genes with low expression in > 50% samples
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Apply TMM normalization
y <- calcNormFactors(y)

# Extract normalized counts
normCounts <- edgeR::cpm(y, log=TRUE)

### PCA

# Run PCA
pcDat  <- prcomp(t(normCounts), scale. = TRUE)

# PC1 vs PC2
autoplot(pcDat, data=group,
         colour="condition", shape="experiment", size=6, frame = F) + theme_classic()

##### calculate distances (was done for first round of revisions)

pca_scores <- data.frame(pcDat$x, Condition = metadata$condition)

# Calculate centers of mass for each condition
centers_of_mass <- aggregate(. ~ Condition, data = pca_scores, mean)
## write as data.frame
write.csv(centers_of_mass, ".../TRAP_centres_of_mass.csv")

# Calculate pairwise distances between centers of mass
distances <- dist(centers_of_mass[,-1]) # Exclude the 'Condition' column before calculating distances

# To view the distances in a more readable format
distance_matrix <- as.matrix(distances)
distance.df <- as.data.frame(distance_matrix)
rownames(distance.df) <- colnames(distance.df) <- centers_of_mass$Condition
write.csv(distance.df, ".../TRAP_geometric_distance_PC1PC2.csv")


### DE analysis between two conditions


# Construct design matrix
design <- model.matrix(~ 0 + conditions)
colnames(design) <- levels(conditions)

# Estimate dispersion
y <- estimateDisp(y, design, robust=TRUE)                  

# Estimate quasi-likelihood (QL) dispersions
fitGlm <- glmQLFit(y, design, robust = TRUE)

# Define contrast of interest
contr <- makeContrasts(HC - CC_s, levels=design)

# Apply QL F-test for differential expression
qlf <- glmQLFTest(fitGlm, contrast=contr)

# Extract HFD vs CHD DE results
res <- topTags(qlf, n = Inf, sort.by = "PValue")$table               

# Gene mapping
ensembl <- useMart('ensembl', dataset="mmusculus_gene_ensembl")

geneAnnotation <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
                        filters="ensembl_gene_id", 
                        values=rownames(res), mart=ensembl)

# Combine DE results with gene names
comb <- merge(geneAnnotation, data.frame(res), 
              by.x="ensembl_gene_id", 
              by.y="row.names")

# Rename columns
colnames(comb)[1:2] <- c("geneID", "geneName")

# Annotate genes based on log2FC and FDR or P adj thresholds 
comb$regulation <- ""
comb[which(comb$logFC > 1 & comb$PValue < 0.01),"regulation"] <- "up"
comb[which(comb$logFC < -1 & comb$PValue < 0.01),"regulation"] <- "down"
comb[which(abs(comb$logFC) <= 1 | comb$PValue >= 0.01),"regulation"] <- "non-significant"

table(comb$regulation)

# Save significant DE genes 
write.table(comb, file = "DEGs.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)



