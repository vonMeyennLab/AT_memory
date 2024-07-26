
# This R script has been used to generate the following figures:
### Figures 4D, 4E

library(MOFA2)
library(data.table)
library(ggplot2)
library(tidyverse)
library(DESeq2)

# Load sample metadata for H3K4me3, H3K27me3, H3K27ac, H3K4me1, ATACseq and TRAPseq
# Metadata table has the following columns:
# ID | omics | replicate | diet | experiment | condition | timePoint | bamFilePath | peakFilePath
metadata <- read.table("../metadata.txt", header = TRUE, sep = "\t")

# Define samples and conditions
condition <- metadata$condition
samples <- metadata$ID

################################################################################

# Following code chunk has been used for each histone mark and ATACseq data
# to identify the most variable 3,000 peak regions

# Load union peak based raw counts (for each histone mark)
# Peak regions are represented in (chr):(start)-(end)_(annotated_closest_gene) format
peakCounts <- readRDS("../filtPeakCounts.rds")

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = as.matrix(peakCounts), 
                              colData = metadata,
                              design = ~condition)

# Remove peak regions with < 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Apply variance stabilizing transformation
vsd <- vst(dds)

# Extract normalized peak counts
normCounts <- assay(vsd)

# Identify most variable 3,000 peak regions
varPeaks <- rownames(normCounts)[head(order(rowSds(as.matrix(normCounts)),
                                            decreasing=TRUE),3000)]

# Select normalized counts of most variable peak regions
varPeaks_normCounts <- normCounts[varPeaks,]

# Transform normalized count matrix into data frame
h3k4me3 <- data.frame(varPeaks_normCounts) # h3k4me3 as an example

# Add NA to missing replicate(s)
h3k4me3$CCC_3 <- NA

# Sort data frame in sample order
h3k4me3 <- h3k4me3[,samples]

################################################################################

# Following code chunk has been used only for TRAPseq data
# to identify the most variable 3,000 genes

# Load raw gene counts
trap <- read.table("../TRAP_rawCounts.txt", header = TRUE, sep = "\t")

# Convert data frame into matrix
trap <- as.matrix(trap)

# Add gene annotations
ensembl <- useMart('ensembl', dataset="mmusculus_gene_ensembl")

geneAnnotation <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
                        filters="ensembl_gene_id", values=rownames(trap), mart=ensembl)

# Combine gene names with raw counts
trap <- merge(geneAnnotation, trap, by.x="ensembl_gene_id", by.y="row.names")

# Exclude irrelevant (noisy) genes from the most variable gene selection
trap <- trap[-grep("^Gm|Rik|Olf", trap$external_gene_name),]

# Update row names by gene names
rownames(trap) <- NULL
trap <- column_to_rownames(trap, "external_gene_name")

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = as.matrix(trap), 
                              colData = metadata,
                              design = ~condition)

# Remove genes with < 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Apply variance stabilizing transformation
vsd <- vst(dds)

# Extract normalized gene counts
normCounts <- assay(vsd)

# Identify most variable 3,000 genes
varGenes <- rownames(normCounts)[head(order(rowSds(as.matrix(normCounts)),
                                            decreasing=TRUE),3000)]

# Select normalized counts of most variable genes
varGenes_normCounts <- normCounts[varGenes,]

# Transform normalized count matrix into data frame
trap <- data.frame(varGenes_normCounts)

# Add NA to missing replicates
trap$CC_l_3 <- NA

# Sort data frame in sample order
trap <- trap[,samples]

################################################################################

# Following code chunk has been used to prepare the input
# and run Multi-Omics Factor Analysis (MOFA)

# Combine most variable features from different omics assays into a list
# Ensure for each assay, samples are in the same order
matList <- list("TRAP" = trap,
                "ATAC" = atac,
                "H3K4me3" = h3k4me3,
                "H3K27me3" = h3k27me3,
                "H3K27ac" = h3k27ac,
                "H3K4me1" = h3k4me1)

# Convert into list of matrices
matList <- lapply(matList, function(x) { x <- as.matrix(x)})

# Create MOFA object
MOFAobject <- create_mofa(matList)

# Check summary of different omics assays
plot_data_overview(MOFAobject)

# Define MOFA options
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE

# Build model with 5 latent factors
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 5

# Train model
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
train_opts$weight_views <- TRUE

MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts)

# Run MOFA
MOFAobject <- run_mofa(MOFAobject)

# Add sample metadata to the model
metadata <- rownames_to_column(metadata, "sample")
samples_metadata(MOFAobject) <- metadata

################################################################################

# Following code chunk has been used to explore and visualize MOFA results

# Check if the factors are largely not correlated
# Indicator of good model fit
plot_factor_cor(MOFAobject)

# Visualize % of variance explained by each factor across each data modality 
p <- plot_variance_explained(MOFAobject, plot_total = FALSE) 
p + scale_color_brewer(palette="Reds")

# Total variance explained per assay
plot_variance_explained(MOFAobject, plot_total = T)[[2]]


# Characterizing  latent factors
# Plot Factor1 vs Factor2
p <- plot_factors(MOFAobject, 
                  factors = c(1,2), 
                  color_by = "condition",
                  shape_by = "experiment",
                  dot_size = 6,
                  show_missing = T
)

p <- p + 
  geom_hline(yintercept=-1, linetype="dashed") +
  geom_vline(xintercept=(-0.5), linetype="dashed")

print(p)

# Visualize latent factor1 scores 
plot_factor(MOFAobject, factors = 1, color_by = "condition", shape_by = "experiment", 
            dodge = TRUE, add_violin = TRUE)
