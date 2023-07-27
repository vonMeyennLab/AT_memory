library(MOFA2)
library(data.table)
library(ggplot2)
library(tidyverse)
library(DESeq2)

# Load CnT metadata to check for outliers
metadata <- read.table("~/public/_Projects/AT_HFD_Memory/Epigenetic_memory/QC/CnT/CnT_metadata.txt", header = TRUE, sep = "\t")
metadata <- subset(metadata, timepoint=="AdipoERCre")
metadata <- subset(metadata, histone!="H3K14pr")
#metadata <- subset(metadata, QC!="No")

# load ATAC and CnT union peak based raw counts
atac <- readRDS("~/public/_Projects/AT_HFD_Memory/Epigenetic_memory/QC/ATAC/ATAC_AdipoERCre_filtPeakCounts.rds")
h3k4me1 <- readRDS("~/public/_Projects/AT_HFD_Memory/Epigenetic_memory/ERCre/CnT_diffPeaks_AdipoERCre/H3K4me1/H3K4me1_AdipoERCre_filtPeakCounts.rds")
h3k4me3 <- readRDS("~/public/_Projects/AT_HFD_Memory/Epigenetic_memory/ERCre/CnT_diffPeaks_AdipoERCre/H3K4me3/H3K4me3_AdipoERCre_filtPeakCounts.rds")
h3k27me3 <- readRDS("~/public/_Projects/AT_HFD_Memory/Epigenetic_memory/ERCre/CnT_diffPeaks_AdipoERCre/H3K27me3/H3K27me3_AdipoERCre_filtPeakCounts.rds")
h3k27ac <- readRDS("~/public/_Projects/AT_HFD_Memory/Epigenetic_memory/ERCre/CnT_diffPeaks_AdipoERCre/H3K27ac/H3K27ac_AdipoERCre_filtPeakCounts.rds")

# Sample renaming
colnames(h3k4me1)
colnames(h3k4me1) <- gsub("_AdipoERCre", "", colnames(h3k4me1))
colnames(h3k4me1) <- gsub("_H3K4me1", "", colnames(h3k4me1))

colnames(h3k4me3)
colnames(h3k4me3) <- gsub("_AdipoERCre", "", colnames(h3k4me3))
colnames(h3k4me3) <- gsub("_H3K4me3", "", colnames(h3k4me3))

colnames(h3k27me3)
colnames(h3k27me3) <- gsub("_AdipoERCre", "", colnames(h3k27me3))
colnames(h3k27me3) <- gsub("_H3K27me3", "", colnames(h3k27me3))

colnames(h3k27ac)
colnames(h3k27ac) <- gsub("_AdipoERCre", "", colnames(h3k27ac))
colnames(h3k27ac) <- gsub("_H3K27ac", "", colnames(h3k27ac))

colnames(atac)
colnames(atac) <- gsub("_AdipoERCre", "", colnames(atac))
colnames(atac) <- gsub("_ATAC", "", colnames(atac))
colnames(atac) <- gsub("S09[0-1]_[0-1][0-9]_", "", colnames(atac))
atac <- atac[,colnames(h3k4me1)]

####################################################################################################################################

# Exclude outliers
h3k27ac_sub <- h3k27ac[,-c(6,9,11)]

# Define sample and conditions
samples <- sort(colnames(h3k27ac_sub))
condition <- as.factor(gsub("_[1-3]", "", samples))

# Add sample information
metadata <- data.frame(condition, row.names = samples)

#trap <- trap[,rownames(metadata)]

# Create dds object
dds <- DESeqDataSetFromMatrix(countData = as.matrix(h3k27ac_sub), 
                              colData = metadata,
                              design = ~condition)

# Remove genes with < 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Apply variance stabilizing transformation
vsd <- vst(dds)

# Extract normalized counts
normCounts <- assay(vsd)

# Identify high variance features
varGenes <- rownames(normCounts)[head(order(rowSds(as.matrix(normCounts)),
                                            decreasing=TRUE),3000)]

varGenes_normCounts <- normCounts[varGenes,]

h3k27ac <- data.frame(varGenes_normCounts)

h3k27ac$CC_short_3 <- NA
h3k27ac$CC_long_3 <- NA
h3k27ac$CCC_long_2 <- NA

h3k27ac <- h3k27ac[,colnames(atac)]

####################################################################################################################################

# Load TRAPseq raw counts
trap <- read.table("~/public/_Projects/AT_HFD_Memory/TRAPseq/ERCre/ERCre_TRAP_rawCounts.txt", header = TRUE, sep = "\t")
#trap$CC_long_3 <- NA
#trap <- trap[,colnames(h3k27ac)]
trap <- as.matrix(trap)

# Add gene names
ensembl <- useMart('ensembl', dataset="mmusculus_gene_ensembl")

geneAnnotation <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
                        filters="ensembl_gene_id", values=rownames(trap), mart=ensembl)

trap <- merge(geneAnnotation, trap, by.x="ensembl_gene_id", by.y="row.names")
trap <- trap[-grep("^Gm|Rik|Olf", trap$external_gene_name),]

rownames(trap) <- NULL
trap <- column_to_rownames(trap, "ensembl_gene_id")
trap <- trap[,2:24]

# Define sample and conditions
samples <- sort(colnames(trap))
condition <- as.factor(gsub("_[1-3]", "", samples))

# Add sample information
metadata <- data.frame(condition, row.names = samples)

trap <- trap[,rownames(metadata)]

# Create dds object
dds <- DESeqDataSetFromMatrix(countData = as.matrix(trap), 
                              colData = metadata,
                              design = ~condition)

# Remove genes with < 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Apply variance stabilizing transformation
vsd <- vst(dds)

# Extract normalized counts
normCounts <- assay(vsd)

# Identify high variance features
varGenes <- rownames(normCounts)[head(order(rowSds(as.matrix(normCounts)),
                                            decreasing=TRUE),3000)]

varGenes_normCounts <- normCounts[varGenes,]

trap <- data.frame(varGenes_normCounts)

trap$CC_long_3 <- NA
trap <- trap[,colnames(h3k4me1)]

geneAnnotation <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
                        filters="ensembl_gene_id", values=rownames(trap), mart=ensembl)

trap <- merge(geneAnnotation, trap, by.x="ensembl_gene_id", by.y="row.names")

rownames(trap) <- NULL
trap <- column_to_rownames(trap, "external_gene_name")
trap <- trap[,2:25]

#####################################################################################################################################

h3k27ac$CC_short_3 <- NA

h3k27ac <- h3k27ac[,colnames(h3k4me1)]

group <- group[colnames(h3k4me1),]

#####################################################################################################################################

rawCnTCounts <- atac

# List unique factors
CnTsamples <- colnames(rawCnTCounts)

# Annotate peaks
rawCnTCounts <- data.frame(rawCnTCounts)
rawCnTCounts$loci <- rownames(rawCnTCounts)

CnTpeakRawCounts <- separate(rawCnTCounts, "loci", c("chr", "pos"), sep=":")
CnTpeakRawCounts <- separate(CnTpeakRawCounts, "pos", c("start", "end"), sep="-")

CnTpeakRawCounts <- makeGRangesFromDataFrame(CnTpeakRawCounts,
                                             keep.extra.columns=TRUE,
                                             ignore.strand=FALSE,
                                             seqinfo=NULL,
                                             seqnames.field=c("chromosome", "chrom",
                                                              "chr", "chromosome_name"),
                                             start.field="start",
                                             end.field=c("end", "stop"),
                                             strand.field="strand",
                                             starts.in.df.are.0based=FALSE)

CnTpeakRawCounts <- renameSeqlevels(CnTpeakRawCounts, mapSeqlevels(seqlevels(CnTpeakRawCounts), "UCSC"))

CnTpeakAnno <- annotatePeak(CnTpeakRawCounts, tssRegion=c(-2000, 2000), 
                            TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")

# Convert the peak annotations into data frame
CnTpeakAnno_df <- data.frame(CnTpeakAnno)
CnTpeakAnno_df$ID <- paste0(CnTpeakAnno_df$seqnames, ":", CnTpeakAnno_df$start, "-", CnTpeakAnno_df$end, "_", CnTpeakAnno_df$SYMBOL)
CnTpeakAnno_df <- column_to_rownames(CnTpeakAnno_df, "ID")
CnTpeakAnno_df <- CnTpeakAnno_df[-grep("predicted gene|microRNA|RIKEN|Riken|pseudogene", CnTpeakAnno_df$GENENAME),]

atac <- CnTpeakAnno_df[,c(6:29)]

#####################################################################################################################################

# Indicate the outliers
h3k4me1[,"CCC_long_2"] <- NA 
h3k4me3[,"CCC_long_2"] <- NA 
h3k27me3[,"CCC_long_2"] <- NA 
h3k27ac[,"CCC_long_2"] <- NA 

h3k4me1[,"CC_long_3"] <- NA 
h3k4me3[,"CC_long_3"] <- NA 
h3k27me3[,"CC_long_3"] <- NA 
h3k27ac[,"CC_long_3"] <- NA 

# Define group info
samples <- colnames(h3k4me1)
condition <- c(rep("C", 3), rep("CC", 5), rep("CCC", 2),
               rep("H", 3), rep("HC", 3), rep("HH", 3), rep("HHC", 3))
experiment <- c(rep("short", 3), rep("long", 2),
                rep("short", 3), rep("long", 2),
                rep("short", 3), rep("short", 3),
                rep("long", 3), rep("long", 3))

samples <- samples[c(1:5,7:10,12:24)]
group <- data.frame(condition, experiment, row.names=samples)
group$cond <- paste0(group$condition, "_", group$experiment)

# Create list of diff omics assays
matList <- list("TRAP"=trap,
                "ATAC" = atac,
                "H3K4me3"=h3k4me3,
                "H3K27me3"=h3k27me3,
                "H3K27ac"=h3k27ac,
                "H3K4me1"=h3k4me1)

matList <- lapply(matList, function(x) { x <- as.matrix(x)})

saveRDS(matList[2:6], "~/public/_Projects/AT_HFD_Memory/Epigenetic_memory/MOFA/MOFA_wo_TRAP/mofa2_input_varFeatures.rds")
matList <- readRDS("~/public/_Projects/AT_HFD_Memory/Epigenetic_memory/MOFA/MOFA_wo_TRAP/mofa2_input_varFeatures.rds")

# Update ATAC and TRAP matrices
matList$TRAP <- matList$TRAP[,samples]
matList$H3K27ac <- matList$H3K27ac[,samples]
matList$H3K27me3 <- matList$H3K27me3[,samples]
matList$H3K4me3 <- matList$H3K4me3[,samples]
matList$H3K4me1 <- matList$H3K4me1[,samples]

trap[,"CCC_long_2"] <- NA 

#####################################################################################################################################

# Create MOFA object
MOFAobject <- create_mofa(matList)

MOFAobject

plot_data_overview(MOFAobject)

# Define MOFA options
data_opts <- get_default_data_options(MOFAobject)
data_opts
data_opts$scale_views <- TRUE

# Build model
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 5

model_opts

# Train model
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
train_opts$weight_views <- TRUE
#train_opts$maxiter <- 1000

MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts)

MOFAobject <- run_mofa(MOFAobject)

slotNames(MOFAobject)
names(MOFAobject@data)
dim(MOFAobject@data$W)
names(MOFAobject@expectations)

####################################################################################################################################

# Add sample metadata to the model
group <- rownames_to_column(group, "sample")

samples_metadata(MOFAobject) <- group

# Check if the factors are largely not correlated
# Indicator of good model fit
plot_factor_cor(MOFAobject)

# % of variance explained by each factor across each data modality 
p <- plot_variance_explained(MOFAobject, plot_total = FALSE) 
p + scale_color_brewer(palette="Reds")

# Total variance explained per view
plot_variance_explained(MOFAobject, plot_total = T)[[2]]

# Cumulative variance explained by factors
r2 <- MOFAobject@cache$variance_explained$r2_per_factor[[1]]

r2.dt <- r2 %>%
  as.data.table %>% .[,factor:=as.factor(1:MOFAobject@dimensions$K)] %>%
  melt(id.vars=c("factor"), variable.name="view", value.name = "r2") %>%
  .[,cum_r2:=cumsum(r2), by="view"]

ggline(r2.dt, x="factor", y="cum_r2", color="view") +
  labs(x="Factor number", y="Cumulative variance explained (%)") +
  theme(
    legend.title = element_blank(), 
    legend.position = "top",
    axis.text = element_text(size=rel(0.8))
  )

####################################################################################################################################

# Characterizing factor

# Plot Factor1 vs Factor2
p <- plot_factors(MOFAobject, 
                  factors = c(1,4), 
                  color_by = "condition",
                  shape_by = "experiment",
                  dot_size = 6,
                  show_missing = T
)

p <- p + 
  geom_hline(yintercept=-1, linetype="dashed") +
  geom_vline(xintercept=(-0.5), linetype="dashed")

print(p)

# Factor 1
plot_weights(MOFAobject,
             view = "H3K4me1",
             factor = 1,
             nfeatures = 20,     # Top number of features to highlight
             scale = T ) 

plot_top_weights(MOFAobject,
                 view = "H3K4me1",
                 factor = 1,
                 nfeatures = 40,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
                 view = "H3K4me3",
                 factor = 1,
                 nfeatures = 40,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
                 view = "H3K27ac",
                 factor = 1,
                 nfeatures = 40,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
                 view = "ATAC",
                 factor = 1,
                 nfeatures = 40,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
                 view = "TRAP",
                 factor = 4,
                 nfeatures = 40,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
                 view = "TRAP",
                 factor = 5,
                 nfeatures = 40,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)

plot_factor(MOFAobject, factors = 1, color_by = "condition", dodge = TRUE, add_violin = TRUE)
plot_factor(MOFAobject, factors = 2, color_by = "condition", dodge = TRUE, add_violin = TRUE)
plot_factor(MOFAobject, factors = 3, color_by = "condition", dodge = TRUE, add_violin = TRUE)
plot_factor(MOFAobject, factors = 4, color_by = "condition", dodge = TRUE, add_violin = TRUE)
plot_factor(MOFAobject, factors = 5, color_by = "condition", dodge = TRUE, add_violin = TRUE)

saveRDS(MOFAobject, "~/public/_Projects/AT_HFD_Memory/Epigenetic_memory/MOFA/MOFA_wo_TRAP/mofa2_out_varFeatures.rds")

MOFAobject <- readRDS("~/public/_Projects/AT_HFD_Memory/Epigenetic_memory/MOFA/mofa2_out_varFeatures.rds")

####################################################################################################################################

# Extract factor values
factors <- MOFAobject@expectations$Z$group1

# Extract feature weights for each factor
trap_wts <- MOFAobject@expectations$W$TRAP

# Extract factor based variance
p <- plot_variance_explained(MOFAobject, plot_total = FALSE)

factorVar_longDF <- data.frame(p$data)

factorVar_wideDF <- spread(factorVar_longDF[,1:3], key = view, value = value)
factorVar_wideDF <- column_to_rownames(factorVar_wideDF, "factor")

pheatmap(as.matrix(factorVar_wideDF), 
         clustering_method = "ward.D2",
         color             = brewer.pal(n=9, name="Blues"),
         scale = "none", 
         cluster_rows = F,
         cluster_cols = F, 
         show_rownames = TRUE,
         fontsize_row = 9, fontsize_col = 9,
         annotation_legend = TRUE,
         fontsize=8,
         cellheight        = 40, 
         cellwidth         = 60,
         border_color      = 'black')
