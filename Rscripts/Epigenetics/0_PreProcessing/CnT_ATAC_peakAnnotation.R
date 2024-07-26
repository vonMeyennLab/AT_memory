
# This R script has been used to generate the following figures:
### Figures S9A, S9B

library(ChIPseeker)
library(ChIPpeakAnno)
library(GenomicRanges)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# Load sample metadata
# Metadata table has the following columns:
# ID | omics | replicate | diet | experiment | condition | timePoint | bamFilePath | peakFilePath
metadata <- read.table("../metadata.txt", header = TRUE, sep = "\t")

# Load blacklist regions
blReg <- read.table("../Mouse/mm10_blacklist_nochr.bed", sep = "\t", header = FALSE)
blRegGR <- GRanges(blReg$V1, IRanges(start = blReg$V2, end = blReg$V3), strand = "*")

# List unique factors
histone <- unique(metadata$omics)
diet <- unique(metadata$diet)
experiment <- unique(metadata$experiment)
replicate <- unique(metadata$replicate)
timePoint <- unique(metadata$timepoint)

################################################################################

## Peak QC

# Initialize data frames
peakN       <- c()    # peak count
peakSummary <- c()    # peak summary

# Calculate peak count and peak width for each sample
for(i in (1:nrow(metadata)))
{
  sampleInfo <- strsplit(metadata[i,"ID"], "_")[[1]]
  peakInfo   <- read.table(metadata[i, "peakFile"], header = FALSE, fill = TRUE)
  peakInfo$ID <- paste0(peakInfo$V1, ":", peakInfo$V2, "-", peakInfo$V3)

  # Calculate peak width
  peakInfo$ID           <- metadata[i,"ID"]
  peakInfo$condition    <- sampleInfo[3]
  peakInfo$exp          <- sampleInfo[4]
  peakInfo$histone      <- sampleInfo[5]
  peakInfo$replicate    <- sampleInfo[6]
  peakInfo$timePoint    <- sampleInfo[7]
  peakInfo$width        <- abs(peakInfo$V3 - peakInfo$V2) + 1

  # Calculate peak % overlapping with blacklist region
  histGR <- GRanges(peakInfo$V1, IRanges(start = peakInfo$V2, end = peakInfo$V3), strand = "*")
  peakInfo$numBLoverlap <- GenomicRanges::countOverlaps(histGR, blRegGR)

  countBLoverlap <- length(which(peakInfo$numBLoverlap>0))

  # Summarize peak statistics
  peakSummary <- rbind(peakSummary, peakInfo)
  peakN       <- data.frame(sampleID = metadata[i,"ID"],
                            condition = sampleInfo[3], exp = sampleInfo[4], histone = sampleInfo[5], replicate = sampleInfo[6], timePoint = sampleInfo[7],
                            peakN = nrow(peakInfo) - countBLoverlap, countBLoverlap = countBLoverlap) %>% rbind(peakN, .)
}

# Remove peaks overlapping with blacklist regions
peakSummary <- subset(peakSummary, numBLoverlap == 0)

saveRDS(list(peakSummary = peakSummary, peakN = peakN), "filtPeaks.rds")

################################################################################

## Peak annotation

# Initialize empty list to store peak annotations
sampleL <- list()
samples <- metadata$ID

# Store filtered peaks for each sample into a list
for(sample in samples)
{

  hist <- subset(peakSummary, ID == sample)
  histGR <- GRanges(hist$V1, IRanges(start = hist$V2, end = hist$V3), strand = "*")

  histGR <- renameSeqlevels(histGR, mapSeqlevels(seqlevels(histGR), "UCSC"))
  sampleL[[sample]] <- histGR
}

# Annotate peaks for all samples
peakAnno <- lapply(sampleL, function(x){annotatePeak(x, tssRegion=c(-2000, 2000),
                                                      TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")})

saveRDS(peakAnno, "filtPeakAnno.rds")

################################################################################

## Peak fold enrichment

# Load mm10 genomic annotation tracks
genome_annotation <- readRDS("../genomeAnnotation_mm10.rds")

# Select only genomic annotation tracks of interest
genome_annotation_subset <-
  genome_annotation %>% .[c("CGI Promoters", "Non-CGI Promoters", "3' UTRs", "5' UTRs",
                            "Exons", "Introns", "Intergenic regions", 
                            "PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3" )] %>% GRangesList()

# Function calculating genome size based on genic and intergenic regions
genome_size <- 
  sum(
    GenomicRanges::reduce(genome_annotation$`Gene bodies`, ignore.strand = TRUE) %>% 
      GenomicRanges::width() %>% sum,
    GenomicRanges::reduce(genome_annotation$`Intergenic regions`, ignore.strand = TRUE) %>%
      GenomicRanges::width() %>% sum
  )

calc_fold_enrichment <- 
  function(gr1, gr2, genome_size){
    
    A <- 
      GenomicRanges::intersect(
        GenomicRanges::reduce(gr2, ignore.strand = TRUE),
        GenomicRanges::reduce(gr1, ignore.strand = TRUE)
      ) %>% GenomicRanges::width() %>% sum
    
    B <- 
      GenomicRanges::reduce(gr1, ignore.strand = TRUE) %>%
      GenomicRanges::width() %>% sum
    
    C <- 
      GenomicRanges::reduce(gr2, ignore.strand = TRUE) %>% 
      GenomicRanges::width() %>% sum
    
    log2((A/B)/(C/genome_size))
  }

# Calculate peak fold enrichment of all samples
sampleFE <- list()

for (sample in samples)
{
  hist <- subset(peakSummary, ID == sample)
  hist$chr <- gsub("chr", "", hist$V1)
  histGR <- GRanges(hist$V1, IRanges(start = hist$V2, end = hist$V3), strand = "*")
  
  histFE <- lapply(genome_annotation_subset, function(x){calc_fold_enrichment(histGR, x, genome_size = genome_size)})
  sampleFE[[sample]] <- histFE
}

# Combine peak fold enrichment values of all samples in table format
sampleFE <- lapply(sampleFE, function(x){x <- do.call("rbind", x)})
comb_FE <- do.call("cbind", sampleFE)
colnames(comb_FE) <- names(sampleFE)

# Group genomic regions into genomic features and cCREs
genomicFeatures <- c("CGI Promoters","Non-CGI Promoters", "Exons", "Introns", "3' UTRs", "Intergenic regions")
cCREs <- c("CTCF-only", "DNase-H3K4me3", "PLS","pELS", "dELS")

# Heat map visualization of peak fold enrichment of all samples
breaksList <- seq(-2, 2, by = 0.1)

pheatmap::pheatmap(t(comb_FE[genomicFeatures,samples]),
                   silent            = F,
                   cutree_rows       = 2,
                   cutree_cols       = 2,
                   scale             = "row",
                   color             = colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(length(breaksList)),
                   fontsize_row      = 10, 
                   fontsize_col      = 10,
                   display_numbers   = FALSE,
                   fontsize_number   = 10, 
                   cluster_cols      = FALSE,
                   cluster_rows      = FALSE,
                   cellheight        = 40, 
                   cellwidth         = 40,
                   border_color      = 'black',
                   breaks            = breaksList)

pheatmap::pheatmap(t(comb_FE[cCREs,samples]),
                   silent            = F,
                   cutree_rows       = 2,
                   cutree_cols       = 2,
                   scale             = "row",
                   color             = colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(length(breaksList)),
                   fontsize_row      = 10, 
                   fontsize_col      = 10,
                   display_numbers   = FALSE,
                   fontsize_number   = 10, 
                   cluster_cols      = FALSE,
                   cluster_rows      = FALSE,
                   cellheight        = 40, 
                   cellwidth         = 40,
                   border_color      = 'black',
                   breaks            = breaksList)