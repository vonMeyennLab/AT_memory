setwd("~/public/_Projects/AT_HFD_Memory/Epigenetic_memory/QC/ATAC")

.libPaths("/scratch/lib/R/4.1.0")

library(chromVAR)
library(GenomicRanges)
library(tidyverse)

## Load sample data
metadata <- read.table("ATAC_metadata.txt", header = TRUE, sep = "\t")

# Select only AdipoERCre samples
metadata <- subset(metadata, timepoint == "AdipoERCre")

metadata <- subset(metadata, condition == "C")

# Update ID 
#metadata$ID <- paste0(metadata$ID, "_", metadata$timepoint)
#metadata$ID <- paste0(metadata$Name, "_", metadata$ID)

# Load blacklist region filtered peaks
filtPeaks   <- readRDS("ATAC_filtPeaks.rds")
peakSummary <- filtPeaks$peakSummary

# List peak files to be used
condPeaks <- unique(metadata$ID)

# Merge all called peaks for all samples
mPeak <- GRanges()

for(cp in condPeaks)
{
  peakRes <- subset(peakSummary, ID == cp)
  mPeak   <- GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
}

# Reduce genomic ranges (only for genome wide quantification)
masterPeak <- GenomicRanges::reduce(mPeak)

# List samples
samples <- metadata$ID

# Initialize count matrix
countMat <- matrix(NA, length(masterPeak), length(samples))
dim(countMat)

# Overlap with bam files to get count table
for(i in (1:nrow(metadata)))
{
  bam <- metadata[i, "bamFile"]
  fragment_counts <- chromVAR::getCounts(bam, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
  countMat[, i]   <- counts(fragment_counts)[,1]
}

colnames(countMat) <- samples

df <- data.frame(masterPeak)
rownames(countMat) <- paste0(df$seqnames, ":", df$start, "-", df$end)

# saveRDS(countMat, "ATAC_AdipoCre_filtPeakCounts.rds")
saveRDS(countMat, "ATAC_AdipoERCre_refC_filtPeakCounts.rds")
