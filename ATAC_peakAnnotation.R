setwd("~/public/_Projects/AT_HFD_Memory/Epigenetic_memory")

.libPaths("/scratch/lib/R/4.1.0")

library(ChIPseeker)
library(ChIPpeakAnno)
library(GenomicRanges)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

## Load sample data
metadata <- read.table("ATAC_metadata.txt", header = TRUE, sep = "\t")

# Update ID 
#metadata$ID <- paste0(metadata$ID, "_", metadata$timepoint)
#metadata$ID <- paste0(metadata$Name, "_", metadata$ID)

# Load blacklist regions
blReg <- read.table("~/public/_Projects/MS_H3K18la/Mouse/mm10_blacklist_nochr.bed", sep = "\t", header = FALSE)
blRegGR <- GRanges(blReg$V1, IRanges(start = blReg$V2, end = blReg$V3), strand = "*")

# List unique factors
histone <- unique(metadata$histone)
condition <- unique(metadata$condition)
experiment <- unique(metadata$experiment)
replicate <- unique(metadata$replicate)
timePoint <- unique(metadata$timepoint)

# Initialize data frames
peakN       <- c()    # peak number
peakSummary <- c()    # peak summary

# Calculate peak number and peak width for each sample
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

# Remove the peaks overlapping the blacklist region
peakSummary <- subset(peakSummary, numBLoverlap == 0)

saveRDS(list(peakSummary = peakSummary, peakN = peakN), "ATAC_filtPeaks.rds")

################################################################################################################################

## Peak annotation

# Initialize empty list to store the peak annotations
sampleL <- list()
samples <- metadata$ID

# Calculate overlapping peaks for each mark
for(sample in samples)
{

  hist <- subset(peakSummary, ID == sample)
  histGR <- GRanges(hist$V1, IRanges(start = hist$V2, end = hist$V3), strand = "*")

  histGR <- renameSeqlevels(histGR, mapSeqlevels(seqlevels(histGR), "UCSC"))
  sampleL[[sample]] <- histGR
}

# Annotate the peaks
peakAnno <- lapply(sampleL, function(x){annotatePeak(x, tssRegion=c(-2000, 2000),
                                                      TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")})

saveRDS(peakAnno, "ATAC_filtPeakAnno.rds")

#promoter <- getPromoters(TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, upstream=3000, downstream=3000)
#tagMatrix <- lapply(sampleL, function(x){getTagMatrix(x, windows=promoter)})
 
#saveRDS(tagMatrix, "ATAC_filtPeakTagMatrix.rds")

################################################################################################################################

# Load mm10 genomic annotation
genome_annotation <- readRDS("~/public/AdhidebGhosh/H3K18la/genomeAnnotation_mm10.rds")

# Adding CGI and Non-CGI promoters to the annotation list
genome_annotation$`CGI Promoters`     <- genome_annotation$Promoters %>% .[.$CGI == TRUE]
genome_annotation$`Non-CGI Promoters` <- genome_annotation$Promoters %>% .[.$CGI == FALSE]

genome_annotation_subset <-
  genome_annotation %>% .[c("CGI Promoters", "Non-CGI Promoters", "3' UTRs", "5' UTRs",
                            "Exons", "Introns", "Intergenic regions", 
                            "PLS", "pELS", "dELS", "CTCF", "DNase_H3K4me3" )] %>% GRangesList()
                            
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

## Calculate peak fold enrichment
sampleFE <- list()

for (sample in samples)
{
  hist <- subset(peakSummary, ID == sample)
  hist$chr <- gsub("chr", "", hist$V1)
  histGR <- GRanges(hist$V1, IRanges(start = hist$V2, end = hist$V3), strand = "*")

  histFE <- lapply(genome_annotation_subset, function(x){calc_fold_enrichment(histGR, x, genome_size = genome_size)})
  sampleFE[[sample]] <- histFE
}

sampleFE <- lapply(sampleFE, function(x){x <- do.call("rbind", x)})

# Peak fold enrichment values of all samples in table format
comb_FE <- do.call("cbind", sampleFE)
colnames(comb_FE) <- names(sampleFE)

saveRDS(comb_FE, "ATAC_filtPeakFE.rds")
