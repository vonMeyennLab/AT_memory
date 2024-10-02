
# This R script has been used to generate results used in 
# Figures 3-6 and associated extended data

library(chromVAR)
library(GenomicRanges)
library(tidyverse)

# Load sample metadata
# Metadata table has the following columns:
# ID | omics | replicate | diet | experiment | condition | timePoint | bamFilePath | peakFilePath
metadata <- read.table("../metadata.txt", header = TRUE, sep = "\t")

# Load blacklist region filtered peaks
filtPeaks   <- readRDS("filtPeaks.rds")
peakSummary <- filtPeaks$peakSummary

# List peak files to be used
condPeaks <- unique(metadata$ID)

# Merge all called peaks for all CnT and ATAC samples
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

saveRDS(countMat, "filtPeakCounts.rds")

################################################################################

## Peak counts in promoters

# Annotate peaks
countMat <- data.frame(countMat)
countMat$loci <- rownames(countMat)

peakRawCounts <- separate(countMat, "loci", c("chr", "pos"), sep=":")
peakRawCounts <- separate(peakRawCounts, "pos", c("start", "end"), sep="-")

peakRawCounts <- makeGRangesFromDataFrame(peakRawCounts,
                                              keep.extra.columns=TRUE,
                                              ignore.strand=FALSE,
                                              seqinfo=NULL,
                                              seqnames.field=c("chromosome", "chrom",
                                                               "chr", "chromosome_name"),
                                              start.field="start",
                                              end.field=c("end", "stop"),
                                              strand.field="strand",
                                              starts.in.df.are.0based=FALSE)

peakRawCounts <- renameSeqlevels(peakRawCounts, mapSeqlevels(seqlevels(peakRawCounts), "UCSC"))

peakAnno <- annotatePeak(peakRawCounts, tssRegion=c(-2000, 2000), 
                             TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")

# Convert peak annotations into data frame
peakAnno_df <- data.frame(peakAnno)

# Summarize peak read counts in promoter
peakPr <- subset(peakAnno_df, abs(distanceToTSS) <= 2000)
peakPr <- peakPr[-grep("predicted gene|microRNA|RIKEN|Riken|pseudogene", peakPr$GENENAME),]

peakPr <- peakPr[,c("geneId", samples)]

peakPr_rawCounts <- peakPr %>% 
  group_by(geneId) %>% 
  summarise(across(everything(), sum))

peakPr_rawCounts <- column_to_rownames(peakPr_rawCounts, "geneId")

saveRDS(countMat, "filtPeakPromCounts.rds")

################################################################################

## Peak counts in enhancers

# Load enhancer regions based on ChromHMM analysis
enhancers <- read.table("../ChromHMM/chromHMM_short_enhancers.bed", header = TRUE, sep = "\t") 

# Make GRanges for enhancer regions
EN_GR <- makeGRangesFromDataFrame(enhancers,
                                  keep.extra.columns=TRUE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field="seqnames",
                                  start.field="start",
                                  end.field="end")

EN_GR <- renameSeqlevels(EN_GR, mapSeqlevels(seqlevels(EN_GR), "UCSC"))

# Find peak regions overlapping with enhancers
table(!is.na(findOverlaps(peakRawCounts, EN_GR, select="arbitrary")))

EN_overlapsGR <- GenomicRanges::findOverlaps(peakRawCounts, EN_GR)
EN_ol <- data.frame(peakRawCounts[unique(EN_overlapsGR@from)])  
EN_ol$loci <- paste0(EN_ol$seqnames, ":", EN_ol$start, "-", EN_ol$end)

EN_ol_GR <- GRanges(seqnames = EN_ol$seqnames, IRanges(start = EN_ol$start, end = EN_ol$end), strand = "*", keep.extra.columns=TRUE)

# Select overlapping peaks outside promoter region
peakEn <- subset(peakAnno_df, abs(distanceToTSS) > 2000)
peakEn <- peakEn[-grep("predicted gene|microRNA|RIKEN|Riken|pseudogene", peakEn$GENENAME),]

peakEn$loci <- paste0(peakEn$seqnames, ":", peakEn$start, "-", peakEn$end)

EN_ol_peakCounts <- subset(EN_ol, loci %in% peakEn$loci)

peakEn <- EN_ol_peakCounts[,c("loci", samples)]
rownames(peakEn) <- NULL
peakEn <- column_to_rownames(peakEn, "loci")

saveRDS(peakEn, "../filtPeakShortEnCounts.rds")
