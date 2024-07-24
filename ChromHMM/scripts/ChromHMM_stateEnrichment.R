
# This R script has been used to generate the following figures:
### Figures S10A, S10B, S10C

library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(tidyr)

# Load emission data for the final model and convert into a matrix
emissions_data <- 
  read_tsv("../model/emissions_8.txt") %>%
  column_to_rownames("State (Emission order)") %>%
  as.matrix()

# Heat map of the emission states
pheatmap::pheatmap(emissions_data,
                   silent            = F,
                   cutree_rows       = 2,
                   cutree_cols       = 2,
                   scale             = "none",
                   color             = brewer.pal(n=9, name="Reds"),
                   fontsize_row      = 10, 
                   fontsize_col      = 10,
                   display_numbers   = FALSE,
                   fontsize_number   = 10, 
                   cluster_cols      = FALSE,
                   cluster_rows      = FALSE,
                   cellheight        = 15, 
                   cellwidth         = 30,
                   border_color      = 'black')


################################################################################

# Load genomic feature enrichment of ChromHMM states
genomFeat_enrichment <-
  read_tsv("../ATmem_genomicFeature_enrichment.txt") %>% 
  dplyr::rename("state" = "State (Emission order)") %>%
  .[1:9,] %>%
  mutate(state = c(1:8, "Base"))

# Remove file extensions
names(genomFeat_enrichment) <- genomFeat_enrichment %>%
  names %>%
  str_remove(".bed.gz") 

# Move states to row names
genomFeat_FE <- genomFeat_enrichment %>%
  filter(state != "Base") %>%
  column_to_rownames("state") %>% 
  as.matrix()

# Define genomic features of interest
genomicFeatures <- c("Genome %", "CGI_Promoters", "Non-CGI_Promoters", "Exons", "Introns", "3Prime_UTRs", "Intergenic")

# Heat map of genomic features enrichment
pheatmap::pheatmap(genomFeat_FE[,genomicFeatures],
                   silent            = F,
                   cutree_rows       = 2,
                   cutree_cols       = 2,
                   scale             = "column",
                   color             = colorRampPalette(rev(brewer.pal(n=11, name="RdYlBu")))(length(seq(-2, 2, by = 0.1))),
                   fontsize_row      = 10, 
                   fontsize_col      = 10,
                   display_numbers   = FALSE,
                   fontsize_number   = 10, 
                   cluster_cols      = FALSE,
                   cluster_rows      = FALSE,
                   cellheight        = 15, 
                   cellwidth         = 30,
                   border_color      = 'black')

################################################################################

# Load cCRE enrichment of ChromHMM states
cCRE_enrichment <-
  read_tsv("../ATmem_cCRE_enrichment.txt") %>% 
  dplyr::rename("state" = "State (Emission order)") %>%
  .[1:9,] %>%
  mutate(state = c(1:8, "Base"))

# Remove file extensions
names(cCRE_enrichment) <- cCRE_enrichment %>%
  names %>%
  str_remove(".bed.gz") 

# Move states to row names
cCRE_FE <- cCRE_enrichment %>%
  filter(state != "Base") %>%
  column_to_rownames("state") %>% 
  as.matrix()

# Heat map of ENCODE cCRE enrichment
pheatmap::pheatmap(cCRE_FE,
                   silent            = F,
                   cutree_rows       = 2,
                   cutree_cols       = 2,
                   scale             = "column",
                   color             = colorRampPalette(rev(brewer.pal(n=11, name="RdYlBu")))(length(seq(-2, 2, by = 0.1))),
                   fontsize_row      = 10, 
                   fontsize_col      = 10,
                   display_numbers   = FALSE,
                   fontsize_number   = 10, 
                   cluster_cols      = FALSE,
                   cluster_rows      = FALSE,
                   cellheight        = 15, 
                   cellwidth         = 30,
                   border_color      = 'black')