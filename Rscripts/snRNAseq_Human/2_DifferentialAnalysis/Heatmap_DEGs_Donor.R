#### This code can be used to generate heatmaps; this specific code was used to generated heatmaps for number of DEGs or number of memory DEGs per cell types across cohorts in Extended Data Figure 1 

### Load libraries
library(tidyverse)
library(dplyr)
library(readr)
library(xlsx)
library(readxl)
### change names of columns
# Function to append dataframe name to specific column names
add_df_name_to_columns <- function(df, df_name, columns_to_modify) {
  for (col in columns_to_modify) {
    if (col %in% names(df)) {
      names(df)[names(df) == col] <- paste(df_name, col, sep = "_")
    }
  }
  return(df)
}

# Function to apply the above function to a list of dataframes
modify_list_of_dataframes <- function(df_list, columns_to_modify) {
  lapply(names(df_list), function(df_name) {
    add_df_name_to_columns(df_list[[df_name]], df_name, columns_to_modify)
  })
}
### read files; these files are DEG counts per cell type per donor
setwd(".../DEGs_t0")
files <- list.files(pattern = "count.csv", ignore.case = TRUE)
xx <- files %>%lapply(read.csv, stringsAsFactors=F) 
names(xx) <- files
xx <- modify_list_of_dataframes(xx, c("Down", "Up"))

#merge all data frames in list
xxx <- xx %>% reduce(inner_join, by='cellType')
xxx <- column_to_rownames(xxx, "cellType")
test <- as.data.frame(xxx)

##write into csv

write.csv(test, ".../DEG_Count_Summary.csv")

### heatmap
conditions <- factor(rep(c("Down", "Up"), 7)) #rep depends on how many donors

mapping <- data.frame(row.names = colnames(test),
                      conditions)
breakslist <- c(-2, -1, 0, 1, 2)
pheatmap(test[,1:14], 
         legend_breaks = breakslist,
         #color=getBlueRedScale(), 
         #clustering_method="ward.D2",
         scale="column", cluster_rows=FALSE,
         cluster_cols=F, 
         #cutree_rows = 2,
        # cutree_cols = 2,
         annotation_col = mapping,
         show_rownames = TRUE,
         show_colnames = TRUE,
         border_color = "black", 
         fontsize_row = 10, fontsize_col = 10,
         annotation_legend = TRUE,
         fontsize=8,
         cellheight = 10,
         cellwidth = 10)

#### where memory is strongest across celltypes and cohorts
### load file; this file contains a table with number of memory DEGs per cell type per cohort split by up and down
Memory_stats_for_heatmap <- read_excel("...//Memory_stats_for_heatmap_n.xlsx")
df<- column_to_rownames(Memory_stats_for_heatmap, "cellType")
Memory_Regulation <- factor(rep(c("Up", "Down"),4))
cohorts <- factor(c("MTSS", "MTSS", "LTSS", "LTSS", "LTSS", "LTSS", "NEFA", "NEFA"))
tissue <- factor(c(rep("omAT", 4), rep("scAT",4)))
mapping <- data.frame(row.names = colnames(df),
                      Memory_Regulation, tissue, cohorts)
breakslist <- c(-2, -1, 0, 1, 2)

pheatmap(df, 
         legend_breaks = breakslist,
         #color=getBlueRedScale(), 
         #clustering_method="ward.D2",
         scale="column", cluster_rows=FALSE,
         cluster_cols=F, 
         #cutree_rows = 2,
         # cutree_cols = 2,
         annotation_col = mapping,
         show_rownames = TRUE,
         show_colnames = TRUE,
         border_color = "black", 
         fontsize_row = 10, fontsize_col = 10,
         annotation_legend = TRUE,
         fontsize=8,
         cellheight = 10,
         cellwidth = 10)

