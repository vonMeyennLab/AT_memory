#### This code can be used to generate stackec barplots used in Figure 1 and associated exteneded data. Note: size of texts etc were modified later. 

### load libraries
library(ggplot2)
library(ggpubr)


## load data; this file contains percentage of restored/non-restored DEGs per cell type for both up and down regulated DEGs and the N number of DEGs per cell type. It was compiled from the output of the transcriptional retention analysis and reformatted
df2 <- read_csv(".../TranscriptionalRetention.csv")


df_up <- subset(df2, (df2$State == "Upregulated at T0; restored at T1" | df2$State == "Upregulated at T0; NOT restored at T1") & df2$tissue == "scAT" & df2$cohort == "LTSS" 
                & (df2$cellType == "Adipo" | df2$cellType == "APCs" | df2$cellType == "EndoCs" | df2$cellType == "Macro" | df2$cellType == "MesoCs"
                | df2$cellType == "Peri" | df2$cellType == "SMCs"))

ggbarplot(df_up, "cellType", "Percent",
          fill = "State", color = "State", palette = c("#F47A60","#316879"),
          label = TRUE, lab.col = "black", lab.pos = "in", lab.size = 5) +
          labs(title = "Upregulated Memory DEGs scAT cohort L2") +
          theme(legend.text = element_text(size=15)) +
          theme(legend.title = element_text(size=15)) + geom_text(aes(label = paste0("n = ", Total_N)))

df_down <- subset(df2, (df2$State == "Downregulated at T0; restored at T1" | df2$State == "Downregulated at T0; NOT restored at T1") & 
                    df2$tissue == "scAT" & df2$cohort == "L2" 
                & (df2$cellType == "Adipo" | df2$cellType == "APCs" | df2$cellType == "EndoCs" | df2$cellType == "Macro" | df2$cellType == "MesoCs"
                   | df2$cellType == "Peri" | df2$cellType == "SMCs"))

ggbarplot(df_down, "cellType", "Percent",
          fill = "State", color = "State", palette = c("#F47A60","#316879"),
          label = TRUE, lab.col = "black", lab.pos = "in", lab.size = 5) +
  labs(title = "Downregulated Memory DEGs scAT cohort L2") +
  theme(legend.text = element_text(size=15)) +
  theme(legend.title = element_text(size=15)) + geom_text(aes(label = paste0("n = ", Total_N)))


