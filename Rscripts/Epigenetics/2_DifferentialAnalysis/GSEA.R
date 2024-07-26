##### GSEA analysis using EnrichR; this tool was used for all GSEA analysis. 

### load libraries
suppressMessages({
  options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")
  options(enrichR.live = TRUE)
  options(modEnrichR.use = TRUE)
  options(enrichR.sites.base.address = "https://maayanlab.cloud/")
  options(enrichR.sites = c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr"))
  
  # Set the search to Human genes.
  enrichR::setEnrichrSite(site = "Enrichr")
  
  websiteLive <- TRUE
  dbs <- enrichR::listEnrichrDbs()
  # Get all the possible databases to query.
  dbs <- sort(dbs$libraryName)
})

# Choose the dataset to query against.Can use multiple or only one. For snRNAseq we used WikiPathways, for epigenetics WikiPathways, GOCC and Reactome 
dbs_use <- c("Reactome_2022",
             "GO_Cellular_Component_2021", 
             "WikiPathways_2019_Mouse" )

### load gene files, these could be genes associated to "memory" promoters or enhancers or "memory" DEGs from snRNAseq
## if you run this after the alluvial plot script you may use this codd=e

genes_up <- subset(Overlap_s, Overlap_s$score_HvC == "up" & Overlap_s$score_HCvCC == "up")$external_gene_name
genes_down <- subset(Overlap_s, Overlap_s$score_HvC == "down" & Overlap_s$score_HCvCC == "down")$external_gene_name

# Retrieve the enriched terms and plot them
enriched_terms <- enrichR::enrichr(genest_up, dbs_use)
write.xlsx(enriched_terms, ".../up.terms.xlsx")

SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms, plot.title = NULL, plot.subtitle = "Upregulated terms ", text_labels_size = 5, colors.use = c("red", "#FFFACD"), nchar_wrap = 20, nterms = 8 )

enriched_terms <- enrichR::enrichr(genes_down, dbs_use)
write.xlsx(enriched_terms, ".../down.terms.xlsx")
SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms, plot.title = NULL, plot.subtitle = "Downregulated terms ", text_labels_size = 5, colors.use = c("red", "#FFFACD"), nchar_wrap = 20, nterms = 8 )
