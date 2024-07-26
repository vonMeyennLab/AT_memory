# :file_folder: Epigenetics

## Folders

- [0_Preprocessing]([put link](https://github.com/vonMeyennLab/AT_memory/tree/main/Rscripts/Epigenetics/0_PreProcessing)): Scripts for peak annotation, filtering, quantification of peaks, promoters and enhancers. And script to generate Fig. 5a.
- [1_MOFA]([link](https://github.com/vonMeyennLab/AT_memory/tree/main/Rscripts/Epigenetics/1_MOFA)): Script to prepare input for MOFA for every modality, run MOFA and generate factor plots.
- [2_DiffernetialAnalysis]([link](https://github.com/vonMeyennLab/AT_memory/tree/main/Rscripts/Epigenetics/2_DifferentialAnalysis)): Differential analysis for promoters and enhancers for a hPTM after quantification. Scripts for alluvial plots and GSEA.


## Workflow
After peakcalling peaks are annotated, filtered against a [blacklist](https://github.com/vonMeyennLab/AT_memory/blob/main/Rscripts/Epigenetics/0_PreProcessing/mm10_blacklist_nochr.bed) and then quantified using the bam files. Peaks overlapping promoters are identified and counts are subset for thos peaks and counts are aggregated at gene/promoter level. For enhancers, enhancer coordinates generated using ChromHMM, are used for subsetting peaks. MOFA is run by selecting the top variable features of each quantified modality. For differential analysis peak counts are used and analysis is performed per hPTM. An example of a [metadata file](https://github.com/vonMeyennLab/AT_memory/blob/main/Rscripts/Epigenetics/0_PreProcessing/CnT_metadata_example.txt) is included in this folder.
