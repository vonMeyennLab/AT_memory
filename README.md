# Adipose Tissue Retains an Epigenetic Memory of Obesity After Weight Loss

This respository contains code and files related to the manuscript [Adipose Tissue Retains an Epigenetic Memory of Obesity After Weight Loss](link).

## Abstract
Reducing body weight to improve metabolic health and other comorbidities is a primary goal in treating obesity. However, maintaining weight loss is a considerable challenge, especially as the body is believed to retain an obesogenic memory that defends against body weight changes. Yet, overcoming this hurdle to long-term effective treatment is difficult because the molecular mechanisms underpinning this phenomenon remain largely unknown. Here, by using single-nuclei RNA-sequencing, we show that both human and mouse adipose tissue retain cellular transcriptional changes after appreciable weight loss. Furthermore, we observed that the mouse adipocyte epigenome continues to bear obesity-induced alterations, negatively affecting adipocyte function. In mice, adipocytes carrying this obesogenic epigenetic memory respond differently to nutritional stimuli, resulting in accelerated rebound weight gain. We find that the epigenetic memory in mice can explain future transcriptional deregulation in response to further high-fat diet feeding. Together, our data suggests the existence of an obesogenic memory in mouse adipocytes, and likely other cells, largely based on stable epigenetic changes. These changes appear to prime cells to respond in a pathological manner to an obesogenic environment and may contribute to the problematic "yo-yo" effect on body weight observed with dieting. Targeting these changes could potentially improve long-term weight management and health outcomes.


## Apps to explore data
[ShinyApp to explore snRNAseq data, cell type specific gene expression analysis and epigenetic analysis ](http://nme.ethz.ch/ATMemory.html) 

## Accession Codes
GEO: GSE236580![image](https://github.com/user-attachments/assets/edab7a47-6b3e-4881-bfbf-3b1660748b6e)


## Contents of this Repository
#### 1. :file_folder: ```Rscripts```</p>
&emsp;&emsp;:file_folder: ```Epigenetics ```</p>
&emsp;&emsp;&emsp;&emsp;:file_folder: ```0_PreProcessing```&ensp;*Peak QC, annotation, and quantification*</p>
&emsp;&emsp;&emsp;&emsp;:file_folder: ```1_MOFA```&ensp;*Cut&Tag and ATACseq based multi-omics factor analysis*</p>
&emsp;&emsp;&emsp;&emsp;:file_folder: ```2_DifferentialAnalysis```&ensp;*Differential analysis of promoters, enhancers, alluvial plots, and GSEA*</p>

&emsp;&emsp;:file_folder: ```TRAP```&ensp;*Contains Rscripts for analysis of TRAP-seq data and correlation with snRNAseq data*</p>

&emsp;&emsp;:file_folder: ```snRNAseq_Mouse ```&ensp;*Contains Rscripts to analyse mouse epiAT snRNA-seq data and generate plots*</p>
&emsp;&emsp;&emsp;&emsp;:file_folder: ```0_PreProcessing```&ensp;* QC and reference mapping*</p>
&emsp;&emsp;&emsp;&emsp;:file_folder: ```1_SampleIntegration```&ensp;*Integration of snRNAseq datasets*</p>
&emsp;&emsp;&emsp;&emsp;:file_folder: ```2_DifferentialAnalysis```&ensp;*Cell-type specific DE and transcriptional retention analysis*</p>

&emsp;&emsp;:file_folder: ```snRNAseq_Human ```&ensp;*Contains Rscripts to analyse human snRNA-seq data and generate plots*</p>
&emsp;&emsp;&emsp;&emsp;:file_folder: ```0_PreProcessing```&ensp;* QC and reference mapping*</p>
&emsp;&emsp;&emsp;&emsp;:file_folder: ```1_SampleIntegration```&ensp;*Integration of snRNAseq datasets*</p>
&emsp;&emsp;&emsp;&emsp;:file_folder: ```2_DifferentialAnalysis```&ensp;*Cell-type specific DE and transcriptional retention analysis*</p>

#### 2. :file_folder: ```ChromHMM```</p>
&emsp;&emsp;:file_folder: ```ChromHMM```&ensp;*Commands and outputs of ChroHMM analysis (enhancers)*</p>
&emsp;&emsp;&emsp;&emsp;:file_folder: ```tracks```&ensp;*Enhancer bed files for adipocytes*</p>
&emsp;&emsp;&emsp;&emsp;:file_folder: ```scripts```&ensp;*Commands for Cut&Tag based ChromHMM analysis*</p>

#### 3. :file_folder: ```cellSNP_Vireo_demux```&ensp;*Contains code to run SNP demultiplexing based on tools cellSNP and vireo*</p>

#### 4. :file_folder: ```DEGs```&ensp;*Contains cell type specific DEGs for human and mouse AT (obese vs lean; weight loss vs lean)*</p>

#### 5. :file_folder: ```SessionInfo```</p>

## Associated Repositories 
1. [NextFlow Pipeline for CUT&Tag](https://github.com/vonMeyennLab/nf_cutntag)
2. [NextFlow Pipeline for RNAseq](https://github.com/vonMeyennLab/nf_rnaseq)


