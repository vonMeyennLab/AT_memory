# :file_folder: cellSNP-lite/Vireo Processing

## Files

- [sc_snp.smk](sc_snp.smk) (*Required*): Snakemake workflow file.
- [sc-snp-demux.yml](sc-snp-demux.yml) (*Required*): The conda environment yml file. Is used by the Snakemake workflow to install all required packages
- [config.yml](config.yml) (*Required*): Input configuration specifying 3 important parameters:
	1. The path to the 10X output folder. Must contain "possorted_genome_bam.bam" and "filtered_feature_bc_matrix"
	2. The path to the directory containing the reference SNPs
	3. The VCF file name. We used 'genome1K.phase3.SNP_AF5e2.chr1toX.hg38' from the (1000 genome project)[https://sourceforge.net/projects/cellsnp/files/SNPlist/].
- [sc_snp_scaling.sh](sc_snp_scaling.sh): Example slurm script calling the Snakemake workflow on some samples

## Workflow description

To perform SNP-calling and demultiplexing on the pooled samples, cellsnp-lite is first used to call SNPs on a cell-level using the 1000 genomes-based reference variant call file for hg38 at a resolution of 7.4 million SNPs. SNPs with <20 counts and a <10% minor allele frequency were filtered out, as per the developer recommendations. Finally, the tool vireo is called on the outputs from cellSNP-lite to demultiplex the pooled data using the cellsnp-lite-derived genotype information.