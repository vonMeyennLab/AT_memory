rule cellsnp:
    input:
        bam=os.path.join(config['SAMPLE_DIR'], "{sample}", "possorted_genome_bam.bam"),
        bc=os.path.join(config['SAMPLE_DIR'], "{sample}", "filtered_feature_bc_matrix", "barcodes.tsv.gz"),
        ref=os.path.join(config['REF_DIR'], config['REF_VCF_NAME'] + ".vcf.gz")
    output:
        directory("{sample}_" + config['REF_VCF_NAME'])
    threads: config["THREADS"] 
    conda: "sc-snp-demux.yml"
    shell:
        """
        cellsnp-lite --samFile {input.bam} -R {input.ref} --barcodeFile {input.bc} --outDir {output} -p {threads} --minMAF 0.1 --minCOUNT 20 --gzip 2> {wildcards.sample}_cellsnp.log
        """

rule vireo:
    input:
        "{sample}_" + config['REF_VCF_NAME']
    output:
        directory("{sample}_{n_samples}N_" + config["REF_VCF_NAME"] + ".demux")
    conda: "sc-snp-demux.yml"
    threads: config["THREADS"]
    shell:
        """
        vireo -c {input} -N {wildcards.n_samples} -o {output} 2> {wildcards.sample}_vireo.log
        """
