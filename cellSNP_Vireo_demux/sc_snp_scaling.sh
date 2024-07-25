#!/bin/bash

#SBATCH --job-name="sc_snp_scaling"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20

set +e

TMP_DIR="."
echo $TMP_DIR
cd $TMP_DIR

. /usr/local/ngseq/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

snakemake -p \
	--cores 20 \
	--use-conda \
	-s sc_snp.smk \
    [Sample1]_[Num Cluster1]N_[Ref VCF Name].demux [Sample2]_[Num Cluster2]N_[Ref VCF Name].demux ... \
	--configfile config.yml

cd ..
