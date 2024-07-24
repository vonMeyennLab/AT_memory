##################################################
Concatenated design
##################################################

# Binarize Bam Files
sbatch --time=4:00:00 --ntasks=4 --mem-per-cpu=6000 --wrap="java -mx20000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
BinarizeBam \
-paired /cluster/work/nme/software/apps/ChromHMM/1.22/CHROMSIZES/mm10_nochr.txt \
/scratch/../ChromHMM_design_concatenated_ATmem.txt \
/scratch/../binary_files" \
/

# Learn Model
sbatch --time=24:00:00 --ntasks=20 --mem-per-cpu=5000 --wrap="seq 1 30 | parallel java -mx4000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
LearnModel \
-s 123 -nobrowser -noimage -nobed -noenrich -noautoopen \
/scratch/../binary_files/ \
/scratch/../model/ \
{} \
GRCm38" \
/

# Compare Models
sbatch --time=4:00:00 --ntasks=10 --mem-per-cpu=2000 --wrap="java -mx20000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
CompareModels \
/scratch/../model/emissions_30.txt \
/scratch/../model/ \
/scratch/../compare_models/concatenated_1_to_30" \
/

# Make Segmentation
sbatch --time=24:00:00 --ntasks=10 --mem-per-cpu=6000 --wrap="ls /scratch/../model/model_8.txt | \
parallel java -mx5000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
MakeSegmentation \
{} \
/scratch/../binary_files \
/scratch/../segmentation" \
/

# Genomic Features Enrichment
sbatch --time=4:00:00 --ntasks=10 -R --mem-per-cpu=2000 --wrap="java -mx4000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
OverlapEnrichment \
/scratch/../segmentation/ATmem_8_segments.bed \
/nfs/nas12.ethz.ch/fs1201/green_groups_nme_public/AdhidebGhosh/ChromHMM/annotatr_mm/ \
/scratch/../annotatr_enrichment/ATmem_genomicFeature_enrichment" \
/

# ENCODE cCRE Enrichment
sbatch --time=4:00:00 --ntasks=10 -R --mem-per-cpu=2000 --wrap="java -mx4000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
OverlapEnrichment \
/scratch/../segmentation/ATmem_8_segments.bed \
/nfs/nas12.ethz.ch/fs1201/green_groups_nme_public/AdhidebGhosh/ChromHMM/ENCODE_mm_collapsed/ \
/scratch/../cCRE_enrichment/ATmem_cCRE_enrichment" \
/