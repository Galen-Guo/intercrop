#!/bin/bash
#$ -S /bin/bash
#$ -m eas -M galen.guo@canada.ca
#$ -o /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -e /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -N anvi18map


source ~/miniconda3/bin/activate
conda activate anvio7

#directories
GALEN=/isilon/ottawa-rdc/users/shared/chenw_lab/galen
WORK=$GALEN/Intercrop/
RAW=$WORK/raw
MAPPING=$WORK/mapping
TRIM=$WORK/trimmed
MEGAHIT=$WORK/assembly

ANVIO=$WORK/anvio/subsample




#Gene calling
anvi-run-hmms -c $ANVIO/subsample.db -T 80

anvi-run-kegg-kofams -c $ANVIO/subsample.db -T 80 --hmmer-program hmmsearch --just-do-it
anvi-run-ncbi-cogs -c$ANVIO/subsample.db --cog-data-dir $GALEN/COG --num-threads 80 --search-with diamond --temporary-dir-path $ANVIO/temp
anvi-run-scg-taxonomy -c $ANVIO/subsample.db -T 80
anvi-estimate-scg-taxonomy -c $ANVIO/subsample.db --metagenome-mode --output-file $ANVIO/subsample_scg_est_output.txt


anvi-get-sequences-for-gene-calls -c $ANVIO/subsample.db --get-aa-sequences -o $ANVIO/subsample_amino-acid-sequences.fa

anvi-display-contigs-stats $ANVIO/subsample.db  --report-as-text  --output-file $ANVIO/subsample_contigs_post-hmm-cogs_kegg.txt



anvi-estimate-metabolism -c $ANVIO/subsample.db -O $ANVIO/metabolism \
                         --matrix-format -T 40

####################################################################
### loading into anvio ### 
####################################################################


for sample in `awk '{print $1}' $MAPPING/NS.1718`;
do
echo "$sample$file_ext"
anvi-profile -i $MAPPING/${sample}_anvi.bam -c $ANVIO/contigs.db --output-dir $ANVIO/profile/$sample --sample-name $sample -T 100 --min-contig-length 1000
done

## merge profile

anvi-merge  $ANVIO/profile/*/PROFILE.db -o $ANVIO/profile_merged -c $ANVIO/contigs.db -S Intercrop_subsample 



