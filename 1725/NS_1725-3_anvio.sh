#!/bin/bash
#$ -S /bin/bash
#$ -m eas -M galen.guo@canada.ca
#$ -o /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -e /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -N anvi28


source ~/miniconda3/bin/activate
conda activate anvio7

#directories

GALEN=/isilon/ottawa-rdc/users/shared/chenw_lab/galen
WORK=$GALEN/Intercrop
RAW=$WORK/raw
MAPPING=$WORK/mapping
TRIM=$WORK/trimmed
MEGAHIT=$WORK/assembly
ANVIO=$WORK/anvio


####################################################################
### loading into anvio ### 
####################################################################

mkdir $ANVIO

mkdir $ANVIO/1725/

anvi-script-reformat-fasta $MEGAHIT/1725/final.contigs.fa -o $MEGAHIT/1725/contigs-fixed.fa --simplify-names

mv $MEGAHIT/1725/contigs-fixed.fa $ANVIO/1725/1725_contigs.fa

###create anvio database



#Gene calling

anvi-gen-contigs-database -f $ANVIO/1725/1725_contigs.fa -o $ANVIO/1725/1725_contigs.db -n GRDI -T 80

anvi-run-hmms -c $ANVIO/1725/1725_contigs.db -T 80

anvi-run-kegg-kofams -c $ANVIO/1725/1725_contigs.db -T 80 --hmmer-program hmmsearch --just-do-it
anvi-run-ncbi-cogs -c$ANVIO/1725/1725_contigs.db --cog-data-dir $GALEN/COG --num-threads 80 --search-with diamond --temporary-dir-path $ANVIO/temp
anvi-run-pfams -c $ANVIO/1725/1725_contigs.db --pfam-data-dir $GALEN/database/pfam --hmmer-program hmmscan --num-threads 80


anvi-run-scg-taxonomy -c $ANVIO/1725/1725_contigs.db -T 80
anvi-estimate-scg-taxonomy -c $ANVIO/1725/1725_contigs.db --metagenome-mode --output-file $ANVIO/1725/1725_scg_est_output.txt


anvi-get-sequences-for-gene-calls -c $ANVIO/1725/1725_contigs.db --get-aa-sequences -o $ANVIO/1725/1725_amino-acid-sequences.fa

anvi-display-contigs-stats $ANVIO/1725/1725_contigs.db  --report-as-text  --output-file $ANVIO/1725/1725_contigs_post-hmm-cogs_kegg.txt



