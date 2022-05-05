#!/bin/bash
#$ -S /bin/bash
#$ -m eas -M galen.guo@canada.ca
#$ -o /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -e /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -N anvi18


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
for i in plot5;
do
echo $i
mkdir $ANVIO
mkdir $ANVIO/1718
mkdir $ANVIO/1718/$i/

anvi-script-reformat-fasta $MEGAHIT/1718/$i/final.contigs.fa -o $ANVIO/1718/$i/${i}_1718_contigs.fa --simplify-names

###create anvio database



#Gene calling

anvi-gen-contigs-database -f $ANVIO/1718/$i/${i}_1718_contigs.fa -o $ANVIO/1718/$i/${i}_1718_contigs.db -n intercrop_1718 -T 24

anvi-run-hmms -c $ANVIO/1718/$i/${i}_1718_contigs.db -T 24

anvi-run-kegg-kofams -c $ANVIO/1718/$i/${i}_1718_contigs.db -T 24 --hmmer-program hmmsearch --just-do-it
anvi-run-ncbi-cogs -c $ANVIO/1718/$i/${i}_1718_contigs.db --cog-data-dir $GALEN/COG --num-threads 24 --search-with diamond 
anvi-run-pfams -c $ANVIO/1718/$i/${i}_1718_contigs.db --pfam-data-dir $GALEN/database/pfam --hmmer-program hmmscan --num-threads 24


anvi-run-scg-taxonomy -c $ANVIO/1718/$i/${i}_1718_contigs.db -T 24
anvi-estimate-scg-taxonomy -c $ANVIO/1718/$i/${i}_1718_contigs.db --metagenome-mode --output-file $ANVIO/1718/$i/${i}_1718_scg_est_output.txt


anvi-get-sequences-for-gene-calls -c $ANVIO/1718/$i/${i}_1718_contigs.db --get-aa-sequences -o $ANVIO/1718/$i/${i}_1718_amino-acid-sequences.fa

anvi-display-contigs-stats -c $ANVIO/1718/$i/${i}_1718_contigs.db --report-as-text  --output-file $ANVIO/1718/$i/${i}_1718_contigs_stats.txt

done
