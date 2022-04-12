#!/bin/bash
#$ -S /bin/bash
#$ -m eas -M galen.guo@canada.ca
#$ -o /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -e /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -N anvi18map


source ~/miniconda3/bin/activate
conda activate humann3_env

#directories

GALEN=/isilon/ottawa-rdc/users/shared/chenw_lab/galen
WORK=$GALEN/Intercrop/
RAW=$WORK/raw
MAPPING=$WORK/mapping
TRIM=$WORK/trimmed
MEGAHIT=$WORK/assembly
mkdir $WORK/anvio/subsample
ANVIO=$WORK/anvio/subsample
HUMANN=$WORK/humann3/subsample



for i in `awk '{print $1}' $WORK/NS.1718`;
do
echo $i
seqtk sample -s100 $WORK/humann3/${i}.fq 1000000 > $WORK/humann3/subsample/${i}_1mil.fq
humann --input $WORK/humann3/subsample/${i}_1mil.fq --output $HUMANN --metaphlan-options '--bowtie2db /isilon/ottawa-rdc/users/shared/chenw_lab/galen/' --threads 40
done



mkdir $HUMANN/summary/


for sample in `awk '{print $1}' $WORK/NS.1718`;
do
echo $i
humann_renorm_table --input $HUMANN/${i}_1mil_genefamilies.tsv --output $HUMANN/cpm/${i}_genefamilies-cpm.tsv --units cpm --update-snames
humann_renorm_table --input $HUMANN/${i}_1mil_pathcoverage.tsv --output $HUMANN/cpm/${i}_pathcoverage-cpm.tsv --units cpm --update-snames
humann_renorm_table --input $HUMANN/${i}_1mil_pathabundance.tsv --output $HUMANN/cpm/${i}_pathabundance-cpm.tsv --units cpm --update-snames
done

humann_join_tables -i $HUMANN/cpm/ -o $HUMANN/summary/all_genefamilies.tsv --file_name genefamilies
humann_join_tables -i $HUMANN/cpm/ -o $HUMANN/summary/all_pathcoverage.tsv --file_name pathcoverage
humann_join_tables -i $HUMANN/cpm/ -o $HUMANN/summary/all_pathabundance.tsv --file_name pathabundance



