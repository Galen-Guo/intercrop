#!/bin/bash
#$ -S /bin/bash
#$ -m eas -M galen.guo@canada.ca
#$ -o /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -e /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -N humann18


source ~/miniconda3/bin/activate
conda activate humann3_env

#directories
GALEN=/isilon/ottawa-rdc/users/shared/chenw_lab/galen
WORK=$GALEN/Intercrop/
RAW=$WORK/raw
MAPPING=$WORK/mapping
TRIM=$WORK/trimmed
MEGAHIT=$WORK/assembly
HUMANN=$WORK/humann3/


#for sample in `awk '{print $1}' $WORK/NS.1718`;
#do

#R1="${sample}.pair1.fq.gz"
#R2="${sample}.pair2.fq.gz"
#seqtk mergepe $TRIM/$R1 $TRIM/$R2 > $WORK/humann3/${sample}.fastq
#done


for i in `awk '{print $1}' $WORK/NS.1718`;
do
echo $i
humann --input $HUMANN/${i}.fastq --output $HUMANN/output --metaphlan-options '--bowtie2db /isilon/ottawa-rdc/users/shared/chenw_lab/galen/' --threads 24
done



mkdir $HUMANN/summary/


for sample in `awk '{print $1}' $WORK/NS.1718`;
do
echo $i
humann_renorm_table --input $HUMANN/output/${i}_genefamilies.tsv --output $HUMANN/output/cpm/${i}_genefamilies-cpm.tsv --units cpm --update-snames
humann_renorm_table --input $HUMANN/output/${i}_pathcoverage.tsv --output $HUMANN/output/cpm/${i}_pathcoverage-cpm.tsv --units cpm --update-snames
humann_renorm_table --input $HUMANN/output/${i}_pathabundance.tsv --output $HUMANN/output/cpm/${i}_pathabundance-cpm.tsv --units cpm --update-snames
done

humann_join_tables -i $HUMANN/output/cpm/ -o $HUMANN/output/summary/all_genefamilies.tsv --file_name genefamilies
humann_join_tables -i $HUMANN/output/cpm/ -o $HUMANN/output/summary/all_pathcoverage.tsv --file_name pathcoverage
humann_join_tables -i $HUMANN/output/cpm/ -o $HUMANN/output/summary/all_pathabundance.tsv --file_name pathabundance



