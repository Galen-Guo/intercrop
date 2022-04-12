#!/bin/bash
#$ -S /bin/bash
#$ -m eas -M galen.guo@canada.ca
#$ -o /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -e /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -N humann25


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

echo concatenating fastq
for sample in `awk '{print $1}' $WORK/NS.1725`;
do

R1="${sample}.pair1.fq.gz"
R2="${sample}.pair2.fq.gz"
seqtk mergepe $TRIM/$R1 $TRIM/$R2 > $HUMANN/1725/${sample}.fastq

done

echo concatenating fastq done

echo running humann3
for i in `awk '{print $1}' $WORK/NS.1725`;
do
echo $i
humann --input $HUMANN/1725/${i}.fastq --output $HUMANN/1725/output --metaphlan-options '--bowtie2db /isilon/ottawa-rdc/users/shared/chenw_lab/galen/' --threads 24
done
echo humann3 done

echo cpm and cpm file


mkdir $HUMANN/1725/summary/


for sample in `awk '{print $1}' $WORK/NS.1725`;
do
echo $i
humann_renorm_table --input $HUMANN/1725/output/${i}_genefamilies.tsv --output $HUMANN/1725/output/cpm/${i}_genefamilies-cpm.tsv --units cpm --update-snames
humann_renorm_table --input $HUMANN/1725/output/${i}_pathcoverage.tsv --output $HUMANN/1725/output/cpm/${i}_pathcoverage-cpm.tsv --units cpm --update-snames
humann_renorm_table --input $HUMANN/1725/output/${i}_pathabundance.tsv --output $HUMANN/1725/output/cpm/${i}_pathabundance-cpm.tsv --units cpm --update-snames
done

humann_join_tables -i $HUMANN/1725/output/cpm/ -o $HUMANN/1725/output/summary/all_genefamilies.tsv --file_name genefamilies
humann_join_tables -i $HUMANN/1725/output/cpm/ -o $HUMANN/1725/output/summary/all_pathcoverage.tsv --file_name pathcoverage
humann_join_tables -i $HUMANN/1725/output/cpm/ -o $HUMANN/1725/output/summary/all_pathabundance.tsv --file_name pathabundance

echo summary done

