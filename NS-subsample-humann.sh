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
HUMANN=$WORK/humann3/subsample
ANVIO=$WORK/anvio/subsample



#random pick a 100000 seq from each fastq file

for sample in `awk '{print $1}' $WORK/NS.1718`;
do

R1="${sample}.pair1.fq.gz"
R2="${sample}.pair2.fq.gz"
seqtk mergepe $TRIM/$R1 $TRIM/$R2 > $WORK/humann3/${sample}.fq.gz

done



for i in `awk '{print $1}' $WORK/NS.1718`;
do
echo $i
seqtk sample -s100 $WORK/humann3/${sample}.fq 100000 > $WORK/humann3/${sample}_100000.fastq
humann --input $WORK/humann3/${i}.fastq --output $HUMANN/ --metaphlan-options '--bowtie2db /isilon/ottawa-rdc/users/shared/chenw_lab/galen/' --threads 50
done


mkdir $HUMANN/cpm/

mkdir $HUMANN/summary/



for i in `awk '{print $1}' $WORK/NS.1718`;
do
echo $i
humann_renorm_table --input $HUMANN/${i}_1mil_genefamilies.tsv --output $HUMANN/cpm/${i}_1mil_genefamilies-cpm.tsv --units cpm --update-snames
humann_renorm_table --input $HUMANN/${i}_1mil_pathcoverage.tsv --output $HUMANN/cpm/${i}_1mil_pathcoverage-cpm.tsv --units cpm --update-snames
humann_renorm_table --input $HUMANN/${i}_1mil_pathabundance.tsv --output $HUMANN/cpm/${i}_1mil_pathabundance-cpm.tsv --units cpm --update-snames
done

humann_join_tables -i $HUMANN/cpm/ -o $HUMANN/summary/all_genefamilies.tsv --file_name genefamilies
humann_join_tables -i $HUMANN/cpm/ -o $HUMANN/summary/all_pathcoverage.tsv --file_name pathcoverage
humann_join_tables -i $HUMANN/cpm/ -o $HUMANN/summary/all_pathabundance.tsv --file_name pathabundance



