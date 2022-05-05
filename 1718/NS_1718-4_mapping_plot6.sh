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
ANVIO=$WORK/anvio


####################################################################
### loading into anvio ###
####################################################################


echo mapping

## create mapping directory


## building mapping files
mkdir $MAPPING
mkdir $MAPPING/1718/
mkdir $MAPPING/1718/plot6
cd $MAPPING/1718/plot6



bowtie2-build $ANVIO/1718/plot6/plot6_1718_contigs.fa $MAPPING/1718/plot6/plot6_contigs --threads 24 --large-index

### mapping read to contigs


for sample in `awk '{print $1}' $WORK/NS.1718`;
do

R1="${sample}.pair1.fq.gz"
R2="${sample}.pair2.fq.gz"
sam="${sample}.sam"
bowtie2 -x $MAPPING/1718/plot6_contigs -1 $TRIM/plot6/$R1 -2 $TRIM/plot6/$R2 -S $MAPPING/1718/plot6/$sam --threads 24
done


conda deactivate


echo convert sam to bam

for sample in `awk '{print $1}' $WORK/NS.1718`;
do
sam=".sam"
bam=".bam"
samtools view -S -b $MAPPING/1718/plot6/${sample}.sam > $MAPPING/1718/plot6/${sample}.bam
done


conda activate anvio7


## bam to anvio_bam profiling

for sample in `awk '{print $1}' $WORK/NS.1718`;
do
done
anvi-init-bam $MAPPING/1718/plot6/${sample}.bam -o $MAPPING/1718/plot6/${sample}_anvi.bam


echo mapping_done
