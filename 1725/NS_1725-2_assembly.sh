#!/bin/bash
#$ -S /bin/bash
#$ -m eas -M galen.guo@canada.ca
#$ -o /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -e /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -N assembly25

source ~/miniconda3/bin/activate
conda activate anvio7

#directories

GALEN=/isilon/ottawa-rdc/users/shared/chenw_lab/galen
WORK=$GALEN/Intercrop
RAW=$WORK/raw
output=$WORK/output
input=$WORK/input
MAPPING=$WORK/mapping
TRIM=$WORK/trimmed
MEGAHIT=$WORK/assembly
ANVIO=$WORK/anvio

####################################################################
### co-assembly megahit
### 
####################################################################

echo co-assembling 1725, by sequencing length


F=$(echo $(ls $TRIM/NS.1725.F8P*.pair1.fq.gz)  | sed "s/ /,/g")
R=$(echo $(ls $TRIM/NS.1725.F8P*.pair2.fq.gz)  | sed "s/ /,/g")

megahit -1 $F -2 $R -o $MEGAHIT/1725 -t 80 --min-contig-len 1000  --k-min 37 --k-step 10 --k-min 127 --no-mercy --continue


#megahit -1 $F -2 $R -o $MEGAHIT/1725 -t 80 --min-contig-len 1000  --k-min 37 --k-step 10 --k-min 127 --no-mercy