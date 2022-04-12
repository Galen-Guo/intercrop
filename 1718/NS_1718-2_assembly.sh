#!/bin/bash
#$ -S /bin/bash
#$ -m eas -M galen.guo@canada.ca
#$ -o /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -e /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -N assembly18


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
### co-assembly megahit
### 
####################################################################

echo co-assembling 1718, by sequencing length


F=$(echo $(ls $TRIM/NS.1718.F8P*.pair1.fq.gz)  | sed "s/ /,/g")
R=$(echo $(ls $TRIM/NS.1718.F8P*.pair2.fq.gz)  | sed "s/ /,/g")

megahit -1 $F -2 $R -o $MEGAHIT/1718/ -t 80 --min-contig-len 1000 --presets meta-large --continue  --k-min 27 --k-max 127 --k-step 10 --kmin-1pass --no-mercy 





