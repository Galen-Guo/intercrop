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
# for i in `awk '{print $1}' $WORK/NS.1718.plot34`;
# do
# echo $i
# mv ${i}.pair1.fq.gz $TRIM/plot4
# mv ${i}.pair2.fq.gz $TRIM/plot4
# done
#
# for i in `awk '{print $1}' $WORK/NS.1718.plot5`;
# do
# echo $i
# mv ${i}.pair1.fq.gz $TRIM/plot5
# mv ${i}.pair2.fq.gz $TRIM/plot5
# done
#
# for i in `awk '{print $1}' $WORK/NS.1718.plot6`;
# do
# echo $i
# mv ${i}.pair1.fq.gz $TRIM/plot6
# mv ${i}.pair2.fq.gz $TRIM/plot6
# done


echo co-assembling 1718, by sequencing length

F=$(echo $(ls $TRIM/plot6/NS.1718.F8P*.pair1.fq.gz)  | sed "s/ /,/g")
R=$(echo $(ls $TRIM/plot6/NS.1718.F8P*.pair2.fq.gz)  | sed "s/ /,/g")

megahit -1 $F -2 $R -o $MEGAHIT/1718/plot6/ -t 12 --min-contig-len 1000 --continue  --k-min 37 --k-max 77 --k-step 10