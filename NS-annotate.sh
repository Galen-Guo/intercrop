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

ANVIO=$WORK/anvio/subsample

krakendb=/isilon/common/reference/databases/kraken2_ncbi_nt

for i in `awk '{print $1}' $WORK/NS.1718`;
do
echo $i
    
    file1=${i}.pair1.fq.gz
    file2=${i}.pair2.fq.gz

kraken2 --db $krakendb  --threads 20 --paired --quick --gzip-compressed \
--output $WORK/subsample_kraken_out/ \
$TRIM/$file1 $TRIM/$file2

done
