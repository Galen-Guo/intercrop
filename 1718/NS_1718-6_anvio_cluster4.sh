#!/bin/bash
#$ -S /bin/bash
#$ -m eas -M galen.guo@canada.ca
#$ -o /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -e /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -N N18_profile


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
### BINNING!
####################################################################

anvi-cluster-contigs -c $ANVIO/1718/plot4/plot4_1718_contigs.db -p $ANVIO/1718/plot4/profile_merged/PROFILE.db -C plot4_concoct --driver concoct -T 24 --just-do-it
anvi-cluster-contigs -c $ANVIO/1718/plot4/plot4_1718_contigs.db -p $ANVIO/1718/plot4/profile_merged/PROFILE.db -C anvio_derep_metabat2 --driver plot4_metabat2 -T 24 --just-do-it
conda activate anvio-7.1
anvi-cluster-contigs -c $ANVIO/1718/plot4/plot4_1718_contigs.db -p $ANVIO/1718/plot4/profile_merged/PROFILE.db -S plot4_concoct,plot4_metabat2 --driver dastool --search-engine diamond -C plot4_dastool -T 24 --just-do-it
conda activate anvio7
anvi-summarize -c -c $ANVIO/1718/plot4/plot4_1718_contigs.db -p $ANVIO/1718/plot4/profile_merged/PROFILE.db -o $ANVIO/1718/plot4/sample_summary_dastool -C plot4_dastool --init-gene-coverages
d
