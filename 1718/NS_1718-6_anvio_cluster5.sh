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

anvi-cluster-contigs -c $ANVIO/1718/plot5/plot5_1718_contigs.db -p $ANVIO/1718/plot5/profile_merged/PROFILE.db -C plot5_concoct --driver concoct -T 24 --just-do-it
anvi-cluster-contigs -c $ANVIO/1718/plot5/plot5_1718_contigs.db -p $ANVIO/1718/plot5/profile_merged/PROFILE.db -C plot5_metabat2 --driver metabat2 -T 24 --just-do-it
conda activate anvio-7.1
anvi-cluster-contigs -c $ANVIO/1718/plot5/plot5_1718_contigs.db -p $ANVIO/1718/plot5/profile_merged/PROFILE.db -S plot5_concoct,plot5_metabat2 --driver dastool --search-engine diamond -C plot5_dastool -T 24 --just-do-it
conda activate anvio7
anvi-summarize -c -c $ANVIO/1718/plot5/plot5_1718_contigs.db -p $ANVIO/1718/plot5/profile_merged/PROFILE.db -o $ANVIO/1718/plot5/sample_summary_dastool -C plot5_dastool --init-gene-coverages
