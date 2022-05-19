#!/bin/bash
#$ -S /bin/bash
#$ -m eas -M galen.guo@canada.ca
#$ -o /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -e /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -N N18_clust6


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

anvi-cluster-contigs -c $ANVIO/1718/plot6/plot6_1718_contigs.db -p $ANVIO/1718/plot6/profile_merged/PROFILE.db -C plot6_concoct --driver concoct -T 24 --just-do-it
anvi-cluster-contigs -c $ANVIO/1718/plot6/plot6_1718_contigs.db -p $ANVIO/1718/plot6/profile_merged/PROFILE.db -C plot6_metabat2 --driver metabat2 -T 24 --just-do-it
conda activate anvio-7.1
anvi-cluster-contigs -c $ANVIO/1718/plot6/plot6_1718_contigs.db -p $ANVIO/1718/plot6/profile_merged/PROFILE.db -C plot6_maxbin2 --driver maxbin2 -T 24 --just-do-it

anvi-cluster-contigs -c $ANVIO/1718/plot6/plot6_1718_contigs.db -p $ANVIO/1718/plot6/profile_merged/PROFILE.db -S plot6_concoct,plot6_metabat2,plot6_maxbin2 --driver dastool --search-engine diamond -C plot6_dastool -T 24 --just-do-it
conda activate anvio7
anvi-summarize -c -c $ANVIO/1718/plot6/plot6_1718_contigs.db -p $ANVIO/1718/plot6/profile_merged/PROFILE.db -o $ANVIO/1718/plot6/sample_summary_dastool -C plot6_dastool --init-gene-coverages
