#!/bin/bash
#$ -S /bin/bash
#$ -m eas -M galen.guo@canada.ca
#$ -o /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -e /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -N N25_profile


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
### loading into anvio ### 
####################################################################
# cd $MAPPING/1725
# rename 's/-/_/g' *_anvi*
# rename 's/NS.1725./NS1725_/g' *_anvi*
 ## ### do everything in one file, using more cpu. (5 days)

# cp $WORK/NS.1725 $MAPPING/1725/
# sed -i 's/-/_/g' $MAPPING/1725/NS.1725
# sed -i 's/NS.1725./NS1725_/g' $MAPPING/1725/NS.1725
# mkdir $ANVIO/1725/profile

for sample in `awk '{print $1}' $MAPPING/1725/NS.1725`;
do
echo $sample 
anvi-profile -i $MAPPING/1725/${sample}_anvi.bam -c $ANVIO/1725/1725_contigs.db --output-dir $ANVIO/1725/profile/$sample --sample-name $sample -T 24 --min-contig-length 1000

done

## merge profile

anvi-merge  $ANVIO/1725/profile/*/PROFILE.db -o $ANVIO/1725/profile_merged -c $ANVIO/1725/1725_contigs.db -S NS_1725

