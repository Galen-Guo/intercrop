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
### loading into anvio ###
####################################################################
cd $MAPPING/1718/plot6
rename 's/-/_/g' *_anvi*
rename 's/NS.1718./NS1718_/g' *_anvi*


cp $WORK/NS.1718 $MAPPING/1718/plot6
sed -i 's/-/_/g' $MAPPING/1718/plot6/NS.1718
sed -i 's/NS.1718./NS1718_/g' $MAPPING/1718/plot6/NS.1718
mkdir $ANVIO/1718/plot6/profile

for sample in `awk '{print $1}' $MAPPING/1718/plot6/NS.1718`;
do
echo $sample
anvi-profile -i $MAPPING/1718/plot6/${sample}_anvi.bam -c $ANVIO/1718/plot6/plot6_1718_contigs.db --output-dir $ANVIO/1718/plot6/profile/$sample --sample-name $sample -T 24 --min-contig-length 1000

done

## merge profile

anvi-merge  $ANVIO/1718/plot6/profile/*/PROFILE.db -o $ANVIO/1718/plot6/profile_merged -c $ANVIO/1718/plot6/plot6_1718_contigs.db -S NS_1718_plot6
