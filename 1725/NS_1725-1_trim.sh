#!/bin/bash
#$ -S /bin/bash
#$ -m eas -M galen.guo@canada.ca
#$ -o /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -e /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -N trimming25


source ~/miniconda3/bin/activate
conda activate anvio7

#directories

GALEN=/isilon/ottawa-rdc/users/shared/chenw_lab/galen
WORK=$GALEN/Intercrop
RAW=$WORK/raw
MAPPING=$WORK/mapping
TRIM=$WORK/trimmed_1725
MEGAHIT=$WORK/assembly
ANVIO=$WORK/anvio_1725


###################################################################
### trimming with trimmomatic, removal of adapter using bbduk
###################################################################

echo Trimming

mkdir $WORK/trimmed
TRIM=$WORK/trimmed

adapter=/home/AAFC-AAC/guog/miniconda3/envs/anvio7/opt/bbmap-38.18/resources
	

for i in `awk '{print $1}' $WORK/NS.1725`;
do
echo $i
    
    file1=${i}_R1.fastq.gz
    file2=${i}_R2.fastq.gz

	
bbduk.sh threads=80 -Xmx1g in1=$RAW/$file1 in2=$RAW/$file2 out1=$TRIM/${i}.clean1.fq out2=$TRIM/${i}.clean2.fq outm1=$TRIM/${i}.adapter1.fq outm2=$TRIM/${i}.adapter2.fq ref=$adapter/adapters.fa k=31 hdist=1 stats=$TRIM/${i}_stats.txt
	
java -jar $GALEN/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 80 $TRIM/${i}.clean1.fq $TRIM/${i}.clean2.fq $TRIM/${i}.pair1.fq.gz $TRIM/${i}.unpair1.fq.gz $TRIM/${i}.pair2.fq.gz $TRIM/${i}.unpair2.fq.gz ILLUMINACLIP:$GALEN/Trimmomatic-0.39/adapters/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

rm $TRIM/${i}.clean1.fq
rm $TRIM/${i}.clean2.fq
rm $TRIM/${i}.adapter1.fq
rm $TRIM/${i}.adapter2.fq
done

echo Trimming_done

####################################################################
### QC of fastqc post-trim
####################################################################
echo QC_of_Trimming

mkdir $WORK/fastqc/trimmed_clean_fq 
cd $WORK/fastqc/trimmed_clean_fq/
fastqc -t 80 $TRIM/*.pair* -o $WORK/fastqc/trimmed_clean_fq 
multiqc $WORK/fastqc/trimmed_clean_fq/ $WORK/fastqc/trimmed_clean_fq  -n NS.1725_trimmed_clean_fastqc_report

echo QC_of_Trimming_done

