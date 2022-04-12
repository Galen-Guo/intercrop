#!/bin/bash
#$ -S /bin/bash


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

#renaming file

for i in *NS.17188P*; do mv -v "$i" "${i/.IDT*.F/}"; done


rename -n 's/NS.17188P/NS.1718.F8P/' NS.17188P*
rename -n 's/NS.1725.0018P/NS.1725.F8P/' NS.1725.0018P*

# do the same with the other prefix
####################################################################
### QC of fastqc pre-trim
####################################################################


mkdir $WORK/fastqc && mkdir $WORK/fastqc/raw_fq 
fastqc $RAW/*fastq.gz -t 10 -o $WORK/fastqc/raw_fq -f fastq
multiqc $WORK/fastqc/raw_fq -o $WORK/fastqc/ -n raw_fastqc_report

###################################################################
### trimming with trimmomatic, removal of adapter using bbduk
###################################################################

echo Trimming

mkdir $WORK/trimmed
TRIM=$WORK/trimmed

adapter=/home/AAFC-AAC/guog/miniconda3/envs/anvio7/opt/bbmap-38.18/resources
	

for i in `awk '{print $1}' $RAW/sample_name`;
do
echo $i
    
    file1=${i}_R1.fastq.gz
    file2=${i}_R2.fastq.gz

	
bbduk.sh threads=80 -Xmx1g in1=$RAW/$file1 in2=$RAW/$file2 out1=$TRIM/${i}.clean1.fq out2=$TRIM/${i}.clean2.fq outm1=$TRIM/${i}.adapter1.fq outm2=$TRIM/${i}.adapter2.fq ref=$adapter/adapters.fa k=31 hdist=1 stats=$TRIM/${i}_stats.txt
	
java -jar $GALEN/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 80 $TRIM/${i}.clean1.fq $TRIM/${i}.clean2.fq $TRIM/${i}.pair1.fq.gz $TRIM/${i}.unpair1.fq.gz $TRIM/${i}.pair2.fq.gz $TRIM/${i}.unpair2.fq.gz ILLUMINACLIP:$GALEN/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:6 CROP:120

done

echo Trimming_done

####################################################################
### QC of fastqc post-trim
####################################################################
echo QC_of_Trimming

mkdir $WORK/fastqc/trimmed_clean_fq 
fastqc -t 80 $TRIM/*.pair* -o $WORK/fastqc/trimmed_clean_fq 
multiqc $WORK/fastqc/trimmed_clean_fq/ $WORK/fastqc/ -n trimmed_clean_fastqc_report

echo QC_of_Trimming_done


####################################################################
### co-assembly megahit
### 
####################################################################

echo co-assembling, by sequencing length

mkdir $MEGAHIT
F=$(echo $(ls $TRIM/NS.1725.F8P*.pair1.fq.gz)  | sed "s/ /,/g")
R=$(echo $(ls $TRIM/NS.1725.F8P*.pair2.fq.gz)  | sed "s/ /,/g")

megahit -1 $F -2 $R -o $MEGAHIT -t 96 --min-contig-len 1000 --presets meta-large --continue  k_min=27 k_max=127 --no-mercy
 

####################################################################
### loading into anvio ### 
####################################################################




anvi-script-reformat-fasta $MEGAHIT/final.contigs.fa -o $MEGAHIT/contigs-fixed.fa -l 0 --simplify-names

mv $MEGAHIT/contigs-fixed.fa $MEGAHIT/contigs.fa

###create anvio database

mkdir $ANVIO

#Gene calling

anvi-gen-contigs-database -f $MEGAHIT/contigs.fa -o $ANVIO/contigs.db -n GRDI -T 80

anvi-run-hmms -c $ANVIO/contigs.db -T 80


anvi-run-kegg-kofams -c contigs.db -T 80 --hmmer-program hmmsearch --keep-all-hits --just-do-it


#anvi-setup-ncbi-cogs $ANVIO/COG --num-threads 80

anvi-run-ncbi-cogs -c contigs.db --cog-data-dir $ANVIO/COG --num-threads 80 --search-with diamond --temporary-dir-path $ANVIO/temp


anvi-get-sequences-for-gene-calls -c $ANVIO/contigs.db \
                                    --get-aa-sequences \
                                    -o $ANVIO/amino-acid-sequences.fa

anvi-display-contigs-stats $ANVIO/contigs.db  --report-as-text  --output-file $ANVIO/contigs_post-hmm-cogs_kegg.txt



####################################################################
### mapping with bowtie2
### 
####################################################################

echo mapping

## create mapping directory


## building mapping files
mkdir $MAPPING
cd $MAPPING

bowtie2-build $MEGAHIT/contigs.fa $MAPPING/trimmomatic_merged_contigs --threads 100

## transfering sample name file to mapping folder.

cp $RAW/sample_name $MAPPING

### mapping read to contigs

### mapping read to contigs

cd $TRIM


for sample in `awk '{print $1}' $TRIM/sample_name`;
do
R1="${sample}.pair1.fq.gz"
R2="${sample}.pair2.fq.gz"
sam="${sample}.sam" 
bowtie2 -x $MAPPING/trimmomatic_merged_contigs -1 $R1 -2 $R2 -S $MAPPING/$sam --threads 100
done


conda deactivate


echo convert sam to bam

for sample in `awk '{print $1}' $MAPPING/sample_name`;
do
sam=".sam"
bam=".bam"
samtools view -S -b $MAPPING/${sample}.sam > $MAPPING/${sample}.bam
done


conda activate anvio7


## bam to anvio_bam profiling

for sample in `awk '{print $1}' $MAPPING/sample_name`;
do
anvi-init-bam $MAPPING/${sample}.bam -o $MAPPING/${sample}_anvi.bam
done


echo mapping_done


####################################################################
### GENE FUNCTION CALLS
### interproscan on amino-acid seq from anvio
####################################################################

## split amino acid file into multiple smaller fasta files (5000 seq per file)

mkdir $ANVIO/iprs

perl /isilon/ottawa-rdc/users/shared/chenw_lab/galen/split_fasta.pl -i $ANVIO/amino-acid-sequences.fa -o $ANVIO/iprs/iprs -n 5000

## Create sample list of all file created do manually

# ls $ANVIO/iprs/ > $ANVIO/iprs/iprs_list

## manually divided into smaller files   using split -l 100 filename prefix

split -l 200 $ANVIO/iprs/iprs_list $ANVIO/iprs/iprs_list


## repeat below for each list.

mkdir $ANVIO/interpro
for sample in `awk '{print $1}' $ANVIO/iprs/iprs_listaa`;
do
/isilon/ottawa-rdc/users/shared/chenw_lab/galen/interproscan-5.39-77.0/interproscan.sh -i $ANVIO/iprs/${sample} -f tsv \
	-d $ANVIO/interpro/ \
	--tempdir $WORK/temp/ \
	--disable-precalc \
	-appl Pfam,PIRSF,SUPERFAMILY,TIGRFAM \
	--iprlookup \
	--goterms \
	--pathways
done

### create new folder to house concatenate of all smaller fxnal annotation output (the script dont like a folder with too much clutter, im guessing)

mkdir $ANVIO/interpro/iprs_output

cat $ANVIO/interpro/*.tsv > $ANVIO/interpro/iprs_output/all_iprs.tsv

### script to clean up and allow import to anvio

## create iprs2anvio.sh file found here: https://github.com/xvazquezc/stuff/blob/master/iprs2anvio.sh

/isilon/ottawa-rdc/users/shared/chenw_lab/galen/interproscan-5.39-77.0/iprs2anvio.sh -i $ANVIO/interpro/iprs_output/all_iprs.tsv -o $ANVIO/interpro/all_iprs -g -p -r


### importing functional annotation to anvio

anvi-import-functions -c $ANVIO/contigs.db -i $ANVIO/interpro/iprs_output/all_iprs_iprs2anvio.tsv


####################################################################
### TAXONOMY CALLS
### centrifuge  on amino-acid seq from anvio
####################################################################
mkdir $ANVIO/centrifuge


centrifuge -f -x $CENTRIFUGE_BASE/p+h+v/p+h+v $ANVIO/amino-acid-sequences.fa -S $ANVIO/centrifuge/centrifuge_hits.tsv

## Make sure there is two files in the work directory ($ANVIO/centrifuge/)

anvi-import-taxonomy-for-genes -c $ANVIO/contigs.db -i $ANVIO/centrifuge/centrifuge_report.tsv $ANVIO/centrifuge/centrifuge_hits.tsv -p centrifuge

### adding Single copy gene into database

# run this first, only once.

anvi-setup-scg-databases


anvi-run-scg-taxonomy -c $ANVIO/contigs.db -T 10

anvi-estimate-genome-taxonomy \ 
				-c $ANVIO/contigs.db \
                              	--metagenome-mode \
				--output-file $ANVIO/scg_est_output.txt

# to run after profiling!

anvi-estimate-genome-taxonomy -c $ANVIO/contigs.db \
                              -p $ANVIO/profile_merged/PROFILE.db \
                              --metagenome-mode \
                              --compute-scg-coverages



####################################################################
### Profiling cont'd very long. bam --> anvio
### PROFILE DONT LIKE SAMPLE NAME WITH "-" 
####################################################################

### change bam file name from "-" to "_"

cd $MAPPING
find . -depth -name '*-*' -exec rename '-' '_' {} +


## ### do everything in one file, using more cpu. (5 days)

mkdir $ANVIO/profile
cp $TRIM/sample_name $ANVIO/profile


file_ext="_anvi.bam"
for sample in `awk '{print $1}' $ANVIO/profile/sample_name`;
do
echo "$sample$file_ext"
anvi-profile -i $MAPPING/"$sample$file_ext" -c $ANVIO/contigs.db --output-dir $ANVIO/profile/$sample --sample-name $sample -T 100 --min-contig-length 2500
done

## merge profile

anvi-merge  $ANVIO/profile/*2016*/PROFILE.db -o $ANVIO/profile_merged_2016 -c $ANVIO/contigs.db -S GRDI_2016 -W
anvi-merge  $ANVIO/profile/*2017*/PROFILE.db -o $ANVIO/profile_merged_2017 -c $ANVIO/contigs.db -S GRDI_2017 -W
anvi-merge  $ANVIO/profile/*2018*/PROFILE.db -o $ANVIO/profile_merged_2018 -c $ANVIO/contigs.db -S GRDI_2018 -W

cp $ANVIO/profile_merged_2016/PROFILE.db $ANVIO/2016_PROFILE.db
cp $ANVIO/profile_merged_2017/PROFILE.db $ANVIO/2017_PROFILE.db
cp $ANVIO/profile_merged_2018/PROFILE.db $ANVIO/2018_PROFILE.db

####################################################################
### BINNING! Finally
### Concoct
####################################################################

anvi-cluster-contigs -p $ANVIO/2016_PROFILE.db -c $ANVIO/contigs.db -C concoct_2016 --driver concoct -T 50 --just-do-it
anvi-cluster-contigs -p $ANVIO/2017_PROFILE.db -c $ANVIO/contigs.db -C concoct_2017 --driver concoct -T 50 --just-do-it
anvi-cluster-contigs -p $ANVIO/2018_PROFILE.db -c $ANVIO/contigs.db -C concoct_2018 --driver concoct -T 50 --just-do-it


####################################################################
### metabat2
####################################################################


anvi-cluster-contigs -p $ANVIO/2016_PROFILE.db -c $ANVIO/contigs.db -C metabat2_2016 --driver metabat2 -T 50 --just-do-it
anvi-cluster-contigs -p $ANVIO/2017_PROFILE.db -c $ANVIO/contigs.db -C metabat2_2017 --driver metabat2 -T 50 --just-do-it
anvi-cluster-contigs -p $ANVIO/2018_PROFILE.db -c $ANVIO/contigs.db -C metabat2_2018 --driver metabat2 -T 50 --just-do-it

####################################################################
### maxbin2
####################################################################

anvi-cluster-contigs -p $ANVIO/2016_PROFILE.db -c $ANVIO/contigs.db -C maxbin2_2016 --driver maxbin2 -T 50 --just-do-it
anvi-cluster-contigs -p $ANVIO/2017_PROFILE.db -c $ANVIO/contigs.db -C maxbin2_2017 --driver maxbin2 -T 50 --just-do-it
anvi-cluster-contigs -p $ANVIO/2018_PROFILE.db -c $ANVIO/contigs.db -C maxbin2_2018 --driver maxbin2 -T 50 --just-do-it

####################################################################
### binsanity
####################################################################

anvi-cluster-contigs -p $ANVIO/2016_PROFILE.db -c $ANVIO/contigs.db -C binsanity_2016 --driver binsanity -T 50 --just-do-it
anvi-cluster-contigs -p $ANVIO/2017_PROFILE.db -c $ANVIO/contigs.db -C binsanity_2017 --driver binsanity -T 50 --just-do-it
anvi-cluster-contigs -p $ANVIO/2018_PROFILE.db -c $ANVIO/contigs.db -C binsanity_2018 --driver binsanity -T 50 --just-do-it


####################################################################
### dastool
####################################################################


anvi-cluster-contigs -p $ANVIO/2016_PROFILE.db -c $ANVIO/contigs.db -S concoct_2016,metabat2_2016,maxbin2_2016,binsanity_2016 --search-engine usearch -C dastool_2016 -T 50 --just-do-it
anvi-cluster-contigs -p $ANVIO/2017_PROFILE.db -c $ANVIO/contigs.db -S concoct_2017,metabat2_2017,maxbin2_2017,binsanity_2017 --search-engine usearch -C dastool_2017 -T 50 --just-do-it
anvi-cluster-contigs -p $ANVIO/2018_PROFILE.db -c $ANVIO/contigs.db -S concoct_2018,metabat2_2018,maxbin2_2018,binsanity_2018 --search-engine usearch -C dastool_2018 -T 50 --just-do-it


