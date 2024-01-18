#!/bin/bash


# USAGE: ./1_align_BSseq.sh <fastq file R1> <fastq file R2>

# MANDATORY: Provide Path to "/Antisense_paper/Align_methylation/"

# This script accepts fastq files of R1 and R2 and performs the following steps on it:
    ## Trimming of fastq files with trim_galore (version 0.6.4)
    ## Mapping with BSMAPz (based on version 2.90, https://github.com/zyndagj/BSMAPz) to mm10
    ## Sorting, mapping quality filtering and indexing of bam files using samtools (version 1.10)
    ## Removal of duplicate reads using picard (version 2.7.1)



#Provide path to folder containing scripts and data [MANDATORY]
path="PATH/TO/Antisense_paper/Align_methylation/"

# Get input file and locations  
fq1=$1
fq2=$2

# Grab base of filename for future naming
base="$(basename $fq1 _R1.fastq.gz)"

# Make directories to store output
mkdir -p temp_trimmed_fastqs			
mkdir -p FastQC_trimmed	
mkdir -p Output
mkdir -p Output/bams

# Set up output filenames and locations
mm10BL=${path}"files/mm10-blacklist.v2.bed"
genome=${path}"files/mm10.fa"
trimmed_fq1="temp_trimmed_fastqs/"$base"_R1_val_1.fq.gz"
trimmed_fq2="temp_trimmed_fastqs/"$base"_R2_val_2.fq.gz"
bam=Output/bams/$base.bam
bam_sort=Output/bams/$base.sorted.bam
bam_filt=Output/bams/$base.sorted.q10.bam
bam_dup=Output/bams/$base.sorted.q10.noDup.bam
metric=Output/bams/$base.picard.metrics.txt









echo "Starting analysis of" $base

# Trim fastq files with trim_galore (version 0.6.4)
	## performs quality and adaptor trimming (default settings) 
	## trims 10 nucleotides from 5' end and 5 nucleotides from 3' end of R1 (--clip_R1 10 --three_prime_clip_R1 5)
	## trims 15 nucleotides from 5' end and 5 nucleotides from 3' end of R2 (--clip_R2 15 --three_prime_clip_R2 5)
	## too short read pairs are discarded (--paired)
	## runs FastQC (version 0.11.9) on trimmed fastq files and stores results in 'QualityControl/FastQC_trimmed' (--fastqc_args "--outdir QualityControl/FastQC_trimmed/")
	## writes trimmed fastq files to 'temp_trimmed_fastqs' (files will be removed later on) (-o temp_trimmed_fastqs/)
	## writes trimming report to 'temp_trimmed_fastqs'
	
trim_galore --fastqc_args "--outdir FastQC_trimmed/" -o temp_trimmed_fastqs/ --clip_R1 10 --three_prime_clip_R1 5 --clip_R2 15 --three_prime_clip_R2 5 --paired --basename $base $fq1 $fq2






# Map trimmed reads using BSMAPz (based on version 2.90, https://github.com/zyndagj/BSMAPz)
	## maps reads to mm10 (-d $scriptpath/Genomes/mm10.fa)
	## quality threshold for trimming set to 20 (-q 20)
	## reports unmapped reads (-u)
	## maximum number of equal best hits to count (-w 100)
	## output sam file converted to bam file using samtools (version 1.10)

bsmapz -a $trimmed_fq1 -b $trimmed_fq2 -d $genome -q 20 -u -w 100 -o $bam

# Clean_up
rm $trimmed_fq1
rm $trimmed_fq2






# Sort, filter (mapping quality <10) index bam file using Samtools (version 1.10)
# Remove duplicated using picard (version 2.7.1)
samtools sort $bam -o $bam_sort
samtools view -q 10 -b $bam_sort > $bam_filt

java -jar picard.jar MarkDuplicates \
	I=$bam_filt \
	O=$bam_dup \
	M=$metric \
	ASSUME_SORTED=TRUE \
	REMOVE_DUPLICATES=TRUE

samtools index $bam_dup

# Clean-up
rm $bam
rm $bam_sort
rm $bam_filt

	
		







