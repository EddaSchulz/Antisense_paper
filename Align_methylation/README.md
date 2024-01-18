# Align methylation

## Description
This folder contains all code that was used to process BSseq data for Mutzel et al., 2024. Output of the script is provided at GSEXXX.

## Software dependencies and operating systems
In order to perform these analyses, the following software has to be installed and available on the command line ($PATH):
- trim_galore (v0.6.4) perl wrapper around FastQC and Cutadapt from "https://github.com/FelixKrueger/TrimGalore"
- BSMAPz (v2.90) short read aligner from "https://github.com/zyndagj/BSMAPz"
- bedtools (v2.29.2) collection of C++ scripts from "https://bedtools.readthedocs.io/en/latest/"
- samtools (v1.10) collection of C scripts from "http://www.htslib.org/"
- Picard (v2.7.1) collection of JAVA tools from "https://broadinstitute.github.io/picard/"
- MethyDackel (v0.3.0) collection of C scripts from "https://github.com/dpryan79/MethylDackel"
- bedGraphToBigWig (v4) file format converter from "https://www.encodeproject.org/software/bedgraphtobigwig/"


## Reproduce analysis
The raw FASTQ files should be retrieved from GEO using "fasterq_dump" (GSEXXX) and renamed to the pattern "BSseq_XX_d(0|2|4)_r(1|2)_R(1|2).fastq.gz" (with "d" = day, "r" = replicate and "R" = read).
The mm10 genome should be downloaded as a FASTA file from "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/" and stored as "Antisense_paper/Align_methylation/files/mm10.fa".
 
The scripts should be run in the indicated order. The output should not be renamed. Both scripts require that the path to "Antisense_paper/Align_methylation/" is added in the script itself.
1_align_BSseq.sh needs to be run for each pair of FASTQ files individually (or started in a loop across all samples). 2_mergeReplicates.sh does not require additional info outside the path. 

- 1_align_BSseq.sh: This script aligns individual pairs of FASTQ files to the mm10 genome and produces BAM files. (USAGE: ./1_align_BSseq.sh $fq_read1 $fq_read2)

- 2_mergeReplicates.sh: This script merges the BAM files produced in 1_align_BSseq.sh and generates BEDGRAPH and BIGWIG files of the CpG information. (USAGE: ./2_mergeReplicates.sh)
