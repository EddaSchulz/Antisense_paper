# Nascent transcriptome generation

## Description
This folder contains all code that was used to generate a nascent transcriptome for Mutzel et al., 2024 

## Software dependencies and operating systems
In order to perform these analyses, the following software has to be installed and available on the command line ($PATH):
- bedtools (v2.29.2) collection of C++ scripts from "https://bedtools.readthedocs.io/en/latest/"
- ChromHMM (v1.19) collection of JAVA scripts from "http://compbio.mit.edu/ChromHMM/"
- Deeptools2 (v3.5.1) collection of PYTHON scripts from "https://deeptools.readthedocs.io/en/develop/"
- Samtools (v1.15.1) collection of C scripts from "http://www.htslib.org/"
- Python3 (v3.8.5) programming language from "https://www.python.org/"

The following Python scripts are required:
- Numpy (v1.22.2)
- Pandas (v1.4.0)

R scripts can be run using R (v.4.2.1) software ("https://cran.r-project.org/"). The following R libraries are required:
- egg (v0.4.5)
- GenomicRanges (v1.48.0)
- gplots (v3.1.0)
- gridExtra (v2.3)
- STAN (v2.24.0)
- tidyverse (v1.3.2)


## Reproduce analysis
The raw NGS data (GSE167358) should be previously processed using code from "/https://github.com/EddaSchulz/Xert_paper/". The script requires TT-seq, ATAC-seq and CUT&Tag data (H3K27ac, H3K4me1 and H3K4me3) for XXdXic cells at days 0, 2 and 4 of differentiation.
The TT-seq data should be split by strand. The BAM files should be merged (e.g. using samtools) and named according to the following patterns:
(ATAC|H3K27ac|H3K4me3|H3K4me1)_XX_d(0|2|4).bam
TT_XX_d(0|2|4)_(plus|minus).bam
 
Data analysis and plotting can be performed using the master script in the "/Nascent_transcriptome/master/" directory. It is necessary to provide a path to "AS_paper/Nascent_transcriptome/" (with -p) and to a directory containing the BAM files (with -b).
The work directory can be specified with -d.


All scripts that are used by the master script are stored in "/Nascent_transcriptome/scripts/". All files that are necessary to run the master script, other than the processed sequencing data, are stored in "/Nascent_transcriptome/files/". Sample output of the master scripts is provided in "/Nascent_transcriptome/output/". 
