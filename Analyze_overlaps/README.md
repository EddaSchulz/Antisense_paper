# Analyze Overlaps

## Description
This folder contains all code that was used to generate figure for the genome-wide antisense transcription analysis in Mutzel et al., 2025 

## Software dependencies and operating systems
In order to perform these analyses, the following software has to be installed and available on the command line ($PATH):
- Deeptools2 (v3.5.1) collection of PYTHON scripts from "https://deeptools.readthedocs.io/en/develop/"

R scripts can be run using R (v.4.2.1) software ("https://cran.r-project.org/"). The following R libraries are required:
- bedr (v1.0.7)
- bsseq (v1.32.0)
- clusterProfiler (v4.16.0)
- egg (v0.4.5)
- EnvStats (v2.8.1)
- genomation (v1.28.0)
- GenomicRanges (v1.48.0)
- ggsankey (v0.0.99999)
- gridExtra (v2.3)
- org.Mm.eg.db (v3.21.0)
- Rsubread (v2.10.5)
- tidyverse (v1.3.2)
- UpSetR (v1.4.0)
- viridis (v0.6.2)


## Reproduce analysis
The raw NGS data (GSE167358) should be previously processed using code from "/https://github.com/EddaSchulz/Xert_paper/". The script requires BAM files of TT-seq and CUT&Tag data (H3K27ac, H3K4me1, H3K9me3, H3K27me3, H2AK119ub, H3K36me3 and H3K4me3) for XXdXic cells at days 0, 2 and 4 of differentiation.
The BAM files should be of individual replicates and named according to the following patterns:
(H3K9me3|H3K27me3|H3K36me3|H3K27ac|H3K4me3|H3K4me1|H2AK119ub)_XX_d(0|2|4)_r(1|2).bam
TT_XX_d(0|2|4)_r(1|2).bam

Furthermore, the script requires BIGWIG files of the XX TT-seq data at days 0, 2 and 4 of differentiation (see "/https://github.com/EddaSchulz/Xert_paper/" for further instructions).
The BIGWIG files should be named by the following pattern:
XX_TT_d(0|2|4)_r(1|2).bigwig

"GENCODE_vM25_plus_Xert.gtf" should be retrieved from https://zenodo.org/records/12822424 and stored in "./files/".

While these scripts depends on the /Antisense_paper/Nascent_transcriptome/ section, all relevant output is already provided in /Antisense_paper/Analyze_overlaps/files/. The output of /Antisense_paper/Align_methylation/ should be either generated as described, or retrieved from GEO (GSE253792). The BEDGRAPH files of merged replicates should be stored in a separate directory and be named according to the pattern:
BSseq_XX_d(0|2|4)_CpG.bedGraph

All data analysis and plotting can be performed using the master script in the "/Analyze_overlaps/master/" directory. It is necessary to provide a path to "Antisense_paper/Analyze_overlaps/" (with -p), to a directory containing the BAM files (with -b), to a directory containing the TT-seq BIGWIG files (with -w) and to a directory containing the BSseq BEDGRAPH files (with -m).
The work directory can be specified with -d.


All scripts that are used by the master script are stored in "/Analyze_overlaps/scripts/". All files that are necessary to run the master script, other than the mentioned processed NGS data, are stored in "/Analyze_overlaps/files/". Sample output of the master script is provided in "/Analyze_overlaps/output/". 
