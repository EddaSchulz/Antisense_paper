#!/bin/bash



# USAGE: ./2_mergeReplicates.sh

# This script merges replicate bam file (generates using 1_align_BSseq.sh) for each sample in a for loop and generates bedGraph and bigwig files containing the per-base CpG metrics for merged replicates:
    ## Merging of replicate bam files using samtools (version 1.10)
    ## Generation of bedgraph file containing the per-base CpG metrics using MethyDackel 'extract' (Version: 0.3.0)
    ## Conversion of bedgraph files to bigwigs using bedGraphToBigWig (version 4)



# Define paths and files [MANDATORY]
path="/project/agsgpa/Till/verena_antisense/github/Align_methylation/"
bampath=$path'Output/bams/'
mkdir -p $path'Output/bedgraphs/'
bedgraphpath=$path'Output/bedgraphs/'
mkdir -p $path'Output/bigwigs/'
bigwigpath=$path'Output/bigwigs/'

mm10BL=${path}"files/mm10-blacklist.v2.bed"
genome=${path}"files/mm10.fa"

bams=$bampath'*R1.sorted.q10.noDup.bam'



# Loop over all samples

for bam in $bams

	do
	
	
	
	base=$(basename $bam _r1.sorted.q10.noDup.bam)
	bam2=$bampath$base'_r2.sorted.q10.noDup.bam'
	merge=$bampath$base'.sorted.q10.noDup.bam'


	
	# Merge replicate bam files using Samtools (version 1.10)
	samtools merge -o $merge $bam $bam2
	samtools index $merge
	
	
	
	# Generate bedgraph files containing the per-base CpG metrics using MethyDackel 'extract' (Version: 0.3.0)
	MethylDackel extract $genome $merge -@ 8 --mergeContext --minDepth 5 -o $bedgraphpath\$base/_CpG_pre



	# Filter out CpGs on mitochondrial DNA (chrM) and unassigned contigs (chrUn) and remove blacklisted regions from begraph files using betools intersect (Version: v2.29.2)
	grep -v chrM $bedgraphpath$base'_CpG_pre.bedGraph' | grep -v _ - | bedtools intersect -a - -b $mm10BL -v > $bedgraphpath$base'_CpG.bedGraph'


	rm $bedgraphpath$base'.mergeContext.minDepth2_CpG.bedGraph'
	
	

	# Convert to bedgraph files to bigwig using bedGraphToBigWig (version 4)
	tail -n +2 $bedgraphpath$base'_CpG.bedGraph' | cut -f 1,2,3,4 - > $bedgraphpath$base'_CpG_edit.bedGraph'			
	bedGraphToBigWig $bedgraphpath$base'_CpG_edit.bedGraph' $scriptpath/Genomes/chrom.sizes.mm10 $bigwigpath$base'_CpG.bw'
	rm $bedgraphpath$base'_CpG_edit.bedGraph'
	
	
	
done
	
	

