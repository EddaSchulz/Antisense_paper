#!/usr/bin/env python3

# Import packages
import numpy as np
import pandas as pd
import sys
import os


# Define input/output files and paths
day, input_dir, output_dir = sys.argv[1], sys.argv[2], sys.argv[3]

# Define functions to apply to df column below
concat = lambda x: ','.join(x)

# Perform steps for the plus and minus strand
for i in ['plus', 'minus']:

	# Define input/output files
	infile = input_dir + day + "_XX_" + i + "_full_annot.bed"
	outfile = output_dir + day + "_XX_" + i + "_transcripts_binned.bed"

	# Read input file
	bed_df = pd.read_csv(infile, sep='\t', header = None)
	bed_df.columns = ['chr', 'start', 'end', 'state', 'strand', 'gene', 'gene_type']


	# For the minus strand, df needs to be reversed as bins are back-to-front compared to plus
	if i == 'minus':
		bed_df = bed_df.iloc[::-1]


	########## Merge bins belonging to the same transcript by sequentially aggregating the dataframe ############

	# Merge rows of identical genomic coordinates but differing gene annotations (note: annotation state will be joined, separated by comma)



	bed_df = (bed_df.groupby((bed_df.start != bed_df.start.shift()).cumsum()).agg({'chr':'first', 'start':min, 'end':max, 'state':'first', 'strand':'first', 'gene':concat, 'gene_type':concat})).reset_index(drop=True)


	# Merge 'tss' bins to directly following transcribed bins (i.e. a 'plus/minus' bin) regardless of annotation (note: state will still be named 'tss' after merging)

#CHECK
	if i == 'plus':
		bed_df = (bed_df.groupby( (((bed_df.state.shift() != 'tss') | (bed_df.state != i) | (bed_df.chr != bed_df.chr.shift())) | (bed_df.start != bed_df.end.shift())).cumsum()).agg({'chr':'first', 'start':min, 'end':max, 'state':'first', 'strand':'first', 'gene':'first', 'gene_type':'first'})).reset_index(drop=True)
	else:
		bed_df = (bed_df.groupby( (((bed_df.state.shift() != 'tss') | (bed_df.state != i) | (bed_df.chr != bed_df.chr.shift())) | (bed_df.end != bed_df.start.shift())).cumsum()).agg({'chr':'first', 'start':min, 'end':max, 'state':'first', 'strand':'first', 'gene':'first', 'gene_type':'first'})).reset_index(drop=True)


	# Remove rows with no gene annotation (note: unannotated transcribed bins previously merged with preceeding tss will not be discarded as their annotation state looks like '.,.')
	bed_df = bed_df[((bed_df.gene != '.') & (bed_df.state == i)) | (bed_df.state == "tss")]

	# Merge consecutive rows which have the same gene annotation within potential transcript group (i.e. between 'tss' states)
	bed_df['tss_group'] = ((bed_df.state == 'tss').cumsum())	# add a column which groups potential bins belonging to the same transcript
	bed_df['gene_group'] = ((bed_df.gene != bed_df.gene.shift()).cumsum())	# add a column which groups bins with same gene annotation
	bed_df = (bed_df.groupby(['tss_group', 'gene_group'] ).agg({'chr':'first', 'start':min, 'end':max, 'state':'first', 'strand':'first', 'gene':'first', 'gene_type':'first'})).reset_index(drop=True)


	# The previous step did not merge transcribed bins with the same gene annotation if a bin has been annotated with multiple genes (i.e. 'gene' = gene1,gene2,etc)
	# Loop over potential transcript groups and merge bins where the gene annotation of the 'tss' bin matches the following bins
	grouped = bed_df.groupby(((bed_df.state == 'tss').cumsum()))
	output_list = []
	index_count = 0
	#print(grouped.head())
	for name, group in grouped:

		#print(group)
		chrom, gene_id, start, end, gene_type = [], [], [], [], []
		for index,row in group.iterrows():
			# If the 'tss' bin is not annotated then it won't be merged with following bins
			if row['state'] == 'tss' and set(row['gene'].split(',')) == {'.'}:
				merged_bin = pd.DataFrame({'chr': row['chr'], 'start': row['start'], 'end': row['end'], 'gene': '.', 'gene_type': '.'}, index=[index_count])
				output_list.append(merged_bin)
				break

			# If the 'tss' bin is annotated, then check if it can be merged with the next bins of the same annotation
			else:
				if row['state'] == 'tss':
					chrom.append(row['chr'])
					start.append(row['start'])
					end.append(row['end'])
					gene_id.append(set(row['gene'].split(',')))
					gene_type.append(set(row['gene_type'].split(',')))
					#gene_id[0].discard('.') # remove '.' from list as it stands for unannotated (since only annotated bins should be merged at this point)
					merged_bin = pd.DataFrame({'chr': chrom[0], 'start': start[0], 'end': end[0], 'gene': ','.join(gene_id[0]), 'gene_type': ','.join(gene_type[0])}, index=[index_count])
					output_list.append(merged_bin)

				else:
					if len(gene_id) == 0: # it's possible that the first group does not contain 'tss' in its first row, in that case stop loop for this group
						index_count -= 1
						break
					elif bool(set(row['gene'].split(',')) & gene_id[0]) == True:
						if i == 'plus':
							end[0]=row['end']
						else:
							start[0]=row['start']
						gene=set(row['gene'].split(',')).intersection(gene_id[0])
						merged_bin = pd.DataFrame({'chr': chrom[0], 'start': start[0], 'end': end[0], 'gene': ','.join(gene), 'gene_type': ','.join(gene_type[0])}, index=[index_count])
						output_list[index_count]=merged_bin
					else:
						break
		index_count += 1

	# Store the output as a dataframe
	output_df = pd.concat(output_list, axis=0)

	# For the minus strand, df needs to be reversed back
	if i == 'minus':
		output_df = output_df.iloc[::-1]

	# Write the output to text file
	output_df.to_csv(outfile, sep='\t', index = False, header = False)
