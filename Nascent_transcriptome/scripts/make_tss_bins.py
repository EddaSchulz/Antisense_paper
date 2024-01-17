#!/usr/bin/env python3

# Import packages
import numpy as np
import pandas as pd
import sys
import os


# Define input/output files and paths
day, input_dir, output_dir = sys.argv[1], sys.argv[2], sys.argv[3]


# Perform steps for the plus and minus strand
for i in ['plus', 'minus']:

	# Define input/output files
	infile = input_dir + day + "_XX_" + i + "_transcripts_binned.bed"
	outfile_transList = output_dir + "XX_" + day + "_" + i + "_transcripts_list.bed"
	outfile_binTSS = output_dir + "XX_" + day + "_" + i + "_tss_bins.bed"
	outfile_binTTS = output_dir + "XX_" + day + "_" + i + "_tts_bins.bed"

	# Read input file
	bed_df = pd.read_csv(infile, sep='\t', header = None)
	bed_df.columns = ['chr', 'start', 'end', 'gene', 'gene_type']

	# Add extra column containing transcript infos
	bed_df['Name'] = ("TRANS_XX" + day + "_" + bed_df['chr'] + "_" + i + "_NR" + (bed_df.groupby('chr').cumcount() + 1).astype(str))

	# Write the output to text file
	bed_df.to_csv(outfile_transList, sep='\t', index = False, header = False)

	#Calculate range for bins as 1/40 of total length
	bed_df["length"]=bed_df["end"] - bed_df["start"]
	bed_df["range"]=bed_df["length"] / 40

	# Bin the start and end of each transcript, i.e. create 10-bp-bins at start and end covering a quarter of the transcript each
	bed_df["binStart_start"] = bed_df.apply(lambda x: [x.start+(i*10) for i in range(int(x.range))], axis = 1)
	bed_df["binStart_end"] = bed_df.apply(lambda x: [x.start+(i*10+9) for i in range(int(x.range))], axis = 1)
	bed_df["binEnd_start"] = bed_df.apply(lambda x: [x.end+(i*10-int(x.range)*10) for i in range(int(x.range))], axis = 1)
	bed_df["binEnd_end"] = bed_df.apply(lambda x: [x.end+(i*10-(int(x.range)*10 - 9)) for i in range(int(x.range))], axis = 1)
	bed_df_binStart = bed_df[["chr", "binStart_start", "binStart_end"]].explode(['binStart_start','binStart_end'])
	bed_df_binEnd = bed_df[["chr", "binEnd_start", "binEnd_end"]].explode(['binEnd_start','binEnd_end'])

	# Write the output to text file
	if i == 'plus':
		bed_df_binStart.to_csv(outfile_binTSS, sep='\t', index = False, header = False)
		bed_df_binEnd.to_csv(outfile_binTTS, sep='\t', index = False, header = False)
	else:
		bed_df_binStart.to_csv(outfile_binTTS, sep='\t', index = False, header = False)
		bed_df_binEnd.to_csv(outfile_binTSS, sep='\t', index = False, header = False)
