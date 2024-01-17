#!/usr/bin/env python3
import numpy as np
import pandas as pd
import sys
import os
import re

sample_list = ['XX_d0', 'XX_d2', 'XX_d4']
data_dir = sys.argv[1]
chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX']


#Read in count table for each sample
os.chdir(data_dir)
for sample in sample_list:
    #Reads count table into pandas dataframe
    counts = pd.read_csv('./' + sample + '_binned_counts_sorted.txt', sep='\t', names = ['chr', 'start', 'end', 'plus', 'minus'],
    dtype = {'chr' : 'str', 'start' : 'int', 'end' : 'int', 'plus' : 'int', 'minus' : 'int'})

    #Removes weird chromosomes and sorts by location
    counts = counts[counts['chr'].isin(chr_list)]


    #Creates columns for TT-seq foldchange (plus + minus)
    counts['fc_plus'] = (pd.to_numeric(counts.plus.shift(-1)) + pd.to_numeric(counts.plus.shift(-2)) + 1) / (pd.to_numeric(counts.plus.shift(1)) + pd.to_numeric(counts.plus.shift(2)) + 1)
    counts['fc_minus'] = (pd.to_numeric(counts.minus.shift(1)) + pd.to_numeric(counts.minus.shift(2)) + 1) / (pd.to_numeric(counts.minus.shift(periods = -1)) + pd.to_numeric(counts.minus.shift(periods = -2)) + 1)

    #Samples potential TSSs for plus and minus strands
    plus_tss_df = counts[pd.to_numeric(counts['fc_plus']) > 10]
    minus_tss_df = counts[pd.to_numeric(counts['fc_minus']) > 10]

    #Writes potential TSSs to BED files
    plus_tss_df.to_csv(data_dir + sample + '_pot_tss_plus.bed', sep = '\t', columns = ['chr', 'start', 'end'], index = False, header = False)
    minus_tss_df.to_csv(data_dir + sample + '_pot_tss_minus.bed', sep = '\t', columns = ['chr', 'start', 'end'], index = False, header = False)
