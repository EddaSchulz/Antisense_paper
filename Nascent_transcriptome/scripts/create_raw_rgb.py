#!/usr/bin/env python3
import numpy as np
import pandas as pd
import sys
import os
import re

sample_list = ['d0_XX', 'd2_XX', 'd4_XX']
input_dir = sys.argv[1]
output_dir = sys.argv[2]
chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX']


#Read in raw tracks and add color info for UCSC
os.chdir(input_dir)
for sample in sample_list:
    #Reads raw bed into pandas dataframe and creates rgb tracks
    plus_bed = pd.read_csv('./' + sample + '_plus_raw.bed', sep='\t', names = ['chr', 'start', 'end', 'type'])
    plus_bed['score'] = 0
    plus_bed['strand'] = '+'
    plus_bed['thickStart'] = plus_bed['start']
    plus_bed['thickEnd'] = plus_bed['end']
    plus_bed['rgb'] = np.select([plus_bed['type'] == 'tss'],['255,000,000'],'250,000,200')

    minus_bed = pd.read_csv('./' + sample + '_minus_raw.bed', sep='\t', names = ['chr', 'start', 'end', 'type'])
    minus_bed['score'] = 0
    minus_bed['strand'] = '-'
    minus_bed['thickStart'] = minus_bed['start']
    minus_bed['thickEnd'] = minus_bed['end']
    minus_bed['rgb'] = np.select([minus_bed['type'] == 'tss'],['000,000,255'],'200,000,250')

    #Writes rgb tracks to BED files
    plus_bed.to_csv(output_dir + sample + '_raw_rgb_plus.bed', sep = '\t', index = False, header = False)
    minus_bed.to_csv(output_dir + sample + '_raw_rgb_minus.bed', sep = '\t', index = False, header = False)
