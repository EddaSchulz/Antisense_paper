#!/usr/bin/env python3
import numpy as np
import pandas as pd
import sys
import os
import re

sample_list = ['d0_XX', 'd2_XX', 'd4_XX']
chrom_dir = sys.argv[1]
output_dir = sys.argv[2]

os.chdir(chrom_dir)

#The CutnTag ChromHMM data is similarly sorted by the total signal and then classified as enhancer or promoters based on the total signal
cnt_emissions = pd.read_csv('./emissions_6.txt', sep='\t')
cnt_emissions.columns = ['old_state', 'H3K4me3', 'H3K27ac', 'H3K4me1', 'ATAC']
cnt_emissions['total_signal'] = cnt_emissions.drop('old_state', axis = 1).sum(axis = 1)
cnt_emissions = cnt_emissions.sort_values(by = 'total_signal', ascending = False)
cnt_prom = cnt_emissions.iloc[0,0]
cnt_enh = cnt_emissions.iloc[1,0]

#Iterates over CnT file and standardizes color/states
for sample in sample_list:
    input = open('./' + sample + '_6_dense.bed', 'r')
    output = open(output_dir + sample + '_chromHMM.bed', 'w')
    for line in input:
    	if re.match('track', line):
    		output.write(line)
    	else:
    		tab = re.split('\t', line)
    		state = tab[3]
    		color = tab[8]
    		if state == str(cnt_prom):
    			color = '000,000,255'
    			state = 1
    		elif state == str(cnt_enh):
    			color = '000,000,125'
    			state = 2
    		else:
    			color = '255,255,255'
    			state = 3
    		output.write(tab[0] + '\t' + tab[1] + '\t' + tab[2] + '\t' + str(state) + '\t' + tab[4] + '\t' + tab[5] + '\t' + tab[6] + '\t' + tab[7] + '\t' + color + '\n')
