clear all;
%Enter path to save summary files in
path = '../../../../simulations/allrep/';
%Enter path containing function function_avg_alleles_200218.m
addpath('../../');
nr_p = 30;
nr_alleles = 100;

% enter directory containing raw simulation output files
%filenames=dir('xxxx');


for j = 0:length(filenames)/10
		%Name of summary file
		outputname = [path, 'summary_1lperset_general_as_model_200319_', num2str(j),'.txt'];
	for i = j*10+0:j*10+9
		%Name of raw simulation output file
		filename = [path,'sim_general_as_model_allrep_200319_',num2str(i),'.txt'];
        sim = dlmread(filename);
		t = (size(sim,2)-nr_p)/16;
		temp = function_avg_alleles_200218(sim, nr_p, t);
		dlmwrite(outputname, temp, '-append');
		clearvars -except j i nr_p nr_alleles path rep_lab outputname;
		i
	end
end
