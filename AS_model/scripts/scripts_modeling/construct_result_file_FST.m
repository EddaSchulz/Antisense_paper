clear all; 
%% All_Rep
mypath = '../../../../simulations/allrep/'; %Path to save FST file in
filenames=dir('../../../../simulations/allrep/summary_1lperset_general_as_model_200319_*.txt'); % Path of summary files previously constructed from raw simulation outputs using script summarize_sim.m
outputname = [mypath, 'FST_general_as_model_allrep_200319.txt'];


nr_p = 30;
nbins = 100;
v = 1/1440;
tmax = 500;

%Indices of simulation output
pa1_count = nr_p+1:nr_p+nbins;
pa1_center = nr_p+nbins+1:nr_p+2*nbins;
pb1_count = nr_p+2*nbins+1:nr_p+3*nbins;
pb1_center = nr_p+3*nbins+1:nr_p+4*nbins;
ra1_count = nr_p+4*nbins+1:nr_p+5*nbins;
ra1_center = nr_p+5*nbins+1:nr_p+6*nbins;
rb1_count = nr_p+6*nbins+1:nr_p+7*nbins;
rb1_center = nr_p+7*nbins+1:nr_p+8*nbins;
pra1_count = nr_p+8*nbins+1:nr_p+9*nbins;
pra1_center = nr_p+9*nbins+1:nr_p+10*nbins;
prb1_count = nr_p+10*nbins+1:nr_p+11*nbins;
prb1_center = nr_p+11*nbins+1:nr_p+12*nbins;

pa2_count = nr_p+12*nbins+1:nr_p+13*nbins;
pa2_center = nr_p+13*nbins+1:nr_p+14*nbins;
pb2_count = nr_p+14*nbins+1:nr_p+15*nbins;
pb2_center = nr_p+15*nbins+1:nr_p+16*nbins;
ra2_count = nr_p+16*nbins+1:nr_p+17*nbins;
ra2_center = nr_p+17*nbins+1:nr_p+18*nbins;
rb2_count = nr_p+18*nbins+1:nr_p+19*nbins;
rb2_center = nr_p+19*nbins+1:nr_p+20*nbins;
pra2_count = nr_p+20*nbins+1:nr_p+21*nbins;
pra2_center = nr_p+21*nbins+1:nr_p+22*nbins;
prb2_count = nr_p+22*nbins+1:nr_p+23*nbins;
prb2_center = nr_p+23*nbins+1:nr_p+24*nbins;

avg_fsw1 = nr_p+24*nbins+6*(tmax+1)+1:nr_p+24*nbins+6*(tmax+1)+3;
avg_fsw2 = nr_p+24*nbins+6*(tmax+1)+4:nr_p+24*nbins+6*(tmax+1)+6;
avg_sw1_mean = nr_p+24*nbins+6*(tmax+1)+7:nr_p+24*nbins+6*(tmax+1)+9;
avg_sw2_mean = nr_p+24*nbins+6*(tmax+1)+10:nr_p+24*nbins+6*(tmax+1)+12;
avg_sw1_std = nr_p+24*nbins+6*(tmax+1)+13:nr_p+24*nbins+6*(tmax+1)+15;
avg_sw2_std = nr_p+24*nbins+6*(tmax+1)+16:nr_p+24*nbins+6*(tmax+1)+18;
avg_sw1_min = nr_p+24*nbins+6*(tmax+1)+19:nr_p+24*nbins+6*(tmax+1)+21;
avg_sw2_min = nr_p+24*nbins+6*(tmax+1)+22:nr_p+24*nbins+6*(tmax+1)+24;
avg_sw1_max = nr_p+24*nbins+6*(tmax+1)+25:nr_p+24*nbins+6*(tmax+1)+27;
avg_sw2_max = nr_p+24*nbins+6*(tmax+1)+28:nr_p+24*nbins+6*(tmax+1)+30;
avg_sw1_n = nr_p+24*nbins+6*(tmax+1)+31:nr_p+24*nbins+6*(tmax+1)+33;
avg_sw2_n = nr_p+24*nbins+6*(tmax+1)+34:nr_p+24*nbins+6*(tmax+1)+36;
mean_sw = nr_p+24*nbins+6*(tmax+1)+37:nr_p+24*nbins+6*(tmax+1)+40;

sim = [];
new_output = [];
for fi = 1:size(filenames,1)
    thisfile = [mypath, filenames(fi).name];
    sim_temp = dlmread(thisfile);
    %%Mean of RNAP, RNA and protein distributions over last 50h and all
    %alleles; 1:ini1 2:ini2
    mean_pa1 = sum(sim_temp(:,pa1_count).*sim_temp(:,pa1_center),2)./sum(sim_temp(:,pa1_count),2);
    mean_pa2 = sum(sim_temp(:,pa2_count).*sim_temp(:,pa2_center),2)./sum(sim_temp(:,pa2_count),2);
    mean_pb1 = sum(sim_temp(:,pb1_count).*sim_temp(:,pb1_center),2)./sum(sim_temp(:,pb1_count),2);
    mean_pb2 = sum(sim_temp(:,pb2_count).*sim_temp(:,pb2_center),2)./sum(sim_temp(:,pb2_count),2);
    mean_ra1 = sum(sim_temp(:,ra1_count).*sim_temp(:,ra1_center),2)./sum(sim_temp(:,ra1_count),2);
    mean_ra2 = sum(sim_temp(:,ra2_count).*sim_temp(:,ra2_center),2)./sum(sim_temp(:,ra2_count),2);
    mean_rb1 = sum(sim_temp(:,rb1_count).*sim_temp(:,rb1_center),2)./sum(sim_temp(:,rb1_count),2);
    mean_rb2 = sum(sim_temp(:,rb2_count).*sim_temp(:,rb2_center),2)./sum(sim_temp(:,rb2_count),2);
    mean_pra1 = sum(sim_temp(:,pra1_count).*sim_temp(:,pra1_center),2)./sum(sim_temp(:,pra1_count),2);
    mean_pra2 = sum(sim_temp(:,pra2_count).*sim_temp(:,pra2_center),2)./sum(sim_temp(:,pra2_count),2);
    mean_prb1 = sum(sim_temp(:,prb1_count).*sim_temp(:,prb1_center),2)./sum(sim_temp(:,prb1_count),2);
    mean_prb2 = sum(sim_temp(:,prb2_count).*sim_temp(:,prb2_center),2)./sum(sim_temp(:,prb2_count),2);
    mean_out = [mean_pa1 mean_pa2 mean_pb1 mean_pb2 mean_ra1 mean_ra2 mean_rb1 mean_rb2 mean_pra1 mean_pra2 mean_prb1 mean_prb2];
    output_temp = [sim_temp(:,1:nr_p) sim_temp(:, avg_fsw1) sim_temp(:, avg_fsw2) mean_out];
    new_output = [new_output; output_temp];
end

dlmwrite(outputname, new_output);
