clear all;
addpath('../../');

k1_range = [5 500];
k2_range = [5 500];
k5_range = [1/3600 10];
k6_range = [1/3600 10];
p26_range = [5, 500];
p21_range = [0, 1];
p24_range = [0, 1];

range_min = [log10(k1_range(1)), log10(k2_range(1)), 1, 1, log10(k5_range(1)), log10(k6_range(1)),0,0,0,0, 100000, 100000, 0, 0, 0, 0, 0, 0, 0, 0,p21_range(1),0,0,p24_range(1),1,p26_range(1),0,0,0,0];
range_max = [log10(k1_range(2)), log10(k2_range(2)), 1, 1, log10(k5_range(2)), log10(k6_range(2)),0,0,0,0, 100000, 100000, 0, 0, 0, 0, 0, 0, 0, 0,p21_range(2),0,0,p24_range(2),1,p26_range(2),0,0,0,0];

%% 1 OFF state & 1-step initiation  & No SDI

k1_range = [5 500];
k2_range = [5 500];

pars = lhsu(range_min,range_max,119009);
pars(:,3:4) = 3.65; %dGFP mRNA halflife estimates from Tigges et al., 2009
pars(:,7:8) = 33.35;
pars(:,9:10) = 4.65; %dGFP protein halflife estimates from Tigges et al., 2009
pars(:,17:20) = 0;
pars(:,[1:2,5:6])=10.^(pars(:,[1:2,5:6]));
pars(:,11:12)=100000;
pars(:,13:14)=0;
pars(:,22)=1;
pars(:,23)=0;
pars(:,26)= floor(pars(:,26));
%dlmwrite('par_as_model_max_simpl.txt', pars);

