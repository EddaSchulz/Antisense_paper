
addpath('../../');

k1_range = [10 1000];
k2_range = [10 1000];
k5_range = [1/3600 10];
k6_range = [1/3600 10];
k11_range = [10 1000];
k12_range = [10 1000];
k13_range = [0.05 0.9];
k14_range = [0.05 0.9];

p25_range = [1, 6];
p26_range = [5, 500];

range_min = [log10(k1_range(1)), log10(k2_range(1)), 1, 1, log10(k5_range(1)), log10(k6_range(1)),0,0,0,0, log10(k11_range(1)), log10(k12_range(1)), log10(k13_range(1)), log10(k14_range(1)), 0, 0, 0, 0, 0, 0,0,0,0,0,p25_range(1),p26_range(1),0,0,0,0];
range_max = [log10(k1_range(2)), log10(k2_range(2)), 1, 1, log10(k5_range(2)), log10(k6_range(2)),0,0,0,0, log10(k11_range(2)), log10(k12_range(2)), log10(k13_range(2)), log10(k14_range(2)), 0, 0, 0, 0, 0, 0,0,0,0,0,p25_range(2),p26_range(2),0,0,0,0];

%% 1 OFF state & 1-step initiation 

k1_range = [5 500];
k2_range = [5 500];
range_min = [log10(k1_range(1)), log10(k2_range(1)), 1, 1, log10(k5_range(1)), log10(k6_range(1)),0,0,0,0, log10(k11_range(1)), log10(k12_range(1)), log10(k13_range(1)), log10(k14_range(1)), 0, 0, 0, 0, 0, 0,0,0,0,0,p25_range(1),p26_range(1),0,0,0,0];
range_max = [log10(k1_range(2)), log10(k2_range(2)), 1, 1, log10(k5_range(2)), log10(k6_range(2)),0,0,0,0, log10(k11_range(2)), log10(k12_range(2)), log10(k13_range(2)), log10(k14_range(2)), 0, 0, 0, 0, 0, 0,0,0,0,0,p25_range(2),p26_range(2),0,0,0,0];

pars = lhsu(range_min,range_max, 35700);
pars(:,3:4) = 3.65; %dGFP mRNA halflife estimates from Tigges et al., 2009
pars(:,7:8) = 33.35;
pars(:,9:10) = 4.65; %dGFP protein halflife estimates from Tigges et al., 2009
pars(:,17:20) = 0;
pars(:,[1:2,5:6,11:14])=10.^(pars(:,[1:2,5:6,11:14]));
pars(:,11:12)=100000;
pars(:,13:14)=0;
pars(:,21) = 0.5;
pars(:,22:24)=1;
pars(:,25)= 1;
pars(:,26)= floor(pars(:,26));
%dlmwrite('par_1OFF_1step_ini.txt', pars);

%% 1 Step-Ini & No SDI

pars = lhsu(range_min,range_max,35700);
pars(:,3:4) = 3.65; %dGFP mRNA halflife estimates from Tigges et al., 2009
pars(:,7:8) = 33.35;
pars(:,9:10) = 4.65; %dGFP protein halflife estimates from Tigges et al., 2009
pars(:,17:20) = 0;
pars(:,[1:2,5:6,11:14])=10.^(pars(:,[1:2,5:6,11:14]));
pars(:,11:12)=100000;
pars(:,13:14)=0;
pars(:,21) = 0.5;
pars(:,22:24)=1;
pars(:,23)=0;
pars(:,25)= floor(pars(:,25));
pars(:,26)= floor(pars(:,26));
%dlmwrite('par_1step_ini_no_SDI.txt', pars);

k1_range = [10 1000];
k2_range = [10 1000];
range_min = [log10(k1_range(1)), log10(k2_range(1)), 1, 1, log10(k5_range(1)), log10(k6_range(1)),0,0,0,0, log10(k11_range(1)), log10(k12_range(1)), log10(k13_range(1)), log10(k14_range(1)), 0, 0, 0, 0, 0, 0,0,0,0,0,p25_range(1),p26_range(1),0,0,0,0];
range_max = [log10(k1_range(2)), log10(k2_range(2)), 1, 1, log10(k5_range(2)), log10(k6_range(2)),0,0,0,0, log10(k11_range(2)), log10(k12_range(2)), log10(k13_range(2)), log10(k14_range(2)), 0, 0, 0, 0, 0, 0,0,0,0,0,p25_range(2),p26_range(2),0,0,0,0];

%% 1 OFF state & No SDI
pars = lhsu(range_min,range_max,35700);
pars(:,3:4) = 3.65; %dGFP mRNA halflife estimates from Tigges et al., 2009
pars(:,7:8) = 33.35;
pars(:,9:10) = 4.65; %dGFP protein halflife estimates from Tigges et al., 2009
pars(:,17:20) = 0;
pars(:,[1:2,5:6,11:14])=10.^(pars(:,[1:2,5:6,11:14]));
pars(:,21) = 0.5;
pars(:,22:24)=1;
pars(:,23)=0;
pars(:,25)= 1;
pars(:,26)= floor(pars(:,26));
%dlmwrite('par_1OFF_no_SDI.txt', pars);


