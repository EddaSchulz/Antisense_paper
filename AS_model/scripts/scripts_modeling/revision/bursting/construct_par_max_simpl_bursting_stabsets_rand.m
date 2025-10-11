clear all;
addpath('../../');

sets_max_simpl = dlmread("../../../../simulations/maximally_simplified_model/FST_max_simpl.txt");

fsw1 = 31;
fsw2 = 34;
mem_long=400;
stable_sets = find(min([sets_max_simpl(:,fsw1),sets_max_simpl(:,fsw2)],[],2)>mem_long); 
selected = sets_max_simpl(stable_sets,:);
selected(:,20) = selected(:,24);

%combine each stable set with different bursting kinetics
k15_range = [0.01,100];
k16_range = [0.01,100]; % tOFF
k17_range = [0.01,100];
k18_range = [0.01,100]; % tON

range_min = [log10(k15_range(1)),log10(k16_range(1)),log10(k17_range(1)),log10(k18_range(1))];
range_max = [log10(k15_range(2)),log10(k16_range(2)),log10(k17_range(2)),log10(k18_range(2))];

par_all = [];
for i = 1:size(selected,1)
    pars_burst = lhsu(range_min,range_max,200);
    pars = ones(size(pars_burst,1),30);
    pars(:,15:18) = 10.^pars_burst;
    pars(:,[1:14,19:30]) = repmat(selected(i,[1:14,19:30]),size(pars_burst,1),1);
    par_all = [par_all;pars];
end

%dlmwrite('par_as_model_max_simpl_bursting_stabsets_rand.txt', par_all)