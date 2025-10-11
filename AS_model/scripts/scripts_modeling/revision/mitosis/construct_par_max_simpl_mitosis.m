clear all;
addpath('../../');

sets_max_simpl = dlmread("../../../../simulations/maximally_simplified_model/FST_max_simpl.txt");

fsw1 = 31;
fsw2 = 34;

mem_long=400;
stable_sets = find(min([sets_max_simpl(:,fsw1),sets_max_simpl(:,fsw2)],[],2)>mem_long); 
selected = sets_max_simpl(stable_sets,:);
selected(:,17:18) = 100000; %avoid basal promoter switches
selected(:,15:16) = 1/100000; %avoid basal OFF switches
selected(:,20) = selected(:,24);


%dlmwrite('par_as_model_max_simpl_mitosis.txt', selected);
