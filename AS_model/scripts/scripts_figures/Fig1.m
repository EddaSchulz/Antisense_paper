clear all; close all;
%% Plotting parameters
% Color scheme initial conditions
cm_a = [27 157 142; 148 195 186]./255;
cm_b = [235 109 30; 245 172 118]./255;
% Color scheme memory timescales
cm2 = [204 204 204; 150 150 150; 82 82 82]./255;
graph_size=[6-6/11 6-6/11];
graph_size1=[6/11 6-6/11];
graph_size2=[6-6/11 1.5];
graph_size3=[1.5 1.5];


lw=1;
fs=6;
fst=10;
pos_x=[2.7:4.2:19.5]; %pos_y=[2 6.5 11 15.5 19.5:3.8:30];
pos_x2=[3.7:4.2:19.5];
pos_y=[0:3.8:11.4]; pos_y=flip(pos_y);
pos = []; indx = 1; indy = 1;
for i = 1:length(pos_x).*length(pos_y)
    pos = [pos; pos_x(indx) pos_y(indy)];
    indx = indx+1;
    if indx>length(pos_x)
        indx = 1;
        indy = indy+1;
    end
end

%% Simulations
oripath = '../../simulations/allrep/';
path_1 = '../../simulations/simplified_models_1/';
path_2 = '../../simulations/simplified_models_2/';
path_3 = '../../simulations/simplified_models_3/';

filename{1} = ([oripath, 'fsw_general_as_model_allrep_200319.txt']);
filename{2} = ([path_1, 'fsw_simpl#1_1OFF_state_as_model_200319.txt']);
filename{3} = ([path_1, 'fsw_simpl#2_no_term_as_model_200319.txt']);
filename{4} = ([path_1, 'fsw_simpl#3_1step_ini_as_model_200319.txt']);
filename{5} = ([path_1, 'fsw_simpl#4_no_SDI_as_model_200319.txt']);
filename{6} = ([path_1, 'fsw_simpl#6_no_PR_as_model_200319.txt']);
filename{7} = ([path_1, 'fsw_simpl#7_no_PC_as_model_200319.txt']);
filename{8} = ([path_2, 'fsw_1OFF_1step_ini_as_model_200319.txt']);
filename{9} = ([path_2, 'fsw_1OFF_no_SDI_as_model_200319.txt']);
filename{10} = ([path_2, 'fsw_1step_ini_no_SDI_as_model_200319.txt']);
filename{11} = ([path_3, 'fsw_1OFF_1step_ini_no_SDI_as_model_200319.txt']);

ext_filename{1} = ([oripath, 'fsw_add_par_as_model_allrep_200319_220125.txt']);
ext_filename{2} = ([path_1, 'fsw_add_par_as_model_simpl#1_1OFF_state_200319_220125.txt']);
ext_filename{3} = ([path_1, 'fsw_add_par_as_model_simpl#2_no_term_200319_220125.txt']);
ext_filename{4} = ([path_1, 'fsw_add_par_as_model_simpl#3_1step_ini_200319_220125.txt']);
ext_filename{5} = ([path_1, 'fsw_add_par_as_model_simpl#4_no_SDI_200319_220125.txt']);
ext_filename{6} = ([path_1, 'fsw_add_par_as_model_simpl#6_no_PR_200319_220125.txt']);
ext_filename{7} = ([path_1, 'fsw_add_par_as_model_simpl#7_no_PC_200319_220125.txt']);
ext_filename{8} = ([path_2, 'fsw_add_par_as_model_simpl_1OFF_1step_ini_200319_220125.txt']);
ext_filename{9} = ([path_2, 'fsw_add_par_as_model_simpl_1OFF_no_SDI_200319_220125.txt']);
ext_filename{10} = ([path_2, 'fsw_add_par_as_model_simpl_1step_ini_no_SDI_200319_220125.txt']);
ext_filename{11} = ([path_3, 'fsw_add_par_as_model_simpl_1OFF_1step_ini_no_SDI_200319_220125.txt']);

mymodel{1} = 'Ori\_Model';
mymodel{2} = '1OFF';
mymodel{3} = 'No\_term';
mymodel{4} = '1step\_ini';
mymodel{5} = 'No\_SDI';
mymodel{6} = 'No\_PR';
mymodel{7} = 'No\_PC';
mymodel{8} = '1OFF\_1step\_ini';
mymodel{9} = '1OFF\_no\_SDI';
mymodel{10} = '1step\_ini\_no\_SDI';
mymodel{11} = '1OFF\_1step\_ini\_no\_SDI';

mymodelnames{1} = 'Ori_Model';
mymodelnames{2} = '1OFF';
mymodelnames{3} = 'No_term';
mymodelnames{4} = '1step_ini';
mymodelnames{5} = 'No_SDI';
mymodelnames{6} = 'No_PR';
mymodelnames{7} = 'No_PC';
mymodelnames{8} = '1OFF_1step_ini';
mymodelnames{9} = '1OFF_no_SDI';
mymodelnames{10} = '1step_ini_no_SDI';
mymodelnames{11} = '1OFF_1step_ini_no_SDI';

% RNAP, RNA or Protein level
mylevel{1} = 'RNAP';
mylevel{2} = 'RNA';
mylevel{3} = 'Protein';
j = 1; %Level to perform analysis on

%% Indices
nr_p = 30;
nbins = 10;
var_par = [1 2 5 6 11 12 13 14 25 26];
var_log = [1 2 5 6 11 12 13 14];
fsw1 = nr_p+1:nr_p+3;
fsw2 = nr_p+4:nr_p+6;
thresh_st = 96; %Stability threshold (before: 100)
thresh_unst = 2;% before 10
fr_stable = NaN(length(filename),3);
%% Read in simulations
for f=1:length(filename)
    temp1 = dlmread(filename{f}); 
    temp2 = dlmread(ext_filename{f}); 
    data{f} = [temp1];
    data{f} = [temp1;temp2];
    data_log{f} = data{f};
    data_log{f}(:,var_log) = log10(data{f}(:,var_log)); 
end
sub=randperm(size(data{1},1),size(data{end},1));data{1}=data{1}(sub,:);data_log{1}=data_log{1}(sub,:); %subsample original sim to same size
%% Boxplot
close all;
f=figure(1);clf;
f.Position(2)=f.Position(2)*0.5;
f.Position(4)=f.Position(4)*1.5;

%% Violin Plot of minFST in different model simplifications
close all; 
mfs_all = []; model = [];
for f=1:length(filename)
    mfs_all = [mfs_all; min([data{f}(:,fsw1(j)),data{f}(:,fsw2(j))],[],2)];
    model = [model;cellstr(repmat(mymodel{f}, size(data{f},1), 1))];
end

%Set minFST=1 to 2 too avoid gap in minFST
mfs_all(find(mfs_all==1))=2;

f=figure(1);clf;
f.Position(2)=f.Position(2)*0.5;
f.Position(4)=f.Position(4)*1.5;
%%% Original model
ax=axes;
%Boxes show interquartile range (Q1 to Q3), whiskers extend to 1.5x interquartile
%range (Q1-1.5*(Q3-Q1) to Q3 + 1.5*(Q3-Q1))?
grouporder={'Ori\_Model','1OFF','No\_term','1step\_ini','No\_SDI','No\_PR','No\_PC','1OFF\_1step\_ini','1OFF\_no\_SDI', ...
    '1step\_ini\_no\_SDI','1OFF\_1step\_ini\_no\_SDI'}
vs = violinplot(mfs_all(strcmp(model,'Ori\_Model')), model(strcmp(model,'Ori\_Model')), 'Width', 0.3, 'EdgeColor', cm2(1,:),'BoxColor',cm2(3,:),'ViolinColor',cm2(1,:),'ViolinAlpha',0.3, 'MarkerSize',2, ...
    'MedianMarkerSize', 10,'MedianColor',[217,95,2]./255, 'ShowMean', true,'ShowMedian', true, 'ShowWhiskers', false,'ShowNotches',false, ...
    'ShowBox',true, 'GroupOrder', grouporder);
hold on;yline(10, ':'); 
hold on;yline(100, ':'); 
ylabel('minFST [h]');
set(gca,'YScale','log','XLim',[0.5, 1.5],'YLim',[1.5 600],'YTick',[10 100 500],'YTickLabel',{'10', '100', '>500'},'TickLength',[0.004 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[1.9 3.8 graph_size1]);
xtickangle(ax,45)
%%% Reduced models
ax=axes; %'BoxColor',cm2(3,:)
vs = violinplot(mfs_all(~strcmp(model,'Ori\_Model')), model(~strcmp(model,'Ori\_Model')), 'Width', 0.3, 'EdgeColor', cm2(1,:),'BoxColor',cm2(3,:),'ViolinColor',cm2(1,:),'ViolinAlpha',0.3, 'MarkerSize',2, ...
    'MedianMarkerSize', 10,'MedianColor',[217,95,2]./255, 'ShowMean', true,'ShowMedian', true, 'ShowWhiskers', false,'ShowNotches',false, ...
    'ShowBox',true, 'GroupOrder', grouporder);
hold on;yline(10, ':','short-term memory'); 
hold on;yline(100, ':','long-term memory'); 
set(gca,'YScale','log','XLim',[0.5, 10.5],'YLim',[1.5 600],'YTickLabel',{},'TickLength',[0.004 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos(11,:) graph_size]);
xtickangle(ax,45)

%% Add BarPlot with fraction of sets generating no, short-term, and long-term memory
figure(1);
X=[];mean_mfs=[];
for i=1:length(filename)
    % Classification into stable and unstable sets
    mem_no{i} = find(min([data_log{i}(:,fsw1(j)),data_log{i}(:,fsw2(j))],[],2)<10); 
    mem_short{i} = find(min([data_log{i}(:,fsw1(j)),data_log{i}(:,fsw2(j))],[],2)>10 & ...
        min([data_log{i}(:,fsw1(j)),data_log{i}(:,fsw2(j))],[],2)<100); 
    mem_long{i} = find(min([data_log{i}(:,fsw1(j)),data_log{i}(:,fsw2(j))],[],2)>100); 
    fr_mem(i,1) = length(mem_no{i})./size(data_log{i},1);
    fr_mem(i,2) = length(mem_short{i})./size(data_log{i},1);
    fr_mem(i,3) = length(mem_long{i})./size(data_log{i},1);
    
    to_plot(i,:)=fr_mem(i,:).*100;
    X = [X; categorical({mymodel{i}})];
end

%%% Original model
ax = axes
%top plot
b=bar(to_plot(1,2:3)); 
b(1).FaceColor=cm2(2,:);
ylabel('Parameter Sets [%]','Fontsize',fs)
set(gca,'YLim',[4.5 10],'YTick',[6 8 10],'XTickLabel',{},'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(1)-0.8 pos_y(1)-1 6/11 0.8])
% bottom plot
ax = axes
b2=bar(to_plot(1,2:3));
b2(1).FaceColor=cm2(2,:);
set(gca,'YLim',[0 2.3],'YTick',[0 1 2],'TickLength',[0.02 0],'XTickLabel',{},'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(1)-0.8 pos_y(1)-1.9 6/11 0.8])

%%% Reduced models
ax = axes
%top plot
b=bar(to_plot(2:end,2:3)); %[pos_x(1) pos_y(1)-1.2 graph_size2]
b(1).FaceColor=cm2(2,:);
b(2).FaceColor=cm2(3,:);
%ylabel('Parameter Sets [%]','Fontsize',fs)
set(gca,'YLim',[4.5 10],'YTickLabel',{},'XTickLabel',{},'TickLength',[0.004 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(1) pos_y(1)-1 6-6/11 0.8])
% bottom plot
ax = axes
b2=bar(to_plot(2:end,2:3));
b2(1).FaceColor=cm2(2,:);
b2(2).FaceColor=cm2(3,:);
set(gca,'YLim',[0 2.3],'YTickLabel',{},'XTickLabel',{},'TickLength',[0.004 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(1) pos_y(1)-1.9 6-6/11 0.8])

set(gcf,'renderer','Painters')
print('../../figures/Fig1fg','-depsc','-loose')



%% Example trajectories for minFST illustration
%% Simulate example sets with no, short-term and long-term memory

long_mem = find(min([data_log{1}(:,fsw1(j)),data_log{1}(:,fsw2(j))],[],2)>100); 
rand_long = randperm(length(long_mem),3);
rand_long = long_mem(rand_long);

short_mem = find(min([data_log{1}(:,fsw1(j)),data_log{1}(:,fsw2(j))],[],2)>10 & min([data_log{1}(:,fsw1(j)),data_log{1}(:,fsw2(j))],[],2)<100); 
rand_short = randperm(length(short_mem),3);
rand_short = short_mem(rand_short);

no_mem = find(min([data_log{1}(:,fsw1(j)),data_log{1}(:,fsw2(j))],[],2)<10); 
rand_no = randperm(length(no_mem),3);
rand_no = no_mem(rand_no);

%dlmwrite(['pars_memory_timescales_MFS.txt'], [data{1}(rand_long, 1:nr_p); data{1}(rand_short, 1:nr_p); data{1}(rand_no, 1:nr_p)]);
inputFilename = ['pars_memory_timescales_MFS.txt'];
outputFilename = ['sim_memory_timescales_MFS.txt'];

% Actual simulation commented out to speed up plotting
%addpath('../scripts_modeling/');
%function_general_as_model_200319_2(inputFilename,outputFilename);

% Instead: Read in simulation file
sim=dlmread(outputFilename);

t=501;
ap1 = [nr_p+1:1:nr_p+t];
bp1 = [nr_p+t+1:1:nr_p+2*t];
ar1 = [nr_p+2*t+1:nr_p+3*t];
br1 = [nr_p+3*t+1:nr_p+4*t];
apr1 = [nr_p+4*t+1:nr_p+5*t];
bpr1 = [nr_p+5*t+1:nr_p+6*t];
%average state of pA per time step
pAo1 = [nr_p+6*t+1:nr_p+7*t];
pBo1 = [nr_p+7*t+1:nr_p+8*t];
ap2 = [nr_p+8*t+1:nr_p+9*t];
bp2 = [nr_p+9*t+1:nr_p+10*t];
ar2 = [nr_p+10*t+1:nr_p+11*t];
br2 = [nr_p+11*t+1:nr_p+12*t];
apr2 = [nr_p+12*t+1:nr_p+13*t];
bpr2 = [nr_p+13*t+1:nr_p+14*t];
pAo2 = [nr_p+14*t+1:nr_p+15*t];
pBo2 = [nr_p+15*t+1:nr_p+16*t];

%% Extract first switching Time for each simulated allele
% Loop over parameter sets and allele
%Calculate Ratio of pol on A vs B normalized to transcription strength
p_rat_ini1 = ((sim(:,ap1)+1)./((1./(1./sim(:,1)+1./sim(:,11))).*(1-sim(:,13))))./((sim(:,bp1)+1)./((1./(1./sim(:,2)+1./sim(:,12))).*(1-sim(:,14))));
p_rat_ini2 = ((sim(:,ap2)+1)./((1./(1./sim(:,1)+1./sim(:,11))).*(1-sim(:,13))))./((sim(:,bp2)+1)./((1./(1./sim(:,2)+1./sim(:,12))).*(1-sim(:,14))));
% Ratio one time step delayed 
p_rat_1_1step_del = [10*ones(size(p_rat_ini1,1),1) p_rat_ini1(:,1:end-1)];
p_rat_2_1step_del = [0.1*ones(size(p_rat_ini2,1),1) p_rat_ini2(:,1:end-1)];
fst1 = t*ones(size(sim,1),1);
fst2 = t*ones(size(sim,1),1);
for h = 1:size(sim,1)
% switches direction 1 ini 1 
sw1_ini1_p = find(p_rat_ini1(h,:)<=1 & p_rat_1_1step_del(h,:)>1);
%switches direction 1 ini2
sw1_ini2_p = find(p_rat_ini2(h,:)>=1 & p_rat_2_1step_del(h,:)<1);
%set to total simulation time if no switch happens
if length(sw1_ini1_p)==0
    sw1_ini1_p = t-1;
end
if length(sw1_ini2_p)==0
    sw1_ini2_p = t-1;
end
fst1(h) = sw1_ini1_p(1);
fst2(h) = sw1_ini2_p(1);
end
mfs = min([fst1 fst2],[],2);

%% select one realisation for each simulated set
% initialize random number generator to make results reproducible
rng(0,'twister');
ind_sets = [];
for i=1:size(sim,1)./100
    a=(i-1)*100;
    r = a+randperm(100,1);
    ind_sets = [ind_sets; r];
end
%% Plot example trajectories for minFST Illustration
f=figure(2);clf;
f.Position(2)=440;
f.Position(4)=530;

%Plotting options: Option1: plot both genes together, separate initial
%conditions; Option2: Plot initial conditions together, separate genes;
% Option3: Plot initial conditions and genes separate
opt_pl = 3;

for i=1:length(ind_sets)
    % I) Plot one example trajectories per stability category for Fig1
    if rem(i,3)==1
        figure(2);
        axes;
    if opt_pl==1
        % Plot both genes together, separate initial conditions
        plot([0:t-1],sim(ind_sets(i),ap1),'Color',cm_a(1,:));hold on; plot([0:t-1],sim(ind_sets(i),bp1),'Color',cm_b(1,:)); 
        hold on; xline(fst1(ind_sets(i)), ':','Linewidth',1.5);
        if i==1
            legend('A','B', 'FST'); title('Ini1');
        end
        title({'minFST=',num2str(min([fst1(ind_sets(i)) fst2(ind_sets(i))]))});
        y_lim = [0-max(max([sim(ind_sets(i),ap1);sim(ind_sets(i),ap2);sim(ind_sets(i),bp1);sim(ind_sets(i),bp2)]))*0.1 max(max([sim(ind_sets(i),ap1);sim(ind_sets(i),ap2);sim(ind_sets(i),bp1);sim(ind_sets(i),bp2)]))*1.1];
        set(gca,'YLim',y_lim, 'XLim',[-50 550],'XTick', [0 200 400],'TickLength',[0.02 0],'TickDir','out','Linewidth',lw,'Fontsize',fs,'Units','Centimeters','Position',[[pos_x(4-ceil(i/3)) pos_y(1)] graph_size3]); 
    elseif opt_pl==2
        % Plot initial conditions together, separate genes
        plot([0:t-1],sim(ind_sets(i),ap1),'Color',cm_a(1,:));hold on; plot([0:t-1],sim(ind_sets(i),ap2),'Color',cm_a(2,:)); 
        hold on; xline(fst1(ind_sets(i)), ':','Linewidth',1.5,'Color',cm_a(1,:));
        hold on; xline(fst2(ind_sets(i)), ':','Linewidth',1.5,'Color',cm_a(2,:));
        if i==1
            legend('Ini1','Ini2', 'FST1', 'FST2'); title('Gene A');
        end
        title({'minFST=',num2str(min([fst1(ind_sets(i)) fst2(ind_sets(i))]))});
        y_lim = [0-max(max([sim(ind_sets(i),ap1);sim(ind_sets(i),ap2);sim(ind_sets(i),bp1);sim(ind_sets(i),bp2)]))*0.1 max(max([sim(ind_sets(i),ap1);sim(ind_sets(i),ap2);sim(ind_sets(i),bp1);sim(ind_sets(i),bp2)]))*1.1];
        set(gca,'YLim',y_lim, 'XLim',[-50 550],'XTick', [0 200 400],'TickLength',[0.02 0],'TickDir','out','Linewidth',lw,'Fontsize',fs,'Units','Centimeters','Position',[[pos_x(4-ceil(i/3)) pos_y(1)] graph_size3]); 
    elseif opt_pl==3
        % Plot initial conditions and genes seperate
        plot([0:t-1],sim(ind_sets(i),ap1),'Color',cm_a(1,:));
        y_lim = [0-max(max([sim(ind_sets(i),ap1);sim(ind_sets(i),ap2);sim(ind_sets(i),bp1);sim(ind_sets(i),bp2)]))*0.1 max(max([sim(ind_sets(i),ap1);sim(ind_sets(i),ap2);sim(ind_sets(i),bp1);sim(ind_sets(i),bp2)]))*1.1];
        set(gca,'YLim',y_lim, 'XLim',[-50 550],'XTick', [0 200 400],'TickLength',[0.02 0],'TickDir','out','Linewidth',lw,'Fontsize',fs,'Units','Centimeters','Position',[[pos_x(4-ceil(i/3)) pos_y(1)] [0.75 1.5]]); 
        title({'minFST=',num2str(min([fst1(ind_sets(i)) fst2(ind_sets(i))]))});
        axes;
        plot([0:t-1],sim(ind_sets(i),bp1),'Color',cm_b(1,:));
        set(gca,'YLim',y_lim, 'XLim',[-50 550],'YTickLabel',[],'XTick', [0 200 400],'TickLength',[0.02 0],'TickDir','out','Linewidth',lw,'Fontsize',fs,'Units','Centimeters','Position',[[pos_x2(4-ceil(i/3)) pos_y(1)] [0.75 1.5]]); 
        %title({'FST1=',num2str(fst1(ind_sets(i)))});
    end
    
    axes;
    if opt_pl==1
    % Plot both genes together, separate initial conditions
    plot([0:t-1],sim(ind_sets(i),ap2),'Color',cm_a(2,:));hold on; plot([0:t-1],sim(ind_sets(i),bp2),'Color',cm_b(2,:)); 
    hold on; xline(fst2(ind_sets(i)), ':','Linewidth',1.5);
    if i==1
        legend('A','B', 'FST'); title('Ini2');
    end
    title({'minFST=',num2str(min([fst1(ind_sets(i)) fst2(ind_sets(i))]))});
    set(gca,'YLim',y_lim,'XLim',[-50 550],'XTick', [0 200 400],'TickLength',[0.02 0],'TickDir','out','Linewidth',lw,'Fontsize',fs,'Units','Centimeters','Position',[[pos_x(4-ceil(i/3)) pos_y(2)] graph_size3]);
    elseif opt_pl==2
        %Plot initial conditions together, separate genes
        plot([0:t-1],sim(ind_sets(i),bp1),'Color',cm_b(1,:));hold on; plot([0:t-1],sim(ind_sets(i),bp2),'Color',cm_b(2,:)); 
        hold on; xline(fst1(ind_sets(i)), ':','Linewidth',1.5,'Color',cm_b(1,:));
        hold on; xline(fst2(ind_sets(i)), ':','Linewidth',1.5,'Color',cm_b(2,:));
        if i==1
            legend('Ini1','Ini2', 'FST1', 'FST2'); title('Gene B');
        end
        title({'minFST=',num2str(min([fst1(ind_sets(i)) fst2(ind_sets(i))]))});
        set(gca,'YLim',y_lim,'XLim',[-50 550],'XTick', [0 200 400],'TickLength',[0.02 0],'TickDir','out','Linewidth',lw,'Fontsize',fs,'Units','Centimeters','Position',[[pos_x(4-ceil(i/3)) pos_y(2)] graph_size3]); 
    elseif opt_pl==3
        %Plot both, initial conditions and genes, seperate
        plot([0:t-1],sim(ind_sets(i),ap2),'Color',cm_a(2,:));
        set(gca,'YLim',y_lim,'XLim',[-50 550],'XTick', [0 200 400],'TickLength',[0.02 0],'TickDir','out','Linewidth',lw,'Fontsize',fs,'Units','Centimeters','Position',[[pos_x(4-ceil(i/3)) pos_y(2)] [0.75 1.5]]); 
        title({'minFST=',num2str(min([fst1(ind_sets(i)) fst2(ind_sets(i))]))});
        axes;
        plot([0:t-1],sim(ind_sets(i),bp2),'Color',cm_b(2,:));
        set(gca,'YLim',y_lim,'XLim',[-50 550],'YTickLabel',[],'XTick', [0 200 400],'TickLength',[0.02 0],'TickDir','out','Linewidth',lw,'Fontsize',fs,'Units','Centimeters','Position',[[pos_x2(4-ceil(i/3)) pos_y(2)] [0.75 1.5]]); 
        %title({'FST2=',num2str(fst2(ind_sets(i)))});
    end
    end
    
end
figure(2);
print(gcf,'../../figures/Fig1e','-depsc','-loose')