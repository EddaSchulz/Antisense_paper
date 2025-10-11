clear all; close all;
%% Plotting parameters
% Color scheme initial conditions
cm_a = [27 157 142; 148 195 186]./255;
cm_b = [235 109 30; 245 172 118]./255;
% Color scheme memory timescales
cm2 = [204 204 204; 150 150 150; 82 82 82]./255;
% color schemes mitosis models
cm3 = [95 60 76; 40 66 81]./255;
graph_size=[2*6/11 6-6/11];
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

oripath = '../../simulations/maximally_simplified_model/';
path_fst = '../../simulations/revision/mitosis/';

filename{1} = ([oripath, 'FST_max_simpl.txt']);
filename{2} = ([path_fst, 'FST_revision_mitosis_noPreact_250301.txt']);
filename{3} = ([path_fst, 'FST_revision_mitosis_250301.txt']);

mymodel{1} = 'no-mitosis';
mymodel{2} = 'mitosis-I';
mymodel{3} = 'mitosis-II';

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
thresh_st = 96; 
thresh_unst = 3;
fr_stable = NaN(length(filename),3);
%% Read in simulations
for f=1:length(filename)
    data{f} = dlmread(filename{f}); 
    data_log{f} = data{f};
    data_log{f}(:,var_log) = log10(data{f}(:,var_log)); 
end
% select only the parameter sets that were simulated with mitosis from the
% original simulation w/o mitosis
rel = find(min([data{1}(:,fsw1(j)),data{1}(:,fsw2(j))],[],2)>400);
data{1} = data{1}(rel,:);
data_log{1} = data_log{1}(rel,:);

%% PLOTS

% Boxplot
close all;
f=figure(1);clf;
f.Position(2)=f.Position(2)*0.5;
f.Position(4)=f.Position(4)*1.5;

% Violin Plot of minFST in different model simplifications
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
grouporder={'no-mitosis','mitosis-I','mitosis-II'}
vs = violinplot(mfs_all(strcmp(model,'no-mitosis')), model(strcmp(model,'no-mitosis')), 'Width', 0.3, 'EdgeColor', cm2(1,:),'BoxColor',cm2(3,:),'ViolinColor',cm2(1,:),'ViolinAlpha',0.3, 'MarkerSize',2, ...
    'MedianMarkerSize', 10,'MedianColor',[217,95,2]./255, 'ShowMean', true,'ShowMedian', true, 'ShowWhiskers', false,'ShowNotches',false, ...
    'ShowBox',true, 'GroupOrder', grouporder);
hold on;yline(10, ':'); 
hold on;yline(100, ':'); 
ylabel('minFST [h]');
set(gca,'YScale','log','XLim',[0.5, 1.5],'YLim',[1.5 600],'YTick',[10 100 500],'YTickLabel',{'10', '100', '>500'},'TickLength',[0.004 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[1.9 3.8 graph_size1]);
xtickangle(ax,45)
%%% Reduced models
ax=axes; %'BoxColor',cm2(3,:)
vs = violinplot(mfs_all(~strcmp(model,'no-mitosis')), model(~strcmp(model,'no-mitosis')), 'Width', 0.3, 'EdgeColor', cm2(1,:),'BoxColor',cm2(3,:),'ViolinColor',cm2(1,:),'ViolinAlpha',0.3, 'MarkerSize',2, ...
    'MedianMarkerSize', 10,'MedianColor',[217,95,2]./255, 'ShowMean', true,'ShowMedian', true, 'ShowWhiskers', false,'ShowNotches',false, ...
    'ShowBox',true, 'GroupOrder', grouporder);
%hold on;yline(10, ':','short-term memory'); 
%hold on;yline(100, ':','long-term memory'); 
set(gca,'YScale','log','XLim',[0.5, 2.5],'YLim',[1.5 600],'YTickLabel',{},'TickLength',[0.004 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos(11,:) graph_size]);
xtickangle(ax,45)

% Add BarPlot with fraction of sets generating no, short-term, and long-term memory
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
b(1).FaceColor=cm2(3,:);
ylabel('Parameter Sets [%]','Fontsize',fs)
set(gca,'YLim',[50 105],'YTick',[60 80 100],'XTickLabel',{},'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(1)-0.8 pos_y(1)-1 6/11 0.8])
% bottom plot
ax = axes
b2=bar(to_plot(1,2:3));
b2(1).FaceColor=cm2(3,:);
set(gca,'YLim',[0 33.5],'YTick',[10 20 30],'TickLength',[0.02 0],'XTickLabel',{},'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(1)-0.8 pos_y(1)-1.9 6/11 0.8])

%%% Reduced models
ax = axes
%top plot
b=bar(to_plot(2:end,2:3)); %[pos_x(1) pos_y(1)-1.2 graph_size2]
b(1).FaceColor=cm2(2,:);
b(2).FaceColor=cm2(3,:);
%ylabel('Parameter Sets [%]','Fontsize',fs)
set(gca,'YLim',[50 105],'YTickLabel',{},'XTickLabel',{},'TickLength',[0.004 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(1) pos_y(1)-1 2*6/11 0.8])
% bottom plot
ax = axes
b2=bar(to_plot(2:end,2:3));
b2(1).FaceColor=cm2(2,:);
b2(2).FaceColor=cm2(3,:);
set(gca,'YLim',[0 33.5],'YTickLabel',{},'XTickLabel',{},'TickLength',[0.004 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(1) pos_y(1)-1.9 2*6/11 0.8])

set(gcf,'renderer','Painters')
print('../../figures/Fig3E','-depsc','-loose')



