% Script analyses and generates plots for simplified model
clear all;
close all;

oripath = '../../simulations/maximally_simplified_model/';
mypath='../../simulations/revision/bursting/';
localpath = '';


%% Plotting parameters 
% Color scheme Locus architectures
cm=[127.5 127.5 127.5; 71 55 136; 191 66 134; 246 166 65]./255;
% Color scheme original par sets
cm_ori=[255,255,204;217,240,163;173,221,142;120,198,121;49,163,84;0,104,55]./255;
cm=[127.5 127.5 127.5; 71 55 136; 191 66 134; 246 166 65]./255;
% Color scheme memory timescales
cm2 = [204 204 204; 150 150 150; 82 82 82]./255;
graph_size=[2.5 2.5];
graph_size3=[1.5 1.5];
lw=1;
fs=8;
fst=10;
pos_x=[0.1:2.5:35]; %pos_y=[2 6.5 11 15.5 19.5:3.8:30];
pos_y=[0:2.2:11.4]; pos_y=flip(pos_y);

%% Simulation files max simpl model
filename{1} = ([mypath, 'FST_revision_bursting_stabsets_rand_250527.txt']);
filename{2} = ([oripath, 'FST_max_simpl.txt']);

% Models
mymodel{1} = 'bursting';
mymodel{2} = 'ori';

% Names for Figures
mymodelnames{1} = 'bursting';
mymodelnames{2} = 'ori';

% RNAP, RNA or Protein level
mylevel{1} = 'RNAP';
mylevel{2} = 'RNA';
mylevel{3} = 'Protein';
j = 1; %Level to perform analysis on

%% Indices
nr_p = 30;
nbins = 10;
var_par = [15 16 17 18];
var_log = [5 6 15 16 17 18];
fsw1 = nr_p+1:nr_p+3;
fsw2 = nr_p+4:nr_p+6;
thresh_st = 96; %Stability threshold
thresh_unst = 3;
fr_stable = NaN(length(filename),3);
sim = [];
%% Read in simulations
for f=1:length(filename)
    data{f} = dlmread(filename{f}); 
    data_log{f} = data{f};
    data_log{f}(:,var_log) = log10(data{f}(:,var_log)); 
end
% select only the parameter sets that were simulated with mitosis from the
% original simulation w/o mitosis
rel = find(min([data{2}(:,fsw1(j)),data{2}(:,fsw2(j))],[],2)>400);
data{2} = data{2}(rel,:);
data_log{2} = data_log{2}(rel,:);


%% Violin plot of MFS 
mfs_all = []; model = [];
for f=1:length(filename)
    mfs_all = [mfs_all; min([data{f}(:,fsw1(j)),data{f}(:,fsw2(j))],[],2)];
    model = [model;cellstr(repmat(mymodel{f}, size(data{f},1), 1))];
end

f=figure(1);clf;
f.Position = [680   388   380   580];
%ori
ax=axes;
%Boxes show interquartile range (Q1 to Q3), whiskers extend to 1.5x interquartile
%range (Q1-1.5*(Q3-Q1) to Q3 + 1.5*(Q3-Q1))
vs = violinplot(mfs_all(strcmp(model,'ori')), model(strcmp(model,'ori')), 'Width', 0.3, 'EdgeColor', cm(1,:),'BoxColor',cm2(3,:),'ViolinColor',cm(4,:),'ViolinAlpha',0.3, 'MarkerSize',2, ...
    'MedianMarkerSize', 10,'MedianColor',[217,95,2]./255, 'ShowMean', true,'ShowMedian', true, 'ShowWhiskers', false,'ShowNotches',false,'ShowBox',true);
hold on;yline(10, ':','short-term memory','Fontsize',6); 
hold on;yline(100, ':','long-term memory','Fontsize',6); 
ylabel('minFST [h]');
set(gca,'YScale','log','XLim',[0.5, 1.5],'YLim',[1.5 600],'YTick',[10 100 500],'YTickLabel',{'10', '100', '>500'},'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',6,'Units','Centimeters','Position',[1.9 3.8 1.63/3 6]);
xtickangle(ax,45)

%bursting
ax=axes; %'BoxColor',cm2(3,:)
vs = violinplot(mfs_all(strcmp(model,'bursting')), model(strcmp(model,'bursting')), 'Width', 0.3, 'EdgeColor', cm2(1,:),'BoxColor',cm2(3,:),'ViolinColor',cm(4,:),'ViolinAlpha',0.3, 'MarkerSize',2, ...
    'MedianMarkerSize', 10,'MedianColor',[217,95,2]./255, 'ShowMean', true,'ShowMedian', true, 'ShowWhiskers', false,'ShowNotches',false, ...
    'ShowBox',true);
set(gca,'YScale','log','XLim',[0.5, 1.5],'YLim',[1.5 600],'YTickLabel',{},'TickLength',[0.004 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[2.7 3.8 1.63/3 6]);
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
b=bar(to_plot(2,2:3)); 
b(1).FaceColor=cm2(2,:);
ylabel('Parameter Sets [%]','Fontsize',fs)
set(gca,'YLim',[50 105],'YTick',[60 80 100],'XTickLabel',{},'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',6,'Units','Centimeters','Position',[1.9 11.3 1.63/3 0.7])
% bottom plot
ax = axes
b2=bar(to_plot(2,2:3));
b2(1).FaceColor=cm2(2,:);
set(gca,'YLim',[0 33.5],'YTick',[10 20 30],'XTickLabel',{},'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',6,'Units','Centimeters','Position',[1.9 11-0.8 1.63/3 1])
%%% bursting
ax = axes
%top plot
b=bar(to_plot(1,2:3)); 
b(1).FaceColor=cm2(2,:);
%b(2).FaceColor=cm2(3,:);
set(gca,'YLim',[50 105],'YTick',[60 80 100],'XTickLabel',{},'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',6,'Units','Centimeters','Position',[2.7 11.3 1.63/3 0.7])
% bottom plot
ax = axes
b2=bar(to_plot(1,2:3));
b2(1).FaceColor=cm2(2,:);
%b2(2).FaceColor=cm2(3,:);
set(gca,'YLim',[0 33.5],'YTick',[10 20 30],'XTickLabel',{},'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',6,'Units','Centimeters','Position',[2.7 11-0.8 1.63/3 1])

set(gcf,'renderer','Painters')
print('../../figures/Fig3B','-depsc','-loose')

%% Parameter distribution in sets of different stability 
var_par = [15 16 17 18];
var_log = [5 6 15 16 17 18];

% Define different stability categories:
thresh_no_mem = 10;
thresh_short_mem = [10 100];
thresh_long_mem = 100;

plot_names = {'k^{A}_{ON} (log10,h^{-1})','k^{B}_{ON} (log10,h^{-1})','1/k^{A}_{OFF} (log10,h)','1/k^{B}_{OFF} (log10,h)', 'k^{A}_{ON} / (k^{A}_{ON} + k^{A}_{OFF})', 'k^{B}_{ON} / (k^{B}_{ON} + k^{B}_{OFF})', 't^{B}_{OFF-AS} / t^{A}_{OFF-basal} (log10)',  't^{A}_{OFF-AS} / t^{B}_{OFF-basal} (log10)'};
plot_titles = {'Burst frequency pA', 'Burst frequency pB','Burst duration pA','Burst duration pB', 'Burst fraction A', 'Burst fraction B', 'Stability PR / Stability OFF A',  'Stability PR / Stability OFF B'};

%% Bursting parameter distribution among all and stable sets
nbins=8; %8 12 => starts looking bimodal??
f=figure(2); clf;
f.Position = [680   517   950   450];
    no_mem = find(min([data_log{1}(:,fsw1(j)),data_log{1}(:,fsw2(j))],[],2)<thresh_no_mem);
    short_mem = find(min([data_log{1}(:,fsw1(j)),data_log{1}(:,fsw2(j))],[],2)>thresh_short_mem(1) & ...
        min([data_log{1}(:,fsw1(j)),data_log{1}(:,fsw2(j))],[],2)<thresh_short_mem(2));
    long_mem = find(min([data_log{1}(:,fsw1(j)),data_log{1}(:,fsw2(j))],[],2)>thresh_long_mem);
    for m = 1:length(var_par)+4
        ax=axes;
        if m<length(var_par)+1 & m>2
            temp1 = data_log{1}(:,var_par(m)); %all sets
            temp2 = data_log{1}((no_mem),var_par(m)); %no memory
            temp3 = data_log{1}((short_mem),var_par(m)); %short-term memory
            temp4 = data_log{1}((long_mem),var_par(m)); %long-term memory
        elseif m<3
            temp1 = log(1)-data_log{1}(:,var_par(m)); %all sets
            temp2 = log(1)-data_log{1}((no_mem),var_par(m)); %no memory
            temp3 = log(1)-data_log{1}((short_mem),var_par(m)); %short-term memory
            temp4 = log(1)-data_log{1}((long_mem),var_par(m)); %long-term memory
        elseif m==length(var_par)+1 
            temp1= (1./data{1}(:,15))./(1./data{1}(:,15)+1./data{1}(:,17));
            temp2= (1./data{1}((no_mem),15))./(1./data{1}((no_mem),15)+1./data{1}((no_mem),17));
            temp3= (1./data{1}((short_mem),15))./(1./data{1}((short_mem),15)+1./data{1}((short_mem),17));
            temp4= (1./data{1}((long_mem),15))./(1./data{1}((long_mem),15)+1./data{1}((long_mem),17));            
        elseif m==length(var_par)+2 
            temp1= (1./data{1}(:,16))./(1./data{1}(:,16)+1./data{1}(:,18));
            temp2= (1./data{1}((no_mem),16))./(1./data{1}((no_mem),16)+1./data{1}((no_mem),18));
            temp3= (1./data{1}((short_mem),16))./(1./data{1}((short_mem),16)+1./data{1}((short_mem),18));
            temp4= (1./data{1}((long_mem),16))./(1./data{1}((long_mem),16)+1./data{1}((long_mem),18));
        elseif m==length(var_par)+3 
            temp1 = (data_log{1}(:,6)-data_log{1}(:,15));
            temp2 =(data_log{1}(no_mem,6)-data_log{1}(no_mem,15));
            temp3 = (data_log{1}(short_mem,6)-data_log{1}(short_mem,15));
            temp4 = (data_log{1}(long_mem,6)-data_log{1}(long_mem,15));
        elseif m==length(var_par)+4 
            temp1 = (data_log{1}(:,5)-data_log{1}(:,16));
            temp2 = (data_log{1}(no_mem,5)-data_log{1}(no_mem,16));
            temp3 = (data_log{1}(short_mem,5)-data_log{1}(short_mem,16));
            temp4 = (data_log{1}(long_mem,5)-data_log{1}(long_mem,16));
        elseif m==length(var_par)+5 
            temp1 = data_log{1}(:,5)-data_log{1}(:,16)-data_log{1}(:,2);
            temp2 = data_log{1}(no_mem,5)-data_log{1}(no_mem,16)-data_log{1}(no_mem,2);
            temp3 = data_log{1}(short_mem,5)-data_log{1}(short_mem,16)-data_log{1}(short_mem,2);
            temp4 = data_log{1}(long_mem,5)-data_log{1}(long_mem,16)-data_log{1}(long_mem,2);
        end
        
        h1 = histogram(temp1,nbins);
        edges = h1.BinEdges;
        counts1 = h1.BinCounts;
        temp = h1.BinEdges + h1.BinWidth/2;
        centers1 = temp(1:end-1);
        h2 = histogram(temp2, edges);
        counts2 = h2.BinCounts;
        temp = h2.BinEdges + h2.BinWidth/2;
        centers2 = temp(1:end-1);
        h3 = histogram(temp3, edges);
        counts3 = h3.BinCounts;
        temp = h3.BinEdges + h3.BinWidth/2;
        centers3 = temp(1:end-1);
        h4 = histogram(temp4, edges);
        counts4 = h4.BinCounts;
        temp = h4.BinEdges + h4.BinWidth/2;
        centers4 = temp(1:end-1);
        lp2=plot((centers2),100*counts2./sum(counts2),'Color',cm2(1,:)); hold on;
        lp3=plot((centers3),100*counts3./sum(counts3),'Color',cm2(2,:)); hold on;
        lp4=plot((centers4),100*counts4./sum(counts4),'Color',cm2(3,:)); hold on;
        set(lp2,'LineWidth',lw);set(lp3,'LineWidth',lw);set(lp4,'LineWidth',lw);
        y_lim=max([100*counts1./sum(counts1) 100*counts2./sum(counts2) 100*counts3./sum(counts3) 100*counts4./sum(counts4)])*1.1;
        set(gca,'YLim',[0 y_lim],'YTick',[],'YTickLabel',[], 'TickLength',[0.02 0],'TickDir','out','Linewidth',lw,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(m) pos_y(3) graph_size3]);
        
        
        title(plot_titles{m}, 'Fontsize',6);
        xlabel(plot_names{m}, 'Fontsize',6);
        if m==1
            ylabel('Relative frequency','Fontsize',6);
            
        end
        
        
    end
 
     print(['../../figures/Fig3C'],'-depsc','-loose')

