% Script analyses and generates plots for simplified model
clear all;
close all;

mypath='../../simulations/maximally_simplified_model/';
localpath = '';


%% Plotting parameters 
% Color scheme Locus architectures
cm=[127.5 127.5 127.5; 71 55 136; 191 66 134; 246 166 65]./255;
% Color scheme memory timescales
cm2 = [204 204 204; 150 150 150; 82 82 82]./255;
graph_size=[2.5 2.5];
graph_size3=[1.5 1.5];
lw=1;
fs=8;
fst=10;
pos_x=[2.7:2.5:20.2]; %pos_y=[2 6.5 11 15.5 19.5:3.8:30];
pos_y=[0:2.2:11.4]; pos_y=flip(pos_y);

%% Simulation files max simpl model
filename{1} = ([mypath, 'fsw_rand_par_var_simplified_as_model_200319.txt']);
filename{2} = ([mypath, 'fsw_rand_add_par_simpmodel_200319_220125.txt']);

% Models
mymodel{1} = 'simplified\_model';

% Names for Figures
mymodelnames{1} = 'simplified_model';

% RNAP, RNA or Protein level
mylevel{1} = 'RNAP';
mylevel{2} = 'RNA';
mylevel{3} = 'Protein';
j = 1; %Level to perform analysis on

%% Indices
nr_p = 30;
nbins = 10;
var_par = [1 2 5 6 21 24 26];
var_log = [1 2 5 6];
fsw1 = nr_p+1:nr_p+3;
fsw2 = nr_p+4:nr_p+6;
thresh_st = 96; %Stability threshold
thresh_unst = 3;
fr_stable = NaN(length(filename),3);
sim = [];
%% Read in simulations
for f=1:length(filename)
	temp = dlmread(filename{f}); 
	sim = [sim;temp];	
end
data{1}= sim; 
data_log{1} = data{1};
data_log{1}(:,var_log) = log10(data{1}(:,var_log)); 
clear sim;
%% Classify Sets according to locus architecture
%Stability of PR
thresh_st_pr=1; thresh_ust_pr=0.0167; %stable PR>1h and unstable PR<1min

% 3' Overlap: both genes only induce unstable PR 
ol_3 = find(data{1}(:,5)<thresh_ust_pr & data{1}(:,6)<thresh_ust_pr);
% Intragenic overlap: One promoter induces stable, the other only unstable PR
ol_intra = find((data{1}(:,5)<thresh_ust_pr&data{1}(:,6)>thresh_st_pr)|(data{1}(:,5)>thresh_st_pr&data{1}(:,6)<thresh_ust_pr));
% Order Intragenic Overlaps by embedded and long gene rather than by S and AS strand
ol_intra_S = find((data{1}(:,5)<thresh_ust_pr&data{1}(:,6)>thresh_st_pr)); % sense = long gene
ol_intra_AS = find((data{1}(:,5)>thresh_st_pr&data{1}(:,6)<thresh_ust_pr)); % antisense = long gene
% Define sense as long, and AS as embedded gene:
temp=data{1}(ol_intra_AS,5);data{1}(ol_intra_AS,5)=data{1}(ol_intra_AS,6);data{1}(ol_intra_AS,6)=temp;
temp=data{1}(ol_intra_AS,1);data{1}(ol_intra_AS,1)=data{1}(ol_intra_AS,2);data{1}(ol_intra_AS,2)=temp;
temp=data_log{1}(ol_intra_AS,5);data_log{1}(ol_intra_AS,5)=data_log{1}(ol_intra_AS,6);data_log{1}(ol_intra_AS,6)=temp;
temp=data_log{1}(ol_intra_AS,1);data_log{1}(ol_intra_AS,1)=data_log{1}(ol_intra_AS,2);data_log{1}(ol_intra_AS,2)=temp;
% Promoter Overlap: PR high and symmetric
ol_pr = find(data{1}(:,5)>thresh_st_pr & data{1}(:,6)>thresh_st_pr);

locus_arch{1} = ol_3;
locus_arch{2} = ol_intra;
locus_arch{3} = ol_pr;

locus_arch_names{1} = '3-Overlap';
locus_arch_names{2} = 'Intragenic-Overlap';
locus_arch_names{3} = 'Promoter-Overlap';

%% Violin plot of MFS for different Locus architectures

mfs_all = []; model = []; n_arch = [];
for i=1:length(locus_arch)
    mfs_all = [mfs_all; min([data{1}(locus_arch{i},fsw1(j)),data{1}(locus_arch{i},fsw2(j))],[],2)];
    model = [model;cellstr(repmat(locus_arch_names{i}, size(data{1}(locus_arch{i},1),1), 1))];
    n_arch = [n_arch, size(data{1}(locus_arch{i},:),1)];
end

%Set MFS=1 to 2 too avoid gap in MFS
mfs_all(find(mfs_all==1))=2;

f=figure(1);clf;
f.Position = [680   388   380   580];
ax=axes;
%Boxes show interquartile range (Q1 to Q3), whiskers extend to 1.5x interquartile
%range (Q1-1.5*(Q3-Q1) to Q3 + 1.5*(Q3-Q1))

grouporder={'3-Overlap','Intragenic-Overlap','Promoter-Overlap'}
vs = violinplot(mfs_all, model, 'Width', 0.3, 'EdgeColor', cm(1,:),'BoxColor',cm2(3,:),'ViolinColor',cm(2:4,:),'ViolinAlpha',0.3, 'MarkerSize',2, ...
    'MedianMarkerSize', 10,'MedianColor',[217,95,2]./255, 'ShowMean', true,'ShowMedian', true, 'ShowWhiskers', false,'ShowNotches',false,'ShowBox',true, 'GroupOrder', grouporder);
hold on;yline(10, ':','short-term memory','Fontsize',6); 
hold on;yline(100, ':','long-term memory','Fontsize',6); 
ylabel('minFST [h]');
set(gca,'YScale','log','XLim',[0.5, 3.5],'YLim',[1.5 600],'YTick',[10 100 500],'YTickLabel',{'10', '100', '>500'},'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',6,'Units','Centimeters','Position',[2.7 3.8 1.63 6]);
xtickangle(ax,45)

% Add BarPlot with fraction of sets generating no, short-term, and long-term memory
figure(1);
X=[];mean_mfs=[];
for i=1:length(locus_arch)
    % Classification into stable and unstable sets
    mem_no{i} = find(min([data_log{1}(locus_arch{i},fsw1(j)),data_log{1}(locus_arch{i},fsw2(j))],[],2)<10); 
    mem_short{i} = find(min([data_log{1}(locus_arch{i},fsw1(j)),data_log{1}(locus_arch{i},fsw2(j))],[],2)>10 & ...
        min([data_log{1}(locus_arch{i},fsw1(j)),data_log{1}(locus_arch{i},fsw2(j))],[],2)<100); 
    mem_long{i} = find(min([data_log{1}(locus_arch{i},fsw1(j)),data_log{1}(locus_arch{i},fsw2(j))],[],2)>100); 
    fr_mem(i,1) = length(mem_no{i})./size(data_log{1}(locus_arch{i},1),1);
    fr_mem(i,2) = length(mem_short{i})./size(data_log{1}(locus_arch{i},1),1);
    fr_mem(i,3) = length(mem_long{i})./size(data_log{1}(locus_arch{i},1),1);
    
    to_plot(i,:)=fr_mem(i,:).*100;
end

figure(1);
ax = axes
%top plot
b=bar(to_plot(:,2:3)); 
b(1).FaceColor=cm2(2,:);
b(2).FaceColor=cm2(3,:);
ylabel('Parameter Sets [%]','Fontsize',fs)
set(gca,'YLim',[33 40],'YTick',[35 40],'XTickLabel',{},'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',6,'Units','Centimeters','Position',[2.7 11.3 1.63 0.7])
% bottom plot
ax = axes
b2=bar(to_plot(:,2:3));
b2(1).FaceColor=cm2(2,:);
b2(2).FaceColor=cm2(3,:);
set(gca,'YLim',[0 10],'YTick',[0 5 10],'XTickLabel',{n_arch},'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',6,'Units','Centimeters','Position',[2.7 11-0.8 1.63 1])

%set(gcf,'renderer','Painters')
%print('../../figures/Fig2d','-depsc','-loose')

%% Parameter distribution in sets of different stability for of each locus architecture
var_par = [1 2 5 6 21 24 26];
var_log = [1 2 5 6 ];

% Define different stability categories:
thresh_no_mem = 10;
thresh_short_mem = [10 100];
thresh_long_mem = 100;

plot_names = {'k^{A}_{ini} (log10,h^{-1})','k^{B}_{ini} (log10,h^{-1})', 't^{A}_{OFF} (log10,h)', 't^{B}_{OFF} (log10,h)', 'p_{coll}', 'p_{PR}', 'L (kb)', 'k^{A}_{ini}/k^{B}_{ini} (log10)'};
plot_titles = {'Initíation rate A', 'Initíation rate B', 'Stability OFF state A', 'Stability OFF state B', 'Collision probability', 'Repression probability','Overlap Length',  'Initiation rate ratio'};


nbins=10;
f=figure(2); clf;
f.Position = [680   517   950   450];
 for i = 2:length(locus_arch)
    no_mem = find(min([data_log{1}(locus_arch{i},fsw1(j)),data_log{1}(locus_arch{i},fsw2(j))],[],2)<thresh_no_mem);
    short_mem = find(min([data_log{1}(locus_arch{i},fsw1(j)),data_log{1}(locus_arch{i},fsw2(j))],[],2)>thresh_short_mem(1) & ...
        min([data_log{1}(locus_arch{i},fsw1(j)),data_log{1}(locus_arch{i},fsw2(j))],[],2)<thresh_short_mem(2));
    long_mem = find(min([data_log{1}(locus_arch{i},fsw1(j)),data_log{1}(locus_arch{i},fsw2(j))],[],2)>thresh_long_mem);
    for m = 1:length(var_par)+1
        ax=axes;
        if m<length(var_par)+1
            temp1 = data_log{1}(locus_arch{i},var_par(m)); %all sets
            temp2 = data_log{1}(locus_arch{i}(no_mem),var_par(m)); %no memory
            temp3 = data_log{1}(locus_arch{i}(short_mem),var_par(m)); %short-term memory
            temp4 = data_log{1}(locus_arch{i}(long_mem),var_par(m)); %long-term memory
        else 
            temp1=data_log{1}(locus_arch{i},1)-data_log{1}(locus_arch{i},2);
            temp2=data_log{1}(locus_arch{i}(no_mem),1)-data_log{1}(locus_arch{i}(no_mem),2);
            temp3=data_log{1}(locus_arch{i}(short_mem),1)-data_log{1}(locus_arch{i}(short_mem),2);
            temp4=data_log{1}(locus_arch{i}(long_mem),1)-data_log{1}(locus_arch{i}(long_mem),2);
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
        set(gca,'YLim',[0 y_lim],'YTick',[],'YTickLabel',[], 'TickLength',[0.02 0],'TickDir','out','Linewidth',lw,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(m) pos_y(i) graph_size3]);
        if m==length(var_par)
            set(gca,'XTick', [0 250 500],'XTickLabel', [0 25 50]); % Covert L from 100bp to kb
        end
        if i==2
            title(plot_titles{m}, 'Fontsize',6);
            if m==1
                ylabel('Relative frequency','Fontsize',6);
            end
            
        else
            xlabel(plot_names{m}, 'Fontsize',6);
        end
        
    end
    %print(['../../figures/Fig2e'],'-depsc','-loose')
 end

