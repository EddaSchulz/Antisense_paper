function [out] = function_avg_alleles_200218(sim, nr_p, t)

[C, ia, ib] = unique(sim(:,1:nr_p), 'rows', 'stable');
nr_alleles = 100;
%nr_alleles = size(sim,1)./length(ia);
%t = (size(sim,2)-nr_par)./nr_x;
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
dima = sim(:,26)+sim(:,27)+sim(:,29);
dimb = sim(:,26)+sim(:,28)+sim(:,30);

%Calculate Ratio of pol on A vs B normalized to transcription strength
p_rat_ini1 = ((sim(:,ap1)+1)./((1./(1./sim(:,1)+1./sim(:,11))).*(1-sim(:,13))))./((sim(:,bp1)+1)./((1./(1./sim(:,2)+1./sim(:,12))).*(1-sim(:,14))));
p_rat_ini2 = ((sim(:,ap2)+1)./((1./(1./sim(:,1)+1./sim(:,11))).*(1-sim(:,13))))./((sim(:,bp2)+1)./((1./(1./sim(:,2)+1./sim(:,12))).*(1-sim(:,14))));
%Calculate Ratio of pol on A vs B normalized to transcription strength and to gene length
%p_rat_ini1 = ((sim(:,ap1)+1)./(dima)*(1/(1/sim(:,1)+1/sim(:,11))*(1-sim(:,13)))))./((sim(:,bp1)+1)./(dima)*(1/(1/sim(:,2)+1/sim(:,12))*(1-sim(:,14)))));
%p_rat_ini2 = ((sim(:,ap2)+1)./(dimb)*(1/(1/sim(:,1)+1/sim(:,11))*(1-sim(:,13)))))./((sim(:,bp2)+1)./(dimb)*(1/(1/sim(:,2)+1/sim(:,12))*(1-sim(:,14)))));
%Calculate Ratio of RNA on A vs B normalized to transcription strength
r_rat_ini1 = ((sim(:,ar1)+1)./((1./(1./sim(:,1)+1./sim(:,11))).*(1-sim(:,13))))./((sim(:,br1)+1)./(1./(1./sim(:,2)+1./sim(:,12)).*(1-sim(:,14))));
r_rat_ini2 = ((sim(:,ar2)+1)./((1./(1./sim(:,1)+1./sim(:,11))).*(1-sim(:,13))))./((sim(:,br2)+1)./(1./(1./sim(:,2)+1./sim(:,12)).*(1-sim(:,14))));
%Calculate Ratio of protein A vs B normalized to transcription strength
pr_rat_ini1 = ((sim(:,apr1)+1)./((1./(1./sim(:,1)+1./sim(:,11))).*(1-sim(:,13))))./((sim(:,bpr1)+1)./(1./(1./sim(:,2)+1./sim(:,12)).*(1-sim(:,14))));
pr_rat_ini2 = ((sim(:,apr2)+1)./((1./(1./sim(:,1)+1./sim(:,11))).*(1-sim(:,13))))./((sim(:,bpr2)+1)./(1./(1./sim(:,2)+1./sim(:,12)).*(1-sim(:,14))));

%Binary definition on Pol Level: ON or OFF: This is a problem when steady state pol level/5 <1 (when k_ini very low)
sim_bi = zeros(size(sim,1), size(sim,2));
%sim_bi(sim>=(1/1440*dima.*sim(:,1)./5))= 1;
sim_bi(sim>=(1/1440.*dima.*(1./(1./sim(:,1)+1./sim(:,11)).*(1-sim(:,13)))./3) & sim>=3)= 1;
% First Off switch of A - Ini1
[r c] = find(sim_bi(:,ap1)==0);
firstIndex_a1 = accumarray(r,c,[size(sim_bi(:,ap1),1),1],@min,t);
% First ON switch of A - Ini2
[r c] = find(sim_bi(:,ap2)==1);
firstIndex_a2 = accumarray(r,c,[size(sim_bi(:,ap2),1),1],@min,t);

sim_bi = zeros(size(sim,1), size(sim,2));
%sim_bi(sim>=(1/1440*dimb.*sim(:,2)./5))= 1;
sim_bi(sim>=(1/1440.*dimb.*sim(:,2)./3) & sim>=3)= 1;
% First ON switch of B - Ini1
[r c] = find(sim_bi(:,bp1)==1);
firstIndex_b1 = accumarray(r,c,[size(sim_bi(:,bp1),1),1],@min,t);
% First Off switch of B - Ini2
[r c] = find(sim_bi(:,bp2)==0);
firstIndex_b2 = accumarray(r,c,[size(sim_bi(:,bp2),1),1],@min,t);
sw_time = [firstIndex_a1 firstIndex_b1 firstIndex_a2 firstIndex_b2];


% Ratio 1 time step delayed 
p_rat_1_1step_del = [10*ones(size(p_rat_ini1,1),1) p_rat_ini1(:,1:end-1)];
r_rat_1_1step_del = [10*ones(size(r_rat_ini1,1),1) r_rat_ini1(:,1:end-1)];
pr_rat_1_1step_del = [10*ones(size(pr_rat_ini1,1),1) pr_rat_ini1(:,1:end-1)];
% Ratio 2 time step delayed 
p_rat_2_1step_del = [0.1*ones(size(p_rat_ini2,1),1) p_rat_ini2(:,1:end-1)];
r_rat_2_1step_del = [0.1*ones(size(r_rat_ini2,1),1) r_rat_ini2(:,1:end-1)];
pr_rat_2_1step_del = [0.1*ones(size(pr_rat_ini2,1),1) pr_rat_ini2(:,1:end-1)];

% Median of RNAP, RNA or Protein Ratio over all alleles for each timepoint for each parameter set
p_med_1 = NaN(length(ia),t);
r_med_1 = NaN(length(ia),t);
pr_med_1 = NaN(length(ia),t);
p_med_2 = NaN(length(ia),t);
r_med_2 = NaN(length(ia),t);
pr_med_2 = NaN(length(ia),t);
% Distribution of RNAPs, RNA or Protein during last 50h on all alleles
nbins = 100;
rnap_a_1_counts = NaN(length(ia), nbins); rnap_a_1_centers = NaN(length(ia), nbins);
rnap_b_1_counts = NaN(length(ia), nbins); rnap_b_1_centers = NaN(length(ia), nbins);
rnap_a_2_counts = NaN(length(ia), nbins); rnap_a_2_centers = NaN(length(ia), nbins);
rnap_b_2_counts = NaN(length(ia), nbins); rnap_b_2_centers = NaN(length(ia), nbins);
rna_a_1_counts = NaN(length(ia), nbins); rna_a_1_centers = NaN(length(ia), nbins);
rna_b_1_counts = NaN(length(ia), nbins); rna_b_1_centers = NaN(length(ia), nbins);
rna_a_2_counts = NaN(length(ia), nbins); rna_a_2_centers = NaN(length(ia), nbins);
rna_b_2_counts = NaN(length(ia), nbins); rna_b_2_centers = NaN(length(ia), nbins);
pr_a_1_counts = NaN(length(ia), nbins); pr_a_1_centers = NaN(length(ia), nbins);
pr_b_1_counts = NaN(length(ia), nbins); pr_b_1_centers = NaN(length(ia), nbins);
pr_a_2_counts = NaN(length(ia), nbins); pr_a_2_centers = NaN(length(ia), nbins);
pr_b_2_counts = NaN(length(ia), nbins); pr_b_2_centers = NaN(length(ia), nbins);
mean_sw = NaN(length(ia),4);
avg_fsw1 = t*ones(length(ia),3);
avg_fsw2 = t*ones(length(ia),3);
avg_sw1_mean = NaN(length(ia),3);
avg_sw2_mean = NaN(length(ia),3);
avg_sw1_std = NaN(length(ia),3);
avg_sw2_std = NaN(length(ia),3);
avg_sw1_min = t*ones(length(ia),3);
avg_sw1_max = t*ones(length(ia),3);
avg_sw2_min = t*ones(length(ia),3);
avg_sw2_max = t*ones(length(ia),3);
avg_sw1_n = zeros(length(ia),3);
avg_sw2_n = zeros(length(ia),3);

%Loop over Parameter sets
for i = 1:length(ia)
	%Loop over Timepoints
    for j = 1:t
		% Median of Pol Ratio for each parameter set
		p_med_1(i,j) = median(p_rat_ini1((i-1)*nr_alleles+1:i*nr_alleles,j));
		r_med_1(i,j) = median(r_rat_ini1((i-1)*nr_alleles+1:i*nr_alleles,j));
		pr_med_1(i,j) = median(pr_rat_ini1((i-1)*nr_alleles+1:i*nr_alleles,j));
		p_med_2(i,j) = median(p_rat_ini2((i-1)*nr_alleles+1:i*nr_alleles,j));
		r_med_2(i,j) = median(r_rat_ini2((i-1)*nr_alleles+1:i*nr_alleles,j));
		pr_med_2(i,j) = median(pr_rat_ini2((i-1)*nr_alleles+1:i*nr_alleles,j));
		if j<5
			%Mean switch times
			mean_sw(i,j) = mean(sw_time((i-1)*nr_alleles+1:i*nr_alleles,j));
		end
    end

    [counts, centers] = hist(reshape(sim((i-1)*nr_alleles+1:i*nr_alleles,ap1(end-50:end)), [length([(i-1)*nr_alleles+1:i*nr_alleles])*length([0:50]) 1]), nbins); 
    rnap_a_1_counts(i,:) = counts; rnap_a_1_centers(i,:) = centers;
    [counts, centers] = hist(reshape(sim((i-1)*nr_alleles+1:i*nr_alleles,bp1(end-50:end)), [length([(i-1)*nr_alleles+1:i*nr_alleles])*length([0:50]) 1]), nbins); 
    rnap_b_1_counts(i,:) = counts; rnap_b_1_centers(i,:) = centers;
    [counts, centers] = hist(reshape(sim((i-1)*nr_alleles+1:i*nr_alleles,ar1(end-50:end)), [length([(i-1)*nr_alleles+1:i*nr_alleles])*length([0:50]) 1]), nbins); 
    rna_a_1_counts(i,:) = counts; rna_a_1_centers(i,:) = centers;
    [counts, centers] = hist(reshape(sim((i-1)*nr_alleles+1:i*nr_alleles,br1(end-50:end)), [length([(i-1)*nr_alleles+1:i*nr_alleles])*length([0:50]) 1]), nbins); 
    rna_b_1_counts(i,:) = counts; rna_b_1_centers(i,:) = centers;
    [counts, centers] = hist(reshape(sim((i-1)*nr_alleles+1:i*nr_alleles,apr1(end-50:end)), [length([(i-1)*nr_alleles+1:i*nr_alleles])*length([0:50]) 1]), nbins); 
    pr_a_1_counts(i,:) = counts; pr_a_1_centers(i,:) = centers;
    [counts, centers] = hist(reshape(sim((i-1)*nr_alleles+1:i*nr_alleles,bpr1(end-50:end)), [length([(i-1)*nr_alleles+1:i*nr_alleles])*length([0:50]) 1]), nbins); 
    pr_b_1_counts(i,:) = counts; pr_b_1_centers(i,:) = centers;
    
    [counts, centers] = hist(reshape(sim((i-1)*nr_alleles+1:i*nr_alleles,ap2(end-50:end)), [length([(i-1)*nr_alleles+1:i*nr_alleles])*length([0:50]) 1]), nbins); 
    rnap_a_2_counts(i,:) = counts; rnap_a_2_centers(i,:) = centers;
    [counts, centers] = hist(reshape(sim((i-1)*nr_alleles+1:i*nr_alleles,bp2(end-50:end)), [length([(i-1)*nr_alleles+1:i*nr_alleles])*length([0:50]) 1]), nbins); 
    rnap_b_2_counts(i,:) = counts; rnap_b_2_centers(i,:) = centers;
    [counts, centers] = hist(reshape(sim((i-1)*nr_alleles+1:i*nr_alleles,ar2(end-50:end)), [length([(i-1)*nr_alleles+1:i*nr_alleles])*length([0:50]) 1]), nbins); 
    rna_a_2_counts(i,:) = counts; rna_a_2_centers(i,:) = centers;
    [counts, centers] = hist(reshape(sim((i-1)*nr_alleles+1:i*nr_alleles,br2(end-50:end)), [length([(i-1)*nr_alleles+1:i*nr_alleles])*length([0:50]) 1]), nbins); 
    rna_b_2_counts(i,:) = counts; rna_b_2_centers(i,:) = centers;
    [counts, centers] = hist(reshape(sim((i-1)*nr_alleles+1:i*nr_alleles,apr2(end-50:end)), [length([(i-1)*nr_alleles+1:i*nr_alleles])*length([0:50]) 1]), nbins); 
    pr_a_2_counts(i,:) = counts; pr_a_2_centers(i,:) = centers;
    [counts, centers] = hist(reshape(sim((i-1)*nr_alleles+1:i*nr_alleles,bpr2(end-50:end)), [length([(i-1)*nr_alleles+1:i*nr_alleles])*length([0:50]) 1]), nbins); 
    pr_b_2_counts(i,:) = counts; pr_b_2_centers(i,:) = centers;
    
	fsw1 = t*ones(nr_alleles,3);
	fsw2 = t*ones(nr_alleles,3);
	sw1_mean = t*ones(nr_alleles,3);
	sw2_mean = t*ones(nr_alleles,3);
	sw1_std = zeros(nr_alleles,3);
	sw2_std = zeros(nr_alleles,3);
	sw1_min = t*ones(nr_alleles,3);
	sw1_max = t*ones(nr_alleles,3);
	sw2_min = t*ones(nr_alleles,3);
	sw2_max = t*ones(nr_alleles,3);
	sw1_n = zeros(nr_alleles,3);
	sw2_n = zeros(nr_alleles,3);
    for h = 1:nr_alleles 
		clear sw1_ini1_p sw1_ini1_r sw2_ini1_p sw2_ini1_r sw1_ini2_p sw1_ini2_r sw2_ini2_p sw2_ini2_r;
        % switches direction 1 ini 1 
		sw1_ini1_p = find(p_rat_ini1((i-1)*nr_alleles+h,:)<=1 & p_rat_1_1step_del((i-1)*nr_alleles+h,:)>1);
		sw1_ini1_r = find(r_rat_ini1((i-1)*nr_alleles+h,:)<=1 & r_rat_1_1step_del((i-1)*nr_alleles+h,:)>1);
		sw1_ini1_pr = find(pr_rat_ini1((i-1)*nr_alleles+h,:)<=1 & pr_rat_1_1step_del((i-1)*nr_alleles+h,:)>1);
		% switches direction 2 ini1
		sw2_ini1_p = find(p_rat_ini1((i-1)*nr_alleles+h,:)>1 & p_rat_1_1step_del((i-1)*nr_alleles+h,:)<=1);
		sw2_ini1_r = find(r_rat_ini1((i-1)*nr_alleles+h,:)>1 & r_rat_1_1step_del((i-1)*nr_alleles+h,:)<=1);
		sw2_ini1_pr = find(pr_rat_ini1((i-1)*nr_alleles+h,:)>1 & pr_rat_1_1step_del((i-1)*nr_alleles+h,:)<=1);
		%switches direction 1 ini2
		sw1_ini2_p = find(p_rat_ini2((i-1)*nr_alleles+h,:)>=1 & p_rat_2_1step_del((i-1)*nr_alleles+h,:)<1);
		sw1_ini2_r = find(r_rat_ini2((i-1)*nr_alleles+h,:)>=1 & r_rat_2_1step_del((i-1)*nr_alleles+h,:)<1);
		sw1_ini2_pr = find(pr_rat_ini2((i-1)*nr_alleles+h,:)>=1 & pr_rat_2_1step_del((i-1)*nr_alleles+h,:)<1);
		%switches direction 2 ini2
		sw2_ini2_p = find(p_rat_ini2((i-1)*nr_alleles+h,:)<1 & p_rat_2_1step_del((i-1)*nr_alleles+h,:)>=1);
		sw2_ini2_r = find(r_rat_ini2((i-1)*nr_alleles+h,:)<1 & r_rat_2_1step_del((i-1)*nr_alleles+h,:)>=1);
		sw2_ini2_pr = find(pr_rat_ini2((i-1)*nr_alleles+h,:)<1 & pr_rat_2_1step_del((i-1)*nr_alleles+h,:)>=1);
		if length(sw1_ini1_p) > 0
			sw2_ini1_p = [0 sw2_ini1_p];
			if length(sw2_ini1_p)>length(sw1_ini1_p)
				% No, because if gene stays on last ON period not yetover, should not be included in calc
				%sw1_ini1_p = [sw1_ini1_p t];
				sw2_ini1_p = sw2_ini1_p(1:end-1);
			end
		    % First switch
		    fsw1(h,1) = sw1_ini1_p(1);
		    % Mean of switching periods
		    sw1_mean(h,1) = mean(sw1_ini1_p - sw2_ini1_p);
		    % Std of switching periods
		    sw1_std(h,1) = std(sw1_ini1_p - sw2_ini1_p);
		    sw1_min(h,1) = min(sw1_ini1_p - sw2_ini1_p);
		    sw1_max(h,1) = max(sw1_ini1_p - sw2_ini1_p);
		    %sw1_std(i) = std(sw1 - [0 sw1(1:end-1)]);
		    % # of switches
		    sw1_n(h,1) = length(sw1_ini1_p);    
	    end
	    if length(sw1_ini1_r) > 0
			sw2_ini1_r = [0 sw2_ini1_r];
			if length(sw2_ini1_r)>length(sw1_ini1_r)
				sw2_ini1_r = sw2_ini1_r(1:end-1);
			end
		    fsw1(h,2) = sw1_ini1_r(1);
		    sw1_mean(h,2) = mean(sw1_ini1_r - sw2_ini1_r);
		    sw1_std(h,2) = std(sw1_ini1_r - sw2_ini1_r);
		    sw1_min(h,2) = min(sw1_ini1_r - sw2_ini1_r);
		    sw1_max(h,2) = max(sw1_ini1_r - sw2_ini1_r);
		    sw1_n(h,2) = length(sw1_ini1_r);    
	    end
	    if length(sw1_ini1_pr) > 0
			sw2_ini1_pr = [0 sw2_ini1_pr];
			if length(sw2_ini1_pr)>length(sw1_ini1_pr)
				sw2_ini1_pr = sw2_ini1_pr(1:end-1);
			end
		    fsw1(h,3) = sw1_ini1_pr(1);
		    sw1_mean(h,3) = mean(sw1_ini1_pr - sw2_ini1_pr);
		    sw1_std(h,3) = std(sw1_ini1_pr - sw2_ini1_pr);
		    sw1_min(h,3) = min(sw1_ini1_pr - sw2_ini1_pr);
		    sw1_max(h,3) = max(sw1_ini1_pr - sw2_ini1_pr);
		    sw1_n(h,3) = length(sw1_ini1_pr);    
	    end
		if length(sw1_ini2_p) > 0
			sw2_ini2_p = [0 sw2_ini2_p];
			if length(sw2_ini2_p)>length(sw1_ini2_p)
				sw2_ini2_p = sw2_ini2_p(1:end-1);
			end
			fsw2(h,1) = sw1_ini2_p(1);
			sw2_mean(h,1) = mean(sw1_ini2_p - sw2_ini2_p);
			sw2_std(h,1) = std(sw1_ini2_p - sw2_ini2_p);
			sw2_min(h,1) = min(sw1_ini2_p - sw2_ini2_p);
			sw2_max(h,1) = max(sw1_ini2_p - sw2_ini2_p);
			sw2_n(h,1) = length(sw1_ini2_p);
	    end
	    if length(sw1_ini2_r) > 0
			sw2_ini2_r = [0 sw2_ini2_r];
			if length(sw2_ini2_r)>length(sw1_ini2_r)
				sw2_ini2_r = sw2_ini2_r(1:end-1);
			end
			fsw2(h,2) = sw1_ini2_r(1);
			sw2_mean(h,2) = mean(sw1_ini2_r - sw2_ini2_r);
			sw2_std(h,2) = std(sw1_ini2_r - sw2_ini2_r);
			sw2_min(h,2) = min(sw1_ini2_r - sw2_ini2_r);
			sw2_max(h,2) = max(sw1_ini2_r - sw2_ini2_r);
			sw2_n(h,2) = length(sw1_ini2_r);
	    end
	    if length(sw1_ini2_pr) > 0
			sw2_ini2_pr = [0 sw2_ini2_pr];
			if length(sw2_ini2_pr)>length(sw1_ini2_pr)
				sw2_ini2_pr = sw2_ini2_pr(1:end-1);
			end
			fsw2(h,3) = sw1_ini2_pr(1);
			sw2_mean(h,3) = mean(sw1_ini2_pr - sw2_ini2_pr);
			sw2_std(h,3) = std(sw1_ini2_pr - sw2_ini2_pr);
			sw2_min(h,3) = min(sw1_ini2_pr - sw2_ini2_pr);
			sw2_max(h,3) = max(sw1_ini2_pr - sw2_ini2_pr);
			sw2_n(h,3) = length(sw1_ini2_pr);
	    end
    end
    avg_fsw1(i,:) = mean(fsw1,1);
    avg_fsw2(i,:) = mean(fsw2,1);
    avg_sw1_mean(i,:) = mean(sw1_mean,1);
    avg_sw2_mean(i,:) = mean(sw2_mean,1);
    avg_sw1_std(i,:) = mean(sw1_std,1);
    avg_sw2_std(i,:) = mean(sw2_std,1);
    %min_sw1_min(i,:) = min(sw1_min,[],1);
    avg_sw1_min(i,:) = mean(sw1_min,1);
    avg_sw2_min(i,:) = mean(sw2_min,1);
    avg_sw1_max(i,:) = mean(sw1_max,1);
    avg_sw2_max(i,:) = mean(sw2_max,1);
    avg_sw1_n(i,:) = mean(sw1_n,1);
    avg_sw2_n(i,:) = mean(sw2_n,1); 
end
distr_rnap_1 = [rnap_a_1_counts rnap_a_1_centers rnap_b_1_counts rnap_b_1_centers];
distr_rna_1 = [rna_a_1_counts rna_a_1_centers rna_b_1_counts rna_b_1_centers];
distr_pr_1 = [pr_a_1_counts pr_a_1_centers pr_b_1_counts pr_b_1_centers];

distr_rnap_2 = [rnap_a_2_counts rnap_a_2_centers rnap_b_2_counts rnap_b_2_centers];
distr_rna_2 = [rna_a_2_counts rna_a_2_centers rna_b_2_counts rna_b_2_centers];
distr_pr_2 = [pr_a_2_counts pr_a_2_centers pr_b_2_counts pr_b_2_centers];
out = [sim(ia,1:nr_p) distr_rnap_1 distr_rna_1 distr_pr_1 distr_rnap_2 distr_rna_2 distr_pr_2 p_med_1 r_med_1 pr_med_1 p_med_2 r_med_2 pr_med_2 avg_fsw1 avg_fsw2 avg_sw1_mean avg_sw2_mean avg_sw1_std avg_sw2_std avg_sw1_min avg_sw2_min avg_sw1_max avg_sw2_max avg_sw1_n avg_sw2_n mean_sw];

end
