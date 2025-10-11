function [] = function_general_as_model_250301_mitosis(inputFilename,outputFilename)
fprintf('scan_par_queue %s -> %s ',inputFilename,outputFilename)

%inputFilename = '.txt';
par=dlmread(inputFilename);
output = [];
v = 1/1440;
t_start=0;
t_max=500;
output_time_step = 1;



nr_alleles = 100;

%Loop over all Parameter sets
for loop=1:size(par,1)
	% kinetic rates
	%k(1): Binding rate of Pol to A, k(2):Binding rate of Pol to B, k(3): Deg RNA A, k(4): Deg RNA B, k(5):average time of pA_AS in OFF state, k(6): average time of pB_AS in OFF state, 
	%k7: Production rate of protein A, k8: Production rate of protein B, k9: Deg rate protein A, k10: Deg rate protein B, k11: rate of release or termination of RNAPs from pA_AS,
	%k12: rate of release or termination of RNAPs from pB_AS, k13: probability for pA_AS bound RNAP to terminate, k14: probability for pB_AS bound RNAP to terminate
	%k15: basal on rate pA_AS, k16: basal on rate pB_AS, k17: basal off rate pA_AS, k18: basal off rate pB_AS
	k = par(loop, 1:19);
	% strength of TI repression / probability for pol removal
	p = par(loop, [21:30,20]);
	% binary parameters determining which repressive mechansims are present
	%b(1): Collisions, b(2): Occlusion, b(3): SDI, b(4): promoter repression
	overlap = par(loop,26);
	%length of A or B before overlap
	dima_prior = par(loop,27);
	dimb_prior = par(loop,28);
	%length of A or B after overlap
	dima_post = par(loop,29);
	dimb_post = par(loop,30);
	
	
	for alleles = 1:nr_alleles
	
		% initial conditions: A ON, B OFF
		A_pol = zeros(dima_prior+overlap+dima_post,1);
		B_pol = zeros(dimb_prior+overlap+dimb_post,1);
		pA_AS = 0;
		pA_b = 0;
		pB_b = 0;
		if (p(11)>0 & p(5)>0 & dimb_prior==0)
			pB_AS = 1;
		else
			pB_AS = 0;
		end
		% Pol occupancy on A
		%net initiation rate w/o influence of TI = 1: 1/(1/k(1)+1/k(11))?
		r = floor(length(A_pol)*v*1/(1/k(1)+1/k(11))*(1-k(13)));
		q = randsample(length(A_pol),r);
		A_pol(q) = 1;
		A_RNA = floor((1/(1/k(1)+1/k(11)))*(1-k(13))/k(3));
		B_RNA = 0;
		A_Pr = floor(k(7)*A_RNA/k(9));
		B_Pr = 0;
		
		const_par = [t_start, t_max, output_time_step, overlap, dima_prior, dimb_prior, dima_post, dimb_post];
		
		[t1,ap1,bp1,ar1,br1,apr1,bpr1,pAo1,pBo1] = reaction_general_AS_model_250301_mitosis(k, p, const_par, ...
		A_pol, B_pol, A_RNA, B_RNA, A_Pr, B_Pr, pA_AS, pB_AS, pA_b, pB_b);
		
		%ar_ini1 = ar1;
		%br_ini1 = br1;
		
		% initial conditions: B ON, A OFF
		A_pol = zeros(dima_prior+overlap+dima_post,1);
		B_pol = zeros(dimb_prior+overlap+dimb_post,1);
		if (p(4)>0 & p(5)>0 & dima_prior==0)
			pA_AS = 1;
		else 
			pA_AS = 0;
		end
		pB_AS = 0;
		% Pol occupancy on B
		r = floor(length(B_pol)*v*1/(1/k(2)+1/k(12))*(1-k(14)));
		q = randsample(length(B_pol),r);
		B_pol(q) = 1;
		A_RNA = 0;
		B_RNA = floor((1/(1/k(2)+1/k(12)))*(1-k(14))/k(4));
		A_Pr = 0;
		B_Pr = floor(k(8)*B_RNA/k(10));
		
		const_par = [t_start, t_max, output_time_step, overlap, dima_prior, dimb_prior, dima_post, dimb_post];		
		[t2,ap2,bp2,ar2,br2,apr2,bpr2,pAo2,pBo2] = reaction_general_AS_model_250301_mitosis(k, p, const_par, ...
		A_pol, B_pol, A_RNA, B_RNA, A_Pr, B_Pr, pA_AS, pB_AS, pA_b, pB_b);
		
		%ar_ini2 = ar2;
		%br_ini2 = br2;
		%bs = zeros(size(ar_ini1,1), nr_alleles);
		%bs(:,alleles) = ((ar_ini1>10)&(ar_ini2<10));
		output = [output; [k, p, ap1', bp1', ar1', br1', apr1', bpr1', pAo1', pBo1', ap2', bp2', ar2', br2', apr2', bpr2', pAo2', pBo2']];
	end %end of loop over alleles
	%bs_sum = sum(bs,2)'/nr_alleles;
	%frac_bs=mean(bs_sum(end-50:end));
	%output = [output; [k, p, frac_bs]];
	
end % end of loop over parameter sets

dlmwrite(outputFilename,output)
fprintf('scan_par_queue %s -> %s done\n',inputFilename,outputFilename)
end
	
		
