#include <iostream>
#include <fstream>
#include <ctime>            
#include <boost/random/mersenne_twister.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <sys/time.h>

#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <queue>

#include <math.h>
#include <matrix.h>
#include <mex.h>

boost::random::mt19937 gen;

using namespace std;

int collisions (int &overlap, int &dima_prior, int &dimb_post, deque<int> &a_gene, deque<int> &b_gene, 
        int &sum_ap, int &sum_bp, double *p, boost::random::uniform_real_distribution<> &dist_real)
{
    int q;
    for (q=0; q<overlap; q++)
            {
            if ((a_gene[q+dima_prior] + b_gene[q+dimb_post])>1) {
            	if (dist_real(gen)<p[0]) {
				// for independent probability of each RNAP to dislode 
                //if (dist_real(gen)<=p[0]) {
                //    a_gene[q+dima_prior]=0;
                //    sum_ap--;
                //} 
                //if (dist_real(gen)<=p[0]) {
                //    b_gene[q+dimb_post]=0;
                //    sum_bp--;
                //}
                
                // for one randomly chosen RNAP to fall off
                if (dist_real(gen)<0.5) {
                    a_gene[q+dima_prior]=0;
                    sum_ap--;
                } 
                else {
                    b_gene[q+dimb_post]=0;
                    sum_bp--;
                }
                }
            }
        }
}


int SDI (bool &apol_paused, bool &bpol_paused, int &overlap, int &dima_prior, int &dimb_prior, int &dimb_post, deque<int> &a_gene, deque<int> &b_gene, 
        int &sum_ap, int &sum_bp, double *p, boost::random::uniform_real_distribution<> &dist_real) {
// if pA is inside overlap
if (dima_prior==0) {
	if (apol_paused+b_gene[dimb_post]>1){
		if (dist_real(gen)<p[2]) {
			apol_paused=0;
			//sum_ap--;
		} //else {
			//b_gene[dimb_post]=0;
			//sum_bp--;
		//}
	}
}
// if pB is inside overlap
if (dimb_prior==0) {
	if (bpol_paused+a_gene[dima_prior+overlap-1]>1) {
		if (dist_real(gen)<p[2]) {
			bpol_paused=0;
			//sum_bp--;
		} //else {
			//a_gene[dima_prior+overlap-1]=0;
			//sum_ap--;
		//}
	}
}
}
int prom_rep (int & overlap, int &dima_prior, int &dimb_prior, int &dimb_post, deque<int> &a_gene, deque<int> &b_gene, 
        int &a_prom_as, int &b_prom_as, int &a_prom_b, int &b_prom_b, double *p, boost::random::uniform_real_distribution<> &dist_real) {
	if (dima_prior==0) {
		if (b_gene[dimb_post]==1){
			//cout << "a_prom_as = " << a_prom_as << endl;
			if (a_prom_as==0) { //we could also leave this loop out, then promoter can be repressed again by passing RNAP even if it has not been completely reactivated yet
				if (dist_real(gen)<p[3]) {
					a_prom_as=p[4];
					a_prom_b=0;
					
				}
			}
		}
	}
	if (dimb_prior==0) {
		if (a_gene[dima_prior+overlap-1]==1) {
			//cout << "b_prom_as = " << b_prom_as << "  ";
			if (b_prom_as==0) {
				if (dist_real(gen)<p[10]) {
					b_prom_as=p[4];
					b_prom_b=0;
				}
			}
		}
	}
}

int elongation (bool &apol_paused, bool &bpol_paused, int &overlap, int &dima_prior, int &dimb_prior, int &dimb_post, deque<int> &a_gene, deque<int> &b_gene, int &sum_ap, int &sum_bp, 
        bool &apol_ini, bool &bpol_ini, int &a_rna, int &b_rna, 
        vector<int> &all_a, int &temp_sum_ar, int &temp_sum_br, int &temp_sum_ap, int &temp_sum_bp, int &temp_sum_pA, int &temp_sum_pB, int &a_prom_as, int &b_prom_as, int &a_prom_b, int &b_prom_b, 
        boost::random::uniform_real_distribution<> &dist_real, int &step, double *p)
{
 	// Pol on A move forward
    a_gene.push_front(apol_ini);
    sum_ap = sum_ap - a_gene.back() + a_gene.front();
    a_rna = a_rna + (int)a_gene.back();
    a_gene.pop_back();
    //cout << "a_prom_as before 1st prom rep= " << a_prom_as << endl;
    if (p[2]>0)  {
		SDI(apol_paused, bpol_paused, overlap,dima_prior,dimb_prior,dimb_post,a_gene,b_gene,sum_ap,sum_bp, p, dist_real);
	}
    if (p[0]>0) {
        collisions (overlap, dima_prior, dimb_post, a_gene, b_gene, sum_ap, sum_bp, p, dist_real);
    }   
    if (p[10]>0) {
		prom_rep(overlap, dima_prior, dimb_prior, dimb_post, a_gene, b_gene, a_prom_as, b_prom_as, a_prom_b, b_prom_b, p, dist_real);
	}
	//cout << "a_prom_as after 1st prom rep= " << a_prom_as << endl;
    // Pol on B move backward
    b_gene.push_back(bpol_ini);
    sum_bp = sum_bp - b_gene.front() + b_gene.back();
    b_rna = b_rna + (int)b_gene.front();
    b_gene.pop_front();
    //cout << "a_prom_as before 2nd prom rep= " << a_prom_as << endl;
    if (p[2]>0)  {
		SDI(apol_paused, bpol_paused, overlap,dima_prior,dimb_prior,dimb_post,a_gene,b_gene,sum_ap,sum_bp, p, dist_real);
	}
    if (p[0]>0) {
        collisions (overlap, dima_prior, dimb_post, a_gene, b_gene, sum_ap, sum_bp, p, dist_real);
    }    
    if (p[3]>0) {
		prom_rep(overlap, dima_prior, dimb_prior, dimb_post, a_gene, b_gene, a_prom_as, b_prom_as, a_prom_b, b_prom_b, p, dist_real);
	}        
	//cout << "a_prom_as after 2nd prom rep= " << a_prom_as << endl;  
    all_a[step] = a_rna;
                      
    temp_sum_ar = temp_sum_ar + a_rna;
    temp_sum_br = temp_sum_br + b_rna;
    temp_sum_ap = temp_sum_ap + sum_ap;
    temp_sum_bp = temp_sum_bp + sum_bp;
    temp_sum_pA = temp_sum_pA + a_prom_as;
    temp_sum_pB = temp_sum_pB + b_prom_as;
    //cout << "temp_sum_pB in elong funct= " << temp_sum_pB << endl;  
}

void mexFunction(int nhls, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    //initiate random number generator
    struct timeval t1;
    gettimeofday(&t1, NULL);
    //use time as seed for random number generator => each round of simulation with different seed => will generate different random numbers
    //If same seed was used in different rounds of simulation, same 'random' numbers would be generated in each round
    gen.seed(static_cast<unsigned int>(t1.tv_usec * t1.tv_sec));
    
	//declare variables // pointers that point at input and output variables
	mxArray *a_pol, *b_pol, *par_k, *par_p, *const_par; 
    mxArray *time_out, *a_pol_out, *b_pol_out, *a_rna_out, *b_rna_out, *a_pr_out, *b_pr_out, *a_prom_out, *b_prom_out; 
	const mwSize *dims;
	double *ap, *bp, *k, *p, *cp, 
            t, t_max;
            //*ar, *br,*apr, *apro, *bpro, 
    double *apo, *bpo, *aro, *bro, *apro, *bpro, *pAo, *pBo,
            *time, v = 1.0/1440.0, out_step;
            //v is the time needed for elongating 1 segment of 100nt in [h]: v = 100nt/(40bp/sec*60*60)
	int dim_ap, dim_bp, a_rna, b_rna, a_pr, b_pr, a_prom_as, b_prom_as, a_prom_b, b_prom_b, a_prom_mit, b_prom_mit, overlap, dima_prior, dimb_prior, dima_post, dimb_post, off_states;
	int n_steps, n_out_steps, n_small_steps, step, sel_rx;
	
	//associate inputs
	//k, p, b, const_par, A_pol, B_pol, A_RNA, B_RNA, pA, pB
    par_k = mxDuplicateArray(prhs[0]);
    par_p = mxDuplicateArray(prhs[1]);
    const_par = mxDuplicateArray(prhs[2]);
	a_pol = mxDuplicateArray(prhs[3]);
    b_pol = mxDuplicateArray(prhs[4]);  
    a_rna = mxGetScalar(prhs[5]);
    b_rna = mxGetScalar(prhs[6]);
    a_pr = mxGetScalar(prhs[7]);
    b_pr = mxGetScalar(prhs[8]);
    //ON:a_prom_as=0, Off1:a_prom_as=1,..., Offn: a_prom_as=n
    a_prom_as = mxGetScalar(prhs[9]); // AS txn induced OFF
    b_prom_as = mxGetScalar(prhs[10]);
    a_prom_b = mxGetScalar(prhs[11]); // basal OFF 
    b_prom_b = mxGetScalar(prhs[12]);
    a_prom_mit = 0; //MItosis OFF
    b_prom_mit = 0;
              
	//associate pointers
    ap = mxGetPr(a_pol);
    bp = mxGetPr(b_pol);
    k = mxGetPr(par_k);
    p = mxGetPr(par_p);
    cp = mxGetPr(const_par);
    // collisions present?
    //b0 = b[0];
    // SDI present?
    //b2 = b[2];
    
    //const_par = [t_start, t_max, output_time_step, overlap, dima_prior, dimb_prior, dima_post, dimb_post]
    t = cp[0];
    t_max =  cp[1];
    out_step = cp[2];
    overlap = cp[3];
    dima_prior = cp[4]; 
    dimb_prior = cp[5]; 
    dima_post = cp[6]; 
    dimb_post = cp[7]; 
    
    // output state of system every 1h
    n_out_steps = 1+ceil(t_max/out_step);
    // each output step (=1h) consists of this many elongation steps
    n_small_steps = round(out_step/v);
    // total # of elongation steps that will occur during total simulation time
    n_steps = (n_out_steps-1)*n_small_steps;
	
	//figure out dimensions and associate output
	dims = mxGetDimensions(prhs[3]);
	dim_ap = (int)dims[0];
    dims = mxGetDimensions(prhs[4]);
	dim_bp = (int)dims[0];
    
    //output
    // [t,ap,bp,ar,br, test]
    time_out = plhs[0] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    a_pol_out = plhs[1] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    b_pol_out = plhs[2] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    a_rna_out = plhs[3] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    b_rna_out = plhs[4] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    a_pr_out = plhs[5] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    b_pr_out = plhs[6] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    a_prom_out = plhs[7] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    b_prom_out = plhs[8] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    apo = mxGetPr(a_pol_out);
    bpo= mxGetPr(b_pol_out);
    aro = mxGetPr(a_rna_out);
    bro = mxGetPr(b_rna_out);
    apro = mxGetPr(a_pr_out);
    bpro = mxGetPr(b_pr_out);
    pAo= mxGetPr(a_prom_out);
    pBo= mxGetPr(b_prom_out);
    time = mxGetPr(time_out);
    //test = plhs[5] = mxCreateDoubleMatrix(n_out_steps, 1, mxREAL);
    //to = mxGetPr(test);
    
    //write pol mxarray in deque and sum up polymerases
    deque<int> a_gene, b_gene;
    int q, sum_ap=0, sum_bp=0;
    
    for (q=0; q<dim_ap; q++) 
        {
            a_gene.push_back (ap[q]);
            sum_ap = sum_ap + ap[q];
        }
    for (q=0; q<dim_bp; q++) 
        {
            b_gene.push_back (bp[q]);
            sum_bp = sum_bp + bp[q];
        }
        
	bool apol_ini=0, bpol_ini=0, apol_paused=0, bpol_paused=0;
	boost::random::uniform_real_distribution<> dist_real(0.0, 1.0);
	
    double t_next, r, sum_rx;
    vector<double> rx(16), cum_rx(16);
    vector<int> all_a(n_steps);
    
    int m = 0, mit_count = 0, l = 0, temp_sum_ar = 0, temp_sum_br = 0, temp_sum_apr = 0, temp_sum_bpr = 0, temp_sum_ap = 0, temp_sum_bp = 0, temp_sum_pA = 0, temp_sum_pB = 0;
    all_a[0] = a_rna;
    aro[l] = a_rna;
    bro[l] = b_rna;
    apro[l] = a_pr;
    bpro[l] = b_pr;
    apo[l] = sum_ap;
    bpo[l] = sum_bp;
    pAo[l] = a_prom_as;
    pBo[l] = b_prom_as;
    
    
    //while (act_time<t_max) -> while less than n_steps elongation events have occured the simulation time is not over
    for (step=0; step<n_steps; step++)
        {
			//if avg time in AS induced off state =0, AS promoter state is always ON
			if (k[4]==0) {a_prom_as=0;}
			if (k[5]==0) {b_prom_as=0;}	
			//if avg time in basal off state =0, basal promoter state is always ON
			if (k[14]==0) {a_prom_b=0;}
			if (k[15]==0) {b_prom_b=0;}			
		////Gillespie reactions
        t_next = 0;
        
        apol_ini=0, bpol_ini=0;  
        // If t_next < next elongation step: perform stochastic Gillespie step
        while (t_next<v) {
           //?? rx[0] = p[0]*(a_prom_as==0)*(b_gene.front()==0); // initiation event on Gene A
           //?? indexing?? rx[1] = p[1]*(b_prom_as==0)*(a_gene.back()==0); ; // initiation event on Gene B
           //RNAP goes from unbound -> paused -> initiation ready -> elongated; while it is in paused state it can spontaneously fall off
            if (dima_prior==0 & dist_real(gen)<p[1]) {
            	rx[0] = k[0]*((a_prom_as+a_prom_b+a_prom_mit)==0)*(apol_ini==0)*(apol_paused==0)*(b_gene[dimb_post]==0);
            }
            else {
            	rx[0] = k[0]*((a_prom_as+a_prom_b+a_prom_mit)==0)*(apol_ini==0)*(apol_paused==0);
            }
            if (dimb_prior==0 & dist_real(gen)<p[1]) {
            	rx[1] = k[1]*((b_prom_as+b_prom_b+b_prom_mit)==0)*(bpol_ini==0)*(bpol_paused==0)*(a_gene[dima_prior+overlap-1]==0);
            } else {
            	rx[1] = k[1]*((b_prom_as+b_prom_b+b_prom_mit)==0)*(bpol_ini==0)*(bpol_paused==0);
            }
            rx[2] = k[2]*a_rna; //A RNA degradation
            rx[3] = k[3]*b_rna; //B RNA degradation
            //if (p[4]>0 & a_prom_as>0 & b_gene[dimb_post]==0) {
			if (p[4]>0 & a_prom_as>0) {
				rx[4] = p[4]/k[4]*(a_prom_as>0)*(b_gene[dimb_post]==0); //Reactivation of pA from Off_i to Off_i+1 or from OFFn back to ON (reversal rate = number OFF states/ average time in off state)
				//cout << "rx prom derep =" << rx[4];
			} else {
				rx[4] = 0;
			}
			if (p[4]>0 & b_prom_as>0) {
				rx[5] = p[4]/k[5]*(b_prom_as>0)*(a_gene[dima_prior+overlap-1]==0); //Reactivation of pB from Off_i to Off_i+1
			} else {
				rx[5] = 0;
			}
			// Protein production
			rx[6] =  k[6]*a_rna;
			rx[7] =  k[7]*b_rna;
			// Protein degradation
			rx[8] = k[8]*a_pr;
			rx[9] = k[9]*b_pr;
			// RNAP pause release or termination from pA and pB
			rx[10] = k[10]*(apol_paused==1);
			rx[11] = k[11]*(bpol_paused==1);
			// basal promoter transitions (busting)
			if ((a_prom_as+a_prom_b+a_prom_mit)==0) {// transition to basal OFF state 
				rx[12] = 1/k[16];
			} else {
				rx[12] = 0;
			}
			if ((b_prom_as+b_prom_b+b_prom_mit)==0) {// transition to basal OFF state 
				rx[13] = 1/k[17];
			} else {
				rx[13] = 0;
			}
			if (a_prom_b>0 & a_prom_mit==0) {// transition to from basal OFF to ON state 
				rx[14] = 1/k[14];
			} else {
				rx[14] = 0;
			}
			if (b_prom_b>0 & b_prom_mit==0) {// transition to from basal OFF to ON state 
				rx[15] = 1/k[15];
			} else {
				rx[15] = 0;
			}
			
			
			
            //To implement elongation as Gillespie reaction: additional reaction with rate 1/v => call function elongation if this reaction occurs, but then also time needs to be differently measured

			sum_rx = 0, sel_rx = 100;
            for (q=0; q<rx.size(); q++) {
                sum_rx = sum_rx + rx[q];
                cum_rx[q] = sum_rx;
                //cout << rx[q] << " "; // << endl;
            }
         	// Calculate timepoint at which next reaction occurs
            if (sum_rx>0) {
                t_next = t_next + (-log(dist_real(gen)))/(sum_rx);
            // if no Gillespie reaction occurs, set t_next to value > v => next reaction will be an elongation
            } else {t_next=1;}
            
            // if a Gillespie reaction occurs, decide which one
            if (t_next<v){
                r = dist_real(gen);
                for (q=0; q<rx.size(); q++) {
                    if (r<(cum_rx[q]/sum_rx)){
                        sel_rx = q;
                        break;
                    }
            	 }
            	//execute reaction
            	switch (sel_rx){
                    case 0:
                        apol_paused = 1; break;
                    case 1:
                        bpol_paused = 1; break;
                    case 2:
                        a_rna--; break;
                    case 3:
                        b_rna--; break;
                    case 4:
                        a_prom_as--; break;
                    case 5:
                        b_prom_as--; break;
                    case 6:
                        a_pr++; break;
                    case 7:
                        b_pr++; break;
                    case 8:
                        a_pr--; break;
                    case 9:
                        b_pr--; break;
                    case 10:
						apol_paused=0;
						// RNAP release or termination
						if (dist_real(gen)>=k[12]) {
						apol_ini=1;
						}
						break;
					case 11:
						bpol_paused=0; 
						// RNAP release or termination
						if (dist_real(gen)>=k[13]) {
						bpol_ini=1; 
						}
						break;
					case 12: a_prom_b = 1; break;
					case 13: b_prom_b = 1; break;
					case 14: a_prom_b = 0; break;
					case 15: b_prom_b = 0; break;
                    default:
                        break;
                }
            }
            }
            //elongation and TI    
            //passes reference variables (no pointers can be passed)   
        	//cout << "a_prom_as before elong = " << a_prom_as << endl;
        	//cout << "temp_sum_pB in main funct before elong= " << temp_sum_pB << endl;
        	elongation (apol_paused, bpol_paused, overlap, dima_prior, dimb_prior, dimb_post, a_gene, b_gene, sum_ap, sum_bp, 
            apol_ini, bpol_ini, a_rna, b_rna, 
            all_a, temp_sum_ar, temp_sum_br, temp_sum_ap, temp_sum_bp, temp_sum_pA, temp_sum_pB, a_prom_as, b_prom_as, a_prom_b, b_prom_b, dist_real, step, p);
            // Update Protein count after each elongation step for smoothing
            temp_sum_apr = temp_sum_apr + a_pr;
			temp_sum_bpr = temp_sum_bpr + b_pr;
            
            
            //cout << "a_prom_as after elong= " << a_prom_as << endl;
            //cout << "temp_sum_pB in main funct after elong= " << temp_sum_pB << endl; 
            
            //write out smoothed variable
        m++; mit_count++;          
        if (m==n_small_steps) {
            l++;
            aro[l] = temp_sum_ar/n_small_steps;
            bro[l] = temp_sum_br/n_small_steps;
            apro[l] = temp_sum_apr/n_small_steps;
            bpro[l] = temp_sum_bpr/n_small_steps;
            apo[l] = temp_sum_ap/n_small_steps;
            bpo[l] = temp_sum_bp/n_small_steps;
            pAo[l] = temp_sum_pA/n_small_steps;
            pBo[l] = temp_sum_pB/n_small_steps;

            time[l] = step*v;
            temp_sum_ar = 0;
            temp_sum_br = 0;
            temp_sum_apr = 0;
            temp_sum_bpr = 0;
            temp_sum_ap = 0;
            temp_sum_bp = 0;
            temp_sum_pA = 0;
            temp_sum_pB = 0;
            a_prom_mit = 0;
			b_prom_mit = 0;
            m = 0;
        }
        if (mit_count==16*n_small_steps) { // set promoter OFF for 1h every 16h
			a_prom_mit = 1;
			b_prom_mit = 1;
			mit_count = 0;	
		}
	}
    if (m>0) {
        l++;
        aro[l] = temp_sum_ar/m;
        bro[l] = temp_sum_br/m;
        apro[l] = temp_sum_apr/m;
        bpro[l] = temp_sum_bpr/m;
        apo[l] = temp_sum_ap/m;
        bpo[l] = temp_sum_bp/m;
        pAo[l] = temp_sum_pA/n_small_steps;
        pBo[l] = temp_sum_pB/n_small_steps;

        time[l] = step*v;
    }
}
        
            
            
        
                
	
