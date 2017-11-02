/*************************************************************************
FindCP.c  -  Change Point Detection Method
-------------------
last modified   : Mon Dec 19 2016 (NS)
email           : hawyang@princeton.edu, nsong@princeton.edu
*************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_erf.h>
#include "changepoint.h"

#define MAX_X 10.0
#define DX    0.01
#define TOL   1.0e-8

// Functions for calculation of critical values
double C_ac();
double Vostrikova();
struct changepoint* AddCPNode();

// Find a change point in traj of timestep delta_t with bounds cpl (left) and cpr (right), 
// exclusive. Type-I error of alpha. Found change point added to tree pointed to by cp. Na is length
// of traj, rc indicates recursive tracing
int FindCP(struct changepoint **cp, double *traj, double delta_t, int cpl, int cpr, double alpha, 
	int * Ncp, int Na, int rc ) {
	int n,i,k,k_max;				//indices to iterate through traj	
	int LB, RB;                 	//left and right bounds of confidence region
	int cp2=0, cp1=0, cp_max;   	//largest change point in absolute index

	//Relevant variables for gaussian mean change point detection
	double kL, nL;												// Double analogs in integers
	double lmean, lvar, rmean, rvar, wvar, wmean; 				// Calculation of log-likelihood
	double critical_region;			// Critical value for specified alpha and N
	double *llrt, llrt_max=0.0;		// Maximum log-likelihood
	double *cumsum, *cumsumSq;		// Arrays to store values
	enum bool dummy_bool = false;	// If applicable, denotes if found CP has been inserted
	double dim = 2;					// Two dimensions for log-likelihood calculation

	// Used in calculation of type of change point
	//double *vlvar, vrvar;
	//double lvar_max, rvar_max, lmean_max, rmean_max, wvar_max; 
	//double mBIC, vBIC, mvBIC;
	//int type = 0;
	
	n = cpr - cpl;
	LB = cpl;
	RB = cpr;
	cp_max = 0;
	if ( n > 1 ) {
		nL = (double) n;

		// Create workspace
		cumsum = (double *) malloc( (n+1) * sizeof(double));
		cumsumSq = (double *) malloc( (n+1) * sizeof(double));
		llrt = (double *) calloc((n+1),sizeof(double));
		vlvar = (double *) malloc( (n+1) * sizeof(double));

		// Calculate the critical region using Horvath's approximation
		// Note that the dimensions are two since two parameters are considered
		critical_region = C_ac(alpha, n, dim); 
			
		cumsum[0] = 0;
		cumsumSq[0] = 0;
		for (i=1; i<n; i++) {
			// Estimate Gaussian means at different k
			cumsum[i] = cumsum[i-1] + traj[cpl+i];
			cumsumSq[i] = cumsumSq[i-1] + (traj[cpl+i]*traj[cpl+i]);
		}

		// Overall mean
		wmean = cumsum[n-1]/(nL-1.0);
		for(i = cpl+1, wvar = 0.0; i < cpr; i++) {
			// Overall variance
			wvar += (traj[i] - wmean) * (traj[i] - wmean)/(nL-2.0)/2.0/delta_t;

			// Variance to be used in variance change point classification
			vlvar[i-cpl] = wvar*(nL-1.0) * 2.0 * delta_t/(i-cpl);
		}

		for (k=2,k_max=0, lvar = 1.0, rvar = 1.0; k<n-2; k++) {
			kL = (double) k;

			lmean = cumsum[k]/k;
			// Estimation of variance
			lvar = (cumsumSq[k]/k - lmean*lmean)/delta_t/2.0;

			// Can also explicitly calculate variance
			// lvar = 0.0;
			// for(i = cpl+1; i <= cpl+k; i++) {
			// 	lvar = lvar + (traj[i] - lmean) * (traj[i] - lmean)/k/delta_t/2.0;
			// }

			rmean = (cumsum[n-1] - cumsum[k])/(nL-kL-1.0);	 
			// Estimation of variance 
			rvar = ((cumsumSq[n-1] - cumsumSq[k])/(nL-kL-1.0) - rmean*rmean)/delta_t/2.0;

			// Can also explicitly calculate variance
			// rvar = 0.0;
			// for(i = cpl+k+1; i < cpr; i++) {
			// 	rvar = rvar + (traj[i] - rmean) * (traj[i] - rmean)/(nL-kL-1.0)/delta_t/2.0;
			// }

			lmean = lmean/delta_t;
			rmean = rmean/delta_t;

			// Traditional Generalized Log-likelihood ratio test
			llrt[k] = 0.5 * ((nL-1.0)*log(wvar) - kL * log(lvar) - (nL-kL-1.0)* log(rvar));

			if ( llrt[k] > llrt_max ) {	
				// Find the maximum of the likelihood functions
				llrt_max = llrt[k];
				k_max = k;
				lvar_max= lvar*delta_t*2.0;
	        	rvar_max = rvar*delta_t*2.0;
	        	lmean_max = lmean;
	        	rmean_max = rmean;
			}
		}

		if ( llrt_max > 0.5*critical_region ) {
			LB = k_max;
			RB = k_max;
			cp_max = cpl + RB;

			// Determine type of change point

			// Calculate right variance with respect to mu_0
			//vrvar = (wvar*delta_t*2.0*(nL-1.0) - vlvar[k_max] * k_max)/(nL-k_max-1.0);
			// Calculate overall variance
	      	//wvar = (lvar_max*k_max + rvar_max*(nL-k_max-1.0))/(nL-1.0);

	      	/* BIC model selection for type of change point. Can find max of -BIC = 2*LL - penalty
			mBIC = (1-nL) * log(2*M_PI*wvart)-(nL-1.0) - 3 * log(nL-1);
			vBIC = -k_max*log(2*M_PI*vlvar[k_max]) -(nL-k_max-1.0)*log(2*M_PI*vrvar) -(nL-1)-3*log(nL-1);
			mvBIC = -k_max*log(2*M_PI*lvar_max) -(nL-k_max-1.0)*log(2*M_PI*rvar_max) -(nL-1)-4*log(nL-1);

			// Change in mean
			if(mBIC > vBIC && mBIC > mvBIC)
				type = 1;
			// Change in variance
			else if(vBIC > mBIC && vBIC > mvBIC)
				type = 2;
			// Change in mean and variance
			else if (mvBIC > mBIC && mvBIC > vBIC)
				type = 3;
			*/

			// Add change point to current tree
			*cp = AddCPNode(k_max+cpl, LB+cpl, RB+cpl, *cp, Ncp, &dummy_bool);

			// No recursion for checking change points
			if(rc==3){
				free(llrt);
				free(cumsum);
				free(cumsumSq);
				free(vlvar);
				// Return index of change point found
				return k_max+cpl;
			}
			// Recursively find change points
			if((rc==1)||(rc==0)){
				// Go to the left branch
				cp1 = FindCP(cp, traj, delta_t, cpl, LB+cpl, alpha, Ncp, Na, 1);
				// Go to the right branch
				cp1 = FindCP(cp, traj, delta_t, RB+cpl, cpr, alpha, Ncp, Na, 1);
				if (cp1 > cp_max ) 
					cp_max = cp1;
			}
		}
		free(llrt);
		free(vlvar);
		free(cumsum);
		free(cumsumSq);
	}
	// Return maximum index of change point found
	return cp_max;
}
