#include <Rcpp.h>
#include "utils.h"

using namespace Rcpp;

//==============================================================================

// [[Rcpp::export]]
double prwicpp (int type, int n, int x, int nc, int jj, 
                    const IntegerVector& cumss, int nmix, 
                    const IntegerVector& w, 
                    const IntegerVector& fi, 
                    const IntegerVector& li, 
                    const NumericVector& openval, 
		    int cc, 
                    const IntegerVector& PIA, 
                    const IntegerVector& PIAJ, 
                    const NumericVector& intervals,
		    int CJSp1) {      // CJSp1 not used for JSSA
    
    double *phij = NULL;
    double *p = NULL;
    double *g = NULL;
    double *beta = NULL;
    double pdt;
    double pbd;
    int j,s;
    int count; 
    bool dead;
    int cjs = 0;    // offset for first primary session (1 for CJS)
    int b,minb,maxb;
    int d,mind,maxd;
    
    // get session-specific real parameter values
    p = (double *) R_alloc (cumss[jj], sizeof(double));
    phij = (double *) R_alloc (jj, sizeof(double));
    g = (double *) R_alloc (jj, sizeof(double));
    getp (n, x, nc, cumss[jj], openval, cc, PIA, p);
    getphij (n, x, nc, jj, openval, cc, PIAJ, intervals, phij);
    getg (type, n, x, nc, jj, openval, cc, PIAJ, g);

    if (type == 1) {
	cjs = 1;
    }
    else {
	beta = (double *) R_alloc (jj, sizeof(double));
	getbeta (type, n, x, nc, jj, openval, cc, PIAJ, intervals, phij, beta);
    }
    // for (int i=0; i<cumss[jj]; i++) Rprintf("%8.5f ", p[i]); Rprintf("\n");
    
    pdt = 0;
    if (type == 1)
	minb = fi[n];
    else 
	minb = 1;
    maxb = fi[n];
    mind = abs(li[n]);
    maxd = jj;
    if (li[n] < 0) maxd = mind; // possible censoring

    // loop over possible birth and death times
    for (b = minb; b <= maxb; b++) {
        for (d = mind; d <= maxd; d++) {
            dead = false;
            if (type == 1) 
		pbd = 1.0;
	    else
		pbd = beta[b-1];
            for (j = b; j < d; j++) {
                pbd *= phij[j-1];
            }
            if (li[n]>0)  // not censored
                pbd *= 1-phij[d-1];
            // pbd now accounts for birth at b (JSSA) and survival to d
            // next multiply by conditional probability of observed CH 
            for (j = b + cjs; j <= d; j++) {
                // detection probability for each secondary session
                // in primary session j
                //notseen = 1;
                for (s = cumss[j-1]; s < cumss[j]; s++) {   
                    count = w[nc * s + n];
                    if (count<0) {count = -count; dead = true; }
                    
                    if (count>0) {
                        // notseen = 0;
                        pbd *= p[s];
                    }
                    else
                        pbd *= 1 - p[s];
                    if (dead) break;
                }
            }   
            pdt += pbd;
        }
    }
    return pdt;
}

//==============================================================================
