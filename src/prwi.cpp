#include "utils.h"
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;
//==============================================================================

// [[Rcpp::export]]
double prwicpp (int type, int n, int x, int nc, int jj, 
                    const IntegerVector cumss, int nmix, 
                    const IntegerVector w, 
                    const IntegerVector fi, 
                    const IntegerVector li, 
                    const NumericMatrix openval, 
		            const IntegerVector PIA, 
                    const IntegerVector PIAJ, 
                    const NumericVector intervals,
		    int CJSp1) {      // CJSp1 not used for JSSA

    const RVector<int> cumssR(cumss); 
    const RVector<int> wR(w); 
    const RVector<int> fiR(fi); 
    const RVector<int> liR(li); 
    const RMatrix<double> openvalR(openval); 
    const RVector<int> PIAR(PIA); 
    const RVector<int> PIAJR(PIAJ); 
    const RVector<double> intervalsR(intervals);
    
    // work vectors for session-specific real parameter values
    std::vector<double> p(cumss[jj]);  // each secondary session
    std::vector<double> phij(jj);      // each primary session
    std::vector<double> beta(jj);      // not used if type==1
    std::vector<double> g(jj);         // each primary session
    double pdt;
    double pbd;
    int j,s;
    int count; 
    bool dead;
    int cjs = 0;    // offset for first primary session (1 for CJS)
    int b,minb,maxb;
    int d,mind,maxd;
    
    getp (n, x, nc, cumssR[jj], openvalR, PIAR, p);
    getphij (n, x, nc, jj, openvalR, PIAJR, intervalsR, phij);
    getg (type, n, x, nc, jj, openvalR, PIAJR, g);

    if (type == 1) {
	cjs = 1;
    }
    else {
	getbeta (type, n, x, nc, jj, openvalR, PIAJR, intervalsR, phij, beta);
    }
    // for (int i=0; i<cumssR[jj]; i++) Rprintf("%8.5f ", p[i]); Rprintf("\n");
    
    pdt = 0;
    if (type == 1)
	minb = fiR[n];
    else 
	minb = 1;
    maxb = fiR[n];
    mind = abs(liR[n]);
    maxd = jj;
    if (liR[n] < 0) maxd = mind; // possible censoring

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
            if (liR[n]>0)  // not censored
                pbd *= 1-phij[d-1];
            // pbd now accounts for birth at b (JSSA) and survival to d
            // next multiply by conditional probability of observed CH 
            for (j = b + cjs; j <= d; j++) {
                // detection probability for each secondary session
                // in primary session j
                //notseen = 1;
                for (s = cumssR[j-1]; s < cumssR[j]; s++) {   
                    count = wR[nc * s + n];
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
