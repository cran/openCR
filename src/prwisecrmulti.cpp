#include <Rcpp.h>
#include "utils.h"
using namespace Rcpp;

// 2018-03-26 dropped separate CJS, JSSA code
//==============================================================================

// [[Rcpp::export]]
List gethcpp (int nc1, int cc, int nmix, int nk, int jj, int mm, 
               const IntegerVector& PIA, 
               const IntegerVector& cumss, 
               const NumericVector& Tsk, 
               const NumericVector& hk) {
    
    // total hazard for animal n on session j  wrt mask point m 
    // construct index 'hindex' to values in 'h'                
    // 'bk' model results in within-trap variation dependent on 
    // n.j, so h must be recalc each time                       
    // c0 -- index of parameters for trap 0, mixture 0          
    // hc0[c0] -- maps c0 to sequential index 'next'            
    // next -- new index of parameters for each n,j             
    // h -- array of computed hazard for [m,next]               
    
    // mixtures are group-specific for full likelihood, and     
    // individual-specific for conditional likelihood           
    
    int fullns = 0;
    int next = 0;
    int c,c0,i,m,n,k,x,gi,hi,s;
    int PIAval0, PIAvalk; 
    double Tski;          
    
    IntegerVector hc0(cc);
    NumericVector h(nc1 * cumss[jj] * mm * nmix);
    IntegerVector hindex(nc1 * cumss[jj]);
    
    // Recognise when not fully specified by n.s, c  
    // this arises when model = bk, Bk               
    // and when detector covariates vary by time     
    
    for (n=0; n < nc1; n++) {
        for (s=0; s<cumss[jj]; s++) {
            PIAval0 = PIA[i4(n,s,0,0, nc1, cumss[jj], nk)];
            for (k=1; k<nk; k++) {
                PIAvalk = PIA[i4(n,s,k,0, nc1, cumss[jj], nk)];
                if (PIAval0 < 0) {
                    PIAval0 = PIAvalk;
                }
                else if (PIAvalk>0) {
                    if (PIAval0 != PIAvalk) {
                        fullns = 1;
                        break;
                    } 
                }              
            } 
            if (fullns) break;
        }
        if (fullns) break;
    }
    
    for (i=0; i<cc; i++) hc0[i] = -1;
    next = 0;        
    for (n=0; n < nc1; n++) {
        for (s=0; s < cumss[jj]; s++) {
            hi = s*nc1 + n;
            // Case 1. within-trap variation 
            if (fullns) {
                for (k=0; k < nk; k++) {
                    Tski = Tsk[s * nk + k]; 
                    for (x = 0; x < nmix; x++) {
                        c = PIA[i4(n,s,k,x, nc1, cumss[jj], nk)]-1; 
                        if (c >= 0) {
                            for (m = 0; m < mm; m++) { 
                                gi = i3(c,k,m,cc,nk);
                                h[i3(x,m,hi,nmix, mm)] += Tski * hk[gi];
                            }
                        }
                    }
                }
                hindex[hi] = hi;   
            }
            // Case 2. no within-trap variation 
            else {
                c0 = PIA[i4(n,s,0,0, nc1, cumss[jj], nk)] - 1;                    
                if (hc0[c0] < 0) {
                    hc0[c0] = next;
                    next ++;
                    for (k=0; k < nk; k++) {
                        Tski = Tsk[s * nk + k];
                        for (x = 0; x < nmix; x++) {
                            c = PIA[i4(n,s,k,x, nc1, cumss[jj], nk)]-1; 
                            if (c >= 0) {
                                for (m=0; m< mm; m++) { 
                                    gi = i3(c,k,m,cc,nk);
                                    h[i3(x, m, hc0[c0], nmix, mm)] += Tski * hk[gi];
                                }
                            }
                        }
                    }
                }
                hindex[hi] = hc0[c0];
            }
        }
    }
    return List::create(hc0,h,hindex);
}
//=============================================================================

// Code to update pjm at each mask point for session j when animal n was alive 
// capture sites are in w for each occasion (secondary session) s              
// used by prwisecrmulti (exclusive detectors)   

void prwimulti (int n, int j, int x, int cc, int nc, int jj, int kk, int mm, int nmix, 
		const IntegerVector& cumss, 
		const IntegerVector& w, 
		const IntegerVector& PIA, 
		const NumericVector& hk, 
		const NumericVector& Tsk, 
		const NumericVector& h, 
		double p0[], 
		const IntegerVector& hindex, 
		double pjm[]) {

    int c, gi, hi, k, m, s, sn, wi, wxi;
    double Tski;
    bool dead = false;

    // over secondary sessions (occasions) in this primary session 
    for (s = cumss[j-1]; s < cumss[j]; s++) {
	sn = s * nc + n;
	wi = w[sn];
	if (wi < 0) dead = true;  
	k = abs(wi)-1;         // trap number 0..kk-1; k = -1 if not caught 

	// Not captured in any trap on occasion s 
	if (k < 0) {
	    for (m=0; m<mm; m++) {
		hi = i3(x, m, hindex[sn], nmix, mm);
		if (h[hi] > fuzz)
		    pjm[m] *= p0[hi];     
	    }
	}
	// Captured in trap k 
	else {
	    Tski = Tsk[s * kk + k];
	    wxi = i4(n, s, k, x, nc, cumss[jj], kk);   
	    c = PIA[wxi] - 1;
	    if (c >= 0) {    // drops unset traps 
		for (m=0; m<mm; m++) {
		    hi = i3(x,m,hindex[sn],nmix, mm);
		    gi  = i3(c, k, m, cc, kk);
		    pjm[m] *=  Tski * (1-p0[hi]) *  hk[gi] / h[hi];
		}
	    }
	}
	if (dead) break;   // out of s loop
    }
}
//-----------------------------------------------------------------------------

void getp0 (int n, int nc, int nmix, int mm, int ss, 
            const IntegerVector& hindex, 
            const NumericVector& h, 
            double p0[]) {
    int x,m,s,hi;
    for (x=0; x<nmix; x++) {
	for (m=0; m<mm; m++) {
	    for (s=0; s<ss; s++) {
		hi = i3(x, m, hindex[s*nc + n], nmix, mm);
		p0[hi] = exp(-h[hi]);
	    }
	}
    }
}

//-----------------------------------------------------------------------------

// [[Rcpp::export]]
double prwisecrmulticpp (int type, int n, int x, int nc, int jj, 
                             const IntegerVector& cumss, 
                             int kk, int mm, int nmix, 
                             const IntegerVector& w, 
                             const IntegerVector& fi, 
                             const IntegerVector& li, 
                             const NumericVector& hk, 
                             const NumericVector& openval, int cc, 
                             const IntegerVector& PIA, 
                             const IntegerVector& PIAJ, 
                             int binomN, 
                             const NumericVector& Tsk, 
                             const NumericVector& intervals, 
                             const IntegerVector& moveargsi,
                             const NumericVector& h,
                             const IntegerVector& hindex, 
                             int CJSp1,                     // CJSp1 not used for JSSA
			    int movemodel,
			    const CharacterVector& usermodel,
			    int kn,
			    const IntegerVector& kernel,
			    const IntegerVector& mqarray,
			    double cellsize) {

    double *phij = NULL;
    double *moveargs = NULL;
    double *kernelp = NULL;
    double *p0 = NULL;
    double *beta = NULL;
    double pdt;
    double pbd;
    int j, m, s;
    int cjs = 0;
    int b, minb, maxb;
    int d, mind, maxd;
    double *pjm;
    double prwi, sump;

    // precompute p0 to save time 
    int maxi = 0;
    for (s=0;s<cumss[jj];s++)
	if (hindex[s* nc + n] > maxi) maxi = hindex[s* nc + n];
    p0  = (double *) S_alloc (nmix * mm * (maxi+1), sizeof(double));
    getp0 (n, nc, nmix, mm, cumss[jj], hindex, h, p0);
    //for (j=0; j<(nmix * mm * (maxi+1)); j++) p0[j] = exp(-h[j]);
    
    phij = (double *) R_alloc (jj, sizeof(double));
    pjm = (double *) R_alloc (mm, sizeof(double));

    getphij (n, x, nc, jj, openval, cc, PIAJ, intervals, phij);

    if (movemodel>1) {
	moveargs = (double *) R_alloc (jj*2, sizeof(double));
	kernelp = (double *) R_alloc (kn*(jj-1), sizeof(double));
	getmoveargs (n, x, nc, jj, openval, cc, PIAJ, moveargsi, moveargs);
	fillkernelp (kn, jj, movemodel-2, kernel, cellsize, 
		     moveargsi, moveargs, usermodel, kernelp);
    }

    if (type == 6) {     // CJSsecr
	minb = fi[n];
	cjs = 1 - CJSp1;
    }
    else {
	minb = 1;
	cjs = 0;
	beta = (double *) R_alloc (jj, sizeof(double));
	getbeta (type, n, x, nc, jj, openval, cc, PIAJ, intervals, phij, beta);
    }
    maxb = fi[n];

    // possible censoring
    mind = abs(li[n]);
    if (li[n] < 0) 
        maxd = mind;
    else 
        maxd = jj;


    pdt = 0;
    
    for (b = minb; b <= maxb; b++) {
	for (d = mind; d <= maxd; d++) {
            // pbd = probability available for detection from b to d 
	    if (type == 6)   // CJSsecr
		pbd = 1;
	    else
		pbd = beta[b-1];
	    for (j = b; j < d; j++) {
		pbd *= phij[j-1];
	    }
	    if (li[n]>0)    // not censored
	        pbd *= 1-phij[d-1];

            // prwi = probability of observed history given available b to d 
	    if ((movemodel==0) || (movemodel>1)) {
		for (m=0; m<mm; m++) pjm[m] = 1;    // uniform
		// over primary sessions in which may have been alive 
		for (j = b + cjs; j <= d; j++) {
		    prwimulti (n, j, x, cc, nc, jj, kk, mm, nmix, cumss, w, 
			       PIA, hk, Tsk, h, p0, hindex, pjm);
		    if ((j<d) && (movemodel>1)) 
			convolvemq(mm, kn, j, kernelp, mqarray, pjm);
		}
		prwi = 0;
		for (m=0; m<mm;m++) prwi += pjm[m]/mm;
	    }
	    else {   // movemodel == 1
		// over primary sessions in which may have been alive 
		// centers allowed to differ between primary sessions 
		prwi = 1.0;
		for (j = b+cjs; j <= d; j++) {
		    for (m=0; m<mm; m++) pjm[m] = 1;   // uniform
		    prwimulti (n, j, x, cc, nc, jj, kk, mm, nmix, cumss, w, 
			       PIA, hk, Tsk, h, p0, hindex, pjm);
		    sump = 0;
		    for (m=0; m<mm; m++) sump += pjm[m]/mm;
		    prwi *= sump;   // product over primary sessions 
		}
	    }

	    pdt += pbd * prwi; 
	}
    }
    return pdt;
}
//==============================================================================

// [[Rcpp::export]]
double allhistmultcpp (int type, int nc, int jj, 
		       int kk, int mm, int nmix,
		       const IntegerVector& cumss, 
		       const IntegerVector& w,
		       const IntegerVector& fi, 
		       const IntegerVector& li,
		       const NumericVector& hk, 
		       const NumericVector& openval, int cc, 
		       const IntegerVector& PIA, 
		       const IntegerVector& PIAJ, 
		       int binomN, 
		       const NumericVector& Tsk, 
		       const NumericVector& intervals, 
		       const NumericVector& h, 
		       const IntegerVector& hindex, 
		       int   CJSp1,                    // CJSp1 not used for JSSA
		       const IntegerVector& moveargsi, 
		       int movemodel,
		       const CharacterVector& usermodel,
		       int kn,
		       const IntegerVector& kernel,
		       const IntegerVector& mqarray,
		       double cellsize,
		       const NumericVector& pmix) {
    
    int n,x;
    double onep;
    double sump = 0.0;
    for (x=0; x < nmix; x++) {
        for (n=0; n < nc; n++) {
	    onep = prwisecrmulticpp (type, n, x, nc, jj, cumss, kk, mm, nmix, w, fi, li,
				     hk, openval, cc, PIA, PIAJ, binomN, Tsk, intervals,
				     moveargsi, h, hindex, CJSp1, movemodel, usermodel,
				     kn, kernel, mqarray, cellsize); 
            sump += log(onep * pmix[n * nmix + x]);
        }
    }
    return sump;
}
//==============================================================================
