#include "utils.h"
#include <string>
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

// SECR models without multithreading in RcppParallel 
// - some code specific to multicatch traps (prwimulti, getp0)
// - some code for proximity/count detectors (prw)
// - then code applicable to both (prwisecrcpp)

//==============================================================================
// multicatch only
// Code to update pjm at each mask point for session j when animal n was alive 
// capture sites are in w for each occasion (secondary session) s              
// used by prwisecrmulti (exclusive detectors)   

void prwimulti (int n, int j, int x, int cc, int nc, int jj, int kk, int mm, int nmix, 
                const RVector<int> cumss, 
                const RVector<int> w, 
                const RVector<int> PIA, 
                const RVector<double> hk, 
                const RMatrix<double> Tsk, 
                const RMatrix<double> h, 
                const RMatrix<double> p0,
                const RMatrix<int> hindex, 
                std::vector<double> &pjm) {
    
    int c, gi, k, m, s, wi, wxi;
    double Tski;
    bool dead = false;
    
    // over secondary sessions (occasions) in this primary session 
    for (s = cumss[j-1]; s < cumss[j]; s++) {
        wi = w[s * nc + n];
        if (wi < 0) dead = true;  
        k = abs(wi)-1;         // trap number 0..kk-1; k = -1 if not caught 
        
        // Not captured in any trap on occasion s 
        if (k < 0) {
            for (m=0; m<mm; m++) {
                if (h(m, hindex(n,s)) > fuzz)
                    pjm[m] *= p0(m, hindex(n,s));     
            }
        }
        // Captured in trap k 
        else {
            Tski = Tsk(k,s);
            wxi = i4(n, s, k, x, nc, cumss[jj], kk);   
            c = PIA[wxi] - 1;
            if (c >= 0) {    // drops unset traps 
                for (m=0; m<mm; m++) {
                    gi  = i3(c, k, m, cc, kk);
                    pjm[m] *=  Tski * (1-p0(m, hindex(n,s))) *  hk[gi] / h(m, hindex(n,s));
                }
            }
        }
        if (dead) break;   // out of s loop
    }
}
//==============================================================================

// not multicatch
void prw (int n, int j, int x, int nc, int jj, int kk, int mm, int cc, int binomN, 
           const RVector<int> cumss, 
           const RVector<int> PIA, 
           const RVector<double> gk,
           const RMatrix<double> Tsk, 
           const RVector<int> w,
           std::vector<double> &pjm) {
    int c, k, m, s, wxi, gi, wi, count;
    bool dead = false;
    for (s = cumss[j-1]; s < cumss[j]; s++) {
        for (k=0; k<kk; k++) {
            wxi =  i4(n, s, k, x, nc, cumss[jj], kk);
            c = PIA[wxi] - 1;
            if (c >= 0) {    // drops unset traps 
                wi = i3(n, s, k, nc, cumss[jj]);
                count = w[wi];
                if (count<0) {count = -count; dead = true; }
                for (m=0; m<mm; m++) {
                    gi  = i3(c, k, m, cc, kk);
                    pjm[m] *= pski(binomN, count, Tsk(k,s), gk[gi]);
                }
            }
        }
        if (dead) break;   // after processing all traps on this occasion
    }
 }
//==============================================================================

// This code serves both multi-catch and proximity/count detectors
// Multi-catch data are indicated by hindex[0] >= 0
// Proximity/count data are indicated by hindex[0] < 0
// The argument 'gk' is understood to contain hazard values if hindex[0]>=0

// [[Rcpp::export]]
double prwisecrcpp (int type, int n, int x, int nc, int jj, 
		       int kk, int mm, int nmix, 
		       const IntegerVector cumss, 
		       const IntegerVector w, 
		       const IntegerVector fi, 
		       const IntegerVector li,
		       const NumericVector gk,       // doubles as hk for prwimulti
		       const NumericMatrix openval, 
		       const IntegerVector PIA,
		       const IntegerVector PIAJ,
		       int binomN, 
		       const NumericMatrix Tsk, 
		       const NumericVector intervals,
		       const IntegerVector moveargsi, 
		       const NumericMatrix h,
		       const IntegerMatrix hindex, 
		       int   CJSp1,
		       int movemodel,
		       const CharacterVector usermodel,
		       const IntegerMatrix kernel,
		       const IntegerMatrix mqarray,
		       double cellsize) {

    // Accessors needed by RcppParallel for safe reference to R objects
    const RVector<int> cumssR(cumss); 
    const RVector<int> wR(w); 
    const RVector<int> fiR(fi); 
    const RVector<int> liR(li);
    const RVector<double> gkR(gk);
    const RMatrix<double> openvalR(openval);
    const RVector<int> PIAR(PIA); 
    const RVector<int> PIAJR(PIAJ); 
    const RMatrix<double> TskR(Tsk);
    const RVector<double> intervalsR(intervals);
    const RVector<int> moveargsiR(moveargsi);
    const RMatrix<double> hR(h);
    const RMatrix<int> hindexR(hindex); 
    const RMatrix<int> kernelR(kernel);
    const RMatrix<int> mqarrayR(mqarray);

    NumericMatrix p0(h.nrow(), h.ncol());
    RMatrix<double> p0R(p0);
    
    //------------------------------------------------------------
    // local variables
    double pdt;
    double pbd;
    int j, m;
    int b, minb, maxb;
    int d, mind, maxd;
    int cjs = 0;  // CJS offset
    double sump, prwi;
    int cc = openvalR.nrow();
    int kn = kernelR.nrow();

    //------------------------------------------------------------
    // work vectors for session-specific real parameter values
    std::vector<double> phij(jj);      // each primary session
    std::vector<double> beta(jj);      // not used if type==1
    std::vector<double> pjm(mm);
    std::vector<double> moveargs(jj*2);
    std::vector<double> kernelp(kn*(jj-1));
    //------------------------------------------------------------
    // precompute p0 to save time  (multicatch only)
    if (hindexR[0]>=0) {
	int hsize = h.nrow() * h.ncol();
	for (int i=0; i < hsize; i++) {
	    p0R[i] = exp(-h[i]);
	}
    }
    //------------------------------------------------------------
    // survival
    getphij (n, x, nc, jj, openvalR, PIAJR, intervalsR, phij);    
    //------------------------------------------------------------
    // movement kernel
    if (movemodel > 1) {
        getmoveargs (n, x, nc, jj, openvalR, PIAJR, moveargsiR, moveargs);
        fillkernelp (kn, jj, movemodel-2, cellsize, 
                     kernelR, moveargsiR, usermodel, 
                     moveargs, kernelp);
    }
    //------------------------------------------------------------
    if (type==6) {     // CJSsecr
	minb = fiR[n];
	cjs = (1-CJSp1);
    }
    else {             // JSSA
	minb = 1;
        // entry probabilities beta
	getbeta (type, n, x, nc, jj, openvalR, PIAJR, intervalsR, phij, beta);
    }
    maxb = fiR[n];
    mind = abs(liR[n]);
    maxd = jj;
    if (liR[n] < 0) maxd = mind; // possible censoring
    //------------------------------------------------------------

    pdt = 0;
    for (b = minb; b <= maxb; b++) {
	for (d = mind; d <= maxd; d++) {
	    if (type == 6)     // CJSsecr
		pbd = 1;
	    else
		pbd = beta[b-1];
	    for (j = b; j < d; j++) {
		pbd *= phij[j-1];
	    }
	    if (liR[n]>0)    // not censored
		pbd *= 1-phij[d-1];

	    if ((movemodel==0) || (movemodel>1)) {
	        for (m=0; m < mm; m++) pjm[m] = 1;   // uniform
	        for (j = b + cjs; j <= d; j++) {
	            if (hindexR[0] >= 0) {
	                prwimulti (n, j, x, cc, nc, jj, kk, mm, nmix, cumssR, wR, 
                            PIAR, gkR, TskR, hR, p0R, hindexR, pjm);
	            }
	            else {
	                prw (n, j, x, nc, jj, kk, mm, cc, binomN, cumssR, PIAR, gkR, 
                      TskR, wR, pjm);
	            }
	            if ((j<d) && (movemodel>1)) 
	                convolvemq(mm, kn, j, mqarrayR, kernelp, pjm);
	        }
	        prwi = 0;
	        for (m=0; m<mm; m++) prwi += pjm[m] / mm;
	    }
	    else {  // movemodel == 1
	        prwi = 1.0;
	        for (j = b + cjs; j <= d; j++) {
	            for (m=0; m<mm; m++) pjm[m] = 1;  // uniform
	            if (hindexR[0] >= 0) {
	                prwimulti (n, j, x, cc, nc, jj, kk, mm, nmix, cumssR, wR, 
                            PIAR, gkR, TskR, hR, p0R, hindexR, pjm);
	            }
	            else {
	                prw(n, j, x, nc, jj, kk, mm, cc, binomN, cumssR, PIAR, gkR, 
                     TskR, wR, pjm);
	            }
	            sump = 0;
	            for (m=0; m < mm; m++) sump += pjm[m]/mm;
	            prwi *= sump;   // product over primary sessions 
	        }
	    }
	    pdt += pbd * prwi;
	}
    }
    return pdt; 
}
//==============================================================================
