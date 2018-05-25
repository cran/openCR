#include <Rcpp.h>
#include <RcppParallel.h>
#include "utils.h"

using namespace Rcpp;
using namespace RcppParallel;

//==============================================================================

struct Somesecrhistories : public Worker {
    
    // input data
    int   x;
    int   type;
    int   mm;
    int   nc;
    int   binomN;
    int   CJSp1;                   
    const RMatrix<double> pmix;
    const RVector<double> intervals;
    const RVector<int>    cumss;
    const RVector<int>    w;
    const RVector<int>    fi;
    const RVector<int>    li;
    const RVector<double> gk; 
    const RMatrix<double> openval;
    const RVector<int>    PIA;
    const RVector<int>    PIAJ;
    const RMatrix<double> Tsk;
    const RMatrix<double> h;
    const RMatrix<int>    hindex;
    int                   movemodel;
    const RVector<int>    moveargsi;
    const RMatrix<int>    kernel;
    const RMatrix<int>    mqarray;
    double                cellsize;    
    int kk, jj, kn, cc;

    // output likelihoods
    RVector<double> output;
    
    // Constructor to initialize an instance of Somehistories 
    // The RMatrix class can be automatically converted to from the Rcpp matrix type
    Somesecrhistories(
        int x, int type, int mm, int nc, int binomN,  int CJSp1,                    
        const NumericMatrix pmix,
        const NumericVector intervals,
        const IntegerVector cumss,
        const IntegerVector w,
        const IntegerVector fi, 
        const IntegerVector li,
        const NumericVector gk, 
        const NumericMatrix openval,
        const IntegerVector PIA,
        const IntegerVector PIAJ, 
        const NumericMatrix Tsk,
        const NumericMatrix h,
        const IntegerMatrix hindex, 
        int   movemodel,
        const IntegerVector moveargsi,
        const IntegerMatrix kernel,
        const IntegerMatrix mqarray,
        double cellsize,
        NumericVector output)    
        : 
        x(x), type(type), mm(mm), nc(nc), binomN(binomN), CJSp1(CJSp1), 
        pmix(pmix), intervals(intervals), cumss(cumss), w(w), fi(fi), li(li), gk(gk), 
        openval(openval), PIA(PIA), PIAJ(PIAJ), Tsk(Tsk), h(h), hindex(hindex), 
        movemodel(movemodel), moveargsi(moveargsi), kernel(kernel), mqarray(mqarray), 
        cellsize(cellsize), output(output) {
        // now can initialise these derived counts
        kk = Tsk.nrow();             // number of detectors
        jj = intervals.size() + 1;   // number of primary sessions
        kn = kernel.nrow();          // number of cells in kernel
        cc = openval.nrow();         // number of parameter combinations

    }

    void prwimulti (int j, int n, std::vector<double> &pjm) {
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
                        // pjm[m] *= p0(m, hindex(n,s));     
                        pjm[m] *= exp(-h(m, hindex(n,s))); 
                }
            }
            // Captured in trap k on occasion s
            else {
                Tski = Tsk(k,s);
                wxi = i4(n, s, k, x, nc, cumss[jj], kk);   
                c = PIA[wxi] - 1;
                if (c >= 0) {    // drops unset traps 
                    for (m=0; m<mm; m++) {
                        gi  = i3(c, k, m, cc, kk);
                        // in this context gk is understood to be hazard hk
                        //pjm[m] *=  Tski * (1-p0(m, hindex(n,s))) *  gk[gi] / h(m, hindex(n,s));
                        pjm[m] *=  Tski * (1-exp(-h(m, hindex(n,s)))) *  gk[gi] / h(m, hindex(n,s));
                    }
                }
            }
            if (dead) break;   // out of s loop
        }
    }

    void prw (int j, int n, std::vector<double> &pjm) {
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
    void pjm1 (std::vector<double> &pjm) {
        for (int m=0; m<mm; m++) pjm[m] = 1;
    }
    double sumpjm (std::vector<double> &pjm) {
        double sump = 0.0;
        for (int m=0; m<mm; m++) sump += pjm[m];
	return sump;
    }
    
    double oneprwisecrcpp (int n) {
        double pdt;
        double pbd;
        int j, m;
        int cjs = 0;          // CJS offset
        int b, minb, maxb;
        int d, mind, maxd;
        double prwi, sump;
        
        // work vectors for session-specific real parameter values etc.
        std::vector<double> phij(jj);      // each primary session
        std::vector<double> beta(jj);      // not used if type==1
        std::vector<double> pjm(mm);
        std::vector<double> moveargs(jj*2);
        std::vector<double> kernelp(kn*(jj-1));

        getphij (n, x, nc, jj, openval, PIAJ, intervals, phij);
            
        if (movemodel > 1) {
            getmoveargs (n, x, nc, jj, openval, PIAJ, moveargsi, moveargs);
            fillkernelparallel (kn, jj, movemodel-2, cellsize, 
                         kernel, moveargsi, moveargs, kernelp);
        }        
        if (type == 6) {     // CJSsecr
            minb = fi[n];
            cjs = 1 - CJSp1;
        }
        else {
            minb = 1;
	    cjs = 0;
            getbeta (type, n, x, nc, jj, openval, PIAJ, intervals, phij, beta);
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
                    pjm1(pjm);     // uniform
                    // over primary sessions in which may have been alive, update pjm
                    for (j = b + cjs; j <= d; j++) {
                        if (hindex[0] >= 0) {
                            prwimulti (j, n, pjm);
                        }
                        else {
                            prw (j, n, pjm);                        
                        }
                        if ((j<d) && (movemodel>1)) 
                        convolvemq(mm, kn, j, mqarray, kernelp, pjm);  
                    }
		    //debug
		    //for (int m=0; m<100; m++) 
		    //Rprintf("n %4d b %4d d %4d m %4d pjm[m] %8.6f \n", n,b,d,m,pjm[m]);

		    prwi = sumpjm(pjm) / mm;
                }
                else {  // movemodel == 1
                    prwi = 1.0;
                    for (j = b + cjs; j <= d; j++) {
                        for (m=0; m<mm; m++) pjm[m] = 1;  // uniform
                        if (hindex[0] >= 0) {
                            prwimulti (j, n, pjm);
                        }
                        else {
                            prw (j, n, pjm);                        
                        }
                        sump = 0;
                        for (m=0; m < mm; m++) sump += pjm[m]/mm;
                        prwi *= sump;   // product over primary sessions 
                    }
                }
                // Rprintf("n %4d b %4d d %4d pbd %8.6f prwi %8.6f \n", n,b,d,pbd,prwi);
                pdt += pbd * prwi;
            }
        }
        return pdt; 
    }

    // function call operator that work for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {
        
        for (std::size_t n = begin; n < end; n++) {
            output[n] = pmix(x,n) * oneprwisecrcpp (n);
        }
    }
};

// [[Rcpp::export]]
NumericVector allhistsecrparallelcpp (int x, int type, int mm, int nc,
                                  int binomN, int CJSp1, int grain,
                                  const NumericMatrix pmix,
                                  const NumericVector intervals, 
                                  const IntegerVector cumss, 
                                  const IntegerVector w,
                                  const IntegerVector fi, 
                                  const IntegerVector li,
                                  const NumericVector gk, 
                                  const NumericMatrix openval,
                                  const IntegerVector PIA, 
                                  const IntegerVector PIAJ, 
                                  const NumericMatrix Tsk, 
                                  const NumericMatrix h,
                                  const IntegerMatrix hindex, 
                                  int   movemodel,
                                  const IntegerVector moveargsi, 
                                  const IntegerMatrix kernel,
                                  const IntegerMatrix mqarray,
                                  double cellsize) {
 
    NumericVector output(nc); 
    
    // Construct and initialise
    Somesecrhistories somehist (x, type, mm, nc, binomN, CJSp1, pmix, intervals,
                            cumss, w, fi, li, gk, openval, PIA, PIAJ, Tsk, 
                            h, hindex,
                            movemodel, moveargsi, kernel, mqarray, 
                            cellsize, output);
    
    // Run operator() on multiple threads
    parallelFor(0, nc, somehist, grain);
   
    // somehist.operator()(0,nc);    // for debugging avoid multithreading to allow R calls

    // Return consolidated result
    return output;
}
//==============================================================================
