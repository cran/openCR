#include <Rcpp.h>
#include <RcppParallel.h>
#include "utils.h"

using namespace Rcpp;
using namespace RcppParallel;

//==============================================================================
// 2018-11-12 deleted pmix argument
// 2019-05-06 revised movement algorithm
// 2019-05-19 combined prw and prwmulti; use binomN == -2 for multi

struct Somesecrhistories : public Worker {
    
    // input data
    int   x;
    int   type;
    int   mm;
    int   nc;
    int   binomN;
    int   CJSp1; 
    int   grain;
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
    int                   movementcode;
    int                   edgecode;
    const String          usermodel;
    const RVector<int>    moveargsi;
    const RMatrix<int>    kernel;
    const RMatrix<int>    mqarray;
    
    double                cellsize;    
    int  kk, jj, kn, cc;
    bool indiv;
    
    int cjs;
    std::vector<double> pdotbd;
    
    // output likelihoods
    RVector<double> output;
    
    // Constructor to initialize an instance of Somehistories 
    // The RMatrix class can be automatically converted to from the Rcpp matrix type
    Somesecrhistories(
        int x, int type, int mm, int nc, int binomN,  int CJSp1, int grain,                    
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
        int   movementcode,
        int   edgecode,
        String usermodel,
        const IntegerVector moveargsi,
        const IntegerMatrix kernel,
        const IntegerMatrix mqarray,
        double cellsize,
        NumericVector output)    
        : 
            x(x), type(type), mm(mm), nc(nc), binomN(binomN), CJSp1(CJSp1),  
            grain(grain),
            intervals(intervals), cumss(cumss), w(w), fi(fi), li(li), gk(gk), 
            openval(openval), PIA(PIA), PIAJ(PIAJ), Tsk(Tsk), h(h), hindex(hindex), 
            movementcode(movementcode), edgecode(edgecode), usermodel(usermodel),  
            moveargsi(moveargsi), kernel(kernel), mqarray(mqarray), 
            cellsize(cellsize), output(output) {
        
        // now can initialise these derived counts
        kk = Tsk.nrow();             // number of detectors
        jj = intervals.size() + 1;   // number of primary sessions
        kn = kernel.nrow();          // number of cells in kernel
        cc = openval.nrow();         // number of parameter combinations
        
        pdotbd.resize(jj*jj);
        std::fill(pdotbd.begin(), pdotbd.end(), 1.0);
        
        if (type == 6) {         // CJSsecr
            cjs = 1 - CJSp1;   // offset 1 session iff true CJS (CJSp1=0)
        }
        else {
            cjs = 0;
            // default initialisation - works when all n individuals the same
        }
        
        indiv = individual();
    }
    //==============================================================================
    
    bool individual() {
        // return true if any variation among animals in parameters affecting turnover
        std::vector<double> phij(jj);      
        std::vector<double> beta(jj);      
        std::vector<double> moveargs(jj*2);
        std::vector<double> phi0(jj);      
        std::vector<double> beta0(jj);      
        std::vector<double> moveargs0(jj*2);
        bool result = false;
        int n,j;
        
        getmoveargs (0, x, nc, jj, openval, PIAJ, moveargsi, moveargs0);
        getphij (0, x, nc, jj, openval, PIAJ, intervals, phi0);
        if (type != 6) {
            getbeta (type, 0, x, nc, jj, openval, PIAJ, intervals, phi0, beta0);
        }
        
        for (n=1; n<nc; n++) {
            getmoveargs (n, x, nc, jj, openval, PIAJ, moveargsi, moveargs);
            getphij (n, x, nc, jj, openval, PIAJ, intervals, phij);
            if (type != 6) {
                getbeta (type, n, x, nc, jj, openval, PIAJ, intervals, phij, beta);   
            }
            for (j = 0; j<jj; j++) {
                if (moveargs[j] != moveargs0[j]) {result = true; break;}
                if (moveargs[j+jj] != moveargs0[j+jj]) {result = true; break;}
                if (phij[j] != phi0[j]) {result = true; break;}
                if (type != 6) {
                    if (beta[j] != beta0[j]) {result = true; break;}
                }
            }
            if (result) break;
        }
        return(result);
    }
    //==============================================================================
    
    void getpdotbd (int n, std::vector<double> &pdotbd) {
        
        // precompute pdotbd 2019-05-24
        std::vector<double> pjmat(jj * mm);
        std::vector<double> alpha(mm);    
        std::vector<double> moveargs(jj*2);
        std::vector<double> kernelp(kn*(jj-1));
        
        int j,b,d,m;
        double prodp;
        //--------------------------------------------------------------------
        // Fill kernel if movement model required
        if (movementcode > 1) {
            getmoveargs (n, x, nc, jj, openval, PIAJ, moveargsi, moveargs);
            if (grain<1) { 
                // single thread allows usermodel
                fillkernelp (kn, jj, movementcode-2, cellsize, kernel, moveargsi, 
                    usermodel, moveargs, kernelp);
            }
            else {
                // cannot call R function from RcppParallel worker ncores > 1, 
                // so use this function with no usermodel option:
                fillkernelparallel (kn, jj, movementcode-2, cellsize, 
                    kernel, moveargsi, moveargs, kernelp);
            }
        }        
        //--------------------------------------------------------------------
        pr0njmx(n, x, cumss, nc, jj, kk, mm, cc, binomN, PIA, gk, Tsk, pjmat);
        std::fill(pdotbd.begin(), pdotbd.end(), 0.0);
        //--------------------------------------------------------------------
        for (b = 1; b <= jj; b++) {
            if (grain<1) {
                Rprintf("\n");
            }
            for (d = b; d <= jj; d++) {
                
                //------------------------------------------------------------
                //  No movement
                if (movementcode == 0) {
                    if ((b+cjs) <= jj) { 
                        for (m=0; m<mm; m++) alpha[m] = pjmat[m * jj + (b+cjs) - 1] / mm;
                        for (j = b + cjs + 1; j <= d; j++) {
                            for (m=0; m<mm; m++) alpha[m] *= pjmat[m * jj + j - 1];
                        }
                        pdotbd[(d-1)*jj + b - 1] = 1 - std::accumulate(alpha.begin(),
                            alpha.end(), 0.0);
                    }
                }
                else if (movementcode == 1) {
                    //------------------------------------------------------------
                    //  Uncorrelated centres
                    prodp = 1.0;
                    for (j = b + cjs; j <= d; j++) {
                        for (m=0; m<mm; m++) alpha[m] = pjmat[m * jj + j - 1] / mm;
                        // product over primary sessions 
                        prodp *= std::accumulate(alpha.begin(), alpha.end(), 0.0);
                    }
                    pdotbd[(d-1)*jj + b - 1] = 1 - prodp;
                }
                else {
                    //------------------------------------------------------------
                    //  Modelled movement > 1
                    if ((b+cjs) <= jj) { 
                        for (m=0; m<mm; m++) alpha[m] = pjmat[m * jj + (b+cjs) - 1] / mm;
                        for (j = b + cjs + 1; j <= d; j++) {
                            convolvemq(mm, kn, j-1, edgecode, mqarray, kernelp, alpha);     
                            for (m = 0; m < mm; m++) alpha[m] *= pjmat[m * jj + j - 1];
                        }
                        pdotbd[(d-1)*jj + b - 1] = 1 - std::accumulate(alpha.begin(),
                            alpha.end(), 0.0);
                    }
                }
                //------------------------------------------------------------
                if (grain<1) {
                    Rprintf("n %d b %d d %d pdotbd %10.8E", n, b, d, pdotbd[(d-1)*jj + b - 1]); 
                }
            }
        }
    }
    //==============================================================================
    
    void prw (int j, int n, std::vector<double> &pjm) {
        int c, gi, k, m, s, wi, wxi, count;
        double Tski, H;
        bool dead = false;
        // multi-catch traps
        if (binomN == -2) {
            // over secondary sessions (occasions) in this primary session 
            for (s = cumss[j-1]; s < cumss[j]; s++) {
                wi = w[s * nc + n]; 
                if (wi < 0) dead = true;  
                k = abs(wi)-1;         // trap number 0..kk-1; k = -1 if not caught 
                
                // Not captured in any trap on occasion s 
                if (k < 0) {
                    for (m=0; m<mm; m++) {
                        H = h(m, hindex(n,s));
                        if (H > fuzz)
                            pjm[m] *= exp(-H);
                    }
                }
                // Captured in trap k on occasion s
                else {
                    wxi = i4(n, s, k, x, nc, cumss[jj], kk);   
                    c = PIA[wxi] - 1;
                    if (c >= 0) {    // ignore unset traps 
                        Tski = Tsk(k,s);
                        for (m=0; m<mm; m++) {
                            H = h(m, hindex(n,s));
                            gi  = i3(c, k, m, cc, kk);
                            // in this context gk is understood to be hazard hk
                            pjm[m] *=  Tski * (1-exp(-H)) *  gk[gi] / H;
                        }
                    }
                }
                if (dead) break;   // out of s loop
            }
            
            // Rprintf("n %4d j %3d pjmsum %10.8E \n", n,j, 
            // std::accumulate(pjm.begin(), pjm.end(), 0.0));   
            
        }
        // all other detectors
        else {
            for (s = cumss[j-1]; s < cumss[j]; s++) {
                for (k=0; k<kk; k++) {
                    wxi =  i4(n, s, k, x, nc, cumss[jj], kk);
                    c = PIA[wxi] - 1;
                    if (c >= 0) {    // ignore unset traps 
                        Tski = Tsk(k,s);
                        wi = i3(n, s, k, nc, cumss[jj]);
                        count = w[wi];
                        if (count<0) {count = -count; dead = true; }
                        for (m=0; m<mm; m++) {
                            gi  = i3(c, k, m, cc, kk);
                            pjm[m] *= pski(binomN, count, Tski, gk[gi]);
                        }
                    }
                }
                if (dead) break;   // after processing all detectors on this occasion
            }
        }
    }
    //==============================================================================
    
    double oneprwisecrcpp (int n) {
        double pdt;
        double pbd;
        int j, m;
        int b, minb, maxb;
        int d, mind, maxd;
        double prwi = 1.0;
        
        // work vectors for session-specific real parameter values etc.
        std::vector<double> phij(jj);       // each primary session
        std::vector<double> beta(jj);       // not used if type==1
        std::vector<double> alpha(mm);      // trial 2019-05-05
        std::vector<double> moveargs(jj*2);
        std::vector<double> kernelp(kn*(jj-1));
        
        getphij (n, x, nc, jj, openval, PIAJ, intervals, phij);
        
        if (movementcode > 1) {
            getmoveargs (n, x, nc, jj, openval, PIAJ, moveargsi, moveargs);
            if (grain<1) { 
                // single thread allows usermodel
                fillkernelp (kn, jj, movementcode-2, cellsize, kernel, moveargsi, 
                    usermodel, moveargs, kernelp);
            }
            else {
                // cannot call R function from RcppParallel worker ncores > 1, 
                // so use this function with no usermodel option:
                fillkernelparallel (kn, jj, movementcode-2, cellsize, 
                    kernel, moveargsi, moveargs, kernelp);
            }
        }        
        if (type == 6) {     // CJSsecr
            minb = fi[n];
        }
        else {
            minb = 1;
            getbeta (type, n, x, nc, jj, openval, PIAJ, intervals, phij, beta);
            // allow individual variation in phi, beta, moveargs
            // by replacing pdotbd with custom value for this n
        }
        
        maxb = fi[n];
        
        // possible censoring
        mind = abs(li[n]);
        if (li[n] < 0) 
            maxd = mind;
        else 
            maxd = jj;
        
        pdt = 0.0;
        for (b = minb; b <= maxb; b++) {
            for (d = mind; d <= maxd; d++) {
                // pbd = probability available for detection from b to d inclusive 
                if (type == 6)   // CJSsecr
                    pbd = 1.0;
                else
                    pbd = beta[b-1];
                for (j = b; j < d; j++) {
                    pbd *= phij[j-1];
                }
                if (li[n]>0)    // not censored
                    pbd *= 1-phij[d-1];
                
                prwi = 1.0;
                if ((b+cjs)<=d) {
                    // prwi = probability of observed history given available b to d 
                    if (movementcode == 0) {
                        for (m=0; m<mm; m++) alpha[m] = 1.0/mm;
                        // initialise alpha to uniform density
                        // loop over primary sessions in which may have been alive
                        for (j = b + cjs; j <= d; j++) {
                            prw (j, n, alpha);                        
                        }
                        prwi = std::accumulate(alpha.begin(), alpha.end(), 0.0);
                    }
                    else if (movementcode == 1) {
                        prwi = 1.0;
                        for (j = b + cjs; j <= d; j++) {
                            for (m=0; m<mm; m++) alpha[m] = 1.0/mm;
                            prw (j, n, alpha);                        
                            // product over primary sessions 
                            prwi *= std::accumulate(alpha.begin(), alpha.end(), 0.0);
                        }
                    }
                    else {  // movementcode > 1
                        for (m=0; m<mm; m++) alpha[m] = 1.0/mm;
                        prw (b+cjs, n, alpha);                        
                        for (j = b + cjs + 1; j <= d; j++) {
                            convolvemq(mm, kn, j-1, edgecode, mqarray, kernelp, alpha);  
                            prw (j, n, alpha);                        
                        }		    
                        prwi = std::accumulate(alpha.begin(), alpha.end(), 0.0);
                    }
                }
                pdt += pbd * prwi / pdotbd[(d-1)*jj + b - 1];
                // Rprintf("x %4d n %4d b %3d d %3d pbd %10.8E prwi %10.8E pdotbd %10.8E \n",
                //	x,n,b,d,pbd, prwi,  pdotbd[(d-1)*jj + b - 1]);       
            }
        }
        // Rprintf("x %4d n %4d pdt %10.8E \n", x,n,pdt);       
        return pdt;    // may be zero 
    }
    
    // function call operator that work for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {        
        for (std::size_t n = begin; n < end; n++) {
            output[n] = oneprwisecrcpp (n);
        }
    }
};

// [[Rcpp::export]]
NumericVector allhistsecrparallelcpp (int x, int type, int mm, int nc,
    int binomN, int CJSp1, int grain,
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
    int   movementcode,
    int   edgecode,
    const String usermodel,
    const IntegerVector moveargsi, 
    const IntegerMatrix kernel,
    const IntegerMatrix mqarray,
    double cellsize) {
    
    NumericVector output(nc); 
    
    // Construct and initialise
    Somesecrhistories somehist (
            x, type, mm, nc, binomN, CJSp1, 
            grain, intervals, cumss, 
            w, fi, li, gk, openval, PIA, PIAJ, Tsk, 
            h, hindex,
            movementcode, edgecode,
            usermodel,
            moveargsi, kernel, mqarray, 
            cellsize, output);
    
    if (grain>0) {
        // Run operator() on multiple threads
        parallelFor(0, nc, somehist, grain);
    }
    else {
        // for debugging avoid multithreading and allow R calls e.g. Rprintf
        somehist.operator()(0,nc);    
    }
    
    // Return consolidated result
    return output;
}
//==============================================================================
