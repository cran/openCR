#include <Rcpp.h>
#include "utils.h"
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

//==============================================================================

// [[Rcpp::export]]
List makegkcpp (
        int cc, 
        int kk, 
        int mm, 
        int detectfn, 
        int sigmai, 
        const NumericMatrix openval, 
        const NumericMatrix traps,
        const NumericMatrix mask
) 
{
    const RMatrix<double> openvalR(openval); 
    const RMatrix<double> trapsR(traps);
    const RMatrix<double> maskR(mask);
    int k, m, c, gi;
    NumericVector gk(cc * kk * mm); 
    NumericVector hk(cc * kk * mm);
    for (k=0; k<kk; k++) {
        for (m=0; m<mm; m++) {
            for (c=0; c<cc; c++) {
                gi = i3(c,k,m,cc, kk);
                hk[gi] = hfn(k, m, c, openvalR, trapsR, maskR, sigmai, detectfn);
                gk[gi] = 1 - exp(-hk[gi]);
            }
        }
    }
    return List::create(gk, hk);
}
//==============================================================================

struct Hckm : public Worker {
    
    // input data
    int sigmai;
    int detectfn;
    const RMatrix<double> openval;
    const RMatrix<double> traps;
    const RMatrix<double> mask;
    
    // output vector to write to
    RVector<double> hk;
    RVector<double> gk;
    
    // initialize from Rcpp input and output matrixes (the RMatrix class
    // can be automatically converted to from the Rcpp matrix type)
    Hckm(int sigmai, int detectfn,
         const NumericMatrix openval, 
         const NumericMatrix traps, 
         const NumericMatrix mask, 
         NumericVector hk,
         NumericVector gk)
        : sigmai(sigmai), detectfn(detectfn), openval(openval), 
          traps(traps), mask(mask), hk(hk), gk(gk) {}

    // function call operator that work for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {
        int cc = openval.nrow();
        int kk = traps.nrow();
        for (std::size_t m = begin; m < end; m++) {
            for (int k=0; k < kk; k++) {
                for (int c=0; c < cc; c++) {
                    int gi = i3(c,k,m, cc, kk);
                    hk[gi] = hfn(k, m, c, openval, traps, mask, sigmai, detectfn);
                    gk[gi] = 1 - exp(-hk[gi]);
                }
            }
            
        }
    }
};

// [[Rcpp::export]]
List makegkParallelcpp (int detectfn, int sigmai, int grain, int ncores,
                        const NumericMatrix& openval, 
                        const NumericMatrix& traps,
                        const NumericMatrix& mask
) 
{
    NumericVector hk(openval.nrow() * traps.nrow() * mask.nrow()); 
    NumericVector gk(openval.nrow() * traps.nrow() * mask.nrow()); 
    
    Hckm hckm (sigmai, detectfn, openval, traps, mask, hk, gk);
    
    if (ncores>1) {
        parallelFor(0, mask.nrow(), hckm, grain, ncores);
    }
    else {
        hckm.operator()(0,mask.nrow());    // for debugging avoid multithreading to allow R calls
    }
    return List::create(gk, hk);
}
//==============================================================================

struct Hckmd : public Worker {
    
    // input data
    int sigmai;
    int detectfn;
    const RMatrix<double> openval;
    const RMatrix<double> distmat;

    // output vector to write to
    RVector<double> hk;
    RVector<double> gk;
    
    // initialize from Rcpp input and output matrixes (the RMatrix class
    // can be automatically converted to from the Rcpp matrix type)
    Hckmd(
        int sigmai, 
        int detectfn,
        const NumericMatrix openval, 
        const NumericMatrix distmat, 
        NumericVector hk,
        NumericVector gk)
        : sigmai(sigmai), detectfn(detectfn), openval(openval), 
            distmat(distmat), hk(hk), gk(gk) {}
    
    // function call operator that work for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {
        int cc = openval.nrow();
        int kk = distmat.nrow();
        for (std::size_t m = begin; m < end; m++) {
            for (int k=0; k < kk; k++) {
                for (int c=0; c < cc; c++) {
                    int gi = i3(c,k,m, cc, kk);
                    hk[gi] = hfnd(k, m, c, openval, distmat, sigmai, detectfn);
                    gk[gi] = 1 - exp(-hk[gi]);
                }
            }
            
        }
    }
};

// [[Rcpp::export]]
List makegkParalleldcpp (
        int detectfn, 
        int sigmai, 
        int grain, 
        int ncores,
        const NumericMatrix& openval, 
        const NumericMatrix& distmat
) 
{
    NumericVector hk(openval.nrow() * distmat.nrow() * distmat.ncol()); 
    NumericVector gk(openval.nrow() * distmat.nrow() * distmat.ncol()); 
    
    Hckmd hckmd (sigmai, detectfn, openval, distmat, hk, gk);
    
    if (ncores>1) {
        parallelFor(0, distmat.ncol(), hckmd, grain, ncores);
    }
    else {
        hckmd.operator()(0,distmat.ncol());    // for debugging avoid multithreading to allow R calls
    }
    return List::create(gk, hk);
}
//==============================================================================
