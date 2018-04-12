#include <Rcpp.h>
#include "utils.h"

using namespace Rcpp;

//==============================================================================

// [[Rcpp::export]]
List makegkcpp (int cc, int kk, int mm, int detectfn, int sigmai, 
             const NumericVector& openval, 
             const NumericVector& traps,
             const NumericVector& mask
             ) 
{
    int k, m, c, gi;
    NumericVector gk(cc * kk * mm); 
    NumericVector hk(cc * kk * mm);
    for (k=0; k<kk; k++) {
        for (m=0; m<mm; m++) {
            for (c=0; c<cc; c++) {
                gi = i3(c,k,m,cc, kk);
                hk[gi] = hfn(k, m, c, openval, cc,
                             traps, mask, kk, mm, sigmai, detectfn);
		gk[gi] = 1 - exp(-hk[gi]);
            }
        }
    }
    return List::create(gk, hk);
}
//==============================================================================
