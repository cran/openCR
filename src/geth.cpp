#include <Rcpp.h>
#include <RcppParallel.h>
#include "utils.h"

using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

//==============================================================================
// multicatch only

// [[Rcpp::export]]
List gethcpp (int nc1, int cc, int nmix, int nk, int jj, int mm, 
              const IntegerVector PIA, 
              const IntegerVector cumss, 
              const NumericVector Tsk, 
              const NumericVector hk) {
    
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
//==============================================================================
