#include <Rcpp.h>
#include "utils.h"

using namespace Rcpp;

//==============================================================================

//-----------------------------------------------------------
// return JSSA probability animal n detected at least once   
// Pledger et al. 2010 Eqn (3) (mixtures outside)            
//-----------------------------------------------------------

// // [[Rcpp::export]]
// NumericVector PCH0cpp (int type, int x, int nc, int jj, 
//                        const IntegerVector& cumss, int nmix,
//                        const NumericVector& openval0, int cc0, 
//                        const IntegerVector& PIA0,  
//                        const IntegerVector& PIAJ,  
//                        const NumericVector& intervals)
// {
//     double *p = NULL;
//     double *phij = NULL;
//     double *g = NULL;
//     double *beta = NULL;
//     NumericVector pdt(nc, 0.0);
//     double pj, ptmp;
//     int j,n,s;
//     int b,d;
//     
//     p = (double *) R_alloc (cumss[jj], sizeof(double));
//     phij = (double *) R_alloc (jj, sizeof(double));
//     g = (double *) R_alloc (jj, sizeof(double));
//     beta = (double *) R_alloc (jj, sizeof(double));
//     
//     for (n = 0; n<nc; n++) {
//         getp (n, x, nc, cumss[jj], openval0, cc0, PIA0, p);
//         getphij (n, x, nc, jj, openval0, cc0, PIAJ, intervals, phij);
//         getg (type, n, x, nc, jj, openval0, cc0, PIAJ, g);
//         getbeta (type, n, x, nc, jj, openval0, cc0, PIAJ, intervals, phij, beta);
//         
//         for (b = 1; b <= jj; b++) {
//             for (d = b; d <= jj; d++) {
//                 pj = beta[b-1];              // entered at b
//                 for (j = b; j < d; j++) {
//                     pj *= phij[j-1];           // survived
//                 }
//                 pj *= 1-phij[d-1];             // departed at d
//                 
//                 for (j = b; j <= d; j++) {
//                     ptmp = 1;
//                     for (s = cumss[j-1]; s < cumss[j]; s++) {
//                         ptmp *= 1 - p[s];       // not detected
//                     }
//                     if (type == 27) {
//                         pj *= g[j] + (1-g[j]) * ptmp;   // absent or present and not detected
//                     }
//                     else {
//                         pj *= ptmp;
//                     }
//                 }
//                 
//                 pdt[n] += pj;
//             }
//         }
//     }
//     return pdt;
// }
// //==============================================================================

// [[Rcpp::export]]
NumericVector PCH1cpp (int type, int x, int nc, int jj, 
                       const IntegerVector& cumss, int nmix,
                       const NumericVector& openval0, int cc0, 
                       const IntegerVector& PIA0,  
                       const IntegerVector& PIAJ,  
                       const NumericVector& intervals)
{
    double *p = NULL;
    double *phij = NULL;
    double *g = NULL;
    double *beta = NULL;
    NumericVector pdt(nc, 0.0);
    double pbd, ptmp;
    int j,n,s;
    int b,d;
    
    p = (double *) R_alloc (cumss[jj], sizeof(double));
    phij = (double *) R_alloc (jj, sizeof(double));
    g = (double *) R_alloc (jj, sizeof(double));
    beta = (double *) R_alloc (jj, sizeof(double));
    
    for (n = 0; n<nc; n++) {
        getp (n, x, nc, cumss[jj], openval0, cc0, PIA0, p);
        getphij (n, x, nc, jj, openval0, cc0, PIAJ, intervals, phij);
        getg (type, n, x, nc, jj, openval0, cc0, PIAJ, g);
        getbeta (type, n, x, nc, jj, openval0, cc0, PIAJ, intervals, phij, beta);
        
        for (b = 1; b <= jj; b++) {
            for (d = b; d <= jj; d++) {
                pbd = beta[b-1];              // entered at b
                for (j = b; j < d; j++) {
                    pbd *= phij[j-1];         // survived
                }
                pbd *= 1-phij[d-1];           // departed at d
                
                ptmp = 1;
                for (j = b; j <= d; j++) {
                    for (s = cumss[j-1]; s < cumss[j]; s++) {
                        ptmp *= 1 - p[s];       // not detected
                    }
                }
                pdt[n] += pbd * (1 - ptmp);
            }
        }
    }
    return pdt;
}
//==============================================================================

void pr0njmx (int n, int x, 
	      const IntegerVector& cumss, 
	      int nc,  int jj, int kk, int mm, int cc0, int binomN,  
              const IntegerVector& PIA0, 
              const NumericVector& gk0, 
              const NumericVector& Tsk, 
              double pjm[]) {

    int c, i, j, k, m, s, ci, gi, pi, ss;
    double Tski;
    double size;
    
    for (i = 0; i < jj*mm; i++) pjm[i] = 1.0;
    ss = cumss[jj];
    for (j = 0; j < jj; j++) {
	for (k = 0; k < kk; k++) {
            // consider occasions from focal session (j)
	    for (s = cumss[j]; s < cumss[j+1]; s++) {
		ci =  i4(n, s, k, x, nc, ss, kk);
		c = PIA0[ci] - 1;
		if (c >= 0) {    
		    Tski = Tsk[s * kk + k];                      // k x s matrix  
		    if (abs(Tski-1) > 1e-10) {                   // effort <> 1.0 
			if (binomN==1) 
			    size = (int) Tski;   
			else 
			    size = binomN;     
		    }
		    else size = 1;

		    for (m = 0; m < mm; m++) {
			gi = i3(c, k, m, cc0, kk);
			pi = m * jj + j;
			if (binomN == 0)
			    pjm[pi] *= exp(-Tski *  -log(1-gk0[gi]));  // Poisson 
			else {
			    if (size == 1)
				pjm[pi] *= 1-gk0[gi];   
			    else
				pjm[pi] *= pow(1-gk0[gi], size);  // Binomial or Bernoulli
			}
		    }
		}
	    }
	}
    }
}
//==============================================================================

// [[Rcpp::export]]
NumericVector PCH1secrcpp (int type, int x, int nc, int jj, 
                           const IntegerVector& cumss, 
                           int kk, int mm, 
                           const NumericVector& openval0, 
                           int cc0, 
                           const IntegerVector& PIA0, 
                           const IntegerVector& PIAJ, 
                           const NumericVector& gk0, 
                           int binomN, 
                           const NumericVector& Tsk, 
                           const NumericVector& intervals,
                           const IntegerVector& moveargsi,
                           int movemodel,
                           const CharacterVector& usermodel,
                           int kn,
                           const IntegerVector& kernel,
                           const IntegerVector& mqarray,
                           double cellsize) {
    
    double *phij = NULL;
    double *beta = NULL;
    double *moveargs = NULL;
    double *kernelp = NULL;
    NumericVector pdt(nc, 0.0);
    double pbd;
    int j, m, n;
    int b,d;
    double *pjm;
    double *pm;
    double prw0, sumpj;  //, prodpj;
    pm = (double *) R_alloc (mm, sizeof(double));
    pjm = (double *) R_alloc (jj * mm, sizeof(double));
    phij = (double *) R_alloc (jj, sizeof(double));
    beta = (double *) R_alloc (jj, sizeof(double));
    if (movemodel > 1) {
        moveargs = (double *) R_alloc (jj*2, sizeof(double));
        kernelp = (double *) R_alloc (kn*(jj-1), sizeof(double));
    }
    
    for (n = 0; n<nc; n++) {
        getphij (n, x, nc, jj, openval0, cc0, PIAJ,  intervals, phij);
        getbeta (type, n, x, nc, jj, openval0, cc0, PIAJ, intervals, phij, beta);
        // primary-session-specific Pr for this animal
        pr0njmx(n, x, cumss, nc, jj, kk, mm, cc0, binomN, PIA0, gk0, Tsk, pjm);
        if (movemodel>1) {
            getmoveargs (n, x, nc, jj, openval0, cc0, PIAJ, moveargsi, moveargs);
            fillkernelp (kn, jj, movemodel-2, kernel, cellsize, 
                         moveargsi, moveargs, usermodel, kernelp);
        }
        for (b = 1; b <= jj; b++) {
            for (d = b; d <= jj; d++) {
                pbd = beta[b-1];
                for (j = b; j < d; j++) {
                    pbd *= phij[j-1];
                }
                pbd *= 1-phij[d-1];
                // static home ranges: take sum over M of product over J
                if ((movemodel==0) || (movemodel>1)) {
                    for (m=0; m<mm; m++) pm[m] = 1;
                    for (j = b; j <= d; j++) {
                        for (m=0; m<mm; m++) pm[m] = pm[m] * pjm[m*jj + j - 1];
                        if ((j<d) && (movemodel>1)) 
                            convolvemq(mm, kn, j, kernelp, mqarray, pm);
                    }
                    prw0 = 0;
                    for (m=0; m < mm; m++) prw0 += pm[m] / mm;
                    
                }
                else {   // movemodel == 1
                    // over primary sessions in which may have been alive 
                    // centers allowed to differ between primary sessions 
                    // uncorrelated home ranges: take product over J of sum over M
                    prw0 = 1.0;
                    for (j = b; j <= d; j++) {
                        sumpj = 0;
                        for (m=0; m<mm; m++) {
                            sumpj += pjm[m * jj + j -1];
                        }
                        prw0 *= sumpj / mm;   
                    }
                }
                pdt[n] += pbd * (1 - prw0);
            }
        }
    }
    return pdt;
}
//==============================================================================

// [[Rcpp::export]]
NumericVector PCH0secrjcpp (int type, int x, int nc, int jj,
			   const IntegerVector& cumss, 
			   int kk, int mm, 
			   const NumericVector& openval0, 
			   int cc0, 
			   const IntegerVector& PIA0, 
			   const NumericVector& gk0, 
			   int binomN, 
			   const NumericVector& Tsk) {
    // by primary session
    NumericVector pdt(nc * jj, 0.0);
    int j, m, n;
    double *pjm;
    pjm = (double *) R_alloc (jj * mm, sizeof(double));
    for (n = 0; n<nc; n++) {
	// primary-session-specific Pr for this animal
	pr0njmx(n, x, cumss, nc, jj, kk, mm, cc0, binomN, PIA0, gk0, Tsk, pjm);
	for (j=0; j<jj; j++) {
	    for (m=0; m<mm; m++) {
		pdt[j * nc + n] += pjm[m * jj + j] / mm;
            }
        }
    }
    return pdt;
}
//==============================================================================

