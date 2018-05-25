#include "utils.h"
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

//==============================================================================

//-----------------------------------------------------------
// return JSSA probability animal n detected at least once   
// Pledger et al. 2010 Eqn (3) (mixtures outside)            
//-----------------------------------------------------------

// [[Rcpp::export]]
NumericVector PCH1cpp (int type, int x, int nc, int jj, 
                       const IntegerVector cumss, int nmix,
                       const NumericMatrix openval0, 
                       const IntegerVector PIA0,  
                       const IntegerVector PIAJ,  
                       const NumericVector intervals)
{
    RVector<int> cumssR(cumss);
    RMatrix<double> openval0R(openval0);
    RVector<int> PIA0R(PIA0);
    RVector<int> PIAJR(PIAJ);
    RVector<double> intervalsR(intervals);
    
    NumericVector pdt(nc, 0.0);
    double pbd, ptmp;
    int j,n,s;
    int b,d;
    std::vector<double> p(cumssR[jj]);
    std::vector<double> phij(jj);
    std::vector<double> g(jj);
    std::vector<double> beta(jj);
    
    for (n = 0; n<nc; n++) {
        getp (n, x, nc, cumssR[jj], openval0R, PIA0R, p);
        getphij (n, x, nc, jj, openval0R, PIAJR, intervalsR, phij);
        getg (type, n, x, nc, jj, openval0R, PIAJR, g);
        getbeta (type, n, x, nc, jj, openval0R, PIAJR, intervalsR, phij, beta);
        
        for (b = 1; b <= jj; b++) {
            for (d = b; d <= jj; d++) {
                pbd = beta[b-1];              // entered at b
                for (j = b; j < d; j++) {
                    pbd *= phij[j-1];         // survived
                }
                pbd *= 1-phij[d-1];           // departed at d
                
                ptmp = 1;
                for (j = b; j <= d; j++) {
                    for (s = cumssR[j-1]; s < cumssR[j]; s++) {
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
              const RVector<int> cumss, 
              int nc,  int jj, int kk, int mm, int cc0, int binomN,  
              const RVector<int> PIA0, 
              const RVector<double> gk0, 
              const RMatrix<double> Tsk, 
              std::vector<double> &pjm) {
    
    int c, i, j, k, m, s, ci, gi, pi, ss;

    for (i = 0; i < jj*mm; i++) pjm[i] = 1.0;
    ss = cumss[jj];
    for (j = 0; j < jj; j++) {
	for (k = 0; k < kk; k++) {
            // consider occasions from focal session (j)
	    for (s = cumss[j]; s < cumss[j+1]; s++) {
		ci =  i4(n, s, k, x, nc, ss, kk);
		c = PIA0[ci] - 1;
		if (c >= 0) {    
		    for (m = 0; m < mm; m++) {
			gi = i3(c, k, m, cc0, kk);
			pi = m * jj + j;
			if (binomN == 0)
			    pjm[pi] *= exp(-Tsk(k,s) *  -log(1-gk0[gi]));  // Poisson 
			else {
			    pjm[pi] *=  pski(binomN, 0, Tsk(k,s), gk0[gi]);
			}
		    }
		}
	    }
	}
    }
}
//==============================================================================

// [[Rcpp::export]]
NumericVector PCH1secrcpp (int type, bool individual, int x, int nc, int jj, 
                           const IntegerVector cumss, 
                           int kk, int mm, 
                           const NumericMatrix openval0, 
                           const IntegerVector PIA0, 
                           const IntegerVector PIAJ, 
                           const NumericVector gk0, 
                           int binomN, 
                           const NumericMatrix Tsk, 
                           const NumericVector intervals,
                           const IntegerVector moveargsi,
                           int movemodel,
                           const CharacterVector usermodel,
                           const IntegerMatrix kernel,
                           const IntegerMatrix mqarray,
                           double cellsize) {
    
    const RVector<int> cumssR(cumss); 
    const RMatrix<double> openval0R(openval0);
    const RVector<int> PIA0R(PIA0); 
    const RVector<int> PIAJR(PIAJ); 
    const RVector<double> gk0R(gk0);
    const RMatrix<double> TskR(Tsk);
    const RVector<double> intervalsR(intervals);
    const RVector<int> moveargsiR(moveargsi);
    const RMatrix<int> kernelR(kernel);
    const RMatrix<int> mqarrayR(mqarray);
    
    NumericVector pdt(nc, 0.0);
    double pbd;
    int j, m, n, nc1;
    int b,d;
    double prw0, sumpj;
    int cc0 = openval0R.nrow();
    int kn = kernelR.nrow();
    
    std::vector<double> p(cumssR[jj]);
    std::vector<double> phij(jj);
    std::vector<double> beta(jj);
    std::vector<double> pm(mm);
    std::vector<double> pjm(jj * mm);
    std::vector<double> moveargs(jj*2);
    std::vector<double> kernelp(kn * (jj-1));

    // allow shortcut when all animals the same
    if (individual) nc1 = nc; else nc1 = 1;
    
    for (n = 0; n<nc1; n++) {
        getphij (n, x, nc, jj, openval0R, PIAJR,  intervalsR, phij);
        getbeta (type, n, x, nc, jj, openval0R, PIAJR, intervalsR, phij, beta);
        // primary-session-specific Pr for this animal
        pr0njmx(n, x, cumssR, nc, jj, kk, mm, cc0, binomN, PIA0R, gk0R, TskR, pjm);
        if (movemodel>1) {
            getmoveargs (n, x, nc, jj, openval0R, PIAJR, moveargsiR, moveargs);
            fillkernelp (kn, jj, movemodel-2, cellsize, 
                         kernelR, moveargsiR, usermodel, 
                         moveargs, kernelp);
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
                            convolvemq(mm, kn, j, mqarrayR, kernelp, pm);
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
    
    // for 'shortcut' case (no individual variation) copy first pdt
    if (!individual)
        for (n=1; n<nc; n++) pdt[n] = pdt[0];
    
    return pdt;
}
//==============================================================================

// [[Rcpp::export]]
NumericVector PCH0secrjcpp (int type, int x, int nc, int jj,
			   const IntegerVector cumss, 
			    int kk, int mm, int cc0,
			   const IntegerVector PIA0, 
			   const NumericVector gk0, 
			   int binomN, 
			   const NumericMatrix Tsk) {
    
    const RVector<int> cumssR(cumss); 
    const RVector<int> PIA0R(PIA0); 
    const RVector<double> gk0R(gk0);
    const RMatrix<double> TskR(Tsk);
    
    // by primary session
    NumericVector pdt(nc * jj, 0.0);
    int j, m, n;
    std::vector<double> pjm(jj*mm);
    for (n = 0; n<nc; n++) {
	// primary-session-specific Pr for this animal
	pr0njmx(n, x, cumssR, nc, jj, kk, mm, cc0, binomN, PIA0R, gk0R, TskR, pjm);
	for (j=0; j<jj; j++) {
	    for (m=0; m<mm; m++) {
		pdt[j * nc + n] += pjm[m * jj + j] / mm;
            }
        }
    }
    return pdt;
}
//==============================================================================

struct pch1struct : public Worker {
    int x;
    int type;
    int jj;
    int kk;
    int mm;
    int nc;
    int kn;
    int cc0;
    const RVector<int> cumss;
    const RMatrix<double> openval0;
    const RVector<int> PIA0;
    const RVector<int> PIAJ;
    const RVector<double> gk0;
    int   binomN;
    const RMatrix<double> Tsk;
    const RVector<double> intervals;
    const RVector<int> moveargsi;
    int   movemodel;
    const RMatrix<int> kernel;
    const RMatrix<int> mqarray;
    double cellsize;
    
    // output likelihoods, one per animal
    RVector<double> output;
    
    // Constructor to initialize an instance of Somehistories 
    // The RMatrix class can be automatically converted to from the Rcpp matrix type
    pch1struct(
        int x, int type, int jj, int mm, int nc, 
        const IntegerVector cumss,
        const NumericMatrix openval0,
        const IntegerVector PIA0,
        const IntegerVector PIAJ,
        const NumericVector gk0,
        int binomN,
        const NumericMatrix Tsk,
        const NumericVector intervals,
        const IntegerVector moveargsi,
        int   movemodel,
        const IntegerMatrix kernel,
        const IntegerMatrix mqarray,
        double cellsize,        
        NumericVector output)    
        : 
        x(x), type(type), jj(jj), mm(mm), nc(nc), 
        cumss(cumss), openval0(openval0), PIA0(PIA0), PIAJ(PIAJ), 
        gk0(gk0), binomN(binomN), Tsk(Tsk), intervals(intervals), 
        moveargsi(moveargsi), movemodel(movemodel), kernel(kernel), mqarray(mqarray), 
        cellsize(cellsize), output(output) 
    {
        kn = kernel.nrow();
        cc0 = openval0.nrow();
        kk = Tsk.nrow();
    }
    
    double onepch1cpp (int n) {
        double pdt = 0.0;
        double pbd;
        int j,m;
        int b,d;
        double prw0, sumpj;
        
        // work vectors for session-specific real parameter values
        std::vector<double> phij(jj);      // each primary session
        std::vector<double> beta(jj);      // not used if type==1
        std::vector<double> pm(mm);
        std::vector<double> pjm(jj * mm);
        std::vector<double> moveargs(jj*2);
        std::vector<double> kernelp(kn * (jj-1));
        
        getphij (n, x, nc, jj, openval0, PIAJ,  intervals, phij);
        getbeta (type, n, x, nc, jj, openval0, PIAJ, intervals, phij, beta);
        // primary-session-specific Pr for this animal
        pr0njmx(n, x, cumss, nc, jj, kk, mm, cc0, binomN, PIA0, gk0, Tsk, pjm);
        if (movemodel>1) {
            getmoveargs (n, x, nc, jj, openval0, PIAJ, moveargsi, moveargs);
            fillkernelparallel (kn, jj, movemodel-2, cellsize, 
                         kernel, moveargsi, moveargs, kernelp);
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
                            convolvemq(mm, kn, j, mqarray, kernelp, pm);
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
                pdt += pbd * (1 - prw0);
            }
        }
	return pdt;
    }

    // function call operator that work for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {
        
        for (std::size_t n = begin; n < end; n++) {
            output[n] = onepch1cpp (n);
        }
    }
};

// [[Rcpp::export]]
NumericVector PCH1secrparallelcpp (int x, int type, int grain,
				   bool individual,
                                   int jj,  int mm, int nc,
                                   const IntegerVector cumss, 
                                   const NumericMatrix openval0, 
                                   const IntegerVector PIA0, 
                                   const IntegerVector PIAJ, 
                                   const NumericVector gk0, 
                                   int   binomN, 
                                   const NumericMatrix Tsk, 
                                   const NumericVector intervals,
                                   const IntegerVector moveargsi,
                                   int   movemodel,
                                   const IntegerMatrix kernel,
                                   const IntegerMatrix mqarray,
                                   double cellsize) {
    
    NumericVector output(nc); 
    
    // Construct and initialise
    pch1struct pch1 (x, type, jj, mm, nc, 
                     cumss, openval0, PIA0, PIAJ, gk0, binomN, Tsk, intervals,
                     moveargsi, movemodel, kernel, mqarray, 
                     cellsize, output);
    
    if (individual) {
	// Run operator() on multiple threads
	parallelFor(0, nc, pch1, grain);
    }
    else {
        // all the same: do once, then copy
	pch1.operator()(0,1);
	for (int i = 1; i < nc; i++) output[i] = output[0];
    }

    //pch1.operator()(0,nc);    // for debugging avoid multithreading to allow R calls
    
    // Return consolidated result
    return output;
}
//==============================================================================
