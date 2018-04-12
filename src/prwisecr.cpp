#include <Rcpp.h>
#include "utils.h"
#include <string>

using namespace Rcpp;

//==============================================================================

void convolvemq (
    int    mm,        /* number of points on mask */
    int    kn,        /* number of points on kernel */
    int    j,         /* session number 1..jj */
    double kernelp[], /* p(move|dx,dy) for points in kernel */
    const IntegerVector&  mqarray, /* input */
    double pjm[]      /* return value */
    )
{
    int m, q, mq;
    double *workpjm = NULL;
    workpjm = (double *) S_alloc (mm, sizeof(double));

    /* convolve movement kernel and pjm... */
    for (m = 0; m < mm; m++) {
	for (q=0; q < kn; q++) {           /* over movement kernel */
	    mq = mqarray[q * mm + m];  
	    if (mq >= 0) {                /* post-dispersal site is within mask */
		if (mq>mm) {
		    Rprintf ("%5d %5d\n", m, q);
		    stop("mq > mm");
		}
		workpjm[mq] += pjm[m] * kernelp[kn * (j-1) + q];   /* probability of this move */
	    }
	}
    }
    for (m = 0; m < mm; m++) 
        pjm[m] = workpjm[m];
}

void fillkernelp (int kn, int jj, 
		  int kerneltype, 
		  const IntegerVector& kernel, 
		  double cellsize,
		  const IntegerVector& moveargsi, 
                  double moveargs[], 
		  const CharacterVector& fnname,
		  double kernelp[]) {
    int j,k;
    double r,r2;
    NumericVector p;
    double *sumj;
    std::string fn;
    sumj = (double *) R_alloc(jj, sizeof(double));
    for (j = 0; j < (jj-1); j++) sumj[j] = 0;
    for (k = 0; k < kn; k++) {
        r2 = (kernel[k]*kernel[k] + kernel[k+kn]*kernel[k+kn]) * cellsize * cellsize;
        r = sqrt(r2);
        for (j = 0; j < (jj-1); j++) {
    // Rprintf(" k %4d j %4d  moveargs[j,] %8.6f %8.6f \n",  k,j,moveargs[j], moveargs[j+jj]);
            if (kerneltype == 0)        /* Gaussian kernel */
                kernelp[j * kn + k] = exp(-r2 / 2 / moveargs[j] / moveargs[j]);
	    else if (kerneltype == 1)   /* Negative exponential kernel */
		kernelp[j * kn + k] = exp(-r / moveargs[j]);
	    else if (kerneltype == 2) {  /* User kernel */
                // call R function from C++
		fn =  as<std::string>(fnname[0]);
		Environment env = Environment::global_env();
		Function f = env[fn];
		if (moveargsi[1]>0)
		    p = f(r, moveargs[j], moveargs[j+jj]);
		else if (moveargs[0]>0)
		    p = f(r, moveargs[j]);
		else 
		    p = f(r);
		kernelp[j * kn + k] = p[0];

    // Rprintf(" k %4d j %4d  kernelp[j * kn + k] %8.6f\n",  k,j,kernelp[j * kn + k]); 
	    }
	    else stop("unrecognised kerneltype");
	    sumj[j] += kernelp[j * kn + k];
        }
    }
    /* normalise */
    for (k = 0; k < kn; k++) 
	for (j = 0; j < (jj-1); j++) {
	    kernelp[j * kn + k] = kernelp[j * kn + k] / sumj[j];
	}
}

void prw (int n, int j, int x, int nc, int jj, int kk, int mm, int cc, int binomN, 
           const IntegerVector& cumss, 
           const IntegerVector& PIA, 
           const NumericVector& gk,
           const NumericVector& Tsk, 
           const IntegerVector& w,
           double pjm[]) {
    int c, k, m, s, wxi, gi, wi, count;
    double Tski;
    double size = 1;
    bool dead = false;
    for (s = cumss[j-1]; s < cumss[j]; s++) {
	for (k=0; k<kk; k++) {
	    wxi =  i4(n, s, k, x, nc, cumss[jj], kk);
	    c = PIA[wxi] - 1;
	    if (c >= 0) {    // drops unset traps 
		Tski = Tsk[s * kk + k];
		if (fabs(Tski-1) > 1e-10) {                  // effort <> 1.0 
		    if (binomN==1) 
			size = (int) Tski;
		    else 
			size = binomN;
		}
		wi = i3(n, s, k, nc, cumss[jj]);
		count = w[wi];
		if (count<0) {count = -count; dead = true; }
		for (m=0; m<mm; m++) {
		    gi  = i3(c, k, m, cc, kk);
		    pjm[m] *= countp (count, size, gk[gi]);
		}
	    }
	}
	if (dead) break;   // after processing all traps on this occasion
    }
 }
//==============================================================================

// [[Rcpp::export]]
double prwisecrcpp (int type, int n, int x, int nc, int jj, 
		       int kk, int mm, int nmix, 
		       const IntegerVector& cumss, 
		       const IntegerVector& w, 
		       const IntegerVector& fi, 
		       const IntegerVector& li,
		       const NumericVector& gk, 
		       const NumericVector& openval, 
		       int cc, 
		       const IntegerVector& PIA,
		       const IntegerVector& PIAJ,
		       int binomN, 
		       const NumericVector& Tsk, 
		       const NumericVector& intervals,
		       const IntegerVector& moveargsi, 
		       int   CJSp1,
		       int movemodel,
		       const CharacterVector& usermodel,
		       int kn,
		       const IntegerVector& kernel,
		       const IntegerVector& mqarray,
		       double cellsize) {

    // merging prwiCJSsecr, prwiJSSAsecr 2018-03-26

    double *beta = NULL;
    double *phij = NULL;
    double *moveargs = NULL;
    double *kernelp = NULL;
    double pdt;
    double pbd;
    int j, m;
    int b, minb, maxb;
    int d, mind, maxd;
    int cjs = 0;  // CJS offset
    double *pjm;
    double sump, prwi;
    pjm = (double *) R_alloc (mm, sizeof(double));
    phij = (double *) R_alloc (jj, sizeof(double));
    getphij (n, x, nc, jj, openval, cc, PIAJ, intervals, phij);
    if (movemodel > 1) {
	moveargs = (double *) R_alloc (jj*2, sizeof(double));
	kernelp = (double *) R_alloc (kn*(jj-1), sizeof(double));
	getmoveargs (n, x, nc, jj, openval, cc, PIAJ, moveargsi, moveargs);
	fillkernelp (kn, jj, movemodel-2, kernel, cellsize,
		     moveargsi, moveargs, usermodel, kernelp);
    }
    
    if (type==6) {     // CJSsecr
	minb = fi[n];
	cjs = (1-CJSp1);
    }
    else {
	minb = 1;
	beta = (double *) R_alloc (jj, sizeof(double));
	getbeta (type, n, x, nc, jj, openval, cc, PIAJ, intervals, phij, beta);
    }
    maxb = fi[n];

    mind = abs(li[n]);
    maxd = jj;
    if (li[n] < 0) maxd = mind; // possible censoring

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
	    if (li[n]>0)    // not censored
		pbd *= 1-phij[d-1];

	    if ((movemodel==0) || (movemodel>1)) {
		for (m=0; m < mm; m++) pjm[m] = 1;   // uniform
		for (j = b + cjs; j <= d; j++) {
		    prw (n, j, x, nc, jj, kk, mm, cc, binomN, cumss, PIA, gk, 
			 Tsk, w, pjm);
		    if ((j<d) && (movemodel>1)) 
			convolvemq(mm, kn, j, kernelp, mqarray, pjm);
		}
		prwi = 0;
		for (m=0; m<mm; m++) prwi += pjm[m] / mm;
	    }
	    else {  // movemodel == 1
		prwi = 1.0;
		for (j = b + cjs; j <= d; j++) {
		    for (m=0; m<mm; m++) pjm[m] = 1;  // uniform
		    prw(n, j, x, nc, jj, kk, mm, cc, binomN, cumss, PIA, gk, 
			Tsk, w, pjm);
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

// [[Rcpp::export]]
double allhistcpp (int type,int nc, int jj, 
		   int kk, int mm, int nmix, 
		   const IntegerVector& cumss, 
		   const IntegerVector& w,
		   const IntegerVector& fi, 
		   const IntegerVector& li,
		   const NumericVector& gk, 
		   const NumericVector& openval, int cc, 
		   const IntegerVector& PIA, 
		   const IntegerVector& PIAJ, 
		   int binomN, 
		   const NumericVector& Tsk, 
		   const NumericVector& intervals, 
		   const IntegerVector& moveargsi, 
		   int   CJSp1,                    // CJSp1 not used for JSSA
		   int movemodel,
		   const CharacterVector& usermodel,
		   int kn,
		   const IntegerVector& kernel,
		   const IntegerVector& mqarray,
		   double cellsize,
		   const NumericVector& pmix) {
    int n,x;
    double onep = 1.0;
    double sump = 0.0;
    for (x=0; x < nmix; x++) {
        for (n=0; n < nc; n++) {
	    onep = prwisecrcpp (type, n, x, nc, jj, kk, mm, nmix, cumss, w, fi, li,
				gk, openval, cc, PIA, PIAJ, binomN, Tsk, intervals,
				moveargsi, CJSp1, movemodel, usermodel, kn, kernel, 
				mqarray, cellsize);           
            sump += log(onep * pmix[n * nmix + x]);
            // Rprintf("n %4d onep %10.8f \n", n, onep);
        }
    }
    return sump;
}

//==============================================================================
