#include <Rcpp.h>
#include "utils.h"

using namespace Rcpp;

//==============================================================================

// [[Rcpp::export]]
double getNcpp ( int type, int nc, int ncf, int jj, int nmix, 
               const NumericVector& pmix, 
               const NumericVector& intervals, 
               const NumericVector& openval, int cc, 
               const IntegerVector& PIAJ) {
    
    // column 4 
    double N = 0;
    double *Bj;
    double *Nj;
    double *phij;
    int j;
    int x;
    if ((type == 2) | (type == 3) | (type == 4) | (type == 21) | (type == 22) | (type == 28) )
        N = openval[cc*3];  // + ncf;  
    else if (type == 18) {
        Bj = (double *) R_alloc (jj, sizeof(double));
        for (x=0; x < nmix; x++) {
            getDj (0, x, nc, jj, openval, cc, PIAJ, Bj);
            for (j=0; j < jj; j++) 
                N += pmix[x] * Bj[j];
        }
    }
    else if (type == 19) {
        Nj = (double *) R_alloc (jj, sizeof(double));
        phij = (double *) R_alloc (jj, sizeof(double));
        for (x=0; x < nmix; x++) {
            getphij (0, x, nc, jj, openval, cc, PIAJ, intervals, phij);
            getDj (0, x, nc, jj, openval, cc, PIAJ, Nj);
            N += pmix[x] * Nj[0];
            for (j=1; j < jj; j++) 
                N += pmix[x] * (Nj[j] - (Nj[j-1]*phij[j-1]));
        }
    }
    else stop("unrecognised type");
    return(N);
}
//--------------------------------------------------------------------

// column 4 
// superpopulation size given f, phi, N1 

// [[Rcpp::export]]
double getNBcpp (int n, int x, int nc, int jj, 
              const NumericVector& openval, int cc, 
              const IntegerVector& PIAJ, 
              const NumericVector& intervals) {
    int j;
    double *d;
    double *fj;
    double *phij = NULL;
    double *beta = NULL;
    double sumbeta;
    double N1 = 0;
    d = (double *) R_alloc (jj, sizeof(double));
    fj = (double *) R_alloc (jj, sizeof(double));
    phij = (double *) R_alloc (jj, sizeof(double));
    beta = (double *) R_alloc (jj, sizeof(double));
    getphij (n, x, nc, jj, openval, cc, PIAJ, intervals, phij);
    getfj (n, x, nc, jj, openval, cc, PIAJ, intervals, phij, fj);
    d[0] = 1;
    for (j = 1; j < jj; j++) {
        d[j] = d[j-1] * (phij[j-1] + fj[j-1]);
    }
    beta[0] = 1;
    sumbeta = beta[0];
    for (j = 1; j < jj; j++) {
        beta[j] = fj[j-1] * d[j-1]; 
        sumbeta += beta[j];
    }
    return(sumbeta * N1);
}

//--------------------------------------------------------------------

// [[Rcpp::export]]
double getDcpp ( int type, int nc, int jj, int nmix, 
               const NumericVector& pmix, 
	       const NumericVector& intervals, 
	       const NumericVector& openval, int cc, 
	       const IntegerVector& PIAJ) {
    int j,x;
    double *B;
    double *D;
    double sumB = 0;
    double *phij;
    phij = (double *) R_alloc (jj, sizeof(double));
    B = (double *) R_alloc (jj, sizeof(double));
    D = (double *) R_alloc (jj, sizeof(double));

    for (x=0; x < nmix; x++) {
        getphij (0, x, nc, jj, openval, cc, PIAJ, intervals, phij);
        if (type == 8) {
            getDj (0, x, nc, jj, openval, cc, PIAJ, D);
            B[0] = D[0];
            sumB += B[0] * pmix[x];
            for (j = 1; j < jj; j++) {
                B[j] = D[j] - D[j-1] * phij[j-1];
                sumB += B[j] * pmix[x];
            }
        }
        else {
            getDj (0, x, nc, jj, openval, cc, PIAJ, B);
            for (j = 0; j < jj; j++) {
                sumB += B[j] * pmix[x];
            }
        }
    }
    return (sumB);
}
//==============================================================================
