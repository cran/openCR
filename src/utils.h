#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
using namespace Rcpp;

#define fuzz 1e-200
#define huge 1e10

int i3 (int i, int j, int k, int ii, int jj);
int i4 (int i, int j, int k, int l, int ii, int jj, int kk);
double gpois (int count, double lambda, int uselog);
double gbinom(int count, int size, double p, int uselog);
double gnbinom (int count, int size, double mu, int uselog);
double gbinomFP (int count, double size, double p, int uselog);
double countp (int count, int binomN, double lambda);
double d2 (int k,
        int m,
        const NumericVector& A1,
        const NumericVector& A2,
        int A1rows,
        int A2rows);
double hfn
     (int k, int m, int c, 
      const NumericVector& openval, int cc, 
      const NumericVector& traps,
      const NumericVector& mask, 
      int kk, int mm, int sigmai, int detectfn);
void getp (int n, int x, int nc, int ss, 
            const NumericVector& openval, int cc, 
            const IntegerVector& PIA, 
            double p[]);

void getpj (int n, int x, int nc, int jj, 
             const NumericVector& openval, int cc, 
             const IntegerVector& PIAJ,
             double pj[]);
void getphij (int n, int x, int nc, int jj, 
             const NumericVector& openval, int cc, 
             const IntegerVector& PIAJ,
             const NumericVector& intervals, 
	      double phij[]);
void getmoveargs (int n, int x, int nc, int jj, 
		const NumericVector& openval, int cc, 
		const IntegerVector& PIAJ,
		const IntegerVector& moveargsi,
		double moveargs[]);
void getgamj (int n, int x, int nc, int jj, 
             const NumericVector& openval, int cc, 
             const IntegerVector& PIAJ,
             const NumericVector& intervals, 
             double phij[]);
void getkapj (int n, int x, int nc, int jj, 
             const NumericVector& openval, int cc, 
             const IntegerVector& PIAJ,
             double kapj[]);
void getg (int type, int n, int x, int nc, int jj, 
           const NumericVector& openval, int cc, 
           const IntegerVector& PIAJ,
           double g[]);
void getfj (int n, int x, int nc, int jj, 
           const NumericVector& openval, int cc, 
           const IntegerVector& PIAJ, 
           const NumericVector& intervals, 
           double phij[], double fj[]);
void getlj (int n, int x, int nc, int jj, 
           const NumericVector& openval, 
           int cc, 
           const IntegerVector& PIAJ, 
           const NumericVector& intervals, 
           double lj[]);
int sumj (int uv[], int j, int k);
void getgaml (int n, int x, int nc, int jj, 
             const NumericVector& openval, int cc, 
             const IntegerVector& PIAJ, 
             const NumericVector& intervals, double gam[]);
void getgam (int n, int x, int nc, int jj, 
             const NumericVector& openval, int cc, 
             const IntegerVector& PIAJ, 
             const NumericVector& intervals, double gam[]);
void getbeta0 (int n, int x, int nc, int jj, 
              const NumericVector& openval, int cc, 
              const IntegerVector& PIAJ, 
              double beta[]);
void gettau (int n, int x, int nc, int jj,
             const NumericVector& openval, int cc, 
             const IntegerVector& PIAJ, 
             double tau[], int M);
void getDj (int n, int x, int nc, int jj, 
            const NumericVector& openval, int cc, 
            const IntegerVector& PIAJ,
            double Dj[]);

void getbeta (int type, int n, int x, int nc, int jj, 
              const NumericVector& openval, int cc, 
              const IntegerVector& PIAJ,
              const NumericVector& intervals,
              double phij[], double beta[]);

void convolvemq (
    int    mm,        /* number of points on mask */
    int    kn,        /* number of points on kernel */
    int    j,         /* session number 1..jj */
    double kernelp[], /* p(move|dx,dy) for points in kernel */
    const IntegerVector&  mqarray, /* input */
    double pjm[]      /* return value */
    );

void fillkernelp (
    int kn, 
    int jj, 
    int kerneltype, 
    const IntegerVector& kernel, 
    double cellsize,
    const IntegerVector& moveargsi,
    double moveargs[], 
    const CharacterVector& usermodel,
    double kernelp[]);

#endif
