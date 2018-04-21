#include <Rcpp.h>
using namespace Rcpp;

#define fuzz 1e-200
#define huge 1e10

// index to vector element corresponding to cell i,j,k in 3D array
// stored in column-major order

int i3 (int i, int j, int k, int ii, int jj) {
    return(ii * (jj * k + j) + i);
}

// index to vector element corresponding to cell i,j,k,l in 4D array
// stored in column-major order

int i4 (int i, int j, int k, int l, int ii, int jj, int kk) {
    return (ii *(jj*(kk*l + k) + j) + i);
}

// probability of count
double gpois (int count, double lambda, int uselog)
{
    if (count == 0) {
        if (uselog)
            return (-lambda);
        else
            return (exp(-lambda));
    }
    else
        return (R::dpois(count, lambda, uselog));
}

double gbinom(int count, int size, double p, int uselog)
{
    double x;
    int i;
    if (count == 0) {
        p = 1 - p;
        x = p;
        for (i=1; i< size; i++) x = x*p;
        if (uselog) x = log(x);
        return (x);   // faster 
    }
    else
        return (R::dbinom (count, size, p, uselog));
}
double gnbinom (int count, int size, double mu, int uselog)
{
    // prob = size / (size + mu) 
    size = abs(size);  // in case negative 'binomN' passed 
    
    if (count == 0) {  // faster - added 2010-10-11
    // for portability cast integer to double before taking log
    // cf Plummer R Journal 2011
    // if (uselog) return( log(size/(size+mu)) * log(size));
    if (uselog) return( log(size/(size+mu)) * log(static_cast<double>(size)));
    else return (pow(size/(size+mu), size));
    }
    else
        return (R::dnbinom (count, size, size/(size+mu), uselog));
}

double gbinomFP (int count, double size, double p, int uselog)
// allow non integers 
{
    return ( lgamma(size+1) - lgamma(size-count+1) - lgamma(count+1) +
             count * log(p) + (size - count) * log (1-p) );
}

double countp (int count, int binomN, double lambda) {
    // Poisson 
    if (binomN == 0) {
        return ( gpois (count, lambda, 0));
    }
    // Bernoulli 
    else if (binomN == 1) {
        if (count == 0)
            return ( 1 - lambda );
        else
            return ( lambda );
    }
    // negative binomial 
    else if (binomN < 0) {
	return ( gnbinom (count, binomN, lambda, 0) );
    }
    // binomial 
    else {
	return ( gbinom (count, binomN, lambda, 0) );
    }
}
//--------------------------------------------------------------------------

double d2 (
        int k,
        int m,
        const NumericVector& A1,
        const NumericVector& A2,
        int A1rows,
        int A2rows)
  
// return squared distance between two points given by row k in A1
// and row m in A2, where A1 and A2 have respectively A1rows and A2rows
      
{
    return(
        (A1[k] - A2[m]) * (A1[k] - A2[m]) +
            (A1[k + A1rows] - A2[m + A2rows]) * (A1[k + A1rows] - A2[m + A2rows])
    );
}
//--------------------------------------------------------------------------

// parameters in openval ordered g0, phi, f, N, sigma, pmix 

// hazard 
double hfn
     (int k, int m, int c, 
      const NumericVector& openval, int cc, 
      const NumericVector& traps,
      const NumericVector& mask, int kk, int mm, int sigmai, int detectfn)
{
    double d;
    double sigma;
    double z;

    sigma =  openval[cc*sigmai + c];
    d = sqrt(d2(k, m, traps, mask, kk, mm));

    // HHN
    if (detectfn == 14) { 
        return (openval[c] * exp(-d*d/2/sigma/sigma));
    }
    // HHR
    else if (detectfn == 15) {
	z =  openval[cc*(sigmai+1) + c];
	return (openval[c] * (1 - exp(-pow(d/sigma, -z))));
    }
    // HEX
    else if (detectfn == 16) {
        return (openval[c] * exp(-d/sigma));
    }
    // HAN
    else if (detectfn == 17) {
	z =  openval[cc*(sigmai+1) + c];
        return (openval[c] * exp(-(d-z)*(d-z) / 2 / sigma / sigma));
    }
    // HCG
    else if (detectfn == 18) {
	z =  openval[cc*(sigmai+1) + c];
        return (openval[c] * R::pgamma(d,z,sigma/z,0,0)); 
    }
    // HVP
    else if (detectfn == 19) {
	z = openval[cc*(sigmai+1) + c];
        return (openval[c] * exp(- pow(d /sigma , z)));
    }
    else stop("detectfn not allowed in openCR");
}

void getp (int n, int x, int nc, int ss, 
            const NumericVector& openval, int cc, 
            const IntegerVector& PIA, 
            double p[]) {
     // column 1 
     int s;
     for (s = 0; s < ss; s++) {
         p[s] = openval[PIA[i3(n, s, x, nc, ss )]-1]; 
     }
 }
//--------------------------------------------------------------------

void getphij (int n, int x, int nc, int jj, 
             const NumericVector& openval, int cc, 
             const IntegerVector& PIAJ,
             const NumericVector& intervals, 
             double phij[]) {
    // column 2 
    int j;
    double phi;
    for (j = 0; j < (jj-1); j++) {
        // jj-1 because one fewer intervals than primary sessions  
        phi = openval[cc + PIAJ[i3(n, j, x, nc, jj)]-1];       
        // adjust for interval duration  
        phij[j] = exp(log(phi) * intervals[j]);  
    }
    phij[jj-1] = 0;
}
//--------------------------------------------------------------------

void getmoveargs (int n, int x, int nc, int jj, 
		const NumericVector& openval, int cc, 
		const IntegerVector& PIAJ,
		const IntegerVector& moveargsi,
		double moveargs[]) {
    // column moveargsi (and maybe moveargsi + 1) 
    int j;
    for (j = 0; j < (jj-1); j++) {
        // jj-1 because one fewer intervals than primary sessions  
        moveargs[j] = openval[cc*moveargsi[0] + PIAJ[i3(n, j, x, nc, jj)]-1];  
        if (moveargsi[1]>0)
            moveargs[j+jj] = openval[cc*moveargsi[1] + PIAJ[i3(n, j, x, nc, jj)]-1];  
    }
    moveargs[jj-1] = 0;
    moveargs[2*jj-1] = 0;
}
//--------------------------------------------------------------------

void getpj (int n, int x, int nc, int jj, 
             const NumericVector& openval, int cc, 
             const IntegerVector& PIAJ,
             double pj[]) {
    // column 2 
    int j;
    for (j = 0; j < jj; j++) {
        pj[j] = openval[PIAJ[i3(n, j, x, nc, jj)]-1];       
    }
}
//--------------------------------------------------------------------

void getg (int type, int n, int x, int nc, int jj, 
             const NumericVector& openval, int cc, 
             const IntegerVector& PIAJ,
             double g[]) {
    // column 4 
    int j;
    for (j = 0; j < jj; j++) {
        if (type != 27)
            g[j] = 0;
        else
            g[j] = openval[cc*3 + PIAJ[i3(n,j, x, nc, jj)]-1];       
    }
}
//--------------------------------------------------------------------

void getfj (int n, int x, int nc, int jj, 
           const NumericVector& openval, int cc, 
           const IntegerVector& PIAJ, 
           const NumericVector& intervals, 
           double phij[], double fj[]) {
    // column 3 
    int j;
    double f,phi;
    for (j = 0; j < (jj-1); j++) {
        // jj-1 because one fewer intervals than primary sessions  
        f = openval[cc*2 + PIAJ[i3(n, j, x, nc, jj )]-1]; 
	phi = exp(log(phij[j]) / intervals[j]);
        // adjust for interval duration  
        fj[j] = exp(log(phi+f) * intervals[j]) - phij[j];  
    }
    fj[jj-1] = 0;
}
//--------------------------------------------------------------------

void getlj (int n, int x, int nc, int jj, 
           const NumericVector& openval, 
           int cc, 
           const IntegerVector& PIAJ, 
           const NumericVector& intervals, 
           double lj[]) {
    // column 3 
    int j;
    double l;
    for (j = 0; j < (jj-1); j++) {
        // jj-1 because one fewer intervals than primary sessions  
        l = openval[cc*2 + PIAJ[i3(n, j, x, nc, jj )]-1]; 
        // adjust for interval duration  
        lj[j] = exp(log(l) * intervals[j]);  
    }
    lj[jj-1] = 0;
}
//--------------------------------------------------------------------

int sumj (int uv[], int j, int k) {
    int i;
    int sum = 0;
    if (j>k)
        return (0);
    else {
        for (i=j; i<=k; i++)
            sum += uv[i];      // use 0:(J-1) indices 
	return(sum);
    }
}
//--------------------------------------------------------------------

void getgaml (int n, int x, int nc, int jj, 
             const NumericVector& openval, int cc, 
             const IntegerVector& PIAJ, 
             const NumericVector& intervals, 
	     double gam[]) {
    // column 3 
    int j;
    double phij;
    double lamj;
    for (j = 0; j < (jj-1); j++) {
        phij = openval[cc + PIAJ[i3(n, j, x, nc, jj)]-1]; 
        phij = exp(log(phij) * intervals[j]);  
        lamj = openval[cc*2 + PIAJ[i3(n, j, x, nc, jj)]-1]; 
        lamj = exp(log(lamj) * intervals[j]);  
        gam[j+1] = phij/lamj;
    }
    gam[0] = 0;
}
//--------------------------------------------------------------------

void getgamj (int n, int x, int nc, int jj, 
             const NumericVector& openval, int cc, 
             const IntegerVector& PIAJ, 
             const NumericVector& intervals, 
	     double gamj[]) {
    // column 3 
    int j;
    double gam;
    for (j = 1; j < jj; j++) {
        gam = openval[cc*2 + PIAJ[i3(n, j, x, nc, jj)]-1]; 
        gamj[j] = exp(log(gam) * intervals[j-1]);  
    }
    gamj[0] = 0;
}
//--------------------------------------------------------------------

void getkapj (int n, int x, int nc, int jj, 
             const NumericVector& openval, int cc, 
             const IntegerVector& PIAJ, 
	     double kapj[]) {
    // column 3 
    int j;
    for (j = 1; j < jj; j++) {
        kapj[j] = openval[cc*2 + PIAJ[i3(n, j, x, nc, jj)]-1];
    }
    kapj[0] = 1;
}
//--------------------------------------------------------------------

// getgam <- function (n, x, openval, PIAJ, intervals) {
//     J2 <- 2:(length(intervals)+1)
//     c(0,openval[PIAJ[n, J2, x],3])
// }

void getbeta0 (int n, int x, int nc, int jj, 
              const NumericVector& openval, int cc, 
              const IntegerVector& PIAJ, 
              double beta[]) {
    // column 3 
    int j;
    double sumbeta = 0;
    for (j = 1; j < jj; j++) {
        beta[j] = openval[cc*2 + PIAJ[i3(n, j, x, nc, jj )]-1]; 
        sumbeta += exp(beta[j]);
    }
    beta[0] = 1;
    for (j = 1; j < jj; j++) {
        beta[j] = exp(beta[j]) / (1 + sumbeta);
        beta[0] -= beta[j];
    }
}

//--------------------------------------------------------------------

void gettau (int n, int x, int nc, int jj,
             double openval[], int cc, int PIAJ[], 
	     double tau[], int M) {
    // column 5 
    int j;
    double sumtau = 0;
    for (j = 0; j < M; j++) {
        tau[j] = openval[cc*4 + PIAJ[i3(n, j, x, nc, jj )]-1];
        sumtau += exp(tau[j]);
    }
    tau[M] = 1;
    for (j = 0; j < M; j++) {
        tau[j] = exp(tau[j]) / (1 + sumtau);
        tau[M] -= tau[j];
    }
    for (j = M+1; j < jj; j++)
        tau[j] = 0;
}
//--------------------------------------------------------------------

void getDj (int n, int x, int nc, int jj, 
            const NumericVector& openval, int cc, 
            const IntegerVector& PIAJ, 
            double Dj[]) {
    // column 3 
    int j;
    for (j = 0; j < jj; j++) {
        Dj[j] = openval[cc*2 + PIAJ[i3(n, j, x, nc, jj )]-1]; 
    }
}
//--------------------------------------------------------------------

// per capita recruitment cf Link & Barker 2005, Schwarz 'Gentle Intro'
void getbetaf (int n, int x, int nc, int jj, 
               const NumericVector& openval, int cc, 
               const IntegerVector& PIAJ, 
               double phij[], 
	       const NumericVector& intervals, double beta[]) {
    int j;
    double *d;
    double *fj;
    double sumbeta = 1;
    d = (double *) R_alloc (jj, sizeof(double));
    fj = (double *) R_alloc (jj, sizeof(double));
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
    for (j = 0; j < jj; j++) {
        beta[j] = beta[j] / sumbeta;
    }
}
//--------------------------------------------------------------------

void getbetal (int n, int x, int nc, int jj, 
               const NumericVector& openval, int cc, 
               const IntegerVector& PIAJ, 
               double phij[], 
	       const NumericVector& intervals, double beta[]) {
    int j;
    double *d;
    double *fj;
    double *lambdaj;
    double sumbeta = 1;
    d = (double *) R_alloc (jj, sizeof(double));
    fj = (double *) R_alloc (jj, sizeof(double));
    lambdaj = (double *) R_alloc (jj, sizeof(double));
    getlj (n, x, nc, jj, openval, cc, PIAJ, intervals, lambdaj);
    for (j=0; j<jj; j++) {
        if (lambdaj[j] < phij[j])  
            fj[j] = 0;    
        else
            fj[j] = lambdaj[j] - phij[j];
    }
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
    for (j = 0; j < jj; j++) {
        beta[j] = beta[j] / sumbeta;
    }
}
//--------------------------------------------------------------------

void getbetag (int n, int x, int nc, int jj, 
               const NumericVector& openval, int cc, 
               const IntegerVector& PIAJ, 
               double phij[], 
	       const NumericVector& intervals, 
	       double beta[]) {
    int j;
    double *d;
    double *fj;
    double *gamj;
    double sumbeta = 1;
    d = (double *) R_alloc (jj, sizeof(double));
    fj = (double *) R_alloc (jj, sizeof(double));
    gamj = (double *) R_alloc (jj, sizeof(double));
    getgamj (n, x, nc, jj, openval, cc, PIAJ, intervals, gamj);

    for (j=1; j<jj; j++) {
        if (gamj[j] <= 0)  
            fj[j-1] = 0;    
        else
            fj[j-1] = phij[j-1] * (1/gamj[j] - 1);   // Pradel 1996 p 708 corrected!
    }
    fj[jj-1] = 0;

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
    for (j = 0; j < jj; j++) {
        beta[j] = beta[j] / sumbeta;
    }
}
//--------------------------------------------------------------------

void getbetak (int n, int x, int nc, int jj, 
               const NumericVector& openval, int cc, 
               const IntegerVector& PIAJ, 
               double phij[], 
	       double beta[]) {
    int i,j;
    double *tau;
    double *pj;
    double *fprod;
    double *fj;
    double *kapj;
    tau = (double *) R_alloc (jj, sizeof(double));
    pj = (double *) R_alloc (jj, sizeof(double));
    fprod = (double *) R_alloc (jj, sizeof(double));
    fj = (double *) R_alloc (jj, sizeof(double));
    kapj = (double *) R_alloc (jj, sizeof(double));

    getkapj (n, x, nc, jj, openval, cc, PIAJ, kapj);
    getpj (n, x, nc, jj, openval, cc, PIAJ, pj);

    tau[0] = 1/pj[0];
    for (j=0; j<(jj-1); j++) {
        fj[j] = (kapj[j+1] - kapj[j]/pj[j] * (1 - pj[j]) * phij[j] * pj[j+1]) / 
	    (tau[j] * pj[j+1]);
        tau[j+1] = tau[0];  // get next tau
	for (i=0; i<(j+1); i++) tau[j+1] *= (phij[i] + fj[i]);
    }
    for (j=1; j<(jj-1); j++) {
	fprod[j] = fj[j];
	for (i=0; i<j; i++) fprod[j] *= (phij[i] + fj[i]);
    }
    beta[0] = 1 + fj[0];
    for (j=1; j<(jj-1); j++) beta[0] += fprod[j];
    beta[0] = 1/beta[0];
    beta[1] = beta[0] * fj[0];
    for (j=1; j<(jj-1); j++) {
	beta[j+1] = beta[0] * fprod[j];
    }
    // for (j=0; j<jj; j++) Rprintf(" %8.6f ", beta[j]);
    // Rprintf("\n");
}
//--------------------------------------------------------------------

// getbetag <- function (n, x, openval, PIAJ, phi, intervals) {
//     J <- length(intervals)+1
//     J1 <- 1:(J-1)
//     gam <- getgam (n, x, openval, PIAJ, intervals)
//     f <- ifelse(gam<=0, 0, 1/gam - 1)[2:J]
//     d <- c(1, cumprod(phi[J1]+f))[J1]
//     beta <- f * d 
//     beta <- c(1, beta)
//     beta/sum(beta)
// }

// return parameterisation cf Pledger et al. 2010 p 885 
void getbetaB (int n, int x, int nc, int jj, 
               const NumericVector& openval, int cc, 
               const IntegerVector& PIAJ,  
               double beta[]) {
    int j;
    double *B;
    double sumB = 0;
    B = (double *) R_alloc (jj, sizeof(double));
    getDj (n, x, nc, jj, openval, cc, PIAJ, B);
    for (j = 0; j < jj; j++) {
        sumB += B[j];
    }
    for (j = 0; j < jj; j++) {
        beta[j] = B[j] / sumB;
    }
}
//--------------------------------------------------------------------

void getbetaD (int n, int x, int nc, int jj, 
               const NumericVector& openval, int cc, 
               const IntegerVector& PIAJ,
               double phij[], double beta[]) {
    int j;
    double *B;
    double *D;
    double sumB;
    B = (double *) R_alloc (jj, sizeof(double));
    D = (double *) R_alloc (jj, sizeof(double));
    getDj (n, x, nc, jj, openval, cc, PIAJ, D);
    
    B[0] = D[0];
    sumB = B[0];
    for (j = 1; j < jj; j++) {
        B[j] = D[j] - D[j-1] * phij[j-1];
        sumB += B[j];
    }
    for (j = 0; j < jj; j++) {
        beta[j] = B[j] / sumB;
    }
}
//--------------------------------------------------------------------

void getbeta (int type, int n, int x, int nc, int jj, 
                  const NumericVector& openval, int cc, 
                  const IntegerVector& PIAJ,
                  const NumericVector& intervals,
                  double phij[], double beta[]) {
    if ((type == 2) || (type == 17) || (type == 11) || (type == 13) || 
             (type == 41) || (type == 43)) 
        getbeta0 (n, x, nc, jj, openval, cc, PIAJ, beta);
    else if ((type == 4) || (type == 15) ||  (type == 27) ||  (type == 7) || (type == 9) ||
	     (type == 37) || (type == 39))
        getbetaf (n, x, nc, jj, openval, cc, PIAJ, phij, intervals, beta);
    else if ((type == 3) || (type == 16) || (type == 10) || (type == 12) || (type == 20) ||
	     (type == 40) || (type == 42))
        getbetal (n, x, nc, jj, openval, cc, PIAJ, phij, intervals, beta);
    else if ((type == 14) || (type == 18))
        getbetaB (n, x, nc, jj, openval, cc, PIAJ, beta);
    else if ((type == 8) || (type == 19) || (type == 38))
        getbetaD (n, x, nc, jj, openval, cc, PIAJ, phij, beta);
    else if ((type == 22) || (type == 23) || (type == 24) || (type == 25) || (type == 26))
        getbetag (n, x, nc, jj, openval, cc, PIAJ, phij, intervals, beta);
    else if ((type == 28) || (type == 29))
        getbetak (n, x, nc, jj, openval, cc, PIAJ, phij, beta);
    else stop("no beta for this type");
}
//--------------------------------------------------------------------
