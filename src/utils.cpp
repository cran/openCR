#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

#define fuzz 1e-200
#define huge 1e10

// negative binomial suppressed 2018-04-17

//--------------------------------------------------------------------------
// index to vector element corresponding to cell i,j,k in 3D array
// stored in column-major order

int i3 (int i, int j, int k, int ii, int jj) {
    return(ii * (jj * k + j) + i);
}

//--------------------------------------------------------------------------
// index to vector element corresponding to cell i,j,k,l in 4D array
// stored in column-major order

int i4 (int i, int j, int k, int l, int ii, int jj, int kk) {
    return (ii *(jj*(kk*l + k) + j) + i);
}

//--------------------------------------------------------------------------

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
//--------------------------------------------------------------------------

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
//--------------------------------------------------------------------------

double gbinomFP (int count, double size, double p, int uselog)
// allow non integer size for Binomial
{
    return ( lgamma(size+1) - lgamma(size-count+1) - lgamma(count+1) +
             count * log(p) + (size - count) * log (1-p) );
}
//--------------------------------------------------------------------------

// probability of count for session s, detector k, animal i
// The argument 'g' is understood to be a cumulative hazard if binomN=0,
// a probability otherwise

double pski ( int binomN,
              int count,
              double Tski,
              double g) {
    
    double result = 1.0;
    
    if (binomN == -1) {                              // binary proximity detectors : Bernoulli
        if (abs(Tski-1) > 1e-10) {                   // effort not unity; adjust g 
            g = 1 - pow(1 - g, Tski);
        }
        if (count>0)                                 
            result = g;  
        else 
            result = 1 - g;
    }
    else if (binomN == 0) {                          // count detectors : Poisson 
        if (count == 0) 
            result = exp(-Tski * g);                 // routinely apply Tsk adjustment to cum. hazard 
        else
            result = R::dpois(count, Tski * g, 0); 
    }
    else if (binomN == 1) {                          // count detectors : Binomial, size from Tsk
        result = gbinom (count, round(Tski), g, 0); 
    }
    else if (binomN > 1) {                           // count detectors : Binomial, specified size 
        if (abs(Tski-1) > 1e-10) {                   // effort not unity, adjust g 
            g = 1 - pow(1 - g, Tski);
        }
        result = gbinom (count, binomN, g, 0);
    }
    else stop("binomN < -1 not allowed");  // code multi -2 separately

    return (result);
}
//--------------------------------------------------------------------------

// distance between two points given by row k in traps and row m in mask
double dkm (int k, int m, RMatrix<double> traps, RMatrix<double> mask)
{
    return(sqrt((traps(k,0) - mask(m,0)) * (traps(k,0) - mask(m,0)) +
		(traps(k,1) - mask(m,1)) * (traps(k,1) - mask(m,1))));
}
//--------------------------------------------------------------------------

// parameters in openval ordered g0, phi, f, N, sigma, pmix 

// hazard detection functions 14-19
double hfn
    (int k, int m, int c, 
     const RMatrix<double> openval,  
     const RMatrix<double> traps,
     const RMatrix<double> mask, 
     int sigmai, 
     int detectfn)
{
    double d;
    double sigma;
    double z = 1;
    sigma =  openval(c, sigmai);
    d = dkm(k, m, traps, mask);
    if (detectfn == 14) { 
        return (openval(c,0) * exp(-d*d/2/sigma/sigma));                   // HHN
    }
    else if (detectfn == 16) {
        return (openval(c,0) * exp(-d/sigma));                             // HEX
    }
    else {
        z =  openval(c,sigmai+1);
        if (detectfn == 15) {
            return (openval(c,0) * (1 - exp(-pow(d/sigma, -z))));          // HHR
        }
        else if (detectfn == 17) {
            return (openval(c,0) * exp(-(d-z)*(d-z) / 2 / sigma / sigma)); // HAN
        }
        else if (detectfn == 18) {
            return (openval(c,0) * R::pgamma(d,z,sigma/z,0,0));            // HCG
        }
        else if (detectfn == 19) {
            return (openval(c,0) * exp(- pow(d /sigma , z)));              // HVP
        }
        else stop("detectfn not allowed in openCR");
    }
}
//--------------------------------------------------------------------------

void convolvemq (
        int    mm,                    // number of points on mask 
        int    kn,                    // number of points on kernel
        int    j,                     // session number 1..jj 
        int    edgecode,              // 0 none, no action; 1 wrapped, no action; 2 normalize truncated kernel
        const  RMatrix<int> &mqarray, // input [& 2020-10-31]
        std::vector<double> &kernelp, // p(move|dx,dy) for points in kernel 
        std::vector<double> &pjm      // return value
)
{
    int m, q, mq;
    double sump;
    std::vector<double> workpjm(mm);
    
    // convolve movement kernel and pjm... 
    for (m = 0; m < mm; m++) {
        if (edgecode == 2) {
            // 2020-10-29 adjust for edge-truncated kernel cf convolvemqold
            sump = 0;
            for (q=0; q < kn; q++) {           // over movement kernel 
                if (mqarray(m,q) >= 0) {       // post-dispersal site is within mask 
                    sump += kernelp[kn * (j-1) + q];
                }
            }
        }
        else {
            sump = 1.0;
        }
        if (sump>0) {
            // over movement kernel 
            for (q=0; q < kn; q++) {           
                mq = mqarray(m,q);  
                // post-dispersal site is within mask 
                if (mq >= 0) {                 
                    // probability of this move 
                    workpjm[mq] += pjm[m] * kernelp[kn * (j-1) + q] / sump;   
                }
            }
        }
    }
    for (m = 0; m < mm; m++) {
        pjm[m] = workpjm[m];
    }
}
//--------------------------------------------------------------------------

void fillkernelp (int kn, 
                  int jj, 
                  int kerneltype, 
                  double cellsize,
                  const RMatrix<int> kernel, 
                  const RVector<int> moveargsi, 
                  //const CharacterVector fnname,
                  const String fnname,
                  const std::vector<double> &moveargs, 
                  std::vector<double> &kernelp) {
    int j,k;
    double r,r2,a,a2,b;
    NumericVector p;
    std::vector<double> sumj(jj);
    for (j = 0; j < (jj-1); j++) sumj[j] = 0;
    for (k = 0; k < kn; k++) {
        r2 = (kernel[k]*kernel[k] + kernel[k+kn]*kernel[k+kn]) * cellsize * cellsize;
        r = sqrt(r2);
        for (j = 0; j < (jj-1); j++) {
            if (kerneltype == 0) {         // Gaussian kernel 
                a2 = moveargs[j] * moveargs[j];
                kernelp[j * kn + k] = exp(-r2 / 2 / a2);
	    }
            else if (kerneltype == 1) {   // Negative exponential kernel 
                a = moveargs[j];
                kernelp[j * kn + k] = exp(-r / a);
	    }
            else if (kerneltype == 3) {   // 2-D t kernel 
                a2 = moveargs[j] * moveargs[j];
                b = moveargs[j+jj] + 1;
                kernelp[j * kn + k] = (b-1) / M_PI / a2 / pow(1 + r*r/a2, b);
            }
            else if (kerneltype == 2) {   // User kernel 
                // call R function from C++
                Environment env = Environment::global_env();
                Function f = env[fnname];
                if (moveargsi[1]>0)
                    p = f(r, moveargs[j], moveargs[j+jj]);
                else if (moveargs[0]>0)
                    p = f(r, moveargs[j]);
                else 
                    p = f(r);
                kernelp[j * kn + k] = p[0];
                // Rprintf(" k %4d j %4d  kernelp[j * kn + k] %8.6f\n",  k,j,kernelp[j * kn + k]); 
            }
            else if (kerneltype == 4) {  // uniform kernel 
                kernelp[j * kn + k] = 1.0 / kn;
            }
            else stop("unrecognised kerneltype");
            sumj[j] += kernelp[j * kn + k];
        }
    }
    // normalise 
    for (k = 0; k < kn; k++) {
        for (j = 0; j < (jj-1); j++) {
            kernelp[j * kn + k] = kernelp[j * kn + k] / sumj[j];
        }
    }
}
//--------------------------------------------------------------------------

// version with no option for user function (R calls prohibited in RcppParallel)
void fillkernelparallel (int kn, 
                  int jj, 
                  int kerneltype, 
                  double cellsize,
                  const RMatrix<int> kernel, 
                  const RVector<int> moveargsi, 
                  const std::vector<double> &moveargs, 
                  std::vector<double> &kernelp) {
    int j,k;
    double r,r2,a,a2,b;
    std::vector<double> p(jj);
    std::vector<double> sumj(jj);
    for (j = 0; j < (jj-1); j++) sumj[j] = 0;
    for (k = 0; k < kn; k++) {
        r2 = (kernel[k]*kernel[k] + kernel[k+kn]*kernel[k+kn]) * cellsize * cellsize;
        r = sqrt(r2);
        for (j = 0; j < (jj-1); j++) {
            if (kerneltype == 0) {        // Gaussian kernel 
                a2 = moveargs[j] * moveargs[j];
                kernelp[j * kn + k] = exp(-r2 / 2 / a2);
	    }
            else if (kerneltype == 1) {   // Negative exponential kernel 
		a = moveargs[j];
                kernelp[j * kn + k] = exp(-r / a);
	    }
            else if (kerneltype == 2) {   // User kernel 
                // cannot call R function from RcppParallel worker
                stop("cannot call R function from RcppParallel worker; try ncores = 1");
            }
            else if (kerneltype == 3) {   // 2-D t kernel 
                a2 = moveargs[j] * moveargs[j];
                b = moveargs[j+jj] + 1;
                kernelp[j * kn + k] = (b-1) / M_PI / a2 / pow(1 + r2/a2, b);
            }
            else if (kerneltype == 4) {   // uniform kernel 
                kernelp[j * kn + k] = 1.0 / kn;
            }
            else stop("unrecognised kerneltype");
            sumj[j] += kernelp[j * kn + k];
        }
    }
    // normalise 
    for (k = 0; k < kn; k++) {
        for (j = 0; j < (jj-1); j++) {
            kernelp[j * kn + k] = kernelp[j * kn + k] / sumj[j];
        }
    }
}
//--------------------------------------------------------------------------



void getp (int n, int x, int nc, int ss, 
           const RMatrix<double> openval,  
           const RVector<int> PIA, 
           std::vector<double> &p) {
    // column 1 
    int s;
    for (s = 0; s < ss; s++) {
        p[s] = openval(PIA[i3(n, s, x, nc, ss )]-1, 0); 
    }
}
//--------------------------------------------------------------------------

void getphij (int n, int x, int nc, int jj, 
              const RMatrix<double> openval,  
              const RVector<int> PIAJ,
              const RVector<double> intervals, 
              std::vector<double> &phij) {
    // column 2 
    int j;
    double phi;
    for (j = 0; j < (jj-1); j++) {
        // jj-1 because one fewer intervals than primary sessions  
        phi = openval(PIAJ[i3(n, j, x, nc, jj)]-1, 1);       
        // adjust for interval duration  
        phij[j] = exp(log(phi) * intervals[j]);  
    }
    phij[jj-1] = 0;
}
//--------------------------------------------------------------------------

void getmoveargs (int n, int x, int nc, int jj, 
		  const RMatrix<double> openval,  
		  const RVector<int> PIAJ,
		  const RVector<int> moveargsi,
		  std::vector<double> &moveargs) {
    // column moveargsi (and maybe moveargsi + 1) 
    int j;
    for (j = 0; j < (jj-1); j++) {
        // jj-1 because one fewer intervals than primary sessions  
        moveargs[j] = openval(PIAJ[i3(n, j, x, nc, jj)]-1, moveargsi[0]);  
        if (moveargsi[1]>0)
            moveargs[j+jj] = openval(PIAJ[i3(n, j, x, nc, jj)]-1, moveargsi[1]);  
    }
    moveargs[jj-1] = 0;
    moveargs[2*jj-1] = 0;
}
//--------------------------------------------------------------------------

void getpj (int n, int x, int nc, int jj, 
            const RMatrix<double> openval,  
            const RVector<int> PIAJ,
            std::vector<double> &pj) {
    // column 2 
    int j;
    for (j = 0; j < jj; j++) {
        pj[j] = openval(PIAJ[i3(n, j, x, nc, jj)]-1, 0);       
    }
}
//--------------------------------------------------------------------------

void getg (int type, int n, int x, int nc, int jj, 
           const RMatrix<double> openval,  
           const RVector<int> PIAJ,
           std::vector<double> &g) {
    // column 4 
    int j;
    for (j = 0; j < jj; j++) {
        if (type != 27)
            g[j] = 0;
        else
            g[j] = openval(PIAJ[i3(n,j, x, nc, jj)]-1, 3);       
    }
}
//--------------------------------------------------------------------------

void getfj (int n, int x, int nc, int jj, 
            const RMatrix<double> openval,  
            const RVector<int> PIAJ,
            const RVector<double> intervals, 
            std::vector<double> &phij,
            std::vector<double> &fj) {
    // column 3 
    int j;
    double f,phi;
    for (j = 0; j < (jj-1); j++) {
        // jj-1 because one fewer intervals than primary sessions  
        f = openval(PIAJ[i3(n, j, x, nc, jj )]-1, 2); 
	phi = exp(log(phij[j]) / intervals[j]);
        // adjust for interval duration  
        fj[j] = exp(log(phi+f) * intervals[j]) - phij[j];  
    }
    fj[jj-1] = 0;
}
//--------------------------------------------------------------------------

void getlj (int n, int x, int nc, int jj, 
            const RMatrix<double> openval,  
            const RVector<int> PIAJ,
            const RVector<double> intervals, 
            std::vector<double> &lj) {
    // column 3 
    int j;
    double l;
    for (j = 0; j < (jj-1); j++) {
        // jj-1 because one fewer intervals than primary sessions  
        l = openval(PIAJ[i3(n, j, x, nc, jj )]-1, 2); 
        // adjust for interval duration  
        lj[j] = exp(log(l) * intervals[j]);  
    }
    lj[jj-1] = 0;
}
//--------------------------------------------------------------------------

void getgaml (int n, int x, int nc, int jj, 
              const RMatrix<double> openval,  
              const RVector<int> PIAJ,
              const RVector<double> intervals, 
              std::vector<double> &gam) {
    // column 3 
    int j;
    double phij;
    double lamj;
    for (j = 0; j < (jj-1); j++) {
        phij = openval(PIAJ[i3(n, j, x, nc, jj)]-1, 1); 
        phij = exp(log(phij) * intervals[j]);  
        lamj = openval(PIAJ[i3(n, j, x, nc, jj)]-1, 2); 
        lamj = exp(log(lamj) * intervals[j]);  
        gam[j+1] = phij/lamj;
    }
    gam[0] = 0;
}
//--------------------------------------------------------------------------

void getgamj (int n, int x, int nc, int jj, 
              const RMatrix<double> openval,  
              const RVector<int> PIAJ,
              const RVector<double> intervals, 
              std::vector<double> &gamj) {
    // column 3 
    int j;
    double gam;
    for (j = 1; j < jj; j++) {
        gam = openval(PIAJ[i3(n, j, x, nc, jj)]-1, 2); 
        gamj[j] = exp(log(gam) * intervals[j-1]);  
    }
    gamj[0] = 0;
}
//--------------------------------------------------------------------------

void getkapj (int n, int x, int nc, int jj, 
              const RMatrix<double> openval,  
              const RVector<int> PIAJ,
              std::vector<double> &kapj) {
    // column 3 
    int j;
    for (j = 1; j < jj; j++) {
        kapj[j] = openval(PIAJ[i3(n, j, x, nc, jj)]-1, 2);
    }
    kapj[0] = 1;
}
//--------------------------------------------------------------------------

// getgam <- function (n, x, openval, PIAJ, intervals) {
//     J2 <- 2:(length(intervals)+1)
//     c(0,openval[PIAJ[n, J2, x],3])
// }

void getbeta0 (int n, int x, int nc, int jj, 
               const RMatrix<double> openval,  
               const RVector<int> PIAJ,
               std::vector<double> &beta) {
    // column 3 
    int j;
    double sumbeta = 0;
    for (j = 1; j < jj; j++) {
        beta[j] = openval(PIAJ[i3(n, j, x, nc, jj )]-1,2); 
        sumbeta += exp(beta[j]);
    }
    beta[0] = 1;
    for (j = 1; j < jj; j++) {
        beta[j] = exp(beta[j]) / (1 + sumbeta);
        beta[0] -= beta[j];
    }
}
//--------------------------------------------------------------------------

void gettau (int n, int x, int nc, int jj,
             const RMatrix<double> openval,  
             const RVector<int> PIAJ,
             std::vector<double> &tau,
             int M) {
    // column 5 
    int j;
    double sumtau = 0;
    for (j = 0; j < M; j++) {
        tau[j] = openval(PIAJ[i3(n, j, x, nc, jj )]-1, 4);
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
//--------------------------------------------------------------------------

void getDj (int n, int x, int nc, int jj, 
            const RMatrix<double> openval,  
            const RVector<int> PIAJ,
            std::vector<double> &Dj) {
    // column 3 
    int j;
    for (j = 0; j < jj; j++) {
        Dj[j] = openval(PIAJ[i3(n, j, x, nc, jj )]-1, 2); 
    }
}
//--------------------------------------------------------------------------

// per capita recruitment cf Link & Barker 2005, Schwarz 'Gentle Intro'
void getbetaf (int n, int x, int nc, int jj, 
               const RMatrix<double> openval,  
               const RVector<int> PIAJ,
               std::vector<double> &phij,
               const RVector<double> intervals, 
               std::vector<double> &beta) {
    int j;
    double sumbeta = 1;
    std::vector<double> d(jj);
    std::vector<double> fj(jj);
    getfj (n, x, nc, jj, openval, PIAJ, intervals, phij, fj);
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
//--------------------------------------------------------------------------

void getbetal (int n, int x, int nc, int jj, 
               const RMatrix<double> openval,  
               const RVector<int> PIAJ,
               std::vector<double> &phij,
               const RVector<double> intervals, 
               std::vector<double> &beta) {
    int j;
    double sumbeta = 1;
    std::vector<double> d(jj);
    std::vector<double> fj(jj);
    std::vector<double> lambdaj(jj);

    getlj (n, x, nc, jj, openval, PIAJ, intervals, lambdaj);
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
//--------------------------------------------------------------------------

void getbetag (int n, int x, int nc, int jj, 
               const RMatrix<double> openval,  
               const RVector<int> PIAJ, 
               std::vector<double> &phij,
               const RVector<double> intervals, 
               std::vector<double> &beta) {
    int j;
    double sumbeta = 1;
    std::vector<double> d(jj);
    std::vector<double> fj(jj);
    std::vector<double> gamj(jj);
    getgamj (n, x, nc, jj, openval, PIAJ, intervals, gamj);

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
//--------------------------------------------------------------------------

void getbetak (int n, int x, int nc, int jj, 
               const RMatrix<double> openval,  
               const RVector<int> PIAJ, 
               std::vector<double> &phij,
               std::vector<double> &beta) {
    int i,j;
    std::vector<double> tau(jj);
    std::vector<double> pj(jj);
    std::vector<double> fprod(jj);
    std::vector<double> fj(jj);
    std::vector<double> kapj(jj);

    getkapj (n, x, nc, jj, openval, PIAJ, kapj);
    getpj (n, x, nc, jj, openval, PIAJ, pj);

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
}
//--------------------------------------------------------------------------

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
               const RMatrix<double> openval,  
               const RVector<int> PIAJ, 
               std::vector<double> &beta) {
    int j;
    double sumB = 0;
    std::vector<double> B(jj);
    getDj (n, x, nc, jj, openval, PIAJ, B);
    for (j = 0; j < jj; j++) {
        sumB += B[j];
    }
    for (j = 0; j < jj; j++) {
        beta[j] = B[j] / sumB;
    }
}
//--------------------------------------------------------------------------

void getbetaD (int n, int x, int nc, int jj, 
               const RMatrix<double> openval,  
               const RVector<int> PIAJ, 
               std::vector<double> &phij,
               std::vector<double> &beta) {
        int j;
    double sumB;
    std::vector<double> B(jj);
    std::vector<double> D(jj);
    getDj (n, x, nc, jj, openval, PIAJ, D);
    
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
//--------------------------------------------------------------------------

void getbeta (int type, int n, int x, int nc, int jj, 
	      const RMatrix<double> openval,  
	      const RVector<int> PIAJ, 
	      const RVector<double> intervals,
	      std::vector<double> &phij,
	      std::vector<double> &beta) {
    if ((type == 2) || (type == 17) || (type == 11) || (type == 13) || 
	(type == 41) || (type == 43) || (type == 30) || (type == 31)) 
        getbeta0 (n, x, nc, jj, openval, PIAJ, beta);
    else if ((type == 4) || (type == 15) ||  (type == 27) ||  (type == 7) || 
             (type == 9) || (type == 37) || (type == 39))
        getbetaf (n, x, nc, jj, openval, PIAJ, phij, intervals, beta);
    else if ((type == 3) || (type == 16) || (type == 10) || (type == 12) || 
             (type == 20) || (type == 40) || (type == 42))
        getbetal (n, x, nc, jj, openval, PIAJ, phij, intervals, beta);
    else if ((type == 14) || (type == 18))
        getbetaB (n, x, nc, jj, openval, PIAJ, beta);
    else if ((type == 8) || (type == 19) || (type == 38))
        getbetaD (n, x, nc, jj, openval, PIAJ, phij, beta);
    else if ((type == 22) || (type == 23) || (type == 24) || (type == 25) ||
             (type == 26))
        getbetag (n, x, nc, jj, openval, PIAJ, phij, intervals, beta);
    else if ((type == 28) || (type == 29))
        getbetak (n, x, nc, jj, openval, PIAJ, phij, beta);
    else stop("no beta for this type");
}
//--------------------------------------------------------------------------
