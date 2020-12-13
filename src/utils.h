#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

#define fuzz 1e-200
#define huge 1e10

int i3 (int i, int j, int k, int ii, int jj);
int i4 (int i, int j, int k, int l, int ii, int jj, int kk);

//--------------------------------------------------------------------------

double gpois (int count, double lambda, int uselog);
double gbinom(int count, int size, double p, int uselog);
double gbinomFP (int count, double size, double p, int uselog);
double pski ( int binomN, int count, double Tski, double g);
    
//--------------------------------------------------------------------------

// distance between two points given by row k in traps and row m in mask
double dkm (int k, int m, RMatrix<double> traps, RMatrix<double> mask);
//--------------------------------------------------------------------------

void convolvemq (
    int    mm,        // number of points on mask 
    int    kn,        // number of points on kernel 
    int    j,         // session number 1..jj 
    int    edgecode,              // adjust for incomplete kernel
    const  RMatrix<int> &mqarray, // input [& 2020-10-31]
    std::vector<double> &kernelp, // p(move|dx,dy) for points in kernel 
    std::vector<double> &pjm      // return value 
    );
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
                  std::vector<double> &kernelp);

//--------------------------------------------------------------------------

void fillkernelparallel (int kn, 
                         int jj, 
                         int kerneltype, 
                         double cellsize,
                         const RMatrix<int> kernel, 
                         const RVector<int> moveargsi, 
                         const std::vector<double> &moveargs, 
                         std::vector<double> &kernelp);
//--------------------------------------------------------------------------

double hfn
    (int k, int m, int c, 
     const RMatrix<double> openval,  
     const RMatrix<double> traps,
     const RMatrix<double> mask, 
     int sigmai, int detectfn);
//--------------------------------------------------------------------------

void getp (int n, int x, int nc, int ss, 
           const RMatrix<double> openval,  
           const RVector<int> PIA, 
           std::vector<double> &p);
//--------------------------------------------------------------------------

void getphij (int n, int x, int nc, int jj, 
              const RMatrix<double> openval,  
              const RVector<int> PIAJ,
              const RVector<double> intervals, 
              std::vector<double> &phij);
//--------------------------------------------------------------------------

void getmoveargs (int n, int x, int nc, int jj, 
                  const RMatrix<double> openval,  
                  const RVector<int> PIAJ,
                  const RVector<int> moveargsi,
                  std::vector<double> &moveargs);
//--------------------------------------------------------------------------

void getpj (int n, int x, int nc, int jj, 
            const RMatrix<double> openval,  
            const RVector<int> PIAJ,
            std::vector<double> &pj);
//--------------------------------------------------------------------------

void getg (int type, int n, int x, int nc, int jj, 
           const RMatrix<double> openval,  
           const RVector<int> PIAJ,
           std::vector<double> &g);
//--------------------------------------------------------------------------

void getfj (int n, int x, int nc, int jj, 
            const RMatrix<double> openval,  
            const RVector<int> PIAJ,
            const RVector<double> intervals, 
            std::vector<double> &phij,
            std::vector<double> &fj);
//--------------------------------------------------------------------------

void getlj (int n, int x, int nc, int jj, 
            const RMatrix<double> openval,  
            const RVector<int> PIAJ,
            const RVector<double> intervals, 
            std::vector<double> &lj);
//--------------------------------------------------------------------------

void getgaml (int n, int x, int nc, int jj, 
              const RMatrix<double> openval,  
              const RVector<int> PIAJ,
              const RVector<double> intervals, 
              std::vector<double> &gam);
//--------------------------------------------------------------------------

void getgamj (int n, int x, int nc, int jj, 
              const RMatrix<double> openval,  
              const RVector<int> PIAJ,
              const RVector<double> intervals, 
              std::vector<double> &gamj);
//--------------------------------------------------------------------------

void getkapj (int n, int x, int nc, int jj, 
              const RMatrix<double> openval,  
              const RVector<int> PIAJ,
              std::vector<double> &kapj);
//--------------------------------------------------------------------------

// getgam <- function (n, x, openval, PIAJ, intervals) {
//     J2 <- 2:(length(intervals)+1)
//     c(0,openval[PIAJ[n, J2, x],3])
// }

void getbeta0 (int n, int x, int nc, int jj, 
               const RMatrix<double> openval,  
               const RVector<int> PIAJ,
               std::vector<double> &beta);
//--------------------------------------------------------------------------

void gettau (int n, int x, int nc, int jj,
             const RMatrix<double> openval,  
             const RVector<int> PIAJ,
             std::vector<double> &tau,
             int M);
//--------------------------------------------------------------------------

void getDj (int n, int x, int nc, int jj, 
            const RMatrix<double> openval,  
            const RVector<int> PIAJ,
            std::vector<double> &Dj);
//--------------------------------------------------------------------------

// per capita recruitment cf Link & Barker 2005, Schwarz 'Gentle Intro'
void getbetaf (int n, int x, int nc, int jj, 
               const RMatrix<double> openval,  
               const RVector<int> PIAJ,
               std::vector<double> &phij,
               const RVector<double> intervals, 
               std::vector<double> &beta);
//--------------------------------------------------------------------------

void getbetal (int n, int x, int nc, int jj, 
               const RMatrix<double> openval,  
               const RVector<int> PIAJ,
               std::vector<double> &phij, 
               const RVector<double> intervals, 
               std::vector<double> &beta);
//--------------------------------------------------------------------------

void getbetag (int n, int x, int nc, int jj, 
               const RMatrix<double> openval,  
               const RVector<int> PIAJ, 
               std::vector<double> &phij, 
               const RVector<double> intervals, 
               std::vector<double> &beta);
//--------------------------------------------------------------------------

void getbetak (int n, int x, int nc, int jj, 
               const RMatrix<double> openval,  
               const RVector<int> PIAJ, 
               std::vector<double> &phij,
               std::vector<double> &beta);
//--------------------------------------------------------------------------
                   
// return parameterisation cf Pledger et al. 2010 p 885 
void getbetaB (int n, int x, int nc, int jj, 
               const RMatrix<double> openval,  
               const RVector<int> PIAJ, 
               std::vector<double> &beta);
//--------------------------------------------------------------------------

void getbetaD (int n, int x, int nc, int jj, 
               const RMatrix<double> openval,  
               const RVector<int> PIAJ, 
               std::vector<double> &phij,
               std::vector<double> &beta);
//--------------------------------------------------------------------------

void getbeta (int type, int n, int x, int nc, int jj, 
              const RMatrix<double> openval,  
              const RVector<int> PIAJ, 
              const RVector<double> intervals,
              std::vector<double> &phij,
              std::vector<double> &beta);
//--------------------------------------------------------------------------

List makelookupcpp (const NumericMatrix x);  
//--------------------------------------------------------------------------

void pr0njmx (int n, int x, 
              const RVector<int> cumss, 
              int nc,  int jj, int kk, int mm, int cc0, int binomN,  
              const RVector<int> PIA0, 
              const RVector<double> gk0, 
              const RMatrix<double> Tsk, 
              std::vector<double> &pjm);
//--------------------------------------------------------------------------
