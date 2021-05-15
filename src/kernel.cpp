// movement model functions
// 2021-02-23 moved from utils.cpp
// 2021-02-24 does not protect against extreme R in annularR
// 2021-05-11 replaced sqrt with std::sqrt to remove ambiguity
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

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

// for variable annulus 

struct ipoint
{
    double x;
    double y;
    double theta;
};

struct arc
{
    double x1;
    double x2;
    double y1;
    double y2;
    double dtheta;
    int    cell;
};

int findcell (
        double x,
        double y, 
        const RMatrix<int> kernel ) 
{
    int mind2 = 100000;
    int kn = kernel.nrow();
    int cell = 0;
    double d2;
    for (int k = 0; k<kn; k++) {
        d2 = (x-kernel(k,0)) * (x-kernel(k,0)) + 
            (y-kernel(k,1)) * (y-kernel(k,1));
        if (d2 < mind2) {
            mind2 = d2;    
            cell = k;
        }
    }
    return cell;
}

// in situ update of pts
void addpoints (
        std::vector<ipoint> &pts, 
        double dx, 
        double dy, 
        double dr2, 
        double D, 
        double R) 
{
    double inc = R * R * dr2 - D * D; 
    int sgndy = 1; if (dy<0) sgndy = -1;
    // from https://mathworld.wolfram.com/Circle-LineIntersection.html
    // assume origin at 0,0
    // if intersects, place coordinates of two intersection points in pts
    if (inc > 0 ) {
        ipoint vertex;
        vertex.x = (D*dy + sgndy * dx * std::sqrt(inc)) / dr2;
        vertex.y = (D*dx + fabs(dy) * std::sqrt(inc)) / dr2;
        pts.push_back(vertex);
        vertex.x = (D*dy - sgndy * dx * std::sqrt(inc)) / dr2;
        vertex.y = (D*dx - fabs(dy) * std::sqrt(inc)) / dr2;
        pts.push_back(vertex);
    }
}

// function used by std::sort to define sort order (here ascending theta)
// cf https://www.cplusplus.com/articles/NhA0RXSz/
bool sortByTheta(const ipoint &lhs, const ipoint &rhs) { return lhs.theta < rhs.theta; }

std::vector<double> annulus (
        double p0, 
        double R, 
        const RMatrix<int> kernel) 
{
    int kn = kernel.nrow();
    std::vector<double> result(kn);  // cell-specific probabilities
    std::vector<ipoint> pts;         // points at which circle intersects cell edges
    std::vector<arc> arcs;           // arcs
    
    int i, k, ni;
    int K2 = 0;
    arc newarc;
    double dx;
    double dy;
    double dr2;
    double D;
    
    for (k=0; k<kn; k++) K2 = std::max(kernel[k], K2);  
    
    // place any intersections in pts
    // ------------------------------
    for (i = 0; i<=K2; i++) {
        dx = 0;                // x2-x1
        dy = 2 * K2 + 1;       // y2-y1
        dr2 = dy * dy;         // dx^2 + dy^2
        D = (i + 0.5) * dy;    // x1*y2 - x2 * y1
        addpoints (pts, dx, dy, dr2, D, R*K2);
        addpoints (pts, dx, dy, dr2, -D, R*K2);
        dx = 2 * K2 + 1;       // x2-x1
        dy = 0;                // y2-y1
        dr2 = dx * dx;         // dx^2 + dy^2
        D = (i+0.5) * dx;      // x1*y2 - x2 * y1
        addpoints (pts, dx, dy, dr2, D, R*K2);
        addpoints (pts, dx, dy, dr2, -D, R*K2);
    }
    ni = pts.size();
    
    // compute ni theta = atan2(y,x)
    for (i=0; i<ni; i++) {
        pts[i].theta = atan2(pts[i].y, pts[i].x) + M_PI;  // range 0..2pi
        // debug
        // Rprintf(" i %4d pts[i].x %8.6g pts[i].y %8.6g pts[i].theta %8.6g \n",
        //    i, pts[i].x, pts[i].y, pts[i].theta); 
    }
    
    // sort intersections by theta 
    sort(pts.begin(), pts.end(), sortByTheta);
    
    // load arcs (ni-1) dtheta, x1, x2, y1, y2
    for (i=0; i<(ni-1); i++) {
        newarc.dtheta = pts[i+1].theta - pts[i].theta;
        newarc.x1 = pts[i].x;
        newarc.y1 = pts[i].y;
        newarc.x2 = pts[i+1].x;
        newarc.y2 = pts[i+1].y;
        arcs.push_back(newarc);
    }
    
    // close last arc dtheta ni (+2*M_PI)
    newarc.dtheta = pts[0].theta - pts[ni-1].theta + 2 * M_PI;
    newarc.x1 = pts[ni-1].x;
    newarc.y1 = pts[ni-1].y;
    newarc.x2 = pts[0].x;
    newarc.y2 = pts[0].y;
    arcs.push_back(newarc);
    
    // find cell for each arc and fill results
    for (i=0; i<ni; i++) {
        double cx = arcs[i].x1 + (arcs[i].x2 - arcs[i].x1)/2;
        double cy = arcs[i].y1 + (arcs[i].y2 - arcs[i].y1)/2;
        arcs[i].cell = findcell(cx,cy,kernel);
        // peripheral arcs sum to 2pi
        result[arcs[i].cell] = (1 - p0) * arcs[i].dtheta / (2 * M_PI);
    }
    result[kn/2] = p0;   // centre assumes symmetrical kernel
    return(result);
}

//--------------------------------------------------------------------------

void fillkernelp (int kn, 
    int jj, 
    int kerneltype, 
    bool sparsekernel,
    double cellsize,
    const RMatrix<int> kernel, 
    const RVector<int> moveargsi, 
    const String fnname,
    const std::vector<double> &moveargs, 
    std::vector<double> &kernelp) 
{
    int j,k;
    int n1 = 0;
    int n2 = 0;
    int K2 = 0;
    double r,r2,a,a2,b;
    double p0 = 1;
    double p1 = 0;
    double R = 0;
    int x,y;
    double diag;
    NumericVector p;
    std::vector<double> sumj(jj);
    std::vector<double> tempkernel(kn);
    for (k=0; k<kn; k++) K2 = std::max(kernel[k], K2);  
    for (j = 0; j < (jj-1); j++) {
        sumj[j] = 0;
        if (kerneltype == 5) {    // annular
            p0 = moveargs[j];
        }
        if (kerneltype == 6) {    // annular2
            p0 = moveargs[j];
            p1 = moveargs[j+jj];
            n1 = 0; n2 = 0;
            for (k = 0; k < kn; k++) {
                x = kernel[k];
                y = kernel[k+kn];
                r = std::sqrt(x*x + y*y);
                if (r > 1e-8) {
                    if (r < (K2-1)) n1++;
                    else n2++;
                }
            }
            Rprintf(" n1 %4d n2 %4d\n",  n1, n2); 
        }
        if (kerneltype == 7) {    // annularR
            p0 = moveargs[j];
            R = moveargs[j+jj];
            // compute arc lengths etc.
            tempkernel = annulus(p0, R, kernel);
            if (j==0) for (k = 0; k < kn; k++) {
                Rprintf(" k %4d tempkernel[k] %8.6f\n",  k, tempkernel[k]); 
            }
        }
        for (k = 0; k < kn; k++) {
            x = kernel[k];
            y = kernel[k+kn];
            if (abs(x) == abs(y) && x!=0) diag = std::sqrt(2); else diag = 1;
            r2 = (x*x + y*y) * cellsize * cellsize;
            r = std::sqrt(r2);
            if (kerneltype == 0) {        // BVN Gaussian kernel 
                a2 = moveargs[j] * moveargs[j];
                kernelp[j * kn + k] = exp(-r2 / 2 / a2);
            }
            else if (kerneltype == 1) {   // BVE Negative exponential kernel 
                a = moveargs[j];
                kernelp[j * kn + k] = exp(-r / a);
            }
            else if (kerneltype == 3) {   // BVT 2-D t kernel 
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
            }
            else if (kerneltype == 4) {  // uniform kernel 
                kernelp[j * kn + k] = 1.0 / kn;
            }
            else if (kerneltype == 5) {  // annular kernel 
                if (r<1e-8)
                    kernelp[j * kn + k] = p0;
                else {
                    kernelp[j * kn + k] = (1-p0)/(kn-1);
                }
            }
            else if (kerneltype == 6) {  // annular2 kernel 
                if (r<1e-8)
                    kernelp[j * kn + k] = p0;
                else if (r < (K2-1)) 
                    kernelp[j * kn + k] = p1/n1;
                else 
                    kernelp[j * kn + k] = (1-p0-p1)/n2;
            }
            else if (kerneltype == 7) {  // annularR kernel 
                // variable annulus precomputed
                kernelp[j * kn + k] = tempkernel[k];
            }
            else stop("unrecognised kerneltype");
            if (sparsekernel) {
                kernelp[j * kn + k] = r * diag * kernelp[j * kn + k];
            }
            sumj[j] += kernelp[j * kn + k];
            // if (j==0) Rprintf(" k %4d j %4d  kernelp[j * kn + k] %8.6f\n",  k,j,kernelp[j * kn + k]); 
        }
    }
    // normalise 
    if (kerneltype != 7) {
        for (k = 0; k < kn; k++) {
            for (j = 0; j < (jj-1); j++) {
                kernelp[j * kn + k] = kernelp[j * kn + k] / sumj[j];
            }
        }
    }
}
//--------------------------------------------------------------------------

// version with no option for user function (R calls prohibited in RcppParallel)
void fillkernelparallel (int kn, 
    int jj, 
    int kerneltype, 
    bool sparsekernel,
    double cellsize,
    const RMatrix<int> kernel, 
    const RVector<int> moveargsi, 
    const std::vector<double> &moveargs, 
    std::vector<double> &kernelp) {
    int j,k;
    int n1 = 0;
    int n2 = 0;
    int K2 = 0;
    double r,r2,a,a2,b;
    double p0 = 1;
    double p1 = 0;
    double R = 0;
    int x,y;
    double diag;
    std::vector<double> p(jj);
    std::vector<double> sumj(jj);
    std::vector<double> tempkernel(kn);
    //
    for (k=0; k<kn; k++) K2 = std::max(kernel[k], K2);  
    if (kerneltype == 6) {    // annular2
        n1 = 0; n2 = 0;
        for (k = 0; k < kn; k++) {
            x = kernel[k];
            y = kernel[k+kn];
            r = std::sqrt(x*x + y*y);
            if (r > 1e-8) {
                if (r < (K2-1)) n1++;
                else n2++;
            }
        }
    }
    for (j = 0; j < (jj-1); j++) {
        sumj[j] = 0;
        if (kerneltype == 5) {    // annular
            p0 = moveargs[j];
        }
        if (kerneltype == 6) {    // annular2
            p0 = moveargs[j];
            p1 = moveargs[j+jj];
        }
        if (kerneltype == 7) {
            p0 = moveargs[j];
            R = moveargs[j+jj];
            // compute arc lengths etc.
            tempkernel = annulus(p0, R, kernel);
        }
        
        for (k = 0; k < kn; k++) {
            x = kernel[k];
            y = kernel[k+kn];
            if (abs(x) == abs(y) && x!=0) diag = std::sqrt(2); else diag = 1;
            r2 = (x*x + y*y) * cellsize * cellsize;
            r = std::sqrt(r2);
            if (kerneltype == 0) {        // BVN Gaussian kernel 
                a2 = moveargs[j] * moveargs[j];
                kernelp[j * kn + k] = exp(-r2 / 2 / a2);
            }
            else if (kerneltype == 1) {   // BVE Negative exponential kernel 
                a = moveargs[j];
                kernelp[j * kn + k] = exp(-r / a);
            }
            else if (kerneltype == 3) {   // BVT 2-D t kernel 
                a2 = moveargs[j] * moveargs[j];
                b = moveargs[j+jj] + 1;
                kernelp[j * kn + k] = (b-1) / M_PI / a2 / pow(1 + r2/a2, b);
            }
            else if (kerneltype == 2) {   // User kernel 
                // cannot call R function from RcppParallel worker
                stop("cannot call R function from RcppParallel worker; try ncores = 1");
            }
            else if (kerneltype == 4) {   // uniform kernel 
                kernelp[j * kn + k] = 1.0 / kn;
            }
            else if (kerneltype == 5) {  // annular kernel 
                if (r<1e-8)
                    kernelp[j * kn + k] = p0;
                else 
                    kernelp[j * kn + k] = (1-p0)/(kn-1);
            }
            else if (kerneltype == 6) {  // annular2 kernel 
                if (r<1e-8)
                    kernelp[j * kn + k] = p0;
                else if (r < (K2-1)) 
                    kernelp[j * kn + k] = p1/n1;
                else 
                    kernelp[j * kn + k] = (1-p0-p1)/n2;
            }
            else if (kerneltype == 7) {  // annularR kernel 
                // variable annulus precomputed
                kernelp[j * kn + k] = tempkernel[k];
            }
            else stop("unrecognised kerneltype");
            if (sparsekernel) {
                kernelp[j * kn + k] = r * diag * kernelp[j * kn + k];
            }
            sumj[j] += kernelp[j * kn + k];
        }
    }
    // normalise 
    if (kerneltype != 7) {
        for (k = 0; k < kn; k++) {
            for (j = 0; j < (jj-1); j++) {
                kernelp[j * kn + k] = kernelp[j * kn + k] / sumj[j];
            }
        }
    }
}
//--------------------------------------------------------------------------

