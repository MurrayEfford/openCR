// movement model functions
// 2021-02-23 moved from utils.cpp
// 2021-02-24 does not protect against extreme R in annularR
// 2021-05-11 replaced sqrt with std::sqrt to remove ambiguity
// 2021-07-18 frE, frG, frL (8,9,10)
// 2021-07-25 Boost library for stat distributions
// 2021-07-27 consolidate fillkernelp and fillkernelparallel in single function
// 2021-07-27 export fillkernelcpp for R function 'make.kernel'
// 2021-07-27 grain = -1 for debug messages
// 2021-07-29 BVNz, BVEz, uniformz, frEz  (11,12,13,14)
// 2021-08-07 String input for usermodel converted to std::string for thread safety
// 2021-09-25 BVN2 (16)
// 2021-10-06 protect Boost call from infinite scale lognormal
// 2021-10-12 strict new codes RDE, RDG, RDL etc. 
// 2021-11-29 settlement option

#include "utils.h"

// [[Rcpp::depends(BH)]]

void convolvemq (
        const int    mm,        // number of points on mask 
        const int    kn,        // number of points on kernel
        const int    j,         // session number 1..jj 
        const int    edgecode,  // 0 none, no action
                                // 1 wrapped, no action
                                // 2 normalize truncated kernel
                                // 3 normalize, weighted by settlement
        const RcppParallel::RMatrix<int> &mqarray, 
        const RcppParallel::RMatrix<double> &settlement, 
        std::vector<double> &kernelp, // p(move|dx,dy) for points in kernel 
        std::vector<double> &pjm      // return value
)
{
    int m, q, mq;
    double sump;
    bool settlemodel = settlement.nrow() == mqarray.nrow(); // true iff settlement model
    std::vector<double> workpjm(mm);
    
    // convolve movement kernel and pjm... 
    for (m = 0; m < mm; m++) {
        if (edgecode >= 2) {
            sump = 0;
            for (q=0; q < kn; q++) {           // over movement kernel 
                mq = mqarray(m,q);
                if (mq >= 0) {       // post-dispersal site is within mask 
                    if (settlemodel) {
                        // weight by relative settlement 
                        sump += kernelp[kn * (j-1) + q] * settlement(mq,j-1);  
                    }
                    else {
                        sump += kernelp[kn * (j-1) + q];
                    }
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
                    if (settlemodel) {
                        // weight by relative settlement 
                        workpjm[mq] += pjm[m] * kernelp[kn * (j-1) + q] * settlement(mq,j-1) / sump;   
                    }
                    else {
                        // probability of this move 
                        workpjm[mq] += pjm[m] * kernelp[kn * (j-1) + q] / sump;   
                    }
                }
            }
        }
    }
    for (m = 0; m < mm; m++) {
        pjm[m] = workpjm[m];
    }
}
//--------------------------------------------------------------------------

// from single point m (used with 'anchored')
void convolvemq1 (
        const int    m,                     // initial location on mask
        const int    j,                     // session number 1..jj 
        const int    edgecode,              // 0 none, no action; 1 wrapped, no action; 2 normalize truncated kernel
        const RcppParallel::RMatrix<int> &mqarray, // input [& 2020-10-31]
        const RcppParallel::RMatrix<double> &settlement, 
        const std::vector<double> &kernelp, // p(move|dx,dy) for points in kernel 
        std::vector<int>    &mj,
        std::vector<double> &pj      // return value
)
{
    int kn = mqarray.ncol();    // number of points on kernel
    int q;
    double sump;
    double settle = 1.0;
    
    bool settlemodel = settlement.nrow() == mqarray.nrow(); // true iff settlement model
    if (edgecode >= 2) {
        // adjust for edge-truncated kernel
        sump = 0;
        for (q=0; q < kn; q++) {           // over movement kernel 
            if (mqarray(m,q) >= 0) {       // post-dispersal site is within mask 
                if (settlemodel) settle = settlement(mqarray(m,q),j-1);
                sump += kernelp[kn * (j-1) + q] * settle;
            }
        }
    }
    else {
        sump = 1.0;
    }
    std::fill(mj.begin(), mj.end(), 0);
    std::fill(pj.begin(), pj.end(), 0);
    if (sump>0) {
        // over movement kernel 
        for (q=0; q < kn; q++) {           
            mj[q] = mqarray(m,q);          // destination; negative if off-mask
            if (mj[q] >= 0) {                 
                if (settlemodel) settle = settlement(mqarray(m,q),j-1);
                pj[q] = kernelp[kn * (j-1) + q] * settle / sump;   // probability of reaching point mj[q]
            }
            else {
                pj[q] = 0;
            }
        }
    }
}
//--------------------------------------------------------------------------

// Interface for use in R
// [[Rcpp::export]]
Rcpp::NumericVector convolvemqcpp (
        const int    j,
        const int    edgecode,
        const Rcpp::NumericMatrix mqarray,
        const Rcpp::NumericMatrix settlement, 
        const Rcpp::NumericVector kernelp,
        const Rcpp::NumericVector pjm) {
    
    int mm = mqarray.nrow();            // number of points on mask
    int kn = mqarray.ncol();            // number of points on kernel
    bool settlemodel = settlement.nrow() == mm; // true iff settlement model
    
    int m, q, mq;
    double sump;
    double settle = 1.0;
    std::vector<double> workpjm(mm);
    
    // convolve movement kernel and pjm... 
    for (m = 0; m < mm; m++) {
        if (edgecode >= 2) {
            // 2020-10-29 adjust for edge-truncated kernel cf convolvemqold
            sump = 0;
            for (q=0; q < kn; q++) {           // over movement kernel 
                mq = mqarray(m,q);  
                if (mq >= 0) {       // post-dispersal site is within mask 
                    if (settlemodel) settle = settlement(mq,j-1);
                    sump += kernelp[kn * (j-1) + q] * settle;
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
                    if (settlemodel) settle = settlement(mq,j-1);
                    workpjm[mq] += pjm[m] * kernelp[kn * (j-1) + q] * settle / sump;   
                }
            }
        }
    }
    return Rcpp::wrap(workpjm);
}
//==============================================================================


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
        const RcppParallel::RMatrix<int> kernel ) 
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
        vertex.y = (D*dx + std::fabs(dy) * std::sqrt(inc)) / dr2;
        pts.push_back(vertex);
        vertex.x = (D*dy - sgndy * dx * std::sqrt(inc)) / dr2;
        vertex.y = (D*dx - std::fabs(dy) * std::sqrt(inc)) / dr2;
        pts.push_back(vertex);
    }
}

// function used by std::sort to define sort order (here ascending theta)
// cf https://www.cplusplus.com/articles/NhA0RXSz/
bool sortByTheta(const ipoint &lhs, const ipoint &rhs) { return lhs.theta < rhs.theta; }

std::vector<double> annulus (
        double p0, 
        double R, 
        const RcppParallel::RMatrix<int> kernel) 
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

// 2021-07-27 argument 'grain' added (and ditch fillkernelparallel for grain>0)
// cannot call R function from RcppParallel worker ncores > 1
// hence fnname ignored when grain > 0 and no calls to warning() stop() or Rprintf()

// kerneltype takes these values

// BVN    = 0 
// BVE    = 1 
// user   = 2 
// BVT    = 3
// UNI    = 4 
// RDE    = 8 
// RDG    = 9 
// RDL    = 10 
// BVNzi  = 11
// BVEzi  = 12
// UNIzi  = 13
// RDEzi  = 14
// BVN2   = 16
// RDLS   = 17
// BVC    = 18
// BVCzi  = 19

void fillkernelp (
        const int    jj,          // number of sessions
        const int    kerneltype, 
        const bool   sparsekernel,
        const double cellsize,
        const double r0,          // effective radius of zero cell as proportion of cell width
        const RcppParallel::RMatrix<int> kernel, 
        const RcppParallel::RVector<int> moveargsi, 
        const std::string fnname,
        const std::vector<double> &moveargs, 
        std::vector<double> &kernelp,
        const bool   normalize,
        const int    grain,
        int    &returncode) 
{
    int j,k,x,y;
    int n1 = 0;
    int n2 = 0;
    double r,r2,a2,b2;
    double a = 1.0;
    double b = 1.0;
    double p0 = 1;
    double p1 = 0;
    double R = 0;
    double diag;
    double cellarea = cellsize * cellsize;  // square metres
    double r0w = r0 * cellsize;             // metres
    int kn = kernel.nrow();
    int centrek = 0;
    bool oneparameter = kerneltype==0 || kerneltype==1 || kerneltype==8 || 
        kerneltype == 13 || kerneltype == 18;
    bool twoparameter = kerneltype==3 || kerneltype==9 || kerneltype==10 ||
        kerneltype==11 || kerneltype==12 || kerneltype == 14 || kerneltype == 16 ||
        kerneltype == 17;
    bool zeroinflated = (kerneltype>=11 && kerneltype<=14) || (kerneltype == 19);
    std::vector<double> p(jj-1);
    std::vector<double> sumj(jj-1);
    std::vector<double> tempkernel(kn);
    
    // determine radius (number of cells)
    int K2 = 0;
    for (k=0; k<kn; k++) K2 = std::max(kernel[k], K2);  
    
    // are parameters for session j a repeat of the previous session?
    std::vector<bool> repeat(jj-1, false);
    for (j = 1; j < (jj-1); j++) {
        if (kerneltype == 4) {
            repeat[j] = true;    // uniform - no parameters
        }
        if (oneparameter || twoparameter) {
            repeat[j] = moveargs[j] == moveargs[j-1];   
        }
        if (twoparameter) {
            repeat[j]  = repeat[j] && (moveargs[j+jj]==moveargs[j-1+jj]);
        }
    }
    
    returncode = 1;
    
    //-------------------------------------------------------------------------
    // iterate over sessions, ignoring last
    
    for (j = 0; j < (jj-1); j++) {
        
        // copy from previous session if parameters the same
        if (repeat[j]) {
            for (k = 0; k < kn; k++) {
                kernelp[j * kn + k] =  kernelp[(j-1) * kn + k];
            }
            sumj[j] = sumj[j-1];
        }
        else
        {
            
            //----------------------------------------------
            // annular
            if (kerneltype == 5) {    
                p0 = moveargs[j];
            }
            // annular2
            if (kerneltype == 6) {    
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
                // debug
                // if (grain<0) Rprintf(" n1 %4d n2 %4d\n",  n1, n2); 
            }
            // annularR
            if (kerneltype == 7) {    
                p0 = moveargs[j];
                R = moveargs[j+jj];
                // compute arc lengths etc.
                tempkernel = annulus(p0, R, kernel);
                // debug
                // if (j==0) for (k = 0; k < kn; k++) {
                // if (grain<0) Rprintf(" k %4d tempkernel[k] %8.6f\n",  k, tempkernel[k]); 
                // }
            }
            //----------------------------------------------
            
            if (oneparameter || twoparameter) {
                a = moveargs[j];   
            }
            if (twoparameter) {
                b = moveargs[j+jj];   
            }
            if (a<=0 || (b<=0 && twoparameter)) {
                // if (grain<1 && a<=0) Rcpp::warning("C++ error: a<=0 when filling kernel");
                // if (grain<1 && (b<=0 && twoparameter)) Rcpp::warning("C++ error: b<=0 when filling kernel");
                fill(kernelp.begin(), kernelp.end(), NAN);
                returncode = -3;    // bad parameter
                return;
            }
            sumj[j] = 0;   // zero accumulator for session j
            a2 = a * a;
            b2 = b * b;
            
            // for each point in discretized kernel
            for (k = 0; k < kn; k++) {
                
                x = kernel[k];
                y = kernel[k+kn];
                // 2021-10-14 adjusted to include first circle of cells
                if (std::abs(x) == std::abs(y) && (std::abs(x)<1.5)) 
                    diag = std::sqrt(2); 
                else 
                    diag = 1;
                r2 = (x*x + y*y) * cellarea;
                r = std::sqrt(r2);
                if (x==0 && y==0) {
                    r = 0;          // to be sure
                    centrek = k;    // for post hoc zero inflation
                }
                
                // BVN or BVNzi Gaussian kernel 
                if (kerneltype == 0 || kerneltype == 11) {    
                    if (r>0 || r0==0) {
                        kernelp[j * kn + k] = std::exp(-r2 / 2 / a2)  / 2 / M_PI / a2 * cellarea;
                    }
                    else  {
                        kernelp[j * kn + k] = 1 - std::exp(-r0w*r0w/2/a2);
                    }
                }
                
                // equal mixture of two BVN kernels 
                else if (kerneltype == 16) {        
                    if (r>0 || r0==0) {
                        kernelp[j * kn + k] = 0.5 * cellarea * 
                        (std::exp(-r2 / 2 / a2)  / 2 / M_PI / a2 + 
                        std::exp(-r2 / 2 / b2)  / 2 / M_PI / b2);
                    }
                    else {
                        kernelp[j * kn + k] = 1 - 0.5 * (std::exp(-r0w*r0w/2/a2) + std::exp(-r0w*r0w/2/b2));
                    }
                }
                
                // BVE or BVEzi Laplace kernel 
                else if (kerneltype == 1 || kerneltype == 12) {   
                    if (r>0 || r0==0) {
                        kernelp[j * kn + k] = std::exp(-r / a) / 2 / M_PI / a2 * cellarea;
                    }
                    else {
                        kernelp[j * kn + k] = 1 - (r0w/a + 1) * std::exp(-r0w/a);
                    }
                }
                // BVC kernel 
                else if (kerneltype == 18 || kerneltype == 19) {   
                    if (r>0) {
                        kernelp[j * kn + k] =  1 / (2 * M_PI) * a / pow(r2 + a2, 1.5) * cellarea;       
                    }
                    else {
                        kernelp[j * kn + k] =   1 - a / std::sqrt(a2 + r0w*r0w); 
                    }
                    // Rprintf("r %8.6g, a %8.6g, a2 %8.6g, r0w %8.6g, p %8.6g \n",
                    //    r, a, a2, r0w, kernelp[j * kn + k]);
                }
                // BVT kernel 
                else if (kerneltype == 3) {   
                    if (r>0 || r0==0) {
                        kernelp[j * kn + k] = b / M_PI / a2 / pow(1 + r2/a2, b+1) * cellarea;
                    }
                    else {
                        kernelp[j * kn + k] = 1 - pow(a2 / (a2 + r0w*r0w), b);
                    }
                }
                
                // UNI or UNIzi kernel 
                else if (kerneltype == 4 || kerneltype == 13) {  
                    kernelp[j * kn + k] = 1.0 / kn;
                }
                
                // RDE or RDEzi kernel 
                else if (kerneltype == 8 || kerneltype == 14) {  
                    if (r>0)
                        kernelp[j * kn + k] =  std::exp(-r/a) / a / 2 / M_PI / r * cellarea;
                    else {
                        kernelp[j * kn + k] =  (1 - std::exp(-r0w/a));
                    }
                }
                
                // RDG kernel 
                else if (kerneltype == 9) {   
                    boost::math::gamma_distribution<> gam(b, a);  // shape, scale
                    if (r>0) {
                        kernelp[j * kn + k] =  boost::math::pdf(gam, r) / 2 / M_PI / r * cellarea;
                    }
                    else {
                        kernelp[j * kn + k] =  boost::math::cdf(gam, r0w);
                    }
                }
                
                // RDL kernel 
                else if (kerneltype == 10) {  
                    double mu = std::log(a); 
                    double s = std::sqrt(std::log(1 + 1/b));
                    // 2021-10-06 catch bad input 
                    if (std::isinf(s)) {
                        kernelp[j * kn + k] = 1/kn;    // uniform
                    }
                    else {
                        boost::math::lognormal_distribution<> ln (mu, s);
                        if (r>0) {
                            kernelp[j * kn + k] = boost::math::pdf(ln,r) / 2 / M_PI / r * cellarea;
                        }
                        else {
                            kernelp[j * kn + k] = boost::math::cdf(ln, r0w);
                        }
                    }
                }
                // RDLS kernel 
                else if (kerneltype == 17) {   
                    if (r>0) {
                        kernelp[j * kn + k] =  
                            2 / (M_PI * r * b) / (pow(r/a, 1/b) + pow(r/a, -1/b)) / 2 / M_PI/ r * cellarea;    
                    }
                    else {
                        double rtmp = cellsize/4;  // ad hoc because I don't know F()
                        kernelp[j * kn + k] =  
                            2 / (M_PI * rtmp * b) / (pow(rtmp/a, 1/b) + pow(rtmp/a, -1/b)) / 2 / M_PI/ rtmp * cellarea;
                    }
                }
                //------------------------------------------
                // oddball kernels
                // User kernel - not available when multithreaded (grain >= 1)
                else if (kerneltype == 2 && grain<1) {
                    // call R function from C++
                    // Rcpp::Environment env = Rcpp::Environment::global_env();
                    // Rcpp::Function f = env[fnname];
                    Rcpp::Function f(fnname);
                    if (moveargsi[1]>0)
                        p = Rcpp::as< std::vector<double> >(f(r, moveargs[j], moveargs[j+jj]));
                    else if (moveargsi[0]>0)
                        p = Rcpp::as< std::vector<double> >(f(r, moveargs[j]));
                    else
                        p = Rcpp::as< std::vector<double> >(f(r));
                    kernelp[j * kn + k] = p[0];
                }
                
                // annular kernel 
                else if (kerneltype == 5) {  
                    if (r<1e-8)
                        kernelp[j * kn + k] = p0;
                    else {
                        kernelp[j * kn + k] = (1-p0)/(kn-1);
                    }
                }
                // annular2 kernel 
                else if (kerneltype == 6) {  
                    if (r<1e-8)
                        kernelp[j * kn + k] = p0;
                    else if (r < (K2-1)) 
                        kernelp[j * kn + k] = p1/n1;
                    else 
                        kernelp[j * kn + k] = (1-p0-p1)/n2;
                }
                // annularR kernel 
                else if (kerneltype == 7) {  
                    // variable annulus precomputed
                    kernelp[j * kn + k] = tempkernel[k];
                }
                
                //------------------------------------------
                
                else {
                    // debug
                    // if (grain<1) Rcpp::stop("unrecognised kerneltype");
                    returncode = -2;
                    fill(kernelp.begin(), kernelp.end(), NAN);
                    return;
                }
                
                //------------------------------------------
                // wrap-up point k, session j
                
                // apply relative weighting for sparsekernel
                if (sparsekernel && r>0) {
                    // kernelp[j * kn + k] = 2 * M_PI * r * diag * kernelp[j * kn + k];
                    // kernelp[j * kn + k] = 2 * M_PI * r * diag * kernelp[j * kn + k] / 8;
                    // 2021-10-14 
                    kernelp[j * kn + k] = 2 * M_PI * r * diag * kernelp[j * kn + k] / 8 / cellsize;
                }
                
                // sum pdf across kernel for session j
                sumj[j] += kernelp[j * kn + k];
                
                //------------------------------------------
                
            }  // end loop over points
            
            if (sumj[j]<=0) {
                // debug
                // if (grain<1) Rcpp::warning("kernel probabilities sum to zero");
                fill(kernelp.begin(), kernelp.end(), NAN);
                returncode = -1;
                return;
            }
        }   // end conditional execution (repeat[j] false)
    }   // end loop over sessions
    
    //-------------------------------------------------------------------------
    
    // almost always normalize 
    if (kerneltype != 7) {
        for (j = 0; j < (jj-1); j++) {
            if (sumj[j]>0) {
                if (!normalize) sumj[j] = 1.0;    // override
                for (k = 0; k < kn; k++) {
                    if (zeroinflated) { // 
                        if (kerneltype == 13)  // uniformz
                            b = moveargs[j];   
                        else
                            b = moveargs[j+jj];   
                        if (k == centrek) {
                            kernelp[j * kn + k] = b + (1-b) * kernelp[j * kn + k] / sumj[j];
                        }
                        else {
                            kernelp[j * kn + k] = (1-b) * kernelp[j * kn + k] / sumj[j];
                        }
                    }
                    else {
                        kernelp[j * kn + k] = kernelp[j * kn + k] / sumj[j];
                    }
                }
            }
            else {
                fill(kernelp.begin(), kernelp.end(), NAN);
                returncode = -1;
                return;
            }
        }
    }
    
}   // end fillkernelp
//--------------------------------------------------------------------------


// fillkernelp wrapper for use from R 2021-07-27

// [[Rcpp::export]]
Rcpp::NumericVector fillkernelcpp (
        const Rcpp::IntegerMatrix kernel, 
        const int    kerneltype, 
        const bool   sparsekernel,
        const double cellsize,
        const double r0,
        const int    jj,
        const  std::string fnname,
        const  Rcpp::IntegerVector moveargsi, 
        const  Rcpp::NumericVector &moveargs,
        const bool   normalize
) 
{
    
    int grain = 0; // no need to change this for single-thread calls from R
    int returncode;
    int kn = kernel.nrow();
    
    // we need RcppParallel::RMatrix and RcppParallel::RVector objects to pass to fillkernelp
    // create these and transfer data
    RcppParallel::RMatrix<int> kernelmat(kernel); 
    RcppParallel::RVector<int> moveargsivec(moveargsi);
    std::vector<double> moveargsvec = Rcpp::as< std::vector<double> >(moveargs);
    std::vector<double> kernelp(kn*(jj-1));
    
    fillkernelp (
            jj, 
            kerneltype, 
            sparsekernel, 
            cellsize, 
            r0,
            kernelmat, 
            moveargsivec, 
            fnname, 
            moveargsvec, 
            kernelp, 
            normalize, 
            grain,
            returncode);
    
    if (returncode<0) kernelp[0] = NAN;
    
    return Rcpp::wrap(kernelp);   // Rcpp::NumericVector of length kn*(jj-1)
    
}

