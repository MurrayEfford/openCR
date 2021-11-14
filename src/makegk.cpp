#include "utils.h"

//==============================================================================

// [[Rcpp::export]]
Rcpp::List makegkcpp (
        const int cc, 
        const int kk, 
        const int mm, 
        const int detectfn, 
        const int sigmai, 
        const Rcpp::NumericMatrix openval, 
        const Rcpp::NumericMatrix traps,
        const Rcpp::NumericMatrix mask
) 
{
    const RcppParallel::RMatrix<double> openvalR(openval); 
    const RcppParallel::RMatrix<double> trapsR(traps);
    const RcppParallel::RMatrix<double> maskR(mask);
    int k, m, c, gi;
    Rcpp::NumericVector gk(cc * kk * mm); 
    Rcpp::NumericVector hk(cc * kk * mm);
    for (k=0; k<kk; k++) {
        for (m=0; m<mm; m++) {
            for (c=0; c<cc; c++) {
                gi = i3(c,k,m,cc, kk);
                hk[gi] = hfn(k, m, c, openvalR, trapsR, maskR, sigmai, detectfn);
                gk[gi] = 1 - std::exp(-hk[gi]);
            }
        }
    }
    return Rcpp::List::create(gk, hk);
}

//==============================================================================

struct Hckm : public RcppParallel::Worker {
    
    // input data
    const int sigmai;
    const int detectfn;
    const RcppParallel::RMatrix<double> openval;
    const RcppParallel::RMatrix<double> traps;
    const RcppParallel::RMatrix<double> mask;
    
    // output vector to write to
    RcppParallel::RVector<double> hk;
    RcppParallel::RVector<double> gk;
    
    // initialize from Rcpp input and output matrixes (the RMatrix class
    // can be automatically converted to from the Rcpp matrix type)
    Hckm(int sigmai, int detectfn,
         const Rcpp::NumericMatrix openval, 
         const Rcpp::NumericMatrix traps, 
         const Rcpp::NumericMatrix mask, 
         Rcpp::NumericVector hk,
         Rcpp::NumericVector gk)
        : sigmai(sigmai), detectfn(detectfn), openval(openval), 
          traps(traps), mask(mask), hk(hk), gk(gk) {}

    // function call operator that work for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {
        int cc = openval.nrow();
        int kk = traps.nrow();
        for (std::size_t m = begin; m < end; m++) {
            for (int k=0; k < kk; k++) {
                for (int c=0; c < cc; c++) {
                    int gi = i3(c,k,m, cc, kk);
                    hk[gi] = hfn(k, m, c, openval, traps, mask, sigmai, detectfn);
                    gk[gi] = 1 - std::exp(-hk[gi]);
                }
            }
            
        }
    }
};

// [[Rcpp::export]]
Rcpp::List makegkParallelcpp (
        const int detectfn, 
        const int sigmai, 
        const int grain, 
        const int ncores,
        const Rcpp::NumericMatrix& openval, 
        const Rcpp::NumericMatrix& traps,
        const Rcpp::NumericMatrix& mask
) 
{
    Rcpp::NumericVector hk(openval.nrow() * traps.nrow() * mask.nrow()); 
    Rcpp::NumericVector gk(openval.nrow() * traps.nrow() * mask.nrow()); 
    
    Hckm hckm (sigmai, detectfn, openval, traps, mask, hk, gk);
    
    if (ncores>1) {
        RcppParallel::parallelFor(0, mask.nrow(), hckm, grain, ncores);
    }
    else {
        hckm.operator()(0,mask.nrow());    // for debugging avoid multithreading to allow R calls
    }
    return Rcpp::List::create(gk, hk);
}
//==============================================================================

struct Hckmd : public RcppParallel::Worker {
    
    // input data
    const int sigmai;
    const int detectfn;
    const RcppParallel::RMatrix<double> openval;
    const RcppParallel::RMatrix<double> distmat;

    // output vector to write to
    RcppParallel::RVector<double> hk;
    RcppParallel::RVector<double> gk;
    
    // initialize from Rcpp input and output matrixes (the RMatrix class
    // can be automatically converted to from the Rcpp matrix type)
    Hckmd(
        const int sigmai, 
        const int detectfn,
        const Rcpp::NumericMatrix openval, 
        const Rcpp::NumericMatrix distmat, 
        Rcpp::NumericVector hk,
        Rcpp::NumericVector gk)
        : sigmai(sigmai), detectfn(detectfn), openval(openval), 
            distmat(distmat), hk(hk), gk(gk) {}
    
    // function call operator that work for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {
        int cc = openval.nrow();
        int kk = distmat.nrow();
        for (std::size_t m = begin; m < end; m++) {
            for (int k=0; k < kk; k++) {
                for (int c=0; c < cc; c++) {
                    int gi = i3(c,k,m, cc, kk);
                    hk[gi] = hfnd(k, m, c, openval, distmat, sigmai, detectfn);
                    gk[gi] = 1 - std::exp(-hk[gi]);
                }
            }
            
        }
    }
};

// [[Rcpp::export]]
Rcpp::List makegkParalleldcpp (
        const int detectfn, 
        const int sigmai, 
        const int grain, 
        const int ncores,
        const Rcpp::NumericMatrix& openval, 
        const Rcpp::NumericMatrix& distmat
) 
{
    Rcpp::NumericVector hk(openval.nrow() * distmat.nrow() * distmat.ncol()); 
    Rcpp::NumericVector gk(openval.nrow() * distmat.nrow() * distmat.ncol()); 
    
    Hckmd hckmd (sigmai, detectfn, openval, distmat, hk, gk);
    
    if (ncores>1) {
        RcppParallel::parallelFor(0, distmat.ncol(), hckmd, grain, ncores);
    }
    else {
        hckmd.operator()(0,distmat.ncol());    // for debugging avoid multithreading to allow R calls
    }
    return Rcpp::List::create(gk, hk);
}
//==============================================================================