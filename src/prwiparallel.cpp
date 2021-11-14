#include "utils.h"

//==============================================================================

struct Somehistories : public RcppParallel::Worker {
    
    // input data
    const int   x;
    const int   type;
    const int   nc;
    const int   CJSp1;                   
    const int   jj;
    const int   cc;
    const RcppParallel::RVector<double> intervals;
    const RcppParallel::RVector<int> cumss;
    const RcppParallel::RVector<int> w;
    const RcppParallel::RVector<int> fi;
    const RcppParallel::RVector<int> li;
    const RcppParallel::RMatrix<double> openval;
    const RcppParallel::RVector<int> PIA;
    const RcppParallel::RVector<int> PIAJ;
    
    // output likelihoods, one per detection history
    RcppParallel::RVector<double> output;
    
    // Constructor to initialize an instance of Somehistories 
    // The RcppParallel::RMatrix class can be automatically converted to from the Rcpp matrix type
    Somehistories(
        const int x, 
        const int type, 
        const int nc, 
        const int CJSp1, 
        const int jj, 
        const int cc,                   
        const Rcpp::NumericVector intervals,
        const Rcpp::IntegerVector cumss,
        const Rcpp::IntegerVector w,
        const Rcpp::IntegerVector fi, 
        const Rcpp::IntegerVector li,
        const Rcpp::NumericMatrix openval,
        const Rcpp::IntegerVector PIA,
        const Rcpp::IntegerVector PIAJ,
        Rcpp::NumericVector output)    
        : 
        x(x), type(type), nc(nc), CJSp1(CJSp1), jj(jj), cc(cc), 
        intervals(intervals), cumss(cumss), w(w), fi(fi), li(li),
        openval(openval), PIA(PIA), PIAJ(PIAJ), output(output) 
    {
    }
    
    double oneprwicpp (int n) {
        double pdt;
        double pbd;
        int j,s;
        int count; 
        bool dead;
        int cjs = 0;    // offset for first primary session (1 for CJS)
        int b,minb,maxb;
        int d,mind,maxd;

        // work vectors for session-specific real parameter values
        std::vector<double> p(cumss[jj]);  // each secondary session
        std::vector<double> phij(jj);      // each primary session
        std::vector<double> beta(jj);      // not used if type==1

        getp (n, x, nc, cumss[jj], openval, PIA, p);
        getphij (n, x, nc, jj, openval, PIAJ, intervals, phij);

        if (type == 1) {
            // cjs = 1;
            // 2018-10-29
            cjs = 1 - CJSp1;
        }
        else {
            getbeta (type, n, x, nc, jj, openval, PIAJ, intervals, phij, beta);
        }

        pdt = 0;
        if (type == 1)
            minb = fi[n];
        else 
            minb = 1;
        maxb = fi[n];
        mind = std::abs(li[n]);
        maxd = jj;
        if (li[n] < 0) maxd = mind; // possible censoring
        
        // loop over possible birth and death times
        for (b = minb; b <= maxb; b++) {
            for (d = mind; d <= maxd; d++) {
                dead = false;
                if (type == 1) 
                    pbd = 1.0;
                else
                    pbd = beta[b-1];
                for (j = b; j < d; j++) {
                    pbd *= phij[j-1];
                }
                if (li[n]>0)  // not censored
                    pbd *= 1-phij[d-1];
                // pbd now accounts for birth at b (JSSA) and survival to d
                // next multiply by conditional probability of observed CH 
                for (j = b + cjs; j <= d; j++) {
                    // detection probability for each secondary session
                    // in primary session j
                    //notseen = 1;
                    for (s = cumss[j-1]; s < cumss[j]; s++) {   
                        count = w[nc * s + n];
                        if (count<0) {count = -count; dead = true; }
                        
                        if (count>0) {
                            // notseen = 0;
                            pbd *= p[s];
                        }
                        else
                            pbd *= 1 - p[s];
                        if (dead) break;
                    }
                }   
                pdt += pbd;
            }
        }
        return pdt;
    }

    // function call operator that work for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {
        
        for (std::size_t n = begin; n < end; n++) {
            output[n] = oneprwicpp (n);
        }
    }
};

// [[Rcpp::export]]
Rcpp::NumericVector allhistparallelcpp (
        const int x, 
        const int type, 
        const int nc,
        const int CJSp1, 
        const int grain, 
        const int ncores,
        const Rcpp::NumericVector intervals, 
        const Rcpp::IntegerVector cumss, 
        const Rcpp::IntegerVector w,
        const Rcpp::IntegerVector fi, 
        const Rcpp::IntegerVector li,
        const Rcpp::NumericMatrix openval,
        const Rcpp::IntegerVector PIA, 
        const Rcpp::IntegerVector PIAJ) {
    
    Rcpp::NumericVector output(nc); 
    int jj = intervals.size()+1;
    int cc = openval.nrow();
    
    // Construct and initialise
    Somehistories somehist (x, type, nc, CJSp1, 
                            jj, cc, intervals,
                            cumss, w, fi, li, openval, PIA, PIAJ, output);
    
    Rcpp::checkUserInterrupt();
    
    if (ncores>1) {    
        // Run operator() on multiple threads
        RcppParallel::parallelFor(0, nc, somehist, grain, ncores);
    }
    else {
        // for debugging avoid multithreading and allow R calls e.g. Rprintf
        somehist.operator()(0,nc);    
    }
    
    // Return consolidated result
    return output;
}
//==============================================================================
