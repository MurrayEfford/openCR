#include "utils.h"

//-------------------------------------------------------------
// return JSSA probability animal n detected at least once   
// Pledger et al. 2010 Eqn (3) (mixtures outside)            
// 2020-10-25 changed 'pm' to 'alpha' so consistent with allhistsecrparallelcpp
// 2020-10-28 new algorithm for movement models
// See overlapping code in getpdotbd() of allhistsecrparallelcpp()
// also R-only version in PCH1.R
//-------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::NumericVector PCH1cpp (
        const int type, 
        const int x, 
        const int nc, 
        const int jj, 
        const Rcpp::IntegerVector cumss, int nmix,
        const Rcpp::NumericMatrix openval0, 
        const Rcpp::IntegerVector PIA0,  
        const Rcpp::IntegerVector PIAJ,  
        const Rcpp::NumericVector intervals)
{
    RcppParallel::RVector<int> cumssR(cumss);
    RcppParallel::RMatrix<double> openval0R(openval0);
    RcppParallel::RVector<int> PIA0R(PIA0);
    RcppParallel::RVector<int> PIAJR(PIAJ);
    RcppParallel::RVector<double> intervalsR(intervals);
    
    std::vector<double> pdt(nc);   
    double pbd, ptmp;
    int j,n,s;
    int b,d;
    std::vector<double> p(cumssR[jj]);
    std::vector<double> phij(jj);
    std::vector<double> g(jj);
    std::vector<double> beta(jj);
    
    for (n = 0; n<nc; n++) {
        getp (n, x, nc, cumssR[jj], openval0R, PIA0R, p);
        getphij (n, x, nc, jj, openval0R, PIAJR, intervalsR, phij);
        getg (type, n, x, nc, jj, openval0R, PIAJR, g);
        getbeta (type, n, x, nc, jj, openval0R, PIAJR, intervalsR, phij, beta);
        
        for (b = 1; b <= jj; b++) {
            for (d = b; d <= jj; d++) {
                pbd = beta[b-1];              // entered at b
                for (j = b; j < d; j++) {
                    pbd *= phij[j-1];         // survived
                }
                pbd *= 1-phij[d-1];           // departed at d
                
                ptmp = 1;
                for (j = b; j <= d; j++) {
                    for (s = cumssR[j-1]; s < cumssR[j]; s++) {
                        ptmp *= 1 - p[s];       // not detected
                    }
                }
                pdt[n] += pbd * (1 - ptmp);
            }
        }
    }
    // return pdt;
    return Rcpp::wrap(pdt);
}
//==============================================================================

// probability of not being detected in primary session j
// each mask point m
void pr0njmx (int n, int x, 
    const RcppParallel::RVector<int> cumss, 
    const int nc,  
    const int jj, 
    const int kk, 
    const int mm, 
    const int cc0, 
    const int binomN,  
    const RcppParallel::RVector<int> PIA0, 
    const RcppParallel::RVector<double> gk0, 
    const RcppParallel::RMatrix<double> Tsk, 
    std::vector<double> &pjm) {
    
    int c, i, j, k, m, s, ci, gi, pi, ss;
    double Tski;
    bool hazard;
    hazard = (binomN == -2) || (binomN == 0);
    
    if (hazard) for (i = 0; i < jj*mm; i++) pjm[i] = 0.0;
    else for (i = 0; i < jj*mm; i++) pjm[i] = 1.0;
    ss = cumss[jj];
    for (j = 0; j < jj; j++) {
        // consider occasions from focal session (j)
        for (s = cumss[j]; s < cumss[j+1]; s++) {
            for (k = 0; k < kk; k++) {
                ci =  i4(n, s, k, x, nc, ss, kk);
                c = PIA0[ci] - 1;
                if (c >= 0) {    
                    Tski = Tsk(k,s);
                    for (m = 0; m < mm; m++) {
                        gi = i3(c, k, m, cc0, kk);
                        pi = m * jj + j;
                        if (hazard) {			   
                            pjm[pi] += gk0[gi] * Tski;  // gk0 is actually hazard 
                        }
                        else {
                            pjm[pi] *=  pski(binomN, 0, Tski, gk0[gi]);
                        }
                    }
                }
            }
        }
    }
    if (hazard) {
        for (i = 0; i < jj*mm; i++) pjm[i] = std::exp(-pjm[i]);
    }
}
//==============================================================================

// [[Rcpp::export]]
Rcpp::NumericVector PCH0secrjcpp (
        const int type, 
        const int x, 
        const int nc, 
        const int jj,
        const Rcpp::IntegerVector cumss, 
        const int kk, 
        const int mm, 
        const int cc0,
        const Rcpp::IntegerVector PIA0, 
        const Rcpp::NumericVector gk0, 
        const int binomN, 
        const Rcpp::NumericMatrix Tsk) {
    
    const RcppParallel::RVector<int> cumssR(cumss); 
    const RcppParallel::RVector<int> PIA0R(PIA0); 
    const RcppParallel::RVector<double> gk0R(gk0);
    const RcppParallel::RMatrix<double> TskR(Tsk);
    
    // by primary session
    // Rcpp::NumericVector pdt(nc * jj, 0.0);
    std::vector<double> pdt(nc * jj);
    
    int j, m, n;
    std::vector<double> pjmat(jj*mm);
    for (n = 0; n<nc; n++) {
        // primary-session-specific Pr for this animal
        pr0njmx(n, x, cumssR, nc, jj, kk, mm, cc0, binomN, PIA0R, gk0R, TskR, pjmat);
        for (j=0; j<jj; j++) {
            for (m=0; m<mm; m++) {
                pdt[j * nc + n] += pjmat[m * jj + j] / mm;
            }
        }
    }
    // return pdt;
    return Rcpp::wrap(pdt);
}
//==============================================================================

struct pch1struct : public RcppParallel::Worker {
    const int x;
    const int type;
    const int grain;
    const int jj;
    const int mm;
    const int nc;
    const RcppParallel::RVector<int> cumss;
    const RcppParallel::RMatrix<double> openval0;
    const RcppParallel::RVector<int> PIA0;
    const RcppParallel::RVector<int> PIAJ;
    const RcppParallel::RVector<double> gk0;
    const int   binomN;
    const RcppParallel::RMatrix<double> Tsk;
    const RcppParallel::RVector<double> intervals;
    const RcppParallel::RVector<int> moveargsi;
    const int   movementcode;
    const bool  sparsekernel;
    const bool  anchored;
    const int   edgecode;
    const std::string usermodel;
    const RcppParallel::RMatrix<int> kernel;
    const RcppParallel::RMatrix<int> mqarray;
    const double cellsize;
    const double r0;
    const RcppParallel::RMatrix<double> settlement;

    int   kk;
    int   kn;
    int   cc0;
    int   fillcode;
    
    // output likelihoods, one per animal
    RcppParallel::RVector<double> output;
    
    // Constructor to initialize an instance of Somehistories 
    // The RMatrix class can be automatically converted to from the Rcpp matrix type
    pch1struct(
        const int x, 
        const int type, 
        const int grain, 
        const int jj, 
        const int mm, 
        const int nc, 
        const Rcpp::IntegerVector cumss,
        const Rcpp::NumericMatrix openval0,
        const Rcpp::IntegerVector PIA0,
        const Rcpp::IntegerVector PIAJ,
        const Rcpp::NumericVector gk0,
        const int binomN,
        const Rcpp::NumericMatrix Tsk,
        const Rcpp::NumericVector intervals,
        const Rcpp::IntegerVector moveargsi,
        const int   movementcode,
        const bool  sparsekernel,
        const bool  anchored,
        const int   edgecode,
        const std::string usermodel,
        const Rcpp::IntegerMatrix kernel,
        const Rcpp::IntegerMatrix mqarray,
        const double cellsize,        
        const double r0,
        const Rcpp::NumericMatrix settlement,
        
        Rcpp::NumericVector output)    
        : 
            x(x), type(type), grain(grain), jj(jj), mm(mm), nc(nc), 
            cumss(cumss), openval0(openval0), PIA0(PIA0), PIAJ(PIAJ), 
            gk0(gk0), binomN(binomN), Tsk(Tsk), intervals(intervals), 
            moveargsi(moveargsi), movementcode(movementcode), 
            sparsekernel(sparsekernel), anchored(anchored), edgecode(edgecode), 
            usermodel(usermodel), kernel(kernel), mqarray(mqarray), 
            cellsize(cellsize), r0(r0), settlement(settlement), 
            output(output) 
    {
        kn = kernel.nrow();
        cc0 = openval0.nrow();
        kk = Tsk.nrow();
        fillcode = movementcode-2;
    }
    
    double onepch1cpp (int n) {
        double pdt = 0.0;
        double pbd;
        int j,m, m2, q;
        int b,d;
        double prw0, sumpj;
        int kernelreturncode = 1;
        double bz;
        double sumalpha;
        double prwm;
        
        // work vectors for session-specific real parameter values
        std::vector<double> phij(jj);      // each primary session
        std::vector<double> beta(jj);      // not used if type==1
        std::vector<double> alpha(mm);
        std::vector<double> pjmat(jj * mm);
        std::vector<double> moveargs(jj*2);
        std::vector<double> kernelp(kn * (jj-1));
        std::vector<int>    mj(kn);
        std::vector<double> pj(kn);
        
        getphij (n, x, nc, jj, openval0, PIAJ,  intervals, phij);
        getbeta (type, n, x, nc, jj, openval0, PIAJ, intervals, phij, beta);
        // primary-session-specific Pr for this animal
        pr0njmx(n, x, cumss, nc, jj, kk, mm, cc0, binomN, PIA0, gk0, Tsk, pjmat);
        if (movementcode>1) {
            getmoveargs (n, x, nc, jj, openval0, PIAJ, moveargsi, moveargs);
            // if (grain<0) for (j=0;j<(jj-1); j++) Rprintf("j %d moveargs[j] %8.6g \n", j, moveargs[j]);
            if (movementcode != 17) {
                fillkernelp (jj, fillcode, sparsekernel, cellsize, r0,
                    kernel, moveargsi, usermodel, moveargs, kernelp, true, grain,
                    kernelreturncode);
            }
            if (kernelreturncode<0) return NAN;
        }
        
        for (b = 1; b <= jj; b++) {
            for (d = b; d <= jj; d++) {
                pbd = beta[b-1];
                for (j = b; j < d; j++) {
                    pbd *= phij[j-1];
                }
                pbd *= 1-phij[d-1];
                
                if (movementcode == 0) {
                    //-------------------------------------------------------------
                    // static home ranges: take sum over M of product over J
                    for (m=0; m<mm; m++) alpha[m] = 1.0/mm;
                    for (j = b; j <= d; j++) {
                        for (m=0; m<mm; m++) alpha[m] *= pjmat[m*jj + j - 1];
                    }
                    prw0 = std::accumulate(alpha.begin(), alpha.end(), 0.0);
                }
                else if (movementcode == 1) {
                    //-------------------------------------------------------------
                    // uncorrelated home ranges: take product over J of sum over M
                    // over primary sessions in which may have been alive 
                    // centers allowed to differ between primary sessions 
                    prw0 = 1.0;
                    for (j = b; j <= d; j++) {
                        sumpj = 0;
                        for (m=0; m<mm; m++) {
                            sumpj += pjmat[m * jj + j -1];
                        }
                        prw0 *= sumpj / mm;   
                    }
                }
                else if (anchored) {
                        std::fill(alpha.begin(), alpha.end(), 1.0/mm);
                        for (m=0; m<mm; m++) {
                            convolvemq1(m, 1, edgecode, mqarray, settlement, kernelp, mj, pj); 
                            for (j = b; j <= d; j++) {
                                prwm = 0;
                                for (q=0; q<kn; q++) {
                                    m2 = mj[q];
                                    if (m2>=0) {
                                        prwm += pj[q] * pjmat[m2*jj + j - 1];                           
                                    }
                                }
                                alpha[m] *= prwm;
                            }
                        }
                        prw0 = std::accumulate(alpha.begin(), alpha.end(), 0.0);
                    }
                else {  // movementcode>1
                    //-------------------------------------------------------------
                    // moving home ranges
                    for (m=0; m<mm; m++) alpha[m] = pjmat[m*jj + b - 1] / mm;
                    for (j = b+1; j <= d; j++) {
                        if (movementcode == 17) {  // uncorrelated zero-inflated
                            // bz static, 1-bz uncorrelated
                            // update alpha
                            bz = moveargs[j-2]; 
                            sumalpha = std::accumulate(alpha.begin(), alpha.end(), 0.0);
                            for (m=0; m<mm; m++) {
                                alpha[m] = bz * alpha[m] + (1-bz) * 1/mm * sumalpha;
                            }
                        } 
                        else {
                            convolvemq(mm, kn, j-1, edgecode, mqarray, settlement, kernelp, alpha);
                        }
                        for (m=0; m<mm; m++) alpha[m] *= pjmat[m*jj + j - 1];
                    }
                    prw0 = std::accumulate(alpha.begin(), alpha.end(), 0.0);
                }
                pdt += pbd * (1 - prw0);
            }
        }
        return pdt;
    }
    
    // function call operator that work for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {
        
        for (std::size_t n = begin; n < end; n++) {
            output[n] = onepch1cpp (n);
        }
    }
};

// [[Rcpp::export]]
Rcpp::NumericVector PCH1secrparallelcpp (
        const int   x, 
        const int   type, 
        const int   grain, 
        const int   ncores,
        const bool  individual,
        const int   jj,  
        const int   mm, 
        const int   nc,
        const Rcpp::IntegerVector cumss, 
        const Rcpp::NumericMatrix openval0, 
        const Rcpp::IntegerVector PIA0, 
        const Rcpp::IntegerVector PIAJ, 
        const Rcpp::NumericVector gk0, 
        const int   binomN, 
        const Rcpp::NumericMatrix Tsk, 
        const Rcpp::NumericVector intervals,
        const Rcpp::IntegerVector moveargsi,
        const int   movementcode,
        const bool  sparsekernel,
        const bool  anchored,
        const int   edgecode,
        const std::string usermodel,
        const Rcpp::IntegerMatrix kernel,
        const Rcpp::IntegerMatrix mqarray,
        const double cellsize,
        const double r0,
        const Rcpp::NumericMatrix settlement) {
    
    Rcpp::NumericVector output(nc); 

    // Construct and initialise
    pch1struct pch1 (x, type, grain, jj, mm, nc, 
        cumss, openval0, PIA0, PIAJ, gk0, binomN, Tsk, intervals,
        moveargsi, movementcode, sparsekernel, anchored, edgecode, 
        usermodel, kernel, mqarray, cellsize, r0, settlement, output);
    
    Rcpp::checkUserInterrupt();
    
    if (individual) {
        if (ncores>1) {
            // Run operator() on multiple threads
            RcppParallel::parallelFor(0, nc, pch1, grain, ncores);
        }
        else {
            pch1.operator()(0,nc);    // for debugging avoid multithreading to allow R calls
        }
    }
    else {
        // all the same: do once, then copy
        pch1.operator()(0,1);
        for (int i = 1; i < nc; i++) output[i] = output[0];
    }
    
    // Return consolidated result
    return output;
}
//==============================================================================
