#include <Rcpp.h>
#include <RcppParallel.h>

// see https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/stat_tut/weg/error_eg.html
// and https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/pol_tutorial/changing_policy_defaults.html

// return NAN for invalid inputs
#define BOOST_MATH_DOMAIN_ERROR_POLICY ignore_error
#include <boost/math/distributions.hpp>       // gamma, normal, poisson, binomial, lognormal distributions

#define fuzz 1e-200
#define huge 1e10

int i3 (int i, int j, int k, int ii, int jj);
int i4 (int i, int j, int k, int l, int ii, int jj, int kk);

//--------------------------------------------------------------------------

double gbinom(int count, int size, double p);
double pski ( int binomN, int count, double Tski, double g);

//--------------------------------------------------------------------------

void convolvemq (
        const int    mm,                            // number of points on mask 
        const int    kn,                            // number of points on kernel 
        const int    j,                             // session number 1..jj 
        const int    edgecode,                      // adjust for incomplete kernel
        const  RcppParallel::RMatrix<int> &mqarray, // input [& 2020-10-31]
        std::vector<double> &kernelp,               // p(move|dx,dy) for points in kernel 
        std::vector<double> &pjm                    // return value 
);
//--------------------------------------------------------------------------

void convolvemq1 (
        const int    m,                     // initial location on mask
        const int    j,                     // session number 1..jj 
        const int    edgecode,              // 0 none, no action; 1 wrapped, no action; 2 normalize truncated kernel
        const  RcppParallel::RMatrix<int> &mqarray, // input [& 2020-10-31]
        const  std::vector<double> &kernelp, // p(move|dx,dy) for points in kernel 
        std::vector<int>    &mj,
        std::vector<double> &pj      // return value
);
//--------------------------------------------------------------------------

void fillkernelp (
        const int                        jj, 
        const int                        kerneltype, 
        const bool                       sparsekernel,
        const double                     cellsize,
        const double                     r0,
        const RcppParallel::RMatrix<int> kernel, 
        const RcppParallel::RVector<int> moveargsi, 
        const std::string                fnname,
        const std::vector<double>        &moveargs, 
        std::vector<double>              &kernelp,
        const bool                       normalize,
        const int                        grain,
        int                              &returncode
);

//--------------------------------------------------------------------------

double hfn (
        const int                           k,                         
        const int                           m, 
        const int                           c, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RMatrix<double> traps,
        const RcppParallel::RMatrix<double> mask, 
        const int                           sigmai, 
        const int                           detectfn);
//--------------------------------------------------------------------------

double hfnd (
        const int k, 
        const int m, 
        const int c, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RMatrix<double> distmat,
        const int sigmai, 
        const int detectfn
);
//--------------------------------------------------------------------------

void getp (
        const int n, 
        const int x, 
        const int nc, 
        const int ss, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIA, 
        std::vector<double> &p
);
//--------------------------------------------------------------------------

void getphij (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ,
        const RcppParallel::RVector<double> intervals, 
        std::vector<double> &phij
);
//--------------------------------------------------------------------------

void getmoveargs (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ,
        const RcppParallel::RVector<int>    moveargsi,
        std::vector<double> &moveargs
);
//--------------------------------------------------------------------------

void getpj (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ,
        std::vector<double> &pj
);
//--------------------------------------------------------------------------

void getg (
        const int type, 
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ,
        std::vector<double> &g
);
//--------------------------------------------------------------------------

void getfj (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ,
        const RcppParallel::RVector<double> intervals, 
        const std::vector<double> phij,
        std::vector<double> &fj
);
//--------------------------------------------------------------------------

void getlj (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ,
        const RcppParallel::RVector<double> intervals, 
        std::vector<double> &lj
);
//--------------------------------------------------------------------------

void getgaml (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ,
        const RcppParallel::RVector<double> intervals, 
        std::vector<double> &gam);
//--------------------------------------------------------------------------

void getgamj (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ,
        const RcppParallel::RVector<double> intervals, 
        std::vector<double> &gamj
);
//--------------------------------------------------------------------------

void getkapj (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ,
        std::vector<double> &kapj
);
//--------------------------------------------------------------------------

void getbeta0 (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ,
        std::vector<double> &beta
);
//--------------------------------------------------------------------------

void gettau (
        const int n, 
        const int x, 
        const int nc, 
        const int jj,
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ,
        std::vector<double> &tau,
        const int M
);
//--------------------------------------------------------------------------

void getDj (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ,
        std::vector<double> &Dj
);
//--------------------------------------------------------------------------

// per capita recruitment cf Link & Barker 2005, Schwarz 'Gentle Intro'
void getbetaf (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ,
        const std::vector<double> phij,
        const RcppParallel::RVector<double> intervals, 
        std::vector<double> &beta
);
//--------------------------------------------------------------------------

void getbetal (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ,
        const std::vector<double> phij, 
        const RcppParallel::RVector<double> intervals, 
        std::vector<double> &beta
);
//--------------------------------------------------------------------------

void getbetag (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ, 
        const std::vector<double> phij, 
        const RcppParallel::RVector<double> intervals, 
        std::vector<double> &beta
);
//--------------------------------------------------------------------------

void getbetak (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ, 
        const std::vector<double> phij,
        std::vector<double> &beta
);
//--------------------------------------------------------------------------

// return parameterisation cf Pledger et al. 2010 p 885 
void getbetaB (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ, 
        std::vector<double> &beta
);
//--------------------------------------------------------------------------

void getbetaD (
        const int n, 
        const int x,
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ, 
        const std::vector<double> phij,
        std::vector<double> &beta
);
//--------------------------------------------------------------------------

void getbeta (
        const int type, 
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int>    PIAJ, 
        const RcppParallel::RVector<double> intervals,
        const std::vector<double> phij,
        std::vector<double> &beta
);
//--------------------------------------------------------------------------

Rcpp::List makelookupcpp (const Rcpp::NumericMatrix x); 

//--------------------------------------------------------------------------

void pr0njmx (
        const int n, 
        const int x, 
        const RcppParallel::RVector<int> cumss, 
        const int nc,  
        const int jj, 
        const int kk, 
        const int mm, 
        const int cc0, 
        const int binomN,  
        const RcppParallel::RVector<int>    PIA0, 
        const RcppParallel::RVector<double> gk0, 
        const RcppParallel::RMatrix<double> Tsk, 
        std::vector<double> &pjm
);
//--------------------------------------------------------------------------

