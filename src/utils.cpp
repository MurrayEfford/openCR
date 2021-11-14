#include "utils.h"

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

// customised dbinom 
double gbinom(int count, int size, double p)
{
    double x, q;
    int i;
    if ((count < 0) || (count > 0 && p <= 0)) {
        x = 0;
    }
    else if (count == 0) {
        if (size== 0)
            x = 1;
        else {
            q = 1 - p;
            x = q;
            for (i=1; i< size; i++) x *= q;
        }
    }
    else {
        boost::math::binomial_distribution<> bin(size, p);
        x = boost::math::pdf(bin, count);
    }
    return (x);   
}

//--------------------------------------------------------------------------

// probability of count for session s, detector k, animal i
// The argument 'g' is understood to be a cumulative hazard if binomN=0,
// a probability otherwise

double pski ( int binomN,
    int count,
    double Tski,
    double g) {
    
    double lambda = 0.0;
    double result = 1.0;
    
    if (binomN == -1) {                              // binary proximity detectors : Bernoulli
        if (std::abs(Tski-1) > 1e-10) {                   // effort not unity; adjust g 
            g = 1 - std::pow(1 - g, Tski);
        }
        if (count>0)                                 
            result = g;  
        else 
            result = 1 - g;
    }
    else if (binomN == 0) {                          // count detectors : Poisson 
        lambda = Tski * g;
        if ((count < 0) || (count>0 && lambda<=0)) {         
            result = 0;
        }
        else if (count == 0) {
            result = std::exp(-lambda);            // routinely apply Tsk adjustment to cum. hazard 
        }
        else {
            //result = R::dpois(count, Tski * g * pI, 0); 
            boost::math::poisson_distribution<> pois(lambda);
            result = boost::math::pdf(pois,count);
        }
    }
    else if (binomN == 1) {                          // count detectors : Binomial, size from Tsk
        result = gbinom (count, std::round(Tski), g); 
    }
    else if (binomN > 1) {                           // count detectors : Binomial, specified size 
        if (std::abs(Tski-1) > 1e-10) {                   // effort not unity, adjust g 
            g = 1 - std::pow(1 - g, Tski);
        }
        result = gbinom (count, binomN, g);
    }
    else result = NAN;   // code multi -2 separately
    // Rcpp::stop("binomN < -1 not allowed");
    
    return (result);
}
//--------------------------------------------------------------------------

// distance between two points given by row k in traps and row m in mask
double dkm (
        const int k, 
        const int m, 
        const RcppParallel::RMatrix<double> &traps, 
        const RcppParallel::RMatrix<double> &mask
)
{
    return(std::sqrt((traps(k,0) - mask(m,0)) * (traps(k,0) - mask(m,0)) +
        (traps(k,1) - mask(m,1)) * (traps(k,1) - mask(m,1))));
}
//--------------------------------------------------------------------------

// rectangular distance between two points given by row k in traps and row m in mask
double dkmrect (
        const int k, 
        const int m, 
        const RcppParallel::RMatrix<double> &traps, 
        const RcppParallel::RMatrix<double> &mask
    )
{
    return(std::max(std::fabs(traps(k,0) - mask(m,0)), 
                    std::fabs(traps(k,1) - mask(m,1))));
}
//--------------------------------------------------------------------------

// parameters in openval ordered g0, phi, f, N, sigma, pmix 

// hazard detection functions 14-20
double hfn (
        const int k, 
        const int m, 
        const int c, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RMatrix<double> traps,
        const RcppParallel::RMatrix<double> mask, 
        const int sigmai, 
        const int detectfn)
{
    double d;
    double sigma;
    double z = 1;
    sigma =  openval(c, sigmai);
    if (detectfn == 20) {
        d = dkmrect(k, m, traps, mask);
    }
    else {
        d = dkm(k, m, traps, mask);
    }
    
    if (detectfn == 14) { 
        return (openval(c,0) * std::exp(-d*d/2/sigma/sigma));                   // HHN
    }
    else if (detectfn == 16) {
        return (openval(c,0) * std::exp(-d/sigma));                             // HEX
    }
    else if (detectfn == 20) {
        if (d<=sigma)
            return (openval(c,0));
        else
            return (0);                                                    // HPX
            // return (1e-5);                                              // arbitrary to avoid NA
    }
    else {
        z =  openval(c,sigmai+1);
        if (detectfn == 15) {
            return (openval(c,0) * (1 - std::exp(-std::pow(d/sigma, -z))));      // HHR
        }
        else if (detectfn == 17) {
            return (openval(c,0) * std::exp(-(d-z)*(d-z) / 2 / sigma / sigma));  // HAN
        }
        else if (detectfn == 18) {
            // return (openval(c,0) * R::pgamma(d,z,sigma/z,0,0));               // HCG
            boost::math::gamma_distribution<> gam(z,sigma/z);
            // use complement for upper tail
            return (openval(c,0) * boost::math::cdf(boost::math::complement(gam, d)));                     // HCG
        }
        else if (detectfn == 19) {
            return (openval(c,0) * std::exp(- std::pow(d /sigma , z)));          // HVP
        }
        else return NAN; // Rcpp::stop("detectfn not allowed in openCR");
    }
}
//--------------------------------------------------------------------------

// hazard detection functions 14-20
double hfnd (
        const int k, 
        const int m, 
        const int c, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RMatrix<double> distmat,
        const int sigmai, 
        const int detectfn
    )
{
    double d;
    double sigma;
    double z = 1;
    sigma =  openval(c, sigmai);
    d = distmat(k,m);
    if (detectfn == 14) { 
        return (openval(c,0) * std::exp(-d*d/2/sigma/sigma));                   // HHN
    }
    else if (detectfn == 16) {
        return (openval(c,0) * std::exp(-d/sigma));                             // HEX
    }
    else if (detectfn == 20) {
        if (d<=sigma)
            return (openval(c,0));
        else
            return (0);                                                    // HPX
        // return (1e-5);                                              // arbitrary to avoid NA
    }
    else {
        z =  openval(c,sigmai+1);
        if (detectfn == 15) {
            return (openval(c,0) * (1 - std::exp(-std::pow(d/sigma, -z))));          // HHR
        }
        else if (detectfn == 17) {
            return (openval(c,0) * std::exp(-(d-z)*(d-z) / 2 / sigma / sigma)); // HAN
        }
        else if (detectfn == 18) {
            return (openval(c,0) * R::pgamma(d,z,sigma/z,0,0));            // HCG
        }
        else if (detectfn == 19) {
            return (openval(c,0) * std::exp(- std::pow(d /sigma , z)));              // HVP
        }
        else return NAN; // Rcpp::stop("detectfn not allowed in openCR");
    }
}
//--------------------------------------------------------------------------



void getp (
        const int n, 
        const int x, 
        const int nc, 
        const int ss, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIA, 
        std::vector<double> &p
) {
    // column 1 
    int s;
    for (s = 0; s < ss; s++) {
        p[s] = openval(PIA[i3(n, s, x, nc, ss )]-1, 0); 
    }
}
//--------------------------------------------------------------------------

void getphij (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ,
        const RcppParallel::RVector<double> intervals, 
        std::vector<double> &phij) {
    // column 2 
    int j;
    double phi;
    for (j = 0; j < (jj-1); j++) {
        // jj-1 because one fewer intervals than primary sessions  
        phi = openval(PIAJ[i3(n, j, x, nc, jj)]-1, 1);       
        // adjust for interval duration  
        phij[j] = std::exp(std::log(phi) * intervals[j]);  
    }
    phij[jj-1] = 0;
}
//--------------------------------------------------------------------------

void getmoveargs (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ,
        const RcppParallel::RVector<int> moveargsi,
        std::vector<double> &moveargs) {
    // column moveargsi (and maybe moveargsi + 1) 
    int j;
    int moveaindex;
    int movebindex;
    // jj-1 because one fewer intervals than primary sessions  
    for (j = 0; j < (jj-1); j++) {
        // 2021-05-11 more explicit to avoid overflow bug
        // PIAJ has dimension (n.j.x)
        // openval has dimension (c, nrealpar)
        moveaindex = moveargsi[0];
        movebindex = moveargsi[1];
        if (moveaindex>=0) {
            moveargs[j] = openval(PIAJ[i3(n, j, x, nc, jj)]-1, moveaindex);  
        }
        if (movebindex>moveaindex) {
            moveargs[j+jj] = openval(PIAJ[i3(n, j, x, nc, jj)]-1, movebindex);  
        }
    }
    moveargs[jj-1] = 0;    // no movement beyond final session (jj)
    moveargs[2*jj-1] = 0;  // no movement beyond final session (jj)
}
//--------------------------------------------------------------------------

void getpj (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ,
        std::vector<double> &pj) 
{
    // column 2 
    int j;
    for (j = 0; j < jj; j++) {
        pj[j] = openval(PIAJ[i3(n, j, x, nc, jj)]-1, 0);       
    }
}
//--------------------------------------------------------------------------

void getg (
        const int type, 
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ,
        std::vector<double> &g) 
{
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

void getfj (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ,
        const RcppParallel::RVector<double> intervals, 
        const std::vector<double> phij,
        std::vector<double> &fj) 
{
    // column 3 
    int j;
    double f,phi;
    for (j = 0; j < (jj-1); j++) {
        // jj-1 because one fewer intervals than primary sessions  
        f = openval(PIAJ[i3(n, j, x, nc, jj )]-1, 2); 
        phi = std::exp(std::log(phij[j]) / intervals[j]);
        // adjust for interval duration  
        fj[j] = std::exp(std::log(phi+f) * intervals[j]) - phij[j];  
    }
    fj[jj-1] = 0;
}
//--------------------------------------------------------------------------

void getlj (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ,
        const RcppParallel::RVector<double> intervals, 
        std::vector<double> &lj) 
{
    // column 3 
    int j;
    double l;
    for (j = 0; j < (jj-1); j++) {
        // jj-1 because one fewer intervals than primary sessions  
        l = openval(PIAJ[i3(n, j, x, nc, jj )]-1, 2); 
        // adjust for interval duration  
        lj[j] = std::exp(std::log(l) * intervals[j]);  
    }
    lj[jj-1] = 0;
}
//--------------------------------------------------------------------------

void getgaml (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ,
        const RcppParallel::RVector<double> intervals, 
        std::vector<double> &gam) 
{
    // column 3 
    int j;
    double phij;
    double lamj;
    for (j = 0; j < (jj-1); j++) {
        phij = openval(PIAJ[i3(n, j, x, nc, jj)]-1, 1); 
        phij = std::exp(std::log(phij) * intervals[j]);  
        lamj = openval(PIAJ[i3(n, j, x, nc, jj)]-1, 2); 
        lamj = std::exp(std::log(lamj) * intervals[j]);  
        gam[j+1] = phij/lamj;
    }
    gam[0] = 0;
}
//--------------------------------------------------------------------------

void getgamj (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ,
        const RcppParallel::RVector<double> intervals, 
        std::vector<double> &gamj) 
{
    // column 3 
    int j;
    double gam;
    for (j = 1; j < jj; j++) {
        gam = openval(PIAJ[i3(n, j, x, nc, jj)]-1, 2); 
        gamj[j] = std::exp(std::log(gam) * intervals[j-1]);  
    }
    gamj[0] = 0;
}
//--------------------------------------------------------------------------

void getkapj (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ,
        std::vector<double> &kapj) 
{   
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

void getbeta0 (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ,
        std::vector<double> &beta) 
{
    // column 3 
    int j;
    double sumbeta = 0;
    for (j = 1; j < jj; j++) {
        beta[j] = openval(PIAJ[i3(n, j, x, nc, jj )]-1,2); 
        sumbeta += std::exp(beta[j]);
    }
    beta[0] = 1;
    for (j = 1; j < jj; j++) {
        beta[j] = std::exp(beta[j]) / (1 + sumbeta);
        beta[0] -= beta[j];
    }
}
//--------------------------------------------------------------------------

void gettau (
        const int n, 
        const int x, 
        const int nc, 
        const int jj,
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ,
        std::vector<double> &tau,
        const int M) 
{
    // column 5 
    int j;
    double sumtau = 0;
    for (j = 0; j < M; j++) {
        tau[j] = openval(PIAJ[i3(n, j, x, nc, jj )]-1, 4);
        sumtau += std::exp(tau[j]);
    }
    tau[M] = 1;
    for (j = 0; j < M; j++) {
        tau[j] = std::exp(tau[j]) / (1 + sumtau);
        tau[M] -= tau[j];
    }
    for (j = M+1; j < jj; j++)
        tau[j] = 0;
}
//--------------------------------------------------------------------------

void getDj (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ,
        std::vector<double> &Dj) 
{
    // column 3 
    int j;
    for (j = 0; j < jj; j++) {
        Dj[j] = openval(PIAJ[i3(n, j, x, nc, jj )]-1, 2); 
    }
}
//--------------------------------------------------------------------------

// per capita recruitment cf Link & Barker 2005, Schwarz 'Gentle Intro'
void getbetaf (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ,
        const std::vector<double> phij,
        const RcppParallel::RVector<double> intervals, 
        std::vector<double> &beta) 
{
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

void getbetal (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ,
        const std::vector<double> phij,
        const RcppParallel::RVector<double> intervals, 
        std::vector<double> &beta) 
{
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

void getbetag (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ, 
        const std::vector<double> phij,
        const RcppParallel::RVector<double> intervals, 
        std::vector<double> &beta) 
{
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

void getbetak (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ, 
        const std::vector<double> phij,
        std::vector<double> &beta) 
{
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
void getbetaB (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ, 
        std::vector<double> &beta) 
{
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

void getbetaD (
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ, 
        const std::vector<double> phij,
        std::vector<double> &beta) 
{
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

void getbeta (
        const int type, 
        const int n, 
        const int x, 
        const int nc, 
        const int jj, 
        const RcppParallel::RMatrix<double> openval,  
        const RcppParallel::RVector<int> PIAJ, 
        const RcppParallel::RVector<double> intervals,
        const std::vector<double> phij,
        std::vector<double> &beta) 
{
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
    else ; // Rcpp::stop("no beta for this type");
}
//--------------------------------------------------------------------------
