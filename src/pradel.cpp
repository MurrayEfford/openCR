#include "utils.h"

//==============================================================================

// Pradel model
// NOT ROBUST DESIGN 
// assumes all primary sessions comprise one secondary session so S=J 

//==============================================================================
int sumj (std::vector<int> &uv, int j, int k) {
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
//--------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::NumericVector pradelloglikcpp (
        const int    type,        // 20 = Pradel, 26 = Pradelg 
        const Rcpp::IntegerVector w,           // counts (1:jj, 1:5) 
        const int    nc,          // number of capture histories 
        const int    jj,          // number of sessions 
        const int    nmix,        // number of mixtures 
        const Rcpp::NumericMatrix openval,     // Parameter values - turnover 
        const Rcpp::IntegerVector PIAJ,         // lookup which parameter combination to use n, J, mix 
        const Rcpp::NumericVector intervals    // vector of jj-1 between-occasion intervals 
)
    
{ 
    const RcppParallel::RMatrix<double> openvalR(openval); /// Parameter values - turnover 
    const RcppParallel::RVector<int> PIAJR(PIAJ);          // lookup which parameter combination to use n, J, mix 
    const RcppParallel::RVector<double> intervalsR(intervals);    // vector of jj-1 between-occasion intervals 

    int    j,k;          // indices &  miscellaneous 
    int    x = 0;
    double prd;
    double BL = 0;
    Rcpp::NumericVector value (2, 0.0);   
    
    // int    c, wxi;
    // double pmix[3];             latent class membership probability 
    // int gpar = 3;
     
    std::vector<double> p(jj);
    std::vector<double> phij(jj);
    std::vector<double> gamj(jj);
    std::vector<double> mu(jj);
    std::vector<double> chi(jj);
    std::vector<double> xi(jj);
    std::vector<int> ni(jj);
    std::vector<int> u(jj);
    std::vector<int> v(jj);
    std::vector<int> d(jj);

    for (j=0; j< jj; j++) {
        ni[j] = w[j];
        u[j] = ni[j] - w[jj * 2 + j];
        v[j] = ni[j] - w[jj * 3 + j];
        d[j] = ni[j] - w[jj + j];
        mu[j] = (double) w[jj + j] / ni[j];
    }
    mu[jj-1] = 1.0;
    
    //---------------------------------------------------------
    // mixture proportions                                     
    
     // gpar = 3;
     // for (i=0; i < nmix; i++) pmix[i] = 1;  default 
     // if (nmix>1) {
     // // one extra real parameter 
     // gpar++;
     // for (x=0; x<nmix; x++) {
     // wxi = i3(0,0,x,1,jj);
     // c = PIAJ[wxi] - 1;
     // pmix[x] = openval(c, gpar-1); 
     // }

    //---------------------------------------------------------
    
    // MIXTURES NOT USED YET 2011-12-15 x defaults to 0 
    getp   (0, x, nc, jj, openvalR, PIAJR, p);
    getphij (0, x, nc, jj, openvalR, PIAJR, intervalsR, phij);
 
    if (type==20)
	getgaml (0, x, nc, jj, openvalR, PIAJR, intervalsR, gamj);
    else
	getgamj (0, x, nc, jj, openvalR, PIAJR, intervalsR, gamj);
    
    chi[jj-1] = 1;
    for (j = jj-2; j>=0; j--) {
        chi[j] = 1- phij[j] * (1 - (1-p[j+1]) * chi[j+1]);
    }
    
    xi[0] = 1;
    for (j = 1; j<jj; j++) {
        xi[j] = (1-gamj[j]) + gamj[j] * (1-p[j-1]) / (1-p[j-1]*(1-mu[j-1]) ) * xi[j-1];
    }

    // for (j=0; j < jj; j++) {
    // 	Rprintf("j %4d p[j] %7.6f phij[j] %7.6f gamj[j] %7.6f mu[j] %7.6f chi[j] %7.6f xi[j] %7.6f\n",
    // 		j, p[j], phij[j], gamj[j], mu[j], chi[j], xi[j]);
    // }

    for (j=0; j<jj; j++) {
        if (xi[j]>0)
            value[0] += u[j] * log (xi[j]);
        if (gamj[j]>0)
            value[0]+= sumj(u,0,j-1) * log (gamj[j]);
        if (p[j]>0)
            value[0] += ni[j] * log (p[j]);
        if (p[j]<1)
            value[0] += (sumj(u,0,j) - sumj(v,0,j-1) - ni[j]) * log (1-p[j]);
        if (j < (jj-1))
            value[0] += sumj(v,j+1,jj-1) * log (phij[j]);
        if (mu[j] < 1)
            value[0] += (ni[j] - d[j]) * log (mu[j]) + d[j] * log (1-mu[j]);
        if (chi[j]>0)
            value[0] += (v[j]-d[j]) * log (chi[j]); 
        value[0] += sumj(u,j+1,jj-1) * log (1 - p[j] * (1 - mu[j]));
    }
    
    for (j=0; j<jj; j++) {
        prd = xi[j];
        if (j > 0)
            for (k=0; k<j; k++) prd *= phij[k] * (1-p[k] * (1-mu[k]));
        if (j < jj)
            for (k=j+1; k<jj; k++) prd *= gamj[k];
        BL += prd * p[j];
    }
    
    value[1] = -sumj(u,0,jj-1) * log (BL);
    
    return value;   // return value log likelihood 
    
}
//==============================================================================
