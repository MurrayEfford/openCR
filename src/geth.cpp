#include "utils.h"   // i3, i4

// [[Rcpp::export]]
Rcpp::List gethcpp (
        const int nc1, 
        const int cc, 
        const int nmix, 
        const int nk, 
        const int ss, 
        const int mm, 
        const Rcpp::IntegerVector PIA, 
        const Rcpp::NumericMatrix Tsk, 
        const Rcpp::NumericVector hk) {
    
    // This function fills a vector h representing a 4-D (x,m,n,s) array with
    // the total hazard (summed across traps) for animal n on occasion s 
    // wrt mask point m and latent class x
    
    // Computation is limited to combinations of n, s with unique parameter combinations 
    // (values in PIA and Tsk) and the returned n x s matrix 'hindex' contains the index for
    // each n, s to the unique total in h (for given x, m).
    
    int c,i,m,n,k,x,gi,hi,s;
    double Tski;          
    
    Rcpp::NumericMatrix xmat (nc1*ss, nk*(nmix+1));
    for (n=0; n<nc1; n++) {
        for (s=0; s<ss; s++) {
            for (k=0; k<nk; k++) {                         
                for(x=0; x<nmix; x++) {
                    xmat(nc1*s+n, x*nk + k) = PIA[i4(n,s,k,x, nc1, ss, nk)];
                }
                xmat(nc1*s+n, nmix*nk + k) = Tsk(k,s);		
            }
        }
    }
    Rcpp::List lookup = makelookupcpp(xmat);
    Rcpp::NumericMatrix ymat = Rcpp::as<Rcpp::NumericMatrix>(lookup["lookup"]);
    Rcpp::IntegerVector index = Rcpp::as<Rcpp::IntegerVector>(lookup["index"]);
    Rcpp::IntegerMatrix hindex(nc1, ss);
    int uniquerows = ymat.nrow();
    for (n=0; n<nc1; n++) {
        for (s=0; s<ss; s++) {
            hindex(n,s) = index[nc1*s+n]-1;
        }
    }
    
    int hlength = uniquerows * mm * nmix;
    Rcpp::NumericVector h(hlength);
    for (i=0; i<hlength; i++) h[i] = 0;
    
    // search hindex for each row index in turn, identifying first n,s with the index
    // fill h[] for this row
    hi = 0;
    for (s=0; s < ss; s++) {    // scan by column (see geth2 in secr)
        for (n=0; n < nc1; n++) {
	    if (hindex(n,s) == hi) {
                for (k=0; k < nk; k++) {
                    Tski = Tsk(k, s);
                    for (x=0; x<nmix; x++) {
                        c = PIA[i4(n,s,k,x, nc1, ss, nk)]-1;
                        // c<0 (PIA=0) implies detector not used on this occasion
                        if (c >= 0) {
                            for (m=0; m<mm; m++) {
                                gi = i3(c,k,m,cc,nk);
                                h[i3(x,m,hi, nmix, mm)] += Tski * hk[gi];
                            }
                        }
                    }
                }
                hi++;
            } 
            if (hi >= uniquerows) break;
        }
        if (hi >= uniquerows) break;
    }
    
    return Rcpp::List::create(Rcpp::Named("h") = h,
                        Rcpp::Named("hindex") = hindex);
    
}
//==============================================================================
