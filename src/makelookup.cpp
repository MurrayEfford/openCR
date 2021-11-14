#include "utils.h"

//==============================================================================

// [[Rcpp::export]]
Rcpp::List makelookupcpp (const Rcpp::NumericMatrix x)  
    
{
    // Create lookup table to the unique rows in a matrix
    // Only the first 'unique' rows of y contain valid data on exit
    // The indices are 1..length(unique(x))=nrow(y')
    // MGE 2018-11-05
    
    int i;
    int j;
    int k;
    int dupl = 0;
    int unique = 0;
    int nrow = x.nrow();
    int ncol = x.ncol();

    // outputs
    Rcpp::NumericMatrix y (nrow, ncol);
    Rcpp::IntegerVector index (nrow);   // output lookup rows of x in y
    int resultcode = 1;
    
    // Avoid sort for now as it's complex to keep order of first occurrence, not needed
    // scan for unique rows of x, copying to y
    // assign unique index to original rows as we go
    
    for (j=0; j<ncol; j++)
	y(0,j) = x(0,j);    // first row 
    index[0] = 1;
    unique = 0;
    
    // Loop over rows of input matrix
    for (i = 1; i < nrow; i++) {
        // Is this row unique? Compare with each previous unique in turn 
        for (k=0; k <= unique; k++) {
            dupl = 1;
            for (j = 0; j < ncol; j++) {
                if (x(i,j) != y(k,j))
		{dupl=0; break;}
	    }
	    if (dupl==1) break;  // found previous instance 
	}
	if (dupl==0) { // add unique row
	    unique ++;
	    k = unique;
	    for (j=0; j<ncol; j++)
	        y(unique,j) = x(i,j);
	}
	index[i] = k+1;
    }
    
    y = y(Rcpp::Range(0,unique), Rcpp::Range(0,ncol-1));
    colnames(y) = colnames(x);
    
    return Rcpp::List::create(
        Rcpp::Named("resultcode") = resultcode,
        Rcpp::Named("lookup") = y,
        Rcpp::Named("index") = index);
        
}
//==============================================================================
