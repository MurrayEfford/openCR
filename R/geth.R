gethR <- function (mm, PIA, Tsk, hk) {
    
    ## This function fills a vector h representing a 4-D (x,m,n,s) array with
    ## the total hazard (summed across traps) for animal n on occasion s 
    ## wrt mask point m and latent class x
    
    ## In many cases the same value of total hazard applies across multiple combinations of n and s.
    ## Then the computation is limited to combinations of n, s with unique parameter combinations 
    ## (values in PIA) and the returned matrix 'hindex' contains the index for each n, s to 
    ## the unique total in h (for given x, m).
    
    ## mixtures are group-specific for full likelihood, and     
    ## individual-specific for conditional likelihood           
    
    nk <- nrow(Tsk)
    ss <- ncol(Tsk)
    nmix <- dim(PIA)[5]
    nc1 <- dim(PIA)[2]
    cc <- max(PIA)
    
    xmat <- matrix(0, nc1*ss, nk*(nmix+1))
    for (n in 1:nc1) {
        for (s in 1:ss) {
            for(x in 1:nmix) {
                for (k in 1:nk) {                         
                    xmat[nc1*(s-1)+n, (x-1)*nk + k] = PIA[1,n,s,k,x]
                    xmat[nc1*(s-1)+n, nmix*nk + k] = Tsk[k,s]
                }
            }
        }
    }
    
    lookup <- makelookupcpp(xmat)
    index <- lookup$index
    uniquerows <- max(index)
    hindex <- matrix(0, nc1, ss)
    for (n in 1:nc1) {
        for (s in 1:ss) {
            hindex[n,s] <- index[nc1*(s-1) + n]
        }
    }
    
    ## search hindex for each row index in turn, identifying first n,s with the index
    ## fill h[] for this row
    h <- array(0, dim = c(nmix, mm, uniquerows))
    hi <- 1
    for (n in 1:nc1) {
        for (s in 1:ss) {
            if (hindex[n,s] == hi) {
                for (k in 1:nk) {
                    Tski <- Tsk[k,s]
                    for (x in 1:nmix) {
                        c <- PIA[1,n,s,k,x];
                        ## c==0 (PIA=0) implies detector not used on this occasion
                        if (c > 0) {
                            for (m in 1:mm) {
                                h[x, m, hi] <- h[x, m, hi] + Tski * hk[c,k,m]
                            }
                        }
                    }
                }
                hi <- hi+1
            } 
            if (hi > uniquerows) break
        }
        if (hi> uniquerows) break
    }
    
    list(h = h, hindex = hindex-1)
    
}
