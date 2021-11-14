###############################################################################
## openCR
## pch0.R
## 2018-02-26 openCR 1.0.0
###############################################################################

##-----------------------------------------------------------
## return JSSA probability animal n detected at least once   
## Pledger et al. 2010 Eqn (3) (mixtures outside)            
##-----------------------------------------------------------
PCH1 <- function (type, x, nc, cumss, nmix, openval0, PIA0, PIAJ, intervals) {
    J <- length(intervals)+1
    one <- function(n) {    
        p <- getp (n, x, openval0, PIA0)
        phij <- getphij (n, x, openval0, PIAJ, intervals)
        beta <- getbeta (type, n, x, openval0, PIAJ, intervals, phij)
        pdt <- 0
        for (b in 1:J) {
            for (d in b:J) {
                pbd <- beta[b]                      ## entered at b
                if (d>b)
                    pbd <- pbd * prod (phij[b:(d-1)]) ## survived
                pbd <- pbd * (1-phij[d]);             ## departed at d
                ptmp <- 1
                for (j in b:d) {
                    s <- (cumss[j]+1) : cumss[j+1]
                    ptmp <- ptmp * prod(1 - p[s])       ## not detected
                }
                pdt = pdt + pbd * (1 - ptmp)
            }
        }
        pdt
    }
    sapply(1:nc, one)   # or rep(one(1), nc) if all the same
}

#==============================================================================
# 
# prepare matrix n x j x m of session-specific Pr(omega_i = 0)
pr0njmx <- function (n, x, cumss, jj, mm, binomN, PIA0, gk0, Tsk) {
    pjm <- array(1, dim = c(jj, mm))
    kk <- dim(Tsk)[1]
    for (j in 1:jj) {
        # s is vector of indices to secondary sessions in this primary session
        s <- (cumss[j]+1) : cumss[j+1]
        S <- length(s)
        csk <- PIA0[1,n,s, ,x, drop = FALSE]
        cski <- rep(as.numeric(csk),mm)
        i <- cbind(cski, rep(rep(1:kk, each = S), mm), rep(1:mm, each = S*kk))
        gsk <- array(0, dim=c(S, kk, mm))
        # warning("bug in pr0njmx: indexing assumes gk0 is array but it isn't")
        gsk[cski>0] <- gk0[i]
        size <- t(Tsk[,s])      
        pjm[j, ] <- if (binomN %in% c(-2,0)) {
            ## for Poisson and multcatch detectors assume gk0 is hazard not probability
            apply(gsk, 3, function(x) exp(-sum(size * x)))  
        }
        else if (all(size==1)) {
            apply(1-gsk,3, prod)
        }
        else {
            apply(1-gsk, 3, function(x) prod (x^size))  ## Binomial or Bernoulli
        }
    }
    pjm
}

PCH1secr <- function (type, individual, x, nc, jj, cumss, kk, mm, openval0, PIA0, PIAJ, 
    gk0, binomN, Tsk, intervals, moveargsi, movementcode, sparsekernel,
    edgecode, usermodel, kernel, mqarray, cellsize, r0) {
    One <- function (n) {
        ## precompute this animal's session-specific Pr for mask points
        pjm <- pr0njmx(n, x, cumss, jj, mm, binomN, PIA0, gk0, Tsk)
        phij <- getphij (n, x, openval0, PIAJ, intervals)
        beta <- getbeta (type, n, x, openval0, PIAJ, intervals, phij)
        if (movementcode>1) {
            moveargsi <- pmax(moveargsi,0)
            moveargs <- getmoveargs (n, x, openval0, PIAJ, intervals, moveargsi)
            kernelp <- fillkernelC ( jj, movementcode-2, sparsekernel, kernel, 
                usermodel, cellsize, r0, moveargsi, moveargs, normalize = TRUE)
        }
        pdt <- 0
        for (b in 1:jj) {
            for (d in b:jj) {
                pbd <- beta[b]                        ## entered at b
                if (d>b)
                    pbd <- pbd * prod (phij[b:(d-1)]) ## survived
                pbd <- pbd * (1-phij[d]);             ## departed at d
                
                # static home ranges: take sum over M of product over J
                if (movementcode==0) {
                    prodpj <- apply(pjm[b:d,, drop = FALSE], 2, prod)
                    prw0 <- sum(prodpj) / mm
                }
                else if (movementcode==1) {
                    ## over primary sessions in which may have been alive
                    ## centers allowed to differ between primary sessions
                    ## mobile home ranges: take product over J of sum over M
                    sumpm <- apply(pjm[,,drop=FALSE], 1, sum )
                    prw0 <- prod(sumpm[b:d]) / mm
                }
                else {  ## normal, exponential etc.
                    pm <-  pjm[b, ] / mm
            
                    if (d>b) {
                        for (j in (b+1):d) {
                            pm <- convolvemq(j-1, kernelp, edgecode, mqarray, pm)
                            pm <- pm * pjm[j, ]
                        }
                    }
                    prw0 <- sum(pm)
                }
                pdt <- pdt + pbd * (1 - prw0)     ## sum over b,d
            }
        }
        pdt
    }
    if (individual) {
        sapply(1:nc, One)
    }
    else {
        rep(One(1), nc)
    }
}

