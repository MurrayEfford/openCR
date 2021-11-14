dbinomraw <- function(x, n, p, q)
{
    lgamma(n+1) - lgamma (x+1) - lgamma(n-x+1) + x*log(p) + (n-x) * log(q)
}

kappaloglik <- function (type, openval,  PIA, PIAJ, stratum, nmix, CJSp1, distrib) {
    
    onehistory <- function (n) {
        prw <- prwi(
            1,   # type = 'CJS'
            1,   # n
            1,
            stratum$J,
            stratum$cumss,
            nmix,
            stratum$capthist[n,, drop = FALSE],
            stratum$fi[n],
            stratum$li[n],
            openval,
            PIA[1,n,,,,drop = FALSE],
            PIAJ[1,n,,,drop = FALSE],
            stratum$primaryintervals,
            CJSp1)
        if (prw<=0) {
            -1e10
        }
        else
            freq[n] * log(prw)
    }
    
    J <- length(stratum$primaryintervals) + 1
    ch <- stratum$capthist
    w <- matrix(stratum$JScounts, nrow=J)
    ni <- w[,1]       # number viewed at i
    u <- ni - w[,3]   # number viewed for first time at i
    n <- sum(u)
    p <- getpj   (1, 1, openval, PIAJ)
    phij <- getphij (1, 1, openval, PIAJ, stratum$primaryintervals)
    kapj <- getkapj (1, 1, openval, PIAJ)
    beta <- getbetak (1, 1, openval, PIAJ, phij, stratum$primaryintervals)  ## needed for beta0
    xi <- kapj/sum(kapj)
    ## f0
    if (type == 28) {
        ## superN <- openval[PIAJ[1,1,1,1], 4]
        pi <- beta[1] * p[1] * sum(kapj)
        ## distribution 0 = Poisson, 1 = binomial
        if (distrib == 1) {
            superN <- openval[PIAJ[1,1,1,1], 4] + n
            ## superN <- openval[PIAJ[1,1,1], 4]
            if (n < superN) {
                f0 <- dbinomraw (n, superN, pi, 1-pi)
                if (is.na(f0)) browser()
            }
            else {
                 f0 <- -1e10
            }
        }
        else {
            superN <- openval[PIAJ[1,1,1,1], 4]
            f0 <- dpois (n, pi * superN, log = TRUE)
        }
    }
    else f0 <- 0   ## CL
    
    ## f1
    f1 <- lgamma(n+1) - sum(lgamma(u+1)) + sum(u * log(xi))

    ## f2
    f2 <- 0   ## ignore losses
    
    ## f3
    freq <- covariates(stratum$capthist)$freq
    f3 <- sum(sapply(1:stratum$nc, onehistory))
    return (c(f0, f1, f2, f3))   # return log likelihood components
    
}
