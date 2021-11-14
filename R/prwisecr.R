## R functions to avoid C for prwi (spatial)

## prwisecr

## 2018-02-12, 2019-05-19,
## 2019-06-19 removed nonspatial to prwi.R
## 2019-06-19 merged multicatch into general prwisecr
## 2020-12-12 moved mqsetup to utility.R

## openval has
## column 1 lambda0
## column 2 phi
## column 3 f, g, l  (recruitment parameter)
## column 4 sigma
## column 5 pmix

#########################################################################################

# Probability of count for session s, detector k, animal i
# The argument 'g' is understood to be a cumulative hazard if binomN=0,
# a probability otherwise

pski <- function(binomN, count, Tski, g) {

    result <- 1.0

    if (binomN == -1) {                              ## binary proximity detectors : Bernoulli
        if (any(abs(Tski-1) > 1e-10)) {              ## effort not unity; adjust g
            g <- 1 - (1 - g)^Tski
        }
        if (count>0)
            result <- g
        else
            result <- 1 - g
    }
    else if (binomN == 0) {                          ## count detectors : Poisson
        if (count == 0)
            result <- exp(-Tski * g)                 ## routinely apply Tsk adjustment to cum. hazard
        else
            result <- dpois(count, Tski * g, FALSE)
    }
    else if (binomN == 1) {                          ## count detectors : Binomial, size from Tsk
        result <- dbinom (count, round(Tski), g, FALSE)
    }
    else if (binomN > 1) {                           ## count detectors : Binomial, specified size
        if (abs(Tski-1) > 1e-10) {                   ## effort not unity, adjust g
            g <- 1 - (1 - g)^Tski
        }
        result <- dbinom (count, binomN, g, FALSE)
    }
    else stop("binomN < -1 not allowed")
    result
}

#########################################################################################

convolvemq <- function (
    j,         ## session number 1..jj
    kernelp,   ## p(move|dx,dy) for points in kernel
    edgecode,  ## 0 = no action, 1 = no action, already wrapped, 2 = truncate
    mqarray,   ## input
    pjm
)
{
    mm <- nrow(mqarray)
    kn <- ncol(mqarray)
    workpjm <- numeric(mm)
    
    ## convolve movement kernel and pjm...
    for (m in 1:mm) {
        if (edgecode == 2) {
            ## warning: unresolved problem with some kernel/mask comb
            q <- mqarray[m,] + 1
            q <- q[q>0]   # vector of indices to inside landing points
            sump <- sum(kernelp[kn * (j-1) + q], na.rm = T)
        }
        else {
            sump <- 1.0
        }
        if (is.na(sump)) browser()   ## debug 2020-12-13
        if (sump>0) {
            for (q in 1:kn) {          ## over movement kernel
                mq <- mqarray[m,q] + 1 ## which new point corresponds to kernel point q relative to current mask point m
                if (mq > 0) {          ## post-dispersal site is within mask
                    if (mq>mm) stop("mq > mm")
                    workpjm[mq] <- workpjm[mq] + pjm[m] * kernelp[q,j]   ## probability of this move CHECK DIM KERNELP
                }
            }
        }
    }
    workpjm
}

###################################################################################
prw <- function (n, j, x, kk, binomN, cumss, w, PIA, Tsk, gk, h, p0, hindex, pjm) {
    # gk assumed to be hazard if multi (binomN=-2) or Poisson count (binomN=0)
    dead <- FALSE
    for (s in (cumss[j]+1):cumss[j+1]) {

        if (binomN == -2) {      ## multi-catch traps, 2-D w
            wi <- w[n, s]
            if (wi < 0) dead <- TRUE
            k <- abs(wi)         ## trap number 1..kk k = 0 if not caught

            ## Not captured in any trap on occasion s
            if (k < 1) {
                OK <- h[x,,hindex[n,s]] > 1e-8
                pjm[OK] <- pjm[OK] * p0[x,OK,hindex[n,s]]
            }
            ## Captured in trap k
            else {
                c <- PIA[1,n, s, k, x]
                if (c >= 1) {    ## drops unset traps
                    pjm <- pjm * Tsk[k,s] * (1-p0[x,,hindex[n,s]]) *  gk[c, k, ] / h[x,,hindex[n,s]]
                }
            }
        }
        else {
            for (k in 1:kk) {
                c <- PIA[1,n,s,k,x]
                if (c >= 1) {    # drops unset traps
                    count <- w[n,s,k]
                    if (count<0) {count <- -count; dead <- TRUE }
                    pjm <- pjm * pski(binomN,count,Tsk[k,s], gk[c,k,])
                }
            }
        }

        if (dead) break;   # after processing all traps on this occasion
    }
    ## message("j ", j, " pjm ", sum(pjm))
    pjm
}
###################################################################################

prwisecr <- function (type, n, x, nc, jj, kk, mm, nmix, cumss, w, fi, li, gk,
    openval, PIA, PIAJ, binomN, Tsk, intervals, h, hindex,
    CJSp1, moveargsi, movementcode, sparsekernel, edgecode,
    usermodel, kernel = NULL, mqarray = NULL, cellsize = NULL, 
    r0) {

    ## precompute p0 to save time (multicatch only)
    p0 <- if (binomN == -2) exp(-h) else 1
    phij <- getphij (n, x, openval, PIAJ, intervals)
    if (movementcode > 1) {
        moveargsi <- pmax(moveargsi,0)
        moveargs <- getmoveargs (n, x, openval, PIAJ, intervals, moveargsi)
        
        kernelp <- fillkernelC ( jj, movementcode-2, sparsekernel, kernel, 
            usermodel, cellsize, r0, moveargsi, moveargs, normalize = TRUE)
        
    }
    if(type==6) {
        minb <- fi[n]
        cjs <- 1 - CJSp1
    }
    else {
        minb <- 1
        cjs <- 0
        beta <- getbeta (type, n, x, openval, PIAJ, intervals, phij)
    }
    maxb <- fi[n]
    mind <- abs(li[n])
    maxd <- jj
    if (li[n] < 0) maxd <- mind     # possible censoring
    pdt <- 0
    pdotbd <- 1.0
    for (b in minb:maxb) {
        for (d in mind:maxd) {
            if (type==6) {     ## CJS
                pbd <- 1
            }
            else {
                pbd <- beta[b]
                pdotbd <- 1
            }
            if (b<d) pbd <- pbd * prod(phij[b:(d-1)])
            if ((li[n]>0) & (d<jj))    # if not censored, died
                pbd <- pbd * (1-phij[d])

            prwi <- 1.0
            if (d >= (b+cjs)) {
                if (movementcode == 0) {

                    alpha <- rep(1.0/mm,mm)
                    for (j in (b+cjs):d) {
                        alpha <- prw(n, j, x, kk, binomN, cumss, w, PIA, Tsk, gk, h, p0, hindex, alpha)
                    }
                    prwi <- sum(alpha)
                }
                else if ( movementcode == 1) { # uncorrelated; product over primary sessions
                    prwi <- 1.0
                    for (j in (b+cjs):d) {
                        alpha <- rep(1.0/mm,mm)
                        alpha <- prw(n, j, x, kk, binomN, cumss, w, PIA, Tsk, gk, h, p0, hindex, alpha)
                        prwi <- prwi * sum(alpha)
                    }
                }
                else { # movementcode>1
                    alpha <- rep(1.0/mm, mm)
                    alpha <- prw (n, b+cjs, x, kk, binomN, cumss, w, PIA, Tsk, gk, h, p0, hindex, alpha)
                    if (d>(b+cjs)) {
                        for (j in (b+cjs+1):d) {
                            alpha <- convolvemq(j-1, kernelp, edgecode, mqarray, alpha)
                            alpha <- prw(n, j, x, kk, binomN, cumss, w, PIA, Tsk, gk, h, p0, hindex, alpha)
                        }
                    }
                    prwi <- sum(alpha)
                }
                # message("n ", n, " b ", b, " d ", d, " pbd ", pbd, " prwi ", prwi, " pdotbd ", pdotbd)
            }
            pdt <- pdt + pbd * prwi / pdotbd
        }
    }
    # message("n ", n, " pdt ", pdt)
    pdt
}

#==============================================================================

