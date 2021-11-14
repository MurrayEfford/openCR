## R functions to avoid C for prwi (nonspatial)

## prwi

## 2018-02-12, 2019-05-19, 2019-06-19 (moved from prwisecr.R)
## 2020-12-12 see prwi CJS variants.R for experimental MTE code CJSmte etc.

## openval has
## column 1 p 
## column 2 phi
## column 3 f, g, l  (recruitment parameter)
## column 4 pmix (non-spatial)


################################################################################

prwi <- function (type, n, x, jj, cumss, nmix, w, fi, li, openval, PIA, PIAJ,
    intervals, CJSp1) {
    # get session-specific real parameter values
    p <- getp (n, x, openval, PIA)
    phij <- getphij (n, x, openval, PIAJ, intervals)
    if (type == 1) {
        minb <- fi[n]
        cjs <- 1-CJSp1
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
    pstar <- tapply(p, rep(1:jj, diff(cumss)), prod)
    
    pdt <- 0
    # loop over possible birth and death times
    for (b in minb:maxb) {
        for (d in mind:maxd) {
            # type 1 = CJS is conditional on release at b
            pbd <- if (type==1) 1 else beta[b]  
            if (b<d) 
                pbd <- pbd * prod(phij[b:(d-1)])
            if ((li[n]>0) & (d<jj))    # not censored
                pbd <- pbd * (1-phij[d])
            # pbd now accounts for birth at b and survival to d
            # next multiply by conditional probability of observed CH
            pjt <- 1
            if ((b+cjs) <= d) {
                for (j in (b+cjs):d) {
                    # p[s] = detection probability for each secondary session
                    # in primary session j
                    s <- (cumss[j]+1) : cumss[j+1]
                    counts <- abs(w[n, s])
                    pjt <- pjt * prod(ifelse(counts>0, p[s], 1-p[s]))
                }
            }
            pdt <- pdt + pbd * pjt
        }
    }
    pdt
}
################################################################################

