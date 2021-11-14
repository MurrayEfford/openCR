################################################################################
## package 'openCR'
## openCR.pdot.R
## 2021-04-20 stratified
################################################################################

openCR.pdot <- function (object, bysession = FALSE, stratum = 1) {

    # Return pdot for open population model
    beta     <- object$fit$par
    parindx  <- object$parindx
    link     <- object$link
    fixed    <- object$fixed
    design0  <- object$design    # object$design0  stopgap
    details  <- object$details
    type     <- object$type
    
    if (!is.null(object$stratified) && object$stratified) {
        capthist <- object$capthist[[stratum]]
    }
    else {
        capthist <- object$capthist
    }
    
    primaryintervals <- primaryintervals(object)[[stratum]]
    cumss    <- getcumss(capthist)
    PIA0  <- design0$PIA [stratum,,,,, drop = FALSE]
    PIAJ0 <- design0$PIAJ[stratum,,,, drop = FALSE]
    
    #--------------------------------------------------------------------
    # Fixed beta
    fb <- details$fixedbeta
    if (!is.null(fb)) {
        fb[is.na(fb)] <- beta
        beta <- fb    ## complete
    }
    #--------------------------------------------------------------------
    # Real parameters
    realparval0 <- makerealparameters (design0, beta, parindx, link, fixed)
    nc  <- nrow(capthist)
    J <- length(primaryintervals) + 1
    
    # type <- switch(type, CJS = 1, JSSAb = 2, JSSAl = 3, JSSAf = 4, JSSAg = 22, JSSAgCL = 23,
    #     JSSAfCL = 15, JSSAlCL = 16, JSSAbCL = 17, JSSAB = 18, JSSAN = 19, JSSARET = 21,
    #     Pradel = 20, Pradelg = 26, JSSAk = 28, JSSAkCL = 29)
    # Use central coding from utility.R 2020-11-01
    type <- typecode(type)
    if (type<0) stop ("model type ", type, " not recognised")
    
    distrib <- switch (object$distribution, poisson = 0, binomial = 1)
    binomN <- details$binomN
    
    ## mixture proportions
    if (details$nmix > 1) {
        pmix <- fillpmix2(nc, details$nmix, PIA0, realparval0)
    }
    else {
        pmix <- matrix(1, nrow = details$nmix, ncol = nc)
    }
    onep <- function (x) {
        if (details$R) {
            pch1 <- PCH1(
                type,
                x,
                nc,
                cumss,
                details$nmix,
                realparval0,
                PIA0,
                PIAJ0,
                primaryintervals)
        }
        else {
            pch1 <-  PCH1cpp(
            as.integer(type),
            as.integer(x-1),
            as.integer(nc),
            as.integer(J),
            as.integer(cumss),
            as.integer(details$nmix),
            as.matrix(realparval0),
            as.integer(PIA0),
            as.integer(PIAJ0),
            as.double(primaryintervals))
        }
        pmix[x,] * pch1
    }
    onepbysession <- function (x) {
        one <- function(n) {    
            p <- getp (n, x, realparval0, PIA0)
            sessp <- function (j) {
                s <- (cumss[j]+1) : cumss[j+1]
                1-prod(1 - p[s])       ## Pr detected
            }
            sapply(1:J, sessp)
        }
        p1 <- sapply(1:nc, one)   # or rep(one(1), nc) if all the same
        sweep(p1, MARGIN=2, STATS=pmix[x,], FUN="*")
    }
    if (bysession) {
        p <- sapply(1:details$nmix, onepbysession, simplify = FALSE)
        apply(abind(p, along=3), 1:2, sum)   ## session * ch, summed over mixture classes
    }
    else {
        p <- sapply(1:details$nmix, onep)
        apply(p,1,sum)
    }
}
############################################################################################


