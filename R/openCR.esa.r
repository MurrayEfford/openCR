################################################################################
## package 'openCR'
## openCR.esa.R
## 2019-05-17
## 2021-06-11 fixed bugs mqarray <- 0; makegkParalleld
################################################################################

# CJSsecr = 6
# JSSAsecrf = 7
# JSSAsecrD = 8
# JSSAsecrfCL = PLBsecrf = 9
# JSSAsecrlCL = PLBsecrl = 10
# JSSAsecrbCL = PLBsecrb = 11
# JSSAsecrl = 12
# JSSAsecrb = 13
# JSSAsecrB = 14
# JSSAsecrg = 24
# JSSAsecrgCL = PLBsecrg = 25

# secrCL = 30
# secrD = 31

openCR.esa <- function (object, bysession = FALSE, stratum = 1) {

    if (!(grepl('secr',object$type) & inherits(object,'openCR')) )
        stop ("requires fitted openCR secr model")
    # Return the esa
    beta     <- object$fit$par
    parindx  <- object$parindx
    link     <- object$link
    fixed    <- object$fixed
    details  <- object$details
    type     <- object$type
    binomN   <- object$binomN
    if (is.null(binomN)) binomN <- 1
    
    if (!is.null(object$stratified) && object$stratified) {
        capthist <- object$capthist[[stratum]]
        mask     <- object$mask[[stratum]]
    }
    else {
        capthist <- object$capthist
        mask     <- object$mask
    }
    
    cumss <- getcumss(capthist)
    PIA0  <- object$design0$PIA[stratum,,,,, drop = FALSE]
    PIAJ0 <- object$design0$PIAJ[stratum,,,, drop = FALSE]
    primaryintervals <- primaryintervals(object)[[stratum]]

    ## openCR >= 1.2.0
    individual <- object$design0$individual
    ## openCR < 1.2.0
    if (is.null(individual)) individual <- individualcovariates(PIA0)

    ###### movement kernel and related
    cellsize <- mqarray <- kernel <- 0    # until set for movement model
    movementcode <- movecode(object$movementmodel)
    sparsekernel <- object$sparsekernel
    anchored <- details$anchored
    if (is.null(anchored)) anchored <- FALSE
    r0 <- details$r0
    if (is.null(r0)) r0 <- 1/sqrt(pi)
    
    edgecode <- edgemethodcode(object$edgemethod)
    # 2021-02-21 modified for annular
    if (object$movementmodel %in% .openCRstuff$kernelmodels) {
        k2 <- object$kernelradius
        # cellsize <- attr(mask,'area')^0.5 * 100   ## metres, equal mask cellsize
        cellsize <- spacing(mask)
        kernel <- expand.grid(x = -k2:k2, y = -k2:k2)
        kernel <- kernel[(kernel$x^2 + kernel$y^2) <= (k2+0.5)^2, ]
        if (object$movementmodel == 'annular') {
            r <- (kernel$x^2 + kernel$y^2)
            kernel <- kernel[(r==0) | (r>(k2-0.5)), ]
        }
        if (object$movementmodel == 'annular2') {
            r <- (kernel$x^2 + kernel$y^2)
            origin <- r==0
            ring1 <- r > (k2/2-0.5) && r<(k2/2+0.5)
            ring2 <- r > (k2-0.5)
            kernel <- kernel[origin | ring1 | ring2, ]
        }
        kernel <- linearisekernel (kernel, mask) 
        mqarray <- mqsetup (mask, kernel, cellsize, edgecode)   
    }

    #--------------------------------------------------------------------
    # Fixed beta
    fb <- details$fixedbeta
    if (!is.null(fb)) {
        fb[is.na(fb)] <- beta
        beta <- fb    ## complete
    }
    #--------------------------------------------------------------------
    # Real parameters

    realparval0 <- makerealparameters (object$design0, beta, parindx, link, fixed)
    nc  <- nrow(capthist)
    J <- length(primaryintervals) + 1

    type <- typecode(type)  # convert to integer
    trps <- traps(capthist)
    k <- nrow(trps)
    m <- nrow(mask)
    if (!is.null(mask)) area <- attr(mask,'area')
    else area <- 0
    distrib <- switch (object$distribution, poisson = 0, binomial = 1)

    detectr <- detector(trps)[1]
    if (detectr == 'count') {
        detectr <- if (object$binomN == 0) "poissoncount" else "binomialcount"
    }

    usge <- usage(trps)
    if (is.null(usge)) usge <- matrix(1, nrow = k, ncol = ncol(capthist))
    distmat <- getdistmat(trps, mask, details$userdist, object$detectfn == 20)
    temp <- makegkParalleldcpp (as.integer(object$detectfn), 
        as.integer(.openCRstuff$sigmai[type]),
        as.integer(details$grain),
        as.integer(setNumThreads()),
        as.matrix(realparval0),
        as.matrix(distmat))
    gk0 <- array(temp[[1]], dim=c(nrow(realparval0), k, m)) # 2020-10-28 as array
    hk0 <- array(temp[[2]], dim=c(nrow(realparval0), k, m)) # 2020-10-28 as array

    settlement <- 1   # dummy value for now
    if (!is.null(details$settlemodel) && details$settlemodel) warning ("esa not ready for settlement")
    
    ## mixture proportions
    if (details$nmix > 1) {
        temp <- fillpmix2(nc, details$nmix, PIA0, realparval0)
        pmix <- matrix(temp, ncol = nc)
    }
    else {
        pmix <- matrix(1, nrow = details$nmix, ncol = nc)
    }
    onea <- function (x) {
        if (details$R) {
            pch1 <-  PCH1secr(
                type,
                as.logical(individual),
                x,
                nc,
                J,
                cumss,
                k,
                m,
                realparval0,
                PIA0,
                PIAJ0,
                if (detectr %in% c("multi", "poissoncount")) hk0 else gk0,
                binomN,
                usge,
                primaryintervals,
                object$moveargsi,
                movementcode,                     # what about edgecode? 2020-12-12
                sparsekernel, 
                object$usermodel,
                kernel,
                mqarray,
                cellsize)
        }
        else {
            pch1 <-  PCH1secrparallelcpp(
                as.integer(x-1),
                as.integer(type),
                as.integer(details$grain),
                as.integer(setNumThreads()),
                as.logical(individual),
                as.integer(J),
                as.integer(m),
                as.integer(nc),
                as.integer(cumss),
                as.matrix (realparval0),
                as.integer(PIA0),
                as.integer(PIAJ0),
                as.double (if (detectr == "poissoncount") hk0 else gk0),  ## 2019-05-08, 17
                as.integer(binomN),
                as.matrix (usge),
                as.double (primaryintervals),
                as.integer(object$moveargsi),
                as.integer(movementcode),
                as.logical(sparsekernel),
                as.logical(anchored),
                as.integer(edgecode),
                as.character(object$usermodel),
                as.matrix(kernel),
                as.matrix(mqarray),
                as.double (cellsize),
                as.double (r0),
                as.matrix (settlement)
            )
        }
        pmix[x,] * pch1
    }
    
    oneabysession <- function (x) {
        a0 <- PCH0secrjcpp (
            as.integer(type),
            as.integer(x-1),
            as.integer(nc),
            as.integer(J),
            as.integer(cumss),
            as.integer(k),
            as.integer(m),
            as.integer(nrow(realparval0)),
            as.integer(PIA0),
            as.double (gk0),
            as.integer(binomN),
            as.matrix (usge)) 
        a0 <- matrix(a0, nrow = nc, ncol = J)
        1-sweep(a0, MARGIN=1, STATS=pmix[x,], FUN="*")
    }
    if (bysession) {
        a <- sapply(1:details$nmix, oneabysession, simplify = FALSE)
        a <- t(apply(abind(a, along=3), 1:2, sum))   ## session * ch, summed over latent classes
        a <- a * area * m
        sess <- primarysessions(intervals(capthist)) ## differs from 'primaryintervals'
        OK <- apply(capthist, 1, by, sess, sum)>0
        freq <- sweep(OK, MARGIN=2, STATS=covariates(capthist)$freq, FUN = "*")
        a <- lapply(1:J, function(j) rep(a[j,], freq[j,]))
    }
    else {
        a <- sapply(1:details$nmix, onea)
        a <- apply(a, 1, sum)    ## sum over latent classes
        a <- rep(a, covariates(capthist)$freq)
        a <- a * area * m
    }
    ## adjusts for covariates(object$capthist)$freq - no need to do this downstream
    a 
}
