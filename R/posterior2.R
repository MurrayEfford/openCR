###############################################################################
# posterior2.R
## 2018-11-19 draft of classMembership as standalone function
##            cf posterior.allocation previously called in openCR.fit
## 2019-05-8 UNTESTED
## 2020-11-07 failed due to bug in secr::alive, now fixed
###############################################################################

matchch <- function (x, sq) {
    df <- as.data.frame(matrix(x, nrow = nrow(x)))
    if (!is.null(covariates(x))) {
        if (nrow(covariates(x))>0)
            df <- cbind(df, covariates(x))
    }
    
    dfsq <- as.data.frame(matrix(sq, nrow = nrow(sq)))
    if (!is.null(covariates(sq))) {
        if (nrow(covariates(sq))>0)
            dfsq <- cbind(dfsq, covariates(sq))
    }
    df$freq <- NULL
    dfsq$freq <- NULL
    
    if (!all(colnames(df) == colnames(dfsq)))
        stop ("sq, x do not match")
    cx <- do.call("paste", c(df[, , drop = FALSE], sep = "\r"))
    csq <- do.call("paste", c(dfsq[, , drop = FALSE], sep = "\r"))
    match(cx, csq, nomatch = NA)
}

classMembership <- function (object, ...) UseMethod("classMembership")

classMembership.openCR <- function (object, fullCH = NULL, ...) {
    
    if (!is.null(object$stratified) && object$stratified) stop ("classMembership not ready for stratified data")
    # Return the probability of membership in latent classes for h2 and h3 models
    
    nmix <- object$details$nmix
    capthist <- object$capthist   # squeezed!
    if (nmix == 1) return (rep(1,nrow(capthist)))
    secr <- grepl('secr', object$type)
    
    #--------------------------------------------------------------------
    # Real parameters
    realparval  <- makerealparameters (object$design, complete.beta(object),
                                       object$parindx, object$link, object$fixed)
    # Parameter Index Array
    PIA <- object$design$PIA
    PIAJ <- object$design$PIAJ
    #--------------------------------------------------------------------
    cumss <- getcumss(capthist)
    J <- length(cumss)-1
    nc <- nrow(capthist)   ## beware freq!
    intervals <- intervals(capthist)
    primarysession <- primarysessions(intervals)
    
    pmix <- fillpmix2(nc, nmix, PIA, realparval)
    
    ##--------------------------------------------------------------
    ## get fi,li, re-form capthist as CH
    lost <- which(apply(capthist,1,min, drop = FALSE)<0)
    twoD <- apply(abs(capthist), 1:2, sum, drop = FALSE)
    CH <- twoD
    if (J==1)
        twoD <- as.matrix(apply(twoD, 1, function(x) tapply(x,primarysession,max)))
    else
        twoD <- t(apply(twoD, 1, function(x) tapply(x,primarysession,max)))  # in terms of primary sessions
    fi <- apply(twoD, 1, function(x) min(which(x>0)))
    li <- apply(twoD, 1, function(x) max(which(x>0)))
    twoD[cbind(lost, li[lost])] <- -1
    li[lost] <- -li[lost]
    covariates(CH) <- covariates(capthist)
    covariates(twoD) <- covariates(capthist)
    JScounts <- unlist(JS.counts(twoD))
    ##--------------------------------------------------------------

    px <- matrix(NA, nrow = nc, ncol = nmix)
    
    if (secr) {
        type <- typecode(object$type)
        if (!type %in% c(6:14,24,25,30,31)) stop ("Invalid likelihood type for posterior allocation")
        trps <- traps(object$capthist)
        k <- nrow(trps)
        m <- nrow(object$mask)
        
        if (!is.null(object$mask)) area <- attr(object$mask,'area')
        else area <- 0
        
        ## use same binomN as logliksecr 2019-05-08
        detectr <- detector(trps)[1]
        if (detectr == 'count') {
            detectr <- if (object$binomN == 0) "poissoncount" else "binomialcount"
        }
        binomN <- switch (detectr, multi = -2, proximity = -1, 
                          poissoncount = 0, binomialcount = object$binomN, -9)
        
        usge <- usage(traps(capthist))
        if (is.null(usge) | object$details$ignoreusage)
            usge <- matrix(1, nrow=k, ncol= cumss[J+1])  # in terms of secondary sessions
        
        ## integer code for movement model
        movementcode <- movecode(object$movementmodel)
        sparsekernel <- object$sparsekernel
        anchored <- object$details$anchored
        if (is.null(anchored)) anchored <- FALSE
        edgecode <- edgemethodcode(object$edgemethod)

        cellsize <- mqarray <- 0
        kernel <- mqarray <- matrix(0,1,2)  ## default
        if (movementcode %in% c(2:5, 7,8)) {
            ## movement kernel
            k2 <- object$kernelradius
            # cellsize <- attr(object$mask,'area')^0.5 * 100   ## metres, equal mask cellsize
            cellsize <- spacing(object$mask)
            kernel <- expand.grid(x = -k2:k2, y = -k2:k2)
            kernel <- kernel[(kernel$x^2 + kernel$y^2) <= (k2+0.5)^2, ]   ## autoclip
            if (sparsekernel) {
                ok <- kernel$x==0 | kernel$y == 0 | kernel$x == kernel$y | kernel$x == -kernel$y
                kernel <- kernel[ok,]
            }
            kernel <- linearisekernel (kernel, object$mask) 
            mqarray <- mqsetup (object$mask, kernel, cellsize, edgecode) # 2020-11-02
        }
        
        usge <- usage(traps(capthist))
        if (is.null(usge) | object$details$ignoreusage) 
            usge <- matrix(1, nrow=k, ncol= cumss[J+1])  # in terms of secondary sessions

        ## 2017-11-26 collapse data from exclusive detectors; modified 2018-01-17
        CH <- capthist
        if (detectr == 'multi') {
            CH <- abs(capthist)
            CH <- apply(CH,1:2, which.max) *  (apply(CH,1:2, max)>0)
            lost <- apply(capthist,1:2, min)<0
            CH[lost] <- -CH[lost]
            class (CH) <- 'capthist'
            traps(CH) <- traps(capthist)
        }
        ##--------------------------------------------------------------
        distmat <- getdistmat(trps, object$mask, object$details$userdist, object$detectfn==20)
        temp <- makegkParalleldcpp (as.integer(object$detectfn),
            as.integer(.openCRstuff$sigmai[type]),
            as.integer(object$details$grain),
            as.integer(setNumThreads()),
            as.matrix(realparval),
            as.matrix(distmat))
        gk <- temp[[1]]
        hk <- temp[[2]]
        pmix <- fillpmix2(nc, nmix, PIA, realparval)
        
        if (detectr=='multi') {
            haztemp <- gethcpp(
                as.integer(nc),
                as.integer(nrow(realparval)),
                as.integer(nmix),
                as.integer(k),
                as.integer(cumss[J+1]),
                as.integer(m),
                as.integer(PIA),
                as.matrix(usge),
                as.double(hk))
            haztemp$h <- array(haztemp$h, dim = c(nmix, m, max(haztemp$hindex)+1))
        }
        else {
            haztemp <- list(h = array(-1, dim=c(nmix,1,1)), hindex = matrix(-1))
        }
        
        if (type != 6) {
            for (x in 1:nmix) {
                hx <- if (detectr == "multi") matrix(haztemp$h[x,,], nrow = m) else -1 ## lookup sum_k (hazard)
                hi <- if (detectr == "multi") haztemp$hindex else -1                        ## index to hx
                px[,x] <-   pmix[x,] * allhistsecrparallelcpp(
                    as.integer(x-1),
                    as.integer(type),
                    as.integer(m),
                    as.integer(nc),
                    as.integer(binomN),
                    as.integer(object$details$CJSp1),
                    as.integer(object$details$grain),
                    as.integer(setNumThreads()),
                    as.double (object$intervals),
                    as.integer(cumss),
                    as.matrix (CH),     
                    as.integer(fi),
                    as.integer(li),
                    as.double (if (detectr %in% c("multi", "poissoncount")) hk else gk), ## precomputed probability or hazard
                    as.matrix (realparval),
                    as.integer(PIA),
                    as.integer(object$design$PIAJ),
                    as.matrix (usge),
                    as.matrix (hx),                
                    as.matrix (hi),      
                    as.integer(movementcode),
                    as.logical(sparsekernel),
                    as.logical(anchored),
                    as.integer(edgecode),
                    as.character(object$usermodel),  ## 2019-05-07
                    as.integer(object$moveargsi),
                    as.matrix (kernel),
                    as.matrix (mqarray),
                    as.double (cellsize))
            }
        }
    }
    
    else {
        ## NON-SPATIAL
        type <- typecode(object$type)
        if (!type %in% c(1:4, 15:23, 26,27)) stop ("Invalid likelihood type")
        if (type %in% c(20,26)) {
            stop("mixtures not expected in Pradel models")
        }
        if (!is.null(fullCH)) {
            fullCH <- reduce(fullCH, outputdetector = 'nonspatial', verify = FALSE, dropunused=FALSE)
            # FAILS TO COMPRESS LAST DIM (bug in reduce.capthist 2018-11-20), SO REPEAT
            fullCH <- reduce(fullCH, verify = FALSE, dropunused=FALSE)
        }
        for (x in 1:nmix) {
            px[,x] <-   pmix[x,] * allhistparallelcpp(
                as.integer(x-1),
                as.integer(type),
                as.integer(nc),
                as.integer(object$details$CJSp1),
                as.integer(object$details$grain),
                as.integer(setNumThreads()),
                as.double (object$intervals),
                as.integer(cumss),
                as.integer(CH),
                as.integer(fi),
                as.integer(li),
                as.matrix (realparval),
                as.integer(PIA),
                as.integer(object$design$PIAJ))
        }
    }
    
    px <- sweep(px, MARGIN = 1, STATS = apply(px,1,sum), FUN = "/")
    out <- data.frame(px)
    names(out) <- paste0('class', 1:nmix)
    out$maxclass <- max.col(out)
    if (is.null(fullCH)) {   ## untested
        fullCH <- unsqueeze(object$capthist)
    }
    else {
        if (ms(fullCH)) fullCH <- join(fullCH, drop.sites = !secr)
        nf <- sum(covariates(capthist)$freq)
        if (nrow(fullCH) != nf)
            stop ("fullCH has differing number of histories (", nc, " vs ", nf, ")")
    }
    i <- matchch(fullCH, capthist)
    out <- out[i,]
    rownames(out) <- rownames(fullCH)
    out
}

######################################################################################################
