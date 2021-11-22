# spatial likelihood

# 2019-04-09 split from loglik.R
# 2019-04-09 explicit treatment of count detector; dropuse of data$multi
# 2019-04-14 R option failed with movement in PCH1secr because argument usermodel omitted
# 2019-04-23 removed type 5 'secr' (unused)
# 2019-05-06 1.4.0
# 2019-06-19 onehistory modified for single call to prwisecr (merged prwimulti)
# 2020-09-01 changed return; to return(1e10) in open.secr.loglikfn component 3
# 2020-10-25 fixed bug in open.secr.loglikfn component 3
# 2021-04-19 stratified
# 2021-08-14 c and [ methods for openCRlist
# 2021-11-08 PIA, PIAJ overwritten for naive animal, rather than always copied  

# types

# CJSsecr = 6

# JSSAsecrf 7
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

#---------------------------------------------------------

open.secr.loglikfn <- function (beta, dig = 3, betaw = 8, oneeval = FALSE, data)

    # Return the negative log likelihood
    # Transformed parameter values are passed in the vector 'beta'
    # details$trace=T sends a one-line report to the screen

{
    
    onestratumll <- function(stratum) {
        
        ########################################################################
        # functions for one history at a time are used mostly for debugging
        ########################################################################
        
        onehistory <- function (n) {
            sump <- 0
            for (x in 1:nrow(pmix)) {
                temp <- prwisecr(
                    type,
                    1,   # n
                    x,
                    1,
                    stratum$J,
                    stratum$k,
                    stratum$m,
                    data$details$nmix,
                    stratum$cumss,
                    if (detectr == "multi") stratum$capthist[n,, drop = FALSE] else stratum$capthist[n,,, drop = FALSE],   ## 3-D CH
                    stratum$fi[n],
                    stratum$li[n],
                    if (detectr %in% c("poissoncount", "multi")) hk else gk, ## precomputed probability or hazard
                    realparval,
                    PIA[1,n,,,,drop = FALSE],
                    PIAJ[1,n,,,drop = FALSE],
                    binomN,
                    stratum$usge,
                    stratum$primaryintervals,
                    haztemp$h,
                    haztemp$hindex[n,,drop = FALSE]+1,
                    data$details$CJSp1,
                    data$moveargsi,
                    data$movementcode,
                    data$sparsekernel,
                    data$edgecode,
                    get(data$usermodel),
                    data$kernel,
                    stratum$mqarray,
                    stratum$cellsize,
                    data$details$r0
                )
                
                sump <- sump + pmix[x,n] * temp
            }
            ## message('n ', n, ' sump ', sump)
            if (sump<=0) NA else freq[n] * log(sump)
        }
        
        ##############################################
        # this is the function that is generally used
        ##############################################
        
        allhistsecrparallel <- function () {
            sump <- numeric(stratum$nc)
            for (x in 1:nrow(pmix)) {
                hx <- if (detectr == "multi") matrix(haztemp$h[x,,], nrow = stratum$m) else -1 ## lookup sum_k (hazard)
                hi <- if (detectr == "multi") haztemp$hindex else -1                        ## index to hx
                temp <-  allhistsecrparallelcpp(
                    as.integer(x-1),
                    as.integer(type),
                    as.integer(stratum$m),
                    as.integer(stratum$nc),
                    as.integer(binomN),
                    as.integer(data$details$CJSp1),
                    as.integer(data$details$grain),
                    as.integer(data$ncores),
                    as.double (stratum$primaryintervals),
                    as.integer(stratum$cumss),
                    as.matrix (stratum$capthist),     
                    as.integer(stratum$fi),
                    as.integer(stratum$li),
                    as.double (if (detectr %in% c("multi", "poissoncount")) hk else gk), ## precomputed probability or hazard
                    as.matrix (realparval),
                    as.integer(PIA),
                    as.integer(PIAJ),
                    as.matrix (stratum$usge),
                    as.matrix (hx),                
                    as.matrix (hi),      
                    as.integer(data$movementcode),
                    as.logical(data$sparsekernel),
                    as.logical(data$anchored),
                    as.integer(data$edgecode),
                    as.character(data$usermodel),  
                    as.integer(data$moveargsi),
                    as.matrix (data$kernel),
                    as.matrix (stratum$mqarray),
                    as.double (stratum$cellsize),
                    as.double (data$details$r0),
                    as.matrix (settlement))
                sump <- sump + pmix[x,] * temp
            }
            if (any(is.na(sump)) || any(sump<=0)) NA else freq * log(sump)
        }
        ########################################################################
        
        freq <- covariates(stratum$capthist)$freq
        if (is.null(freq)) freq <- rep(1, stratum$nc)
        if (length(freq) == 1) freq <- rep(freq, stratum$nc)
        ncf <- sum(freq)

        trps <- traps(stratum$capthist)
        if (!is.null(stratum$mask)) area <- attr(stratum$mask,'area')
        else area <- 0
        
        detectr <- detector(trps)[1]
        if (detectr == 'count') {
            detectr <- if (data$binomN == 0) "poissoncount" else "binomialcount"
        }
        binomN <- switch (detectr, 
            multi = -2, proximity = -1, 
            poissoncount = 0, binomialcount = data$binomN, 
            polygon = -1, polygonX = -2, -9)
        if (binomN < -2)
            stop("open-population secr requires multi, proximity or count detector type")
        
        if (data$details$debug>1) browser()

        nc1 <- max(stratum$nc,1)
        S <- stratum$cumss[stratum$J+1]
        PIA  <- data$design$PIA [stratum$i, 1:nc1, 1:S, 1:stratum$k, , drop = FALSE]
        PIAJ <- data$design$PIAJ[stratum$i, 1:nc1, 1:stratum$J, , drop = FALSE]
        
        #-----------------------------------------
        
        # optionally model settlement as function of mask covariates (and others)
        if (!is.null(settle)) {
            settlement <- getmaskpar (
                settle, 
                m = stratum$m,
                stratumi = stratum$i 
                )
        }
        else {
            settlement <- 1
        }
        #-----------------------------------------
        
        ## number of threads was set in openCR.fit
        temp <- makegkParalleldcpp (as.integer(data$detectfn),
            as.integer(.openCRstuff$sigmai[type]),
            as.integer(data$details$grain),
            as.integer(data$ncores),
            as.matrix(realparval),
            as.matrix(stratum$distmat))
        gk <- array(temp[[1]], dim=c(nrow(realparval), stratum$k, stratum$m))  # array form for R use
        hk <- array(temp[[2]], dim=c(nrow(realparval), stratum$k, stratum$m))  # array form for R use
        sumhk <- sum(hk)
        # additional checks 2021-10-09
        if (is.na(sumhk) || !is.finite(sumhk) || sumhk==0) {
            return(NA)   # changed from 1e10 2021-06-01
        }
        
        if (data$details$debug>0) message ("sum(gk) = ", sum(gk))
        
        pmix <- fillpmix2(stratum$nc, data$details$nmix, PIA, realparval)
        S <- ncol(stratum$capthist)
        #-----------------------------------------
        
        if (detectr == "multi") {
            ## R alternative
            if (data$details$R) {
                haztemp <- gethR(stratum$m, PIA, stratum$usge, hk)
            }
            else {
                haztemp <- gethcpp(
                    as.integer(stratum$nc),
                    as.integer(nrow(realparval)),
                    as.integer(data$details$nmix),
                    as.integer(stratum$k),
                    as.integer(stratum$cumss[stratum$J+1]),
                    as.integer(stratum$m),
                    as.integer(PIA),
                    as.matrix(stratum$usge),
                    as.double(hk))
            }
            haztemp$h <- array(haztemp$h, dim = c(data$details$nmix, stratum$m, max(haztemp$hindex)+1))
        }
        else {
            haztemp <- list(h = array(-1, dim=c(data$details$nmix,1,1)), hindex = matrix(-1))
        }
        if (data$details$debug>0) message ("sum(haztemp$h) = ", sum(haztemp$h))
        
        #####################################################################
        ## Vector to store components of log likelihood
        comp <- numeric(4)
        
        #####################################################################
        # Component 1: Probability of observed histories
        if (!data$details$R) {
            ## this is the streamlined option; always uses C code
            prwi <- allhistsecrparallel()
        }
        else {
            ## clunky option using R for debugging
            prwi <- sapply(1:stratum$nc, onehistory, USE.NAMES = FALSE )
        }
        comp[1] <- sum(prwi)
        if (data$details$debug>0) message ("comp[1] = ", comp[1])
        
        #####################################################################
        # Component 2: Probability of unobserved histories a^{-n} in likelihood for uniform-D

        if ((type %in% c(9,10,11,25,30)) ## CL and require global pdot for component 3
            | (type %in% c(7,8,12,13,14,24,31))) {                   ## all other
            LR <- data$learnedresponse
            if (LR) {   ## overwrite gk with model for naive animal
                temp <- makegkParalleldcpp (as.integer(data$detectfn),
                    as.integer(.openCRstuff$sigmai[type]),
                    as.integer(data$details$grain),
                    as.integer(data$ncores),
                    as.matrix(realparval0),
                    as.matrix(stratum$distmat))
                gk <- array(temp[[1]], dim=c(nrow(realparval), stratum$k, stratum$m))  # array form for R use
                hk <- array(temp[[2]], dim=c(nrow(realparval), stratum$k, stratum$m))  # array form for R use
                PIA0 <- data$design0$PIA[stratum$i, 1:nc1, 1:S, 1:stratum$k, , drop = FALSE]
                PIAJ0 <- data$design0$PIAJ[stratum$i, 1:nc1, 1:stratum$J, , drop = FALSE] 
            }
            ## else use gk as gk0, hk as hk0
            pdot <- rep(0, stratum$nc)  # vector: one value for each unique observed history
            for (x in 1:data$details$nmix) {   # loop over latent classes
                if (data$details$R) {
                    pch1 <-  PCH1secr(
                        type,
                        as.logical(data$design0$individual),
                        x,
                        stratum$nc,
                        stratum$J,
                        stratum$cumss,
                        stratum$k,
                        stratum$m,
                        realparval0,
                        if (LR) PIA0 else PIA,
                        if (LR) PIAJ0 else PIAJ,
                        if (detectr %in% c("multi", "poissoncount")) hk else gk,  ## 2019-05-19
                        binomN,
                        stratum$usge,
                        stratum$primaryintervals,
                        data$moveargsi,
                        data$movementcode,
                        data$sparsekernel,
                        data$edgecode,
                        get(data$usermodel),   # bug fixed 2019-04-14, 2019-05-07
                        data$kernel,
                        stratum$mqarray,
                        stratum$cellsize,
                        data$details$r0
                    )
                }
                else {
                    pch1 <-  PCH1secrparallelcpp(
                        as.integer(x-1),
                        as.integer(type),
                        as.integer(data$details$grain),
                        as.integer(data$ncores),
                        as.logical(data$design0$individual),
                        as.integer(stratum$J),
                        as.integer(stratum$m),
                        as.integer(stratum$nc),
                        as.integer(stratum$cumss),
                        as.matrix (realparval0),
                        as.integer(if (LR) PIA0 else PIA),
                        as.integer(if (LR) PIAJ0 else PIAJ),
                        as.double (if (detectr %in% c("multi", "poissoncount")) hk else gk),  ## 2019-05-19
                        as.integer(binomN),
                        as.matrix (stratum$usge),
                        as.double (stratum$primaryintervals),
                        as.integer(data$moveargsi),
                        as.integer(data$movementcode),
                        as.logical(data$sparsekernel),
                        as.logical(data$anchored),
                        as.integer(data$edgecode),
                        as.character(data$usermodel),
                        as.matrix(data$kernel),
                        as.matrix(stratum$mqarray),
                        as.double (stratum$cellsize),
                        as.double (data$details$r0),
                        as.matrix (settlement)
                    )
                }
                pdot <- pdot + pmix[x] * pch1
            }
            pdot <- rep(pdot, freq)
            comp[2] <- - sum(log(pdot))    ## log(1 / a_i)
        }

        #####################################################################
        # Component 3: Probability of observing nc animals (non-CL types)
        if (type %in% c(7,8,12,13,14,24, 31)) {
            if (type %in% c(7, 12, 13, 24, 31))
                Dsuper <- realparval[nrow(realparval)*3 + stratum$i] # Dsuper direct
            else  {        # type %in% c(8, 14)) D or B parameterisation
                Dsuper <- getD(type, stratum$J, data$details$nmix, pmix,
                    realparval, PIAJ,
                    stratum$primaryintervals)
            }
            A <- maskarea(stratum$mask)
            N <- Dsuper * A
            # impose constraint: return with invalid result code if not possible
            if (N < ncf) return (1e10);
            meanpdot <- ncf / sum(1/pdot)
            ## cf CLmeanesa in 'secr'
            
            comp[3] <- switch (data$distrib+1,
                dpois(ncf, N * meanpdot, log = TRUE),
                lnbinomial (ncf, N, meanpdot),
                NA)
        }
        
        #####################################################################
        
        comp
    }   # end of onestratumll
    
    ############################################################################
    # main line
    #------------------------------------------------------------
    # Fixed beta
    fb <- data$details$fixedbeta
    if (!is.null(fb)) {
        fb[is.na(fb)] <- beta
        beta <- fb    ## complete
    }
    #------------------------------------------------------------
    
    if (is.null(data$design$designMatrices$settle)) {
        settle <- 0
    }
    else {
        settle <- getsettle (
            data$design$designMatrices$settle, 
            beta, 
            data$parindx, 
            data$link, 
            data$fixed,
            nmask = max(sapply(data$stratumdata, '[[', 'm')),
            nstrata = length(data$stratumdata),
            nsessions = max(sapply(data$stratumdata, '[[', 'J')),
            parameter = 'settle') 
        data$design$designMatrices$settle <- NULL ## not to confuse makerealparameters
    }
    
    #------------------------------------------------------------
    # Real parameters
    realparval  <- makerealparameters (data$design, beta, data$parindx, data$link, data$fixed)
    if (data$learnedresponse)
        realparval0 <- makerealparameters (data$design0, beta, data$parindx, data$link, data$fixed)
    else realparval0 <- realparval
    if (data$details$debug>0) print(realparval)    

    #------------------------------------------------------------
    # check valid parameter values
    if (!all(is.finite(realparval))) {
        return (1e10)
    }
    
    #------------------------------------------------------------
    type <- typecode(data$type)
    if (type < 0) stop ("Invalid likelihood type")
    
    if (data$details$debug>2) browser()
    
    #------------------------------------------------------------
    ## call onestratumll
    compbystratum <- lapply(data$stratumdata, onestratumll)
    compbystratum <- matrix(unlist(compbystratum), ncol = 4, byrow = TRUE)
    #------------------------------------------------------------
    ## optional multinomial term (not if CJS)
    if (data$details$multinom & !(type %in% c(6))) {
        compbystratum[,4] <- data$logmult   ## precalculated 2021-03-30
    }
    #------------------------------------------------------------
    ## log-likelihood as sum of components
    loglik <- sum(compbystratum)
    #------------------------------------------------------------
    compbystratum <- cbind(compbystratum, apply(compbystratum,1,sum))
    colnames(compbystratum) <- c('Comp1', 'Comp2', 'Comp3', 
        'logmultinom', 'Total')
    if (nrow(compbystratum)>1) {
        compbystratum <- rbind(compbystratum, apply(compbystratum,2,sum))
        rownames(compbystratum) <- c(paste0('Stratum', 
            1:(nrow(compbystratum)-1)), 'Total')
    }
    else {
        rownames(compbystratum) <- ''
    }
    #------------------------------------------------------------
    
    ## debug
    if (data$details$debug>=1) {
        for (r in 1:nrow(compbystratum)) {
            message("Likelihood components, stratum ", r, " ", 
                paste(format(compbystratum[r,], digits=10), collapse = ' '))
        }
        message("Total ", format(loglik, digits = 10))
        if (data$details$debug>1) browser()
    }
    ## optionally display message for this iteration on console or log file 
    .openCRstuff$iter <- .openCRstuff$iter + 1
    fb <- data$details$fixedbeta  # afresh
    if (!is.null(fb)) beta <- beta[is.na(fb)]
    progressstring <- paste(
        c(format(.openCRstuff$iter, width=4), " ",
            formatC(round(loglik,dig), format='f', digits=dig, width=betaw+2), " ",
            formatC(beta, format='f', digits=dig+1, width=betaw)), collapse = " ")
    if (data$details$trace) {
        if (((.openCRstuff$iter-1) %% data$details$trace) == 0) {
            message(progressstring)
            flush.console()
        }
    }
    logfilename <- data$details$log
    if (logfilename != "" && is.character(logfilename)) {
        progressstring <- paste(
            c(format(.openCRstuff$iter, width=4), "   ",
                formatC(round(loglik,dig), format='f', digits=dig+2, width=betaw+3), " ",
                formatC(beta, format='f', digits=dig+3, width=betaw+3)), collapse = " ")
        cat(progressstring, file = logfilename, sep="\n", append=TRUE)
    }
    
    if (oneeval) {
        out <- c(loglik, beta)
        attr(out, 'components') <- compbystratum
        out
    }
    else {
        if (is.finite(loglik)) -loglik   # return the negative loglikelihood
        else 1e10
    }
    
}
############################################################################################

