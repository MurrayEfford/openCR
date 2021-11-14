# non-spatial likelihood

# 2017-05-17 not ready for mixtures
# 2017-05-17 pass just one animal to onehistory
# 2017-05-17 PCH0 called with single representative history
# 2017-05-17 capthist 2D for non-spatial
# 2017-05-22 use data argument (an environment)
# 2018-02-06 drop 'par.' prefix from function names
# 2018-03-26 switch pch0 to pch1
# 2019-04-09 ejected open.secr.loglikfn to logliksecr.R
# 2019-04-09 explicit treatment of count detector; dropuse of data$multi
# 2020-12-07 CJSmte (Markovian Temporary Emigration) trial - not completed
# 2021-04-19 stratified
# 2021-11-08 PIAJ recalc to fix bug in nonspatial learned response models and valgrind issue
# 2021-11-08 PIA, PIAJ overwritten for naive animal, rather than always copied  

# types

# CJS = 1
# JSSAb = 2
# JSSAl = 3
# JSSAf = 4
# CJSmte = 5
# JSSAfCL = PLBf = 15
# JSSAlCL = PLBl = 16
# JSSAbCL = PLBb = 17
# JSSAB = 18
# JSSAN = 19
# Pradel = 20
# JSSARET = 21
# JSSAg = 22
# JSSAgCL = PLBg = 23
# Pradelg = 26
# JSSAfgCL = 27
# JSSAk = 28
# JSSAkCL = PLBk = 29

#---------------------------------------------------------

open.loglikfn <- function (beta, dig = 3, betaw = 8, oneeval = FALSE, data)

# Return the negative log likelihood
# Transformed parameter values are passed in the vector 'beta'
# details$trace=T sends a one-line report to the screen

{
    
    onestratumll <- function(stratum) {
        freq <- covariates(stratum$capthist)$freq
        if (is.null(freq)) freq <- rep(1, stratum$nc)
        if (length(freq) == 1) freq <- rep(freq, stratum$nc)
        ncf <- sum(freq)
        comp <- numeric(4)
        
        nc1 <- max(stratum$nc,1)
        S <- stratum$cumss[stratum$J+1]
        PIA  <- data$design$PIA [stratum$i, 1:nc1, 1:S, , , drop = FALSE]
        PIAJ <- data$design$PIAJ[stratum$i, 1:nc1, 1:stratum$J, , drop = FALSE]

        if (data$details$debug>1) {
            message('Stratum ', stratum$i)
            message('Type    ', type)
            message('J       ', stratum$J)
            message('nmix    ', data$details$nmix)
            message('realparval')
            print(realparval)
            message('table(PIA)')
            print (table(PIA))
            message('intervals')
            print(stratum$primaryintervals)
            flush.console()
            browser()
        }
        
        if (type %in% c(20,26)) {
            # Pradel model
            if (data$details$R) {
                comp <- pradelloglik(type, stratum$JScounts, realparval,  PIAJ, 
                    stratum$primaryintervals)
            }
            else {
                comp[1:2] <- pradelloglikcpp(
                    as.integer(type),
                    as.integer(stratum$JScounts),
                    as.integer(stratum$nc),               ## needed for nrows of PIAJ
                    as.integer(stratum$J),
                    as.integer(data$details$nmix),
                    as.matrix(realparval),
                    as.integer(PIAJ),                     ## index of nc,S,mix to rows
                    as.double(stratum$primaryintervals))  ## number of interval == J-1
            }
        }
        else if (data$details$R & (type %in% c(28,29))) {
            comp[1:4] <- kappaloglik (type, realparval,  PIA, PIAJ, stratum, 
                data$details$nmix, 
                data$details$CJSp1, 
                data$distrib)
        }
        else {
            onehistory <- function (n, pmix) {
                sump <- 0
                for (x in 1:nrow(pmix)) {
                    temp <- prwi(
                        type,
                        1,   # n
                        x,
                        stratum$J,
                        stratum$cumss,
                        data$details$nmix,
                        stratum$capthist[n,, drop = FALSE],
                        stratum$fi[n],
                        stratum$li[n],
                        realparval,
                        PIA [stratum$i, n,,,, drop = FALSE],
                        PIAJ[stratum$i, n,,, drop = FALSE],
                        stratum$primaryintervals,
                        data$details$CJSp1
                        # , data$moveargsi
                    )
                    sump <- sump + pmix[x,n] * temp
                }
                if (any(sump<=0)) {
                    -1e10
                }
                else
                    freq[n] * log(sump)
            }
            allhistparallel <- function () {
                sump <- numeric(stratum$nc)
                for (x in 1:nrow(pmix)) {
                    temp <-  allhistparallelcpp(
                        as.integer(x-1),
                        as.integer(type),
                        as.integer(stratum$nc),
                        as.integer(data$details$CJSp1),
                        as.integer(data$details$grain),
                        as.integer(data$ncores),
                        as.double (stratum$primaryintervals),
                        as.integer(stratum$cumss),
                        as.integer(stratum$capthist),
                        as.integer(stratum$fi),
                        as.integer(stratum$li),
                        as.matrix (realparval),
                        as.integer(PIA),
                        as.integer(PIAJ))
                    sump <- sump + pmix[x,] * temp
                }
                freq * log(sump)  ## return vector of individual LL contributions
            }
            
            pmix <- fillpmix2(stratum$nc, data$details$nmix, PIA, realparval)
            
            #####################################################################
            # Component 1: Probability of observed histories - all models
            if (data$details$R)
                temp <- sapply(1:stratum$nc, onehistory, pmix = pmix)
            else
                temp <- allhistparallel()
            comp[1] <- sum(temp)

            #####################################################################
            # Component 2: Probability of missed animals (all-zero histories)
            # not CJS, CJSmte
            if (type %in% c(2:4,15:19, 21, 22, 23, 27, 28, 29)) {
                LR <- data$learnedresponse
                if (LR) {
                    PIA0 <- data$design0$PIA[stratum$i, 1:nc1, 1:S, , , drop = FALSE] 
                    PIAJ0 <- data$design0$PIAJ[stratum$i, 1:nc1, 1:stratum$J, , drop = FALSE]
                }
                pdot <- rep(0, stratum$nc)
                for (x in 1:data$details$nmix) {   # loop over latent classes
                    if (data$details$R) {
                        pch1 <- PCH1(
                            type,
                            x,
                            stratum$nc,
                            stratum$cumss,
                            data$details$nmix,
                            realparval0,
                            if (LR) PIA0 else PIA,
                            if (LR) PIAJ0 else PIAJ,
                            stratum$primaryintervals)
                    }
                    else {
                        pch1 <-  PCH1cpp(
                            as.integer(type),
                            as.integer(x-1),
                            as.integer(stratum$nc),
                            as.integer(stratum$J),
                            as.integer(stratum$cumss),
                            as.integer(data$details$nmix),
                            as.matrix(realparval0),
                            as.integer(if (LR) PIA0 else PIA),
                            as.integer(if (LR) PIAJ0 else PIAJ),
                            as.double(stratum$primaryintervals))
                    }
                    pdot <- pdot + pmix[x] * pch1
                }
                comp[2] <- - sum(freq * log(pdot))
            }
            
            #####################################################################
            # Component 3: Probability of observing nc animals
            # not CJS, CJSmte
            if (type %in% c(2:4,18,19,21, 22, 28)) {
                if (type %in% c(2,3,4,22,28)) {
                    superN <- realparval[nrow(realparval)*3 + stratum$i] # stratum i Nsuper direct
                }
                else {
                    superN <- getN(type, ncf, stratum$J, data$details$nmix, pmix, 
                        realparval, PIAJ, stratum$primaryintervals)
                }
                meanpdot <- ncf / sum(1/rep(pdot,freq))  ## cf CLmeanesa in 'secr'
                comp[3] <- switch (data$distrib+1,
                    dpois(ncf, superN * meanpdot, log = TRUE),
                    ## lnbinomial (ncf, superN, meanpdot),
                    lnbinomial (ncf, superN + ncf, meanpdot),
                    NA)
            }
            
            #####################################################################
            
            
        }
        comp
    }  # end of onestratumll

    #####################################################################
    # main line
    
    #--------------------------------------------------------------------
    # Fixed beta
    fb <- data$details$fixedbeta
    if (!is.null(fb)) {
        fb[is.na(fb)] <- beta
        beta <- fb    ## complete
    }
    if (data$details$debug>0) {
        print(beta)
    }
    #--------------------------------------------------------------------
    # Real parameters
    realparval  <- makerealparameters (data$design, beta, data$parindx, data$link, data$fixed)
    if (data$learnedresponse)
        realparval0 <- makerealparameters (data$design0, beta, data$parindx, data$link, data$fixed)
    else realparval0 <- realparval
    #-----------------------------------------
    # check valid parameter values
    if (!all(is.finite(realparval))) {
        message ('beta vector : ', paste(beta, collapse=', '))
        message ('real vector : ', paste(realparval, collapse=','))
        warning ("extreme 'beta' in 'openloglikfn' ",
            "(try smaller stepmax in nlm Newton-Raphson?)")
        return (1e10)
    }
    type <- typecode(data$type)
    if (type<0) stop ("Invalid likelihood type")
    if (type %in% 28:29 & any(unlist(data$stratumdata$primaryintervals) !=1)) 
        stop ("kappa parameterisation available only if all intervals = 1")

    if (data$details$debug>2) browser()
    
    #####################################################################
    compbystratum <- lapply(data$stratumdata, onestratumll)
    compbystratum <- matrix(unlist(compbystratum), ncol = 4, byrow = TRUE)
    #####################################################################
    ## optional multinomial term
    if (data$details$multinom & (type %in% c(2,3,4,15,16,17,18,19,21,22,23,27))) {
        compbystratum[,4] <- data$logmult    ## precalculated 2021-03-30
    }
    #####################################################################
    ## log-likelihood as sum of components
    loglik <- sum(compbystratum)
    #####################################################################
    compbystratum <- cbind(compbystratum, apply(compbystratum,1,sum))
    colnames(compbystratum) <- c('Comp1', 'Comp2', 'Comp3', 'logmultinom', 'Total')
    if (nrow(compbystratum)>1) {
        compbystratum <- rbind(compbystratum, apply(compbystratum,2,sum))
        rownames(compbystratum) <- c(paste0('Stratum', 1:(nrow(compbystratum)-1)), 'Total')
    }
    else {
        rownames(compbystratum) <- ''
    }
    #####################################################################
    
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
        c(format(.openCRstuff$iter, width=4), "   ",
        formatC(round(loglik,dig), format='f', digits=dig, width=betaw+2),  " ",
        formatC(beta, format='f', digits=dig+1, width=betaw+1)), collapse = " ")
    if (data$details$trace) {
        if (((.openCRstuff$iter-1) %% data$details$trace) == 0) {
            message(progressstring)
            flush.console()
        }
    }
    logfilename <- data$details$log
    if (logfilename != "" && is.character(logfilename)) {
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

