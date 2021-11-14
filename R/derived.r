## package 'openCR'
## derived.R
## 2011-12-30, 2013-01-20, 2017-11-21, 2017-12-21, 2018-02-15, 2018-05-25, 
## 2018-10-31 fixed bug in f, lambda
## 2018-10-31 print.derivedopenCR applies Dscale to session-specific D

################################################################################

derived.openCRlist <- function (object, newdata = NULL, all.levels = FALSE,
    Dscale = 1,  HTbysession = FALSE, ...) {
    lapply(object, derived, newdata = newdata, all.levels = all.levels,
        Dscale = Dscale, HTbysession = HTbysession, ...)
}

derived.openCR <- function (object, newdata = NULL, all.levels = FALSE,
    Dscale = 1, HTbysession = FALSE, ...) {
    allvars <- unlist(sapply(object$model, all.vars))
    if ('h2' %in% allvars | 'h3' %in% allvars)
        warning ("derived.openCR does not handle finite mixtures")
    if (is.null(newdata)) {
        # newdata <- openCR.make.newdata(object, all.levels = all.levels)
        newdata <- makeNewData(object, all.levels = all.levels)
    }
    onestratum <- function (newdata, stratumi) {
        nnew <- nrow(newdata)
        J <- length(primaryintervals(object)[[stratumi]]) + 1   ## number of primary sessions
        if (!nnew %% J == 0)
            stop ("rows in newdata should be multiple of number of sessions")
        if (is.null(newdata$session))
            stop ("newdata should contain column 'session'")
        ## optional recursive call for multiple levels etc.
        if (nnew > J) {
            ngrp <- nnew %/% J
            newdata <- split(newdata, rep(1:ngrp, each = J))
            out <- vector('list', ngrp)
            for (n in 1:ngrp) {
                out[[n]] <- derived(object, newdata[[n]],
                    Dscale = Dscale, HTbysession = HTbysession,
                    ...)
            }
            class (out) <- c("derivedopenCR","list")
        }
        else {
            
            beta    <- object$fit$par
            parindx <- object$parindx
            link    <- object$link
            fixed   <- object$fixed
            design0 <- object$design0
            details <- object$details
            type    <- object$type
            if (!is.null(object$stratified) && object$stratified) {
                capthist <- object$capthist[[stratumi]]
                mask <- object$mask[[stratumi]]
            }
            else {
                capthist <- object$capthist
                mask <- object$mask
            }
            primaryintervals <- primaryintervals(object)[[stratumi]]
            stratumn <- sum(covariates(capthist)$freq)

            getreal <- function (par) {
                j <- if (grepl("super",par)) 1.0 else J
                fxd <- object$fixed
                if (!is.null(fxd[[par]]))
                    rep(fxd[[par]], j)
                else {
                    predict(object, newdata)[[par]][1:J,'estimate']
                }
            }
            getfbeta <- function (beta, phij) {
                # return fj for inputs phij and beta
                # phij is per session phi (adjusted for intervals)
                d <- beta  # for beta[1]
                for (j in J1) d[j+1] <- d[j] * phij[j] + beta[j+1]
                c(beta[-1]/d[-J], NA)
            }
            getkappa <- function(p, phi, f) {
                J <- length(p)
                tau <- kappa <- numeric(J)
                tau[1] <- 1/p[1]
                kappa[1] <- 1
                for (j in 1:(J-1)) {
                    kappa[j+1] <- kappa[j]/p[j] * (1-p[j]) * phi[j] * p[j+1] + tau[j] * f[j] * p[j+1]
                    tau[j+1] <- tau[1] * prod((phi+f)[1:j])
                }
                kappa[1] <- NA
                kappa
            }
            
            ## alternatively: makerealparameters (design0, beta, parindx, link, fixed)
            J1 <- 1:(J-1)
            nc <- nrow(capthist)
            persession <- function (x) c(x[-J]^primaryintervals, NA)
            stdrate <- function (x) c(x[-J]^(1/primaryintervals), NA)
            # note b, B always on per interval basis - no scaling
            
            phi <- getreal('phi')
            phij <- persession(phi)
            
            fj <- NULL
            out <- NULL
            
            if (type %in% c('CJS', 'CJSsecr')) {
                warning ("derived.openCR is not implemented for type ", type)
                return(NULL)
            }
            
            if (type %in% "JSSAN") {
                Nj <- getreal('N')
                superN <- Nj[1] + sum (Nj[-1] - (Nj * phij)[-J])
                B <- c(Nj[1], Nj[-1] - (Nj*phij)[-J])
                b <- B/sum(B)
            }
            else {
                if (type %in% "JSSAsecrD") {
                    Dj <- getreal('D')
                    B <- c(Dj[1], Dj[-1] - (Dj*phij)[-J])
                    superD <- sum(B)
                    b <- B/sum(B)
                    fj <- c(B[-1]/Dj[-J], NA)
                }
                ## otherwise extract session-specific b
                else if (type %in% c('JSSAb','JSSAbCL', 'PLBb', 
                    'JSSAsecrb','JSSAsecrbCL', 'PLBsecrb', 
                    'JSSARET')) {
                    b <- getreal('b') # predicted$b[,'estimate']
                }
                else if (type %in% c(
                    'JSSAf', 'JSSAfCL', 'PLBf',
                    'JSSAg', 'JSSAgCL', 'PLBg',
                    'JSSAk', 'JSSAkCL', 'PLBk',
                    'JSSAl', 'JSSAlCL', 'PLBl',
                    'JSSAsecrf', 'JSSAsecrfCL', 'PLBsecrf',
                    'JSSAsecrg', 'JSSAsecrgCL', 'PLBsecrg', 
                    'JSSAsecrl', 'JSSAsecrlCL', 'PLBsecrl', 
                    'Pradelg', 'Pradel')) {
                    
                    if (type %in% c('JSSAl','JSSAlCL', 'PLBl', 
                        'JSSAsecrl','JSSAsecrlCL', 'PLBsecrl',
                        'Pradel')) {
                        lambda <- getreal('lambda')
                        f <- lambda-phi
                    }
                    else if (type %in% c('JSSAg','JSSAgCL', 'PLBg',
                        'JSSAsecrg','JSSAsecrgCL', 'PLBsecrg',
                        'Pradelg')) {
                        gamma <- getreal('gamma')
                        gam1 <- c(gamma[-1],NA)
                        f <- phi * (1/gam1 - 1)
                    }
                    else if (type %in% c('JSSAk','JSSAkCL', 'PLBk')) {
                        kappa <- getreal('kappa')
                        kappa <- c(1, kappa[-1])
                        f <- tau <- numeric(J)
                        p <- getreal('p')
                        tau[1] <- 1/p[1]
                        for (j in 1:(J-1)) {
                            f[j] <- (kappa[j+1] - kappa[j]/p[j] * (1 - p[j]) * phi[j] * p[j+1]) / (tau[j] * p[j+1])
                            tau[j+1] <- tau[1] * prod(phi[1:j] + f[1:j])
                        }
                    }
                    else {
                        f <- getreal('f')
                    }
                    
                    lambdaj <- persession(phi + f)
                    ## fj <- persession(f)
                    ## bug fixed 1.3.6 2019-04-02
                    fj <- lambdaj - persession(phi)
                    
                    d <- numeric(J)
                    d[1] <- 1;
                    for (j in 2:J) d[j] <- d[j-1] * lambdaj[j-1]
                    
                    b <- numeric(J)
                    b[1] = 1;
                    for (j in 2:J) b[j] <- fj[j-1] * d[j-1]
                    b <- b / sum(b)
                    
                }
                else if (type %in% "JSSAB") {
                    BN <- getreal('BN')
                    b <-  BN / sum(BN)
                }
                else if (type %in% "JSSAsecrB") {
                    BD <- getreal('BD')
                    b <-  BD / sum(BD)
                }
                else stop("unrecognized model type")
            }
            
            if (is.null(fj)) fj <- getfbeta(b, phij)
            
            df <- data.frame(newdata, JS.counts(capthist),
                time = cumsum(c(0,primaryintervals)))
            parnames <- names(predict(object))
            for (i in parnames) {
                df[,i] <- getreal(i)
            }
            ## 2020-11-01 allow for fixed parameters
            for (i in names(object$fixed)) {
                df[,i] <- object$fixed[[i]]
            }
            #
            # if (grepl('secr', type)) {
            #     df$lambda0 <- getreal('lambda0')
            #     df$sigma <- getreal('sigma')
            # }
            # else {
            #     df$p <- getreal('p')
            # }
            # df$phi <- phi
            
            df$b <- b
            
            ## WARNING: UNRESOLVED INCONSISTENCY IN ADJUSTMENT FOR INTERVAL phi, f
            ## bug fix 2018-10-31
            ## df$lambda <- stdrate(phij+ fj)
            df$lambda <- stdrate(phij) + stdrate(fj)
            df$f <- df$lambda - df$phi
            df$gamma <- c(NA, (df$phi / (df$f + df$phi))[-J])
            
            ## Nonspatial
            ## superN and time-specific Nj
            if (type %in% c('JSSAN',
                'JSSAf', 'JSSAfCL', 'PLBf',
                'JSSAl', 'JSSAlCL', 'PLBl',
                'JSSAg', 'JSSAgCL', 'PLBg',
                'JSSAk', 'JSSAkCL', 'PLBk',
                'JSSAb', 'JSSAbCL', 'PLBb',
                'Pradel','Pradelg',
                'JSSARET', 'JSSAB')) {
                
                df$kappa <- getkappa (df$p, df$phi, df$f)  ## moved from outside nonspatial 2020-12-04
                
                if (type %in% "JSSAN") {
                    ## as before...
                    Nj <- getreal('N')
                    superN <- Nj[1] + sum (Nj[-1] - (Nj * phij)[-J])
                }
                else {
                    if (type %in% c('JSSAf','JSSAb','JSSAl','JSSAk','JSSAg','JSSARET'))
                        superN <- getreal('superN')[1]
                    else if (type %in% "JSSAB" )
                        superN <- sum(BN)
                    else {
                        p <- openCR.pdot(object, bysession = FALSE)
                        superN <- sum(1 / rep(p, covariates(capthist)$freq))
                    }
                    # assume no mixture for now but see Pledger et al 2010 p885
                    Nj <- numeric(J)
                    if (HTbysession) {
                        p <- openCR.pdot(object, bysession = TRUE)  ## session x ch
                        sess <- primarysessions(intervals(capthist)) ## differs from 'intervals'
                        OK <- apply(capthist,1, by, sess, sum)>0
                        freq <- sweep(OK, MARGIN=2, STATS=covariates(capthist)$freq, FUN = "*")
                        p1 <- lapply(1:J, function(j) 1/ rep(p[j,], freq[j,]))
                        Nj <- sapply(p1, sum)
                    }
                    else {
                        Nj[1] <- superN * b[1]
                        for (j in 2:J) Nj[j] <- Nj[j-1] * phij[j-1] + superN*b[j]
                    }
                }
                df$N <- Nj
                out <- list(totalobserved = stratumn, parameters = parnames, superN = superN, estimates = df)
            }
            
            ## Spatial
            ## superD and time-specific Dj
            
            if (type %in% c('JSSAsecrf','JSSAsecrfCL', 'PLBsecrf',
                'JSSAsecrl','JSSAsecrlCL', 'PLBsecrl',
                'JSSAsecrg','JSSAsecrgCL', 'PLBsecrg',
                'JSSAsecrb','JSSAsecrbCL', 'PLBsecrb',
                'JSSAsecrB',
                'JSSAsecrD')) {
                if (type %in% c('JSSAsecrf','JSSAsecrb','JSSAsecrl','JSSAsecrg'))
                    superD <- getreal('superD')[1]
                else if (type %in% "JSSAsecrB" )
                    superD <- sum(BD)
                else {
                    # esa across all sessions (i.e. using probability animal 
                    # initially at x is detected in at least one session)
                    a <- openCR.esa (object, bysession = FALSE)
                    superD <- sum(1/ a)
                }
                # assume no mixture for now but see Pledger et al 2010 p885
                if (type != 'JSSAsecrD') {  # 2020-10-28 do not recalc Dj
                    Dj <- numeric(J)
                    if (HTbysession) {
                        a <- openCR.esa(object, bysession = TRUE)  ## session x ch
                        # 2020-10-28 freq expansion (unsqueezing) shifted to openCR.esa
                        Dj <- sapply(a, function(x) sum(1/x))
                    }
                    else {
                        Dj[1] <- superD * b[1]
                        for (j in 2:J) Dj[j] <- Dj[j-1] * phij[j-1] + superD*b[j]
                    }
                }
                df$D <- Dj
                
                out <- list(totalobserved = stratumn, parameters = parnames, superD = superD,
                    estimates = df, Dscale = Dscale)
            }
            class (out) <- c("derivedopenCR","list")
        }
        out
    }
    newdata <- split(newdata, newdata$stratum)
    out <- mapply(onestratum, newdata, 1:length(newdata), SIMPLIFY = FALSE)
    names(out) <- paste('stratum', names(newdata))
    if (length(out)==1) out[[1]] else out
}

print.derivedopenCR <- function (x, Dscale = NULL, legend = FALSE, ...) {
    
    if (inherits(x[[1]], "derivedopenCR")) {
        x1 <- x[[1]]  ## multi-level derived
    }
    else {
        x1 <- x
        x <- list(x)
    }
    args <- list(...)
    ngrp <- length(x)
    if (is.null(Dscale)) Dscale <- x1$Dscale
    D.units <- switch(paste0('D',Dscale),
        D1 = 'per ha',
        D100 = 'per km^2',
        D10000 = 'per 100 km^2',
        D100000 = 'per 1000 km^2',
        paste0('per ', Dscale, ' ha'))
    
    ##--------------
    ## Header
    
    cat ("Total number observed", x1$totalobserved, "\n")
    cat ("Parameters in model", paste(x1$parameters, collapse=', '), "\n")
    if (!is.null(x1$superN))
        cat ("Superpopulation size", x1$superN, "\n")
    if (!is.null(x1$superD)) {
        cat ("Superpopulation density", x1$superD  * Dscale, D.units, "\n")
    }
    
    ##--------------
    ## Estimates
    
    if (!"row.names" %in% names(args))
        args$row.names <- FALSE
    if (!"digits" %in% names(args))
        args$digits <- 4
    cat ("Session-specific counts and estimates:\n\n")
    if (ngrp>1) {
        cat ("Warning: multi-group newdata; counts and abundance estimates refer to all animals\n\n")
    }
    
    for (n in 1:ngrp) {
        args$x <- x[[n]]$estimates
        for (parm in c('D','BD')) {
            if (parm %in% names(args$x)) 
                args$x[,parm] <- args$x[,parm] * Dscale
        }
        do.call(print.data.frame, args)
        cat ("\n")
    }
    
    ##--------------
    ## Legend
    
    fields <- c(
        '-------@-----------------------------------------',
        'stratum@independent stratum',
        'session@primary session',
        't@primary session',
        'n@number observed',
        'R@number released',
        'm@number already marked',
        'r@number recaptured in later session',
        'z@number known alive but not caught',
        'time@accumulated time since start',
        'p@detection probability per secondary session',
        'lambda0@intercept of detection function per secondary session',
        'sigma@spatial scale of detection function in metres',
        'p@detection probability per secondary session',
        'phi@apparent survival per unit time',
        'b@entry probabilities',
        'lambda@population growth rate per unit time',
        'f@per capita recruitment per unit time',
        'gamma@seniority (cf reverse-time phi)',
        'kappa@recruitment parameter of Link and Barker (2005)',
        'BN@number of recruits',
        'BD@density of recruits',
        'N@population size',
        paste0('D@density ', D.units))
    leg <- data.frame(do.call(rbind, strsplit(fields, '@')))
    names(leg) <- c('Field', 'Definition')
    
    if (legend) {
        fields <- names(x1$estimates)
        ord <- match(fields, leg[,'Field'])
        leg <- leg[c(1,ord),]
        OK <- leg$Field %in% names(x1$estimates) | leg$Field == '-------'
        print(leg[OK ,,drop=FALSE], right=FALSE, row.names=FALSE)
        cat("\n")
    }
    
}
