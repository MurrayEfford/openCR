###############################################################################
## package 'openCR'
## methods.R
## 2011-12-09, 2015-02-06, 2017-05-12, 2017-05-15, 2017-12
## 2018-01-24 AIC moved to own file
###############################################################################

####################################################
# Generic methods for extracting attributes etc
####################################################

####################################################
## Class : openCR
## open population capture-recapture model fit
####################################################

#---------------------------------------------------
# 2021-04-22

strata <- function (object, ...) UseMethod("strata")

strata.default <- function (object, ...)     {
    ## bypass strata attribute for multi-session objects 
    ## use names(object) for lists i.e. multi-session objects
    if (ms(object)) {
        temp <- names(object)
    }
    else {
        temp <- attr(object,'stratum', exact = TRUE)
        if (is.null(temp)) temp <- 1   
    }
    names(temp) <- NULL
    as.character(temp)     
}

'strata<-' <- function (object, value) {
    if (ms(object)) {
        if (length(value) != length(object))
            stop ("invalid replacement value")
        for (i in 1:length(object)) strata(object[[i]]) <- value[i] 
        structure (object, names = as.character(value))
    }
    else {
        if (length(value) > 1)
            stop ("requires only one name per stratum")
        structure (object, stratum = as.character(value))
    }
}
#---------------------------------------------------

coef.openCR <- function (object, alpha=0.05, ...) {
    beta   <- object$fit$par
    if (!is.null(object$beta.vcv))
        sebeta <- suppressWarnings(sqrt(diag(object$beta.vcv)))
    else sebeta <- rep(NA, length(beta))
    z <- abs(qnorm(1-alpha/2))
    temp <- data.frame(
        row.names = object$betanames,
        beta    = beta,
        SE.beta = sebeta,
        lcl = beta - z*sebeta,
        ucl = beta + z*sebeta
        )
    attr(temp, 'alpha') <- alpha
    temp
}

############################################################################################

print.openCR <- function (x, newdata = NULL, alpha = 0.05, svtol = 1e-5, ...) {

    secrmodel <- grepl("secr", x$type)
    det <- detector(traps(x$capthist))
        
    cat ('\n')

    cl <- paste(names(x$call)[-1],x$call[-1], sep=' = ', collapse=', ' )
    cl <- paste('openCR.fit(', cl, ')')

    cat(strwrap(cl, getOption('width')), sep='\n  ')
    cat ('openCR ', x$version, ', ', x$starttime, '\n', sep='')
    cat ('elapsed time ', round(x$proctime/60,3), ' minutes\n', sep='')
    cat ('\n')

    ###################
    ## Data description

    onestratum <- function (capthist, i) {
        freq <- covariates(capthist)$freq 
        if (is.null(freq)) freq <- rep(1,nrow(x$capthist))
        n  <- sum(freq)                 # number caught
        ncapt <- sum(freq * apply(abs(capthist),1,sum))
        intervals <- intervals(capthist)
        nprimary <- sum(intervals>0)+1
        nsecondary <- ncol(capthist)
        temp <- ''
        cat ('\n')
        if (nstrata>1) cat ('Stratum         : ', stratanames[i], temp, '\n')
        
        cat ('N animals       : ', n, temp, '\n')
        cat ('N detections    : ', ncapt, '\n')
        cat ('N sessions      : ', nprimary)
        if (nsecondary>nprimary) cat (paste0(' (secondary ', nsecondary, ')'))
        cat('\n')
        cat ('Intervals       : ', paste(intervals, collapse = ' '), '\n')
        
    }
    stratanames <- strata(x$capthist)
    nstrata <- length(stratanames)
    if (nstrata>1) {
        cat ('Stratified model, ', nstrata, ' strata \n')
        print (summary(x$capthist, terse = TRUE))
        q <- sapply(x$capthist, function(y) attr(y,'q'))
        mapply(onestratum, x$capthist, 1:nstrata)
    }
    else {
        onestratum(x$capthist, 1)
    }
    ####################
    ## Model description
    
    Npar <- max(unlist(x$parindx))
    if (!is.null(x$stratified) && x$stratified)
        n <- sum(unlist(sapply(covariates(x$capthist), '[[', 'freq')))
    else 
        n <- sum(covariates(x$capthist)$freq)
    
    if (!is.null(x$details$fixedbeta)) {
        Npar <- Npar - sum(!is.na(x$details$fixedbeta))
    }
    AICval <- 2*(x$fit$value + Npar)
    AICcval <- ifelse ((n - Npar - 1) > 0,
        2*(x$fit$value + Npar) + 2 * Npar * (Npar+1) / (n - Npar - 1),
        NA)
    cat ('\n')
    cat ('Analysis type   : ', x$type, '\n')
    cat ('Model           : ', model.string(x$model, x$details$userDfn), '\n')
    cat ('Fixed (real)    : ', fixed.string(x$fixed), '\n')
    
    if (secrmodel) {
        if (any(det %in% .openCRstuff$countdetectors)) {
            cat ('Count model     :  ')
            if (x$binomN == 0) cat ('Poisson \n')
            else if (x$binomN == 1) cat ('Binomial, size from usage\n')
            else if (x$binomN > 1) cat('Binomial', x$binomN, '\n')
        }
        cat ('Movement model  : ', x$movementmodel, '\n')
    }
    
    cat ('N parameters    : ', Npar, '\n')
    cat ('Log likelihood  : ', -x$fit$value, '\n')
    cat ('AIC             : ', AICval, '\n')
    cat ('AICc            : ', AICcval, '\n')

    cat ('\n')
    cat(str_pad('Parameter',9,"right"), 'Link', '\n')
    for (i in 1:length(x$link)) {
        cat(str_pad(names(x$link)[i],9,"right"), x$link[[i]], '\n')
    }
    
    cat ('\n')
    cat ('Beta parameters (coefficients)', '\n')
    print(coef(x), ...)

    if (!is.null(x$fit$hessian)) {
        cat ('\n')
        cat ("Eigenvalues : ", round(x$eigH, -log10(svtol)), '\n')
        rankH <- length(which(x$eigH > svtol))
        cat ("Numerical rank of Hessian :", rankH, ' ( svtol =', svtol, ')\n')
        nreal <- sum(apply(x$design$parameterTable,2,function(x) length(unique(x))))
        if (rankH < nreal) cat ("Warning: at least one real parameter is not identifiable\n")
        cat ('\nVariance-covariance matrix of beta parameters', '\n')
        print (x$beta.vcv, ...)
    }
    else
        eigH <- NA

    cat ('\n')
    cat ('Fitted (real) parameters evaluated at base levels of covariates', '\n')
    if (!is.null(x$realpar))
        print( x$realpar )
    else {
        temp <- predict (x, newdata, alpha)
        nd <- length(temp)
        if (is.data.frame(temp)) print(temp, row.names = FALSE, ...)
        else for (new in 1:nd) {
                cat('\n', names(temp)[new],'\n')
                print(temp[[new]], row.names = FALSE, ...)
            }
    }
    cat ('\n')
}
############################################################################################

vcov.openCR <- function (object, realnames = NULL, newdata = NULL, byrow = FALSE, ...) {
## return either the beta-parameter variance-covariance matrix
## or vcv each real parameters between points given by newdata (byrow = FALSE)
## or vcv for real parameters at points given by newdata (byrow = TRUE)

    ## include object$details$fixedbeta! 2018-03-09
    
    if (is.null(object$beta.vcv)) {
        return(NULL)
    }
    else {
        if (is.null(dimnames(object$beta.vcv)))
            dimnames(object$beta.vcv) <- list(object$betanames, object$betanames)
    }

    if (is.null(realnames))
        ## average beta parameters
        return( object$beta.vcv )
    else {
        ## average real parameters
        ## vcv among multiple rows
        nreal <- length(realnames)
        ## allow for fixed beta parameters 
        beta <- complete.beta(object)
        beta.vcv <- complete.beta.vcv(object)
        nbeta <- length(beta)
        validlevels <- object$design$validlevels

        if (byrow) {
            ## need delta-method variance of reals given object$beta.vcv & newdata
            if (is.null(newdata)) {
                # newdata <- openCR.make.newdata (object)
                newdata <- makeNewData (object)
            }
            rowi <- function (i) {
                grad <- matrix(0, nrow = nreal, ncol = nbeta)
                dimnames(grad) <- list(realnames, names(beta))
                for (rn in realnames) {
                    ## grad[rn,] <- fdHess (pars = object$fit$par, fun = reali, rn = rn)$gradient
                    par.rn <- object$parindx[[rn]]
                    dframe <- newdata[i,,drop=FALSE]
                    
                    # dframe <- adjustlevels(rn, dframe, validlevels)
                    # mat <- model.matrix(object$model[[rn]], data = dframe,
                    #                     contrasts.arg = object$details$contrasts)
                    # new wrapper function 2021-07-23
                    
                    mat <- get.model.matrix(object$model[[rn]], rn, dframe, validlevels, object$details$contrasts)
                    
                    ## drop pmix beta0 column from design matrix (always zero)
                    if (rn == 'pmix') mat <- mat[,-1,drop=FALSE]
                    # UNCERTAIN HOW TO TREAT b... 2018-03-07
                    # if (object$link[[rn]] == 'mlogit') mat <- mat[,-1,drop=FALSE]
                    lp <- mat %*% matrix(beta[par.rn], ncol = 1)
                    real <- as.vector(untransform (lp, object$link[[rn]]))
                    ## from Jeff Laake's 'compute.real' in RMark...
                    gr <- switch(object$link[[rn]],
                                 logit    = mat * real * (1-real),
                                 mlogit   = mat * real * (1-real),
                                 log      = mat * real,
                                 loglog   = -mat * real * log(real),
                                 identity = mat,
                                 sin      = mat * cos(asin(2*real-1))/2)
                    grad[rn,par.rn] <- as.vector(gr)
                }
                grad %*% beta.vcv %*% t(grad)
            }

            vcvlist <- lapply(1:nrow(newdata), rowi)
            if (length(vcvlist) == 1) vcvlist <- vcvlist[[1]]
            return(vcvlist)
        }
        else {
            newdata <- as.data.frame(newdata)
            if (nrow(newdata)==0) stop ("vcov.openCR requires newdata when byrow = FALSE")
            rownames <- apply(newdata, 1, function(x) paste(names(newdata), '=', x, sep='',
                                                            collapse=','))
            vcvlist <- list()
            for (rn in realnames) {
                par.rn <- object$parindx[[rn]]
                
                # dframe <- adjustlevels(rn, newdata, validlevels)
                # mat <- model.matrix(object$model[[rn]], data = dframe, 
                #                     contrasts.arg = object$details$contrasts)
                # new wrapper function 2021-07-23
                
                mat <- get.model.matrix(object$model[[rn]], rn, newdata, validlevels, object$details$contrasts)
                
                ## drop pmix beta0 column from design matrix (always zero)
                if (rn == 'pmix') mat <- mat[,-1,drop=FALSE]
                lp <- mat %*% matrix(object$fit$par[par.rn], ncol = 1)
                real <- untransform (lp, object$link[[rn]])
                real <- as.vector(real)
                if (!(object$link[[rn]] %in% c('logit','log','loglog','sin','identity')))
                    stop("not working yet for mlogit etc")
                grad <- switch(object$link[[rn]],
                               logit    = mat * real * (1-real),
                               mlogit   = mat * real * (1-real),
                               log      = mat * real,
                               loglog   = -mat * real * log(real),
                               identity = mat,
                               sin      = mat * cos(asin(2*real-1))/2)
                vcvlist[[rn]] <- grad %*% beta.vcv[par.rn, par.rn] %*% t(grad)
                dimnames(vcvlist[[rn]]) <- list(rownames, rownames)
            }
            names (vcvlist) <- realnames
            return (vcvlist)
        }
    }
}
############################################################################################

## 2017-05-12
openCRlist <- function(...) {
    dots <- match.call(expand.dots = FALSE)$...
    allargs <- list(...)
    if (length(allargs)==1 && inherits(allargs[[1]], 'openCRlist')) {
        return (allargs[[1]])
    }
    else {
        if (is.null(names(allargs))) {
            dots2 <- substitute(list(...))[-1]
            names(allargs) <- sapply(dots2, deparse)
        }
        allargs <- lapply(allargs, function(x) if (inherits(x, 'openCR')) list(x) else x)
        temp <- do.call(c, allargs)
        if (is.null(names(temp)))
            names(temp) <- paste("openCR", 1:length(temp), sep="")
        if (!all(sapply(temp, function(x) inherits(x, 'openCR'))))
            stop ("objects must be of class 'openCR' or 'openCRlist'")
        class(temp) <- 'openCRlist'
        temp
    }
}

############################################################################################

## 2021-08-14
## expects all ... to be openCRlist
c.openCRlist <- function(..., recursive = FALSE) {
    slist <- as.list(...)
    result <- NextMethod('c', slist, recursive = FALSE)
    class(result) <- 'openCRlist'
    # names?
    result
}
############################################################################################

## 2021-08-14
## extract method for openCRlist objects
## retains class
'[.openCRlist' <- function(x, i) {
    y <- NextMethod("[")
    class(y) <- class(x)
    y
}
############################################################################################

# 2018-05-13 timevaryingcov handled in secr
# 'timevaryingcov<-' <- function (object, value) {
#     if (is.null(value))
#         structure (object, timevaryingcov = NULL)
#     else {
#         if (ms(object)) {
#             if (inherits(object, 'capthist'))
#                 stop ("cannot set timevaryingcov on multi-session capthist")
#             temp <- object
#             if (!(is.list(value[[1]]) & (length(value) == length(object))))
#                 value <- list(value)
#             for (i in 1:length(object))
#                 timevaryingcov(temp[[i]]) <- value[[i]]
#             temp
#         }
#         else {
#             if (!is.list(value) | is.null(names(value)))
#                 stop("value should be a list of one or more named vectors")
#             if (!is.null(usage(object))) {
#                 ## traps object with usage
#                 OK <- sapply(value, function(x) ncol(usage(object)) == length(x))
#                 if (any(!OK))
#                     warning ("mismatch between number of occasions in usage and timevaryingcov")
#             }
#             else if (inherits(object, 'capthist')) {
#                 ## capthist object
#                 nsessions <- length(unique(primarysessions(intervals(object))))
#                 OK <- sapply(value, function(x) nsessions == length(x))
#                 if (any(!OK))
#                     warning ("mismatch between number of primary sessions in object and timevaryingcov")
#             }
#             
#             if (is.character(value))
#                 if (!all(value %in% names(covariates(object))))
#                     warning ("mismatch between character vector and covariate names")
#             structure (object, timevaryingcov = value)
#         }
#     }
# }
