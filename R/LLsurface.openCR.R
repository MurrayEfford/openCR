############################################################################################
## package 'openCR'
## LLsurface.openCR.R
## evaluate and plot log likelihood surface for two named beta parameters
## last changed 2020-11-02 - handle ncores as in LLsurface.secr; no parallel clusters
############################################################################################

LLsurface.openCR <- function (object, betapar = c('phi', 'sigma'), xval = NULL, yval = NULL,
    centre = NULL, realscale = TRUE, plot = TRUE, plotfitted = TRUE, ncores = NULL, ...) {

    if (inherits(object, 'list')) {
        temp <- list()
        nopenCR <- length(object)
        for (i in 1:nopenCR) {
            temp[[i]] <- LLsurface.openCR (object[[i]], betapar = betapar, xval=xval, yval=yval,
                centre = centre, realscale = realscale, plot = plot, plotfitted = plotfitted, ...)
        }
        invisible(temp)
    }
    else {
        if (('superN' %in% betapar) & (object$distribution=='binomial') & 
            is.null(xval) & is.null(yval)) {
            warning("specify superN values explicitly if distribution = 'binomial'")
        }
        if (is.null(centre))
             centre <- t(coef(object))['beta',]  ## retains names
        else {
            if (is.null(names(centre)))
                names(centre) <- object$betanames
            else
                if (any(names(centre) != object$betanames))
                    stop ("names of 'centre' do not match 'object$betanames'")
        }
        betaindices <- match(betapar, names(centre))
        if ((length(betapar) != 2) | (any(is.na(betaindices))))
            stop ("requires two named beta parameters")
        if (realscale & any(is.na(match(betapar, names(object$link)))))
            stop ("link function not found - see Notes in help")
        linkx <- ifelse(realscale, object$link[[betapar[1]]], 'identity')
        linky <- ifelse(realscale, object$link[[betapar[2]]], 'identity')
        if (is.null(xval)) {
            betax0 <- centre[betaindices[1]]
            realx0 <- untransform(betax0, linkx)
            xval <- transform (seq(0.8,1.2,0.04) * realx0, linkx)
            xval <- sort(xval)
        }
        else if (realscale) xval <- transform(xval, linkx)
        if (is.null(yval)) {
            betay0 <- centre[betaindices[2]]
            realy0 <- untransform(betay0, linky)
            yval <- transform (seq(0.8,1.2,0.04) * realy0, linky)
            yval <- sort(yval) ## to reverse in case of negative realy0
        }
        else if (realscale) yval <- transform(yval, linky)
        varying <- list(xval,yval)
        names(varying) <- betapar

        grid <- expand.grid(c(as.list(centre[-betaindices]), varying))
        ## restore original order
        grid <- grid[, object$betanames]

        ## drop unnecessary options
        details <- replace(object$details, 'hessian', FALSE)
        details$trace <- FALSE
        details$LLonly <- TRUE

        LL <- function (start) {
            suppressWarnings(
                openCR.fit(capthist = object$capthist, 
                           model = object$model,
                           mask = object$mask, 
                           type = object$type, 
                           detectfn = object$detectfn, 
                           binomN = object$binomN, 
                           movementmodel = object$movementmodel, 
                           start = start, 
                           link = object$link, 
                           fixed = object$fixed, 
                           timecov = object$timecov, 
                           sessioncov = object$sessioncov, 
                           agecov = object$agecov,
                           dframe = object$dframe, 
                           details = details, 
                           method = object$fit$method, 
                           ncores = ncores)[1]
            )
        }
        
        message ('Evaluating log likelihood across grid of ', nrow(grid), ' points...')
        flush.console()

        ## 2021-02-26
        ## temp <- apply (grid, 1, LL)
        temp <- apply (as.matrix(grid), 1, LL)
        temp <- matrix(temp, nrow = length(xval))
        
        if (realscale) {
            xval <- round(untransform(xval, linkx),4)
            yval <- round(untransform(yval, linky),4)
            centre[betapar[1]] <- untransform(centre[betapar[1]], linkx)
            centre[betapar[2]] <- untransform(centre[betapar[2]], linky)
        }

        dimnames(temp) <- list(xval, yval)
        if (plot) {
            contour(x=xval, y=yval, z=temp, xlab=betapar[1], ylab=betapar[2], ...)
            if (plotfitted) {
                points(centre[betapar[1]], centre[betapar[2]], pch = 3)
            }
        }
        invisible(temp)
    }
}

## LLsurface.openCR(fit)
## prompt(LLsurface.openCR, file = 'd:/open populations/openCR/man/LLsurface.Rd')
