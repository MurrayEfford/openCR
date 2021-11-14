################################################################################
## package 'openCR'
## kernel.R
## 2019-06-09, 11
## 2020-10-06 tweak to ed, bR
## 2021-02-23 annular etc.
## 2021-02-25 functions not exported: empirical
## 2021-03-01 functions exported pkernel, matchscale
## 2021-03-16 annular2
## 2021-07-18 many changes - gkernel(), frE, frG, frL
## 2021-07-29 BVNzi, BVEzi, uniformzi, frEzi  (11,12,13,14)
## 2021-09-25 BVN2 (18)
## 2021-10-06 IND INDzi, RDE, RDG, RDL
## 2021-10-12 stdmovement() applied before messing with a movementmodel 
##            to reduce multiple cases
## 2021-10-19 RDLS log-sech (log hyperbolic secant) 
##            cf Halley and Inchausti 2002; Van Houtan et al 2010; Chadoeuf et al. 2018
## 2021-10-25 BVC
## 2021-10-28 BVCzi
################################################################################

## fillkernelp is used in prwisecr.R
## otherwise, these functions are used for exploration of kernel characteristics
## kernel is coded independently within likelihood functions

circleintersectsline <- function (x1, x2, y1, y2, R) {
    # from https://mathworld.wolfram.com/Circle-LineIntersection.html
    # assume origin at 0,0
    dx <- x2-x1
    dy <- y2-y1
    dr2 <- dx^2 + dy^2
    D <- x1*y2 - x2 * y1
    inc <- R^2*dr2 - D^2
    sgndy <- if(dy<0) -1 else 1
    if (inc > 0 ) {
        data.frame(
            x = c(
                (D*dy+sgndy*dx*sqrt(inc)) / dr2,
                (D*dy-sgndy*dx*sqrt(inc)) / dr2
            ),
            y = c(
                (D*dx + abs(dy)*sqrt(inc)) / dr2,
                (D*dx - abs(dy)*sqrt(inc)) / dr2
            )
            
        )
    }
    else {
        NULL
    }
}
#-----------------------------------------------------------------------
## fill cells of kernel with probability of movement from central point

fillkernelC <- function (
    J,              ## number of sessions (J-1 potential movements)
    kerneltype,     ## 0 BVN, 1 BVE, 2 usermodel, 3 BVT, etc.
    sparsekernel,   ## TRUE iff sparse
    kernel,
    usermodel = "",
    cellsize,
    r0,
    moveargsi,
    moveargs,       ## J x npar matrix
    normalize)
{
    kernel <- as.matrix(kernel)   # convert from mask dataframe to matrix
    kernel <- sweep(kernel, FUN="-", STATS=apply(kernel,2,mean), MARGIN=2)
    kernel[,] <- as.integer(round(kernel / cellsize))
    kernelp <- fillkernelcpp (
        as.matrix   (kernel), 
        as.integer  (kerneltype), 
        as.logical  (sparsekernel),
        as.double   (cellsize),
        as.double   (r0),
        as.integer  (J),
        as.character(usermodel),
        as.integer  (moveargsi), 
        as.double   (moveargs),
        as.logical  (normalize)    
    )   
    matrix(kernelp, ncol = J-1)
}

make.kernel <- function (
    movementmodel = c('BVN', 'BVE', 'BVC', 'BVT','RDE', 'RDG', 'RDL', 'UNI'),
    kernelradius = 10, 
    spacing, 
    move.a, 
    move.b, 
    sparsekernel = FALSE, 
    clip = FALSE,
    normalize = TRUE,
    stat = c('estimate','lcl', 'ucl'),
    session = 1,
    r0 = 1/sqrt(pi),
    ...) 
{
    if (inherits(movementmodel, 'openCR')) {
        fit <- movementmodel
        if (fit$version < '2.0.0') stop ("model fitted with openCR < 2.0.0")
        if (fit$movementmodel %in% .openCRstuff$movementmodels) {
            r0 <- fit$details$r0
            if (is.null(r0)) r0 <- 1/sqrt(pi)
            stat <- match.arg(stat)
            pred <- predict(fit, ...)
            kernel <- make.kernel(
                movementmodel = fit$movementmodel, 
                kernelradius  = fit$kernelradius, 
                spacing       = secr::spacing(fit$mask), 
                move.a        = pred$move.a[session, stat],   # gather lcl, ucl as well 2021-08-08
                move.b        = pred$move.b[session, 'estimate'], # gather lcl, ucl as well 2021-08-08
                sparsekernel  = fit$sparsekernel, 
                clip          = TRUE,
                normalize     = TRUE,
                r0            = r0)
            # avoid $ subscripting as kernel, kernelradius ambiguous 2021-09-10
            if (!is.null(fit[['kernel']]) && nrow(kernel) != nrow (fit[['kernel']])) 
                stop ('bad match to kernel')
            kernel
        }
        else {
            NULL
        }
    }
    else {
        if (is.function (movementmodel)) {
            moveargs <- formalArgs(movementmodel)
            usermodel <- as.character(substitute(movementmodel))
            movementmodel <- "user"
        }
        else {
            usermodel <- ""   ## changed from NULL 2021-07-27
            movementmodel <- match.arg(movementmodel[1],
                choices = .openCRstuff$movementmodels) 
            movementmodel <- stdmovement(movementmodel)   ## 2021-10-13
        }
        
        if (missing(move.a) && 
                !(movementmodel %in% c("UNI", "IND")) && 
                (movementmodel %in% .openCRstuff$movementmodels)) { 
            stop ("move.a required for movementmodel ", movementmodel)
        }
        if (missing(move.b) && movementmodel %in% c('annular2','annularR', 
            'BVT', 'RDG', 'RDL', 'BVNzi', 'BVEzi', 'BVCzi', 'RDEzi', 'BVN2')) {
            stop ("move.b required for movementmodel ", movementmodel)
        }
        movementcode <- movecode(movementmodel)
        moveargsi <- c(0,0)
        if (missing(move.a) | (movementmodel %in% c('UNI'))) {
            pars <- move.a <- move.b <- NULL
        }
        else {
            pars <- move.a[1]
            moveargsi[1] <- 1
        }
        if (missing(move.b)) {
            move.b <- NULL
        }
        else {
            pars <- c(pars, move.b[1])
            moveargsi[2] <- 2
        }
        if (is.null(pars)) pars <- c(0,0)
        moveargs <- matrix(pars, nrow = 1)   # J-row matrix for fillkernelp
        k2 <- kernelradius
        kernel <- make.mask(type = 'rectangular',
            spacing = spacing,
            buffer = 0, nx = 2 * k2+1, ny = 2 * k2+1)
        ## centre
        kernel[,] <- sweep(kernel, MARGIN=2, FUN = "-", STATS = rep((k2+0.5)*spacing,2))
        
        if (movementmodel %in% c('annular', 'annular2')) {
            r <- (kernel$x^2 + kernel$y^2)^0.5
            origin <- r==0
            ring2 <- (r >= (k2-0.5) * spacing) & (r<(k2+0.5) * spacing)
            if (movementmodel == 'annular') {
                ok <- (origin | ring2)
            }
            else {
                ring1 <- (r >= (k2/2-0.5) * spacing) & (r<(k2/2+0.5) * spacing)
                ok <- (origin | ring1 | ring2)
            }
            kernel <- subset(kernel, ok)
        }
        
        if (sparsekernel) {
            ok <- kernel$x==0 | kernel$y == 0 | kernel$x == kernel$y | kernel$x == -kernel$y
            kernel <- kernel[ok,]
        }
        
        ## optional clipping 
        if (clip) {
            outside <- (kernel$x^2 + kernel$y^2) > ((k2+0.5)*spacing)^2
            kernel <- subset(kernel, !outside)
        }
        # call wrapper for C++ function fillkernelcpp 
        kernelp <- fillkernelC (
            2,              ## number of sessions+1
            movementcode-2, ## kerneltype 0 BVN, 1 BVE, 2 usermodel, 3 BVT, 4 UNI, 
                            ## 5 annular, 6 annular2, 7 annularR,
                            ## 8 RDE, 9 RDG, 10 RDL,
                            ## 11 BVNzi, 12 BVEzi, 13 UNIzi, 14 RDEzi, 16 BVN2 
                            ## 17 RDLS, 18 BVC, 19 BVCzi
            sparsekernel,   ## TRUE iff sparse
            kernel,         ## coordinates  
            usermodel,      ## name of user function (grain = 0 only) 
            cellsize = spacing,  
            r0       = r0,
            moveargsi,      ## which column has move.a, move.b
            t(matrix(moveargs,2,2)),
            normalize = normalize   ## delay normalization?
        )[,1]
        
        covariates(kernel) <- data.frame(kernelp = kernelp)
        
        ## 2021-03-01 cumulative distribution function
        r <- apply(kernel^2,1,sum)^0.5
        p <- covariates(kernel)$kernelp[order(r)]
        distribution <- data.frame(r = sort(r), cumprob = cumsum(p))
        
        attr(kernel, 'distribution')  <- distribution
        attr(kernel, 'movementmodel') <- movementmodel
        attr(kernel, 'sparsekernel')  <- sparsekernel
        attr(kernel, 'k2')            <- k2
        attr(kernel, 'move.a')        <- move.a
        attr(kernel, 'move.b')        <- move.b
        attr(kernel, 'r0')            <- r0
        class(kernel)                 <- c('kernel','mask','data.frame')
        kernel
    }
}

getpstring <- function (movementmodel, move.a, move.b, sep = ' ') {
    npar <- nparmove(movementmodel)
    pstring <- ''
    if (npar>0) pstring <- paste0('move.a = ', move.a)
    if (npar>1) pstring <- paste(c(pstring, paste0(' move.b = ', move.b)), collapse = sep)
    pstring
}

plot.kernel <- function (x, type = 'kernel', contour = FALSE, 
    levels = NULL, text = FALSE, title = NULL, add = FALSE, xscale = 1, ...) {
    if (length(type) > 1) {
        temp <- lapply(type, plot, x = x, contour = contour, levels = levels, text = text, 
            title = title, add = add, xscale = xscale, ...)
        invisible(data.frame(x, kernelp = covariates(x)$kernelp))
    }
    else {
        spacing <- spacing(x)
        k2 <- attr(x, 'k2')
        move.a <- attr(x, 'move.a')
        move.b <- attr(x, 'move.b')
        movementmodel <- attr(x, 'movementmodel')
        type <- match.arg(type, choices = c('kernel', 'gr', 'fr', 'Fr'))
        if (type == 'gr') {
            r2 <- seq(-(k2+0.5)*spacing, (k2+0.5)*spacing, length.out = 301)    
            r2s <- r2 * xscale
            rugs <- seq(-(k2)*spacing, (k2)*spacing, spacing) * xscale
            gr <- suppressWarnings(gkernel(abs(r2), movementmodel, move.a, move.b, 
                truncate = k2*spacing))
            if (!add) {
                plot(0,0, 
                    xlim = c(min(r2s, na.rm=TRUE), max(r2s, na.rm=TRUE)),
                    ylim = c(0, max(gr, na.rm=TRUE)), type='n',
                    xlab = expression(paste('Distance along diameter ', ~italic(r))),
                    ylab = expression(paste('Probability density ', ~italic(g(r))))
                )
            }
            lines(r2s, gr, ...)
            abline(v=0)
            rug(rugs, ticksize = 0.02)
        }
        else if (type == 'fr') {
            r <- seq(0, (k2+0.1)*spacing, length.out = 301)    
            rs <- r * xscale
            dr <- suppressWarnings(dkernel(r, movementmodel, move.a, move.b, 
                truncate = k2*spacing))
            if (!add) {
                plot(0,0, 
                    xlim = c(0,max(rs, na.rm=TRUE)), 
                    ylim = c(0,max(dr, na.rm=TRUE)), type = 'n',
                    xlab = expression(paste('Distance moved ', ~italic(r))),
                    ylab = expression(paste('Probability density ', ~italic(f(r))))
                )
            }
            lines (rs, dr, ...)
        }
        else if (type == 'Fr') {
            r <- seq(0, (k2)*spacing, length.out = 301)    
            rs <- r * xscale
            pr <- suppressWarnings(pkernel(r, movementmodel, move.a, move.b, 
                truncate = k2*spacing))
            if (!add) {
                plot(0,0, xlim = c(0, max(rs, na.rm=TRUE)), ylim = c(0,1), type = 'n',
                    xlab = expression(paste('Distance moved ', ~italic(r))),
                    ylab = expression(paste('Cumulative probability ', ~italic(F(r)))),
                )
            }
            lines(rs, pr, ...)
        }
        else {
            
            npar <- nparmove(movementmodel)
            pstring <- getpstring(movementmodel, round(move.a,3), round(move.b,3))
            kernelp <- covariates(x)$kernelp
            dots <- list(...)
            if ('border' %in% names(dots)) {
                border <- dots$border
                dots$border <- NULL
            }
            else border <- 1
            if ('meshcol' %in% names(dots)) {
                meshcol <- dots$meshcol
                dots$meshcol <- NULL
            }
            else meshcol <- 'white'
            msk <- x
            class(msk) <- c('mask','data.frame')
            
            if (sum(!is.na(kernelp))==0) {
                warning ("no non-missing kernel probabilities")
                plot(msk, dots = FALSE, meshcol = meshcol, border = border, ...)
            }
            else {
                args <- list (x = msk, dots = FALSE, 
                    meshcol = meshcol, covariate = 'kernelp', border = border)
                args <- c(args, dots)
                do.call(plot, args) 
            }
            
            centrecell <- subset(msk, (x$x^2 + x$y^2)<1e-6)
            # 2021-03-17 modify plot call to allow add in ...
            dots$x <- centrecell
            dots$dots <- FALSE
            dots$meshcol <- 'black'
            dots$col <- NA
            dots$add = TRUE
            do.call(plot, dots)
            
            if (movementmodel %in% c('annular','annular2')) {
                rad <- k2 * spacing
                symbols(0, 0, circles = rad, inches = FALSE, add = TRUE)
                if (movementmodel == 'annular2') {
                    rad <- k2 * spacing
                    symbols(0, 0, circles = rad, inches = FALSE, add = TRUE)
                }
            }
            
            if (movementmodel %in% c('annularR')) {
                rad <- move.b * k2 * spacing
                symbols(0, 0, circles = rad, inches = FALSE, add = TRUE)
            }
            
            if (contour) {
                kp <- kernelp
                z <- matrix(nrow = 2*k2+1, ncol = 2*k2+1)
                kxy <- as.matrix(x) / spacing + k2 + 1
                z[kxy] <- kp
                if (is.null(levels))
                    levels <- pretty(c(0, max(kp, na.rm = TRUE)), 10)
                contour(add = TRUE, (-k2:k2)*spacing, (-k2:k2) * spacing, z, levels = levels)
            }
            if (text) text(x$x, x$y, round(kernelp,3), cex=0.6)
            if (is.null(title)) {
                title <-  paste0('kernel = ', movementmodel, ', spacing = ', spacing, 
                    ', kernelradius = ', k2, ', ', pstring, ', ncells = ', nrow(x))
            }
            mtext (side=3, line = 0.6, title, cex = par()$cex.main)
            invisible(data.frame(x, kernelp=kernelp))
        }
    }
}
################################################################################

expected.d <- function(movementmodel, move.a, move.b, truncate = Inf, 
    mask = NULL, min.d = 1e-4, ...) {
    # Input may be :
    #    fitted openCR model
    #    kernel object
    #    user kernel function g(r)
    #    character name of kernel model

    # build kernel from fitted openCR model
    if (inherits(movementmodel, 'openCR')) {
        if (is.null(movementmodel$movementmodel))  {
            movementmodel <- 'static'
        }
        else if (movementmodel$movementmodel %in% 'static')  {
            movementmodel <- movementmodel$movementmodel
        }
        else if (movementmodel$movementmodel %in% c('IND','INDzi')) {
            if (movementmodel$movementmodel %in% c('INDzi')) {
                k <- make.kernel(movementmodel, ...)
                move.a <- attr(k, 'move.a')
            }
            mask <- movementmodel$mask
            movementmodel <- movementmodel$movementmodel
        }
        else {
            movementmodel <- make.kernel(movementmodel, ...)
        }
    }
    # extract parameters
    if (inherits(movementmodel, 'kernel')) {
        k <- movementmodel
        move.a <- attr(k, 'move.a')
        move.b <- attr(k, 'move.b')
        k2 <- attr(k, 'k2')
        movementmodel <- attr(k, 'movementmodel')
    }
    
    # mean not available
    if (!is.function(movementmodel)) {
        if (movementmodel %in% c('user')){ 
            # warning ('expected distance not available for model : ', movementmodel)
            return(NA)
        }
        
        if (movementmodel %in% 'static') {
            return(0)
        }
        
        # check remaining movement models recognised 
        if (!(movementmodel %in% .openCRstuff$movementmodels)) {
            warning ("unrecognised movement model")
            return(NA)
        }
        
        movementmodel <- stdmovement(movementmodel)
        if (movementmodel %in% c('IND','INDzi')) {
            if (missing(move.a)) move.a <- 0    # zero inflation
            return(mean(as.matrix(dist(mask))) * (1-move.a))
        }
    }
    #--------------------------------------------------------------------------
    # force answer by integration over (0,truncate)
    if (is.finite(truncate) || is.function(movementmodel)) {
        
        if (is.function(movementmodel)) {
            # user function is g(r), so inflate by 2 pi r
            #----------------------------------------------------------------
            # adjustment for truncation added 2021-11-10
            integrand0 <- function(r) 2 * pi * r * movementmodel(r, move.a, move.b)
            if (is.finite(integrand0(0)))
                ptrunc <- try(integrate(integrand0, 0, truncate), silent = TRUE)
            else 
                ptrunc <- try(integrate(integrand0, min.d, truncate), silent = TRUE)
            if (inherits(ptrunc, 'try-error'))
                return(NA)
            else
                ptrunc <- ptrunc$value
            #----------------------------------------------------------------
            integrand <- function(r) 2 * pi * r^2 * movementmodel(r, move.a, move.b)/ptrunc
            zeroinflated <- FALSE
        }
        else {
            zeroinflated <- grepl('zi', movementmodel)
            if (zeroinflated) movementmodel <- gsub("zi","",movementmodel)
            integrand <- function(r) r * dkernel(r, movementmodel, move.a, move.b, 
                truncate = truncate)
        }
        if (is.finite(integrand(0)))
            integ <- try(integrate(integrand, 0, truncate), silent = TRUE)
        else 
            integ <- try(integrate(integrand, min.d, truncate), silent = TRUE)
        if (inherits(integ, 'try-error'))
            return(NA)
        else {

            if (zeroinflated) {
                bz <- if (movementmodel %in% c('IND','UNI')) move.a else move.b
                if (bz>1) stop ("requested zero-inflation outside range 0-1")
            }
            else bz <- 0
            return(integ$value * (1-bz))
        }
    }
    #--------------------------------------------------------------------------
    # otherwise use formula for untruncated expected distance
    if (movementmodel %in% c('BVN'))
        move.a * (pi/2)^0.5         ## Nathan et al 2012 Table 15.1
    else if (movementmodel %in% c('BVE'))
        move.a * 2              ## Nathan et al 2012 Table 15.1
    else if (movementmodel %in% c('BVT')) {
        if (missing(move.b)) stop ("BVT model has two parameters")
        a <- move.a
        b <- move.b + 1
        if (b<3/2)
            Inf
        else
            a * pi^0.5 / 2 * exp(lgamma(b-3/2) - lgamma(b-1)) 
    }
    else if (movementmodel == 'annular') {
        (1-move.a) * truncate
    }
    else if (movementmodel == 'annular2') {
        move.b * truncate/2 + (1-move.a-move.b) * truncate
    }
    else if (movementmodel == 'annularR') {
        (1-move.a) * move.b
    }
    else if (movementmodel %in% c('RDE')) {
        move.a
    }
    else if (movementmodel %in% c('RDG')) {
        move.a * move.b
    }
    else if (movementmodel %in% c('RDL')) {
        mu <- log(move.a)
        s <- sqrt(log(1 + 1/move.b))
        exp(mu + s^2/2)         ## Cousens et al 2008 Table 5.1; 
        ## cf Nathan et al 2012 Table 15.1 use a = exp(mu)
    }
    else if (movementmodel == 'BVNzi') {
        move.a * (pi/2)^0.5 * (1-move.b)
    }
    else if (movementmodel == 'BVEzi') {
        move.a * 2 * (1-move.b)
    }
    else if (movementmodel == 'BVN2') {
        0.5 * (move.a + move.b) * (pi/2)^0.5    # equal parts
    }
    else if (movementmodel %in% c('RDLS')) {
        NA                                     # log-sech
    }
    else if (movementmodel %in% c('BVC', 'BVCzi')) {
        NA                                     # bivariate Cauchy
    }
    else if (movementmodel %in% c('UNI')) {
        NA
    }
    else if (movementmodel %in% c('UNIzi')) {
        NA * (1-move.a)
    }
    else if (movementmodel %in% c('RDEzi')) {
        move.a * (1-move.b)
    }
    else NA
}
################################################################################

bR <- function(R, movementmodel, move.a, move.b) {
    if (movementmodel=='user') return (NA)  # 2020-10-06
    1 - pkernel(R, movementmodel, move.a, move.b)
}

summary.kernel <- function (object, ...) {
    spacing <- spacing(object)
    k2 <- attr(object, 'k2')
    r <- sqrt(apply(object^2,1,sum))
    movementmodel <- attr(object, 'movementmodel')
    move.a <- attr(object, 'move.a')
    move.b <- attr(object, 'move.b')
    kernelp <- covariates(object)$kernelp
    kernelp <- kernelp / sum(kernelp)   # force normalization for E(r)
    result <- list(
        k2 = k2,
        spacing = spacing,
        ncells = nrow(object),
        movementmodel = movementmodel,
        move.a = move.a[1],
        move.b = move.b[1],
        mu = if (movementmodel=='RDL') log(move.a[1]) else NA,
        s  = if (movementmodel=='RDL') sqrt(log(1/move.b[1] + 1)) else NA,
        expectedmove = expected.d(movementmodel, move.a, move.b),
        expectedmovetr = expected.d(movementmodel, move.a, move.b, k2*spacing), 
        expectedmoveemp = sum(r * kernelp),
        ptruncated = bR((k2+0.5)*spacing, movementmodel, move.a, move.b),
        expectedq50 = qkernel(0.5, movementmodel, move.a, move.b),
        expectedq90 = qkernel(0.9, movementmodel, move.a, move.b),
        expectedq50tr = qkernel(0.5, movementmodel, move.a, move.b, truncate = k2*spacing),
        expectedq90tr = qkernel(0.9, movementmodel, move.a, move.b, truncate = k2*spacing)
    )
    class(result) <- 'summary.kernel'
    result
}

print.summary.kernel <- function(x,...) {
    cat('Kernel radius (cells)     : ', x$k2, '\n')
    cat('Spacing (side of cell)    : ', x$spacing, ' (m)\n')
    cat('Number of cells           : ', x$ncells, '\n')
    if (is.function(x$movementmodel))
        cat('Movement model            : user function\n')
    else
        cat('Movement model            : ', x$movementmodel, '\n')
    cat('Parameter(s)              : ', getpstring(x$movementmodel, x$move.a, x$move.b, sep=','), '\n')
    cat('Proportion truncated      : ', x$ptruncated, '\n')
    cat('Movement as truncated at edge of kernel\n')
    cat('Empirical mean distance   : ', x$expectedmoveemp, ' (m)\n')
    cat('Expected distance         : ', x$expectedmovetr, ' (m)\n')
    cat('50th percentile (median)  : ', x$expectedq50tr, ' (m)\n')
    cat('90th percentile           : ', x$expectedq90tr, ' (m)\n')
    cat('Movement, untruncated kernel\n')
    cat('Expected distance         : ', x$expectedmove, ' (m)\n')
    cat('50th percentile (median)  : ', x$expectedq50, ' (m)\n')
    cat('90th percentile           : ', x$expectedq90, ' (m)\n')
}

# useful functions not exported 2021-02-25
# see kerneltest.R

empirical <- function (kernel) {
    r <- apply(kernel^2,1,sum)^0.5
    p <- covariates(kernel)$kernelp
    p <- p[order(r)]
    r <- sort(r)
    out <- cbind(r, cumsum(p))
}

################################################################################
# kernel function (2-D)

gkernel <- function (r, movementmodel = c('BVN', 'BVE', 'BVC', 'BVT','RDE', 'RDG', 'RDL'),
    move.a, move.b, truncate = Inf) {
    movementmodel <- try(match.arg(movementmodel, choices = 
            .openCRstuff$kernelmodels), silent = TRUE)
    if (inherits(movementmodel, 'try-error')) {
        warning("gkernel defined only for movement models ", .openCRstuff$stdmodels)
        rep(NA, length(r))
    }
    else {
        movementmodel <- stdmovement(movementmodel)   # use standard codes
        if (movementmodel %in% c('BVNzi', 'BVEzi','BVCzi','UNIzi','RDEzi')) {
            warning("gkernel not defined for zero-inflated models BVNzi, BVEzi, BVCzi, UNIzi, RDEzi")
            rep(NA, length(r))
        }
        else {
            if (movementmodel %in% c('BVN')) {
                g <- 1/move.a^2/2/pi * exp(-r^2/2/move.a^2)
            }    
            else if (movementmodel %in% c('BVE')) {
                g <- 1/move.a^2/2/pi * exp(-r/move.a)
            }
            else if (movementmodel %in% c('BVC')) {  # a>0
                g <- 1 / (2 * pi) * move.a / (r^2 + move.a^2)^(3/2)       
            }
            else if (movementmodel %in% c('BVT')) {
                g <- move.b / pi/ move.a^2 * (1 + r^2/move.a^2)^-(move.b+1)
            }
            else if (movementmodel %in% c('RDE')) {
                g <- exp(-r/move.a) / move.a / pi / 2 / r
            }
            else if (movementmodel %in% c('RDG')) {
                g <- exp(-r/move.a) / gamma(move.b) / move.a^move.b / pi / 2 * r^(move.b-2)
            }
            else if (movementmodel %in% c('RDL')) {
                mu <- log(move.a)
                s <- sqrt(log(1 + 1/move.b))
                g <- dlnorm(r, mu, s) / 2 / pi / r 
            }
            else if (movementmodel %in% c('UNI','UNIzi')) {
                g <- rep(1,length(r)) / truncate^2 /pi
            }
            else if (movementmodel %in% c('BVN2')) {
                g <- 0.5/move.a^2/2/pi * exp(-r^2/2/move.a^2) +
                    0.5/move.b^2/2/pi * exp(-r^2/2/move.b^2)
            }    
            else if (movementmodel %in% c('RDLS')) {  # a>0, b>0, r>0
                g <- 2 / (pi * r * move.b) / ((r/move.a)^(1/move.b) + (r/move.a)^(-1/move.b)) / 2 / pi / r       
            }
            else {
                g <- rep(NA, length(r))
                warning('pdf not available for movement kernel "', 
                    movementmodel, '"')
            }
            if (truncate != Inf) {
                if (movementmodel %in% c('UNI','UNIzi'))
                    ptrunc <- 1
                else
                    ptrunc <- pkernel(truncate, movementmodel, move.a, move.b, Inf, TRUE)
                g <- g / ptrunc
                g[r>truncate] <- 0
            }
            g
        }
    }
}
################################################################################

# Probability density of distance moved
dkernel <- function (r, movementmodel = c('BVN', 'BVE', 'BVC', 'BVT', 'RDE', 'RDG', 'RDL'),
    move.a, move.b, truncate = Inf) {
    movementmodel <- try(match.arg(movementmodel, choices = 
            .openCRstuff$kernelmodels), silent = TRUE)
    if (inherits(movementmodel, 'try-error')) {
        warning("dkernel defined only for movement models ", .openCRstuff$stdmodels)
        rep(NA, length(r))
    }
    else {
        movementmodel <- stdmovement(movementmodel)   # use standard codes
        if (movementmodel %in% c('BVNzi', 'BVEzi','BVCzi','RDEzi','UNIzi')) {
            warning("dkernel not defined for zero-inflated models BVNzi, BVEzi, BVCzi, UNIzi, RDEzi")
            rep(NA, length(r))
        }
        else {
            if (movementmodel %in% c('BVN')) {
                d <- r/move.a^2 * exp(-r^2/2/move.a^2)
            }    
            else if (movementmodel %in% c('BVE')) {
                d <- r/move.a^2 * exp(-r/move.a)
            }
            else if (movementmodel %in% c('BVT')) {
                d <- 2 * move.b * r / move.a^2 * (1 + r^2/move.a^2)^-(move.b+1)
            }
            else if (movementmodel %in% c('RDE')) {
                d <- exp(-r/move.a) / move.a 
            }
            else if (movementmodel %in% c('RDG')) {
                d <- exp(-r/move.a) / gamma(move.b) / move.a^move.b * r^(move.b-1)
            }
            else if (movementmodel %in% c('RDL')) {
                mu <- log(move.a)
                s <- sqrt(log(1 + 1/move.b))
                d <- dlnorm(r, mu, s) 
            }
            else if (movementmodel %in% c('UNI')) {
                A <- pi * truncate^2
                d <- 2 * pi * r / A
            }
            else if (movementmodel %in% c('BVN2')) {
                d <- 0.5 * r/move.a^2 * exp(-r^2/2/move.a^2) +
                    0.5 * r/move.b^2 * exp(-r^2/2/move.b^2)
            }    
            else if (movementmodel %in% c('RDLS')) {
                a <- log(move.a)
                b <- log(move.b)
                #d <- 2 / (pi * r * move.b) / ((r/move.a)^(1/move.b) + (r/move.a)^(-1/move.b))
                d <- 2 / (pi * r * b) / ((r*exp(-a))^(1/b) + (r*exp(-a))^(-1/b))
            }
            else if (movementmodel %in% c('BVC')) {
                # d <- 2 * r * move.a / (r^2 + move.a^2)^(3/2)
                d <- r * move.a / (r^2 + move.a^2)^(3/2)
            }
            else {
                d <- rep(NA, length(r))
                warning('pdf not available for movement kernel "', 
                    movementmodel, '"')
            }
            if (truncate != Inf) {
                if (movementmodel %in% c('UNI'))
                    ptrunc <- 1
                else
                    # using lower tail to normalise truncated distribution
                    ptrunc <- pkernel(truncate, movementmodel, move.a, move.b, Inf, TRUE)
                if (ptrunc == 0) {
                    warning(movementmodel, " numerical error: mass of truncated distribution = 0")
                    ptrunc <- NA
                }
                d <- d / ptrunc
                d[r>truncate] <- 0
            }
            d
        }
    }
}
################################################################################

# Distribution function for distance moved
pkernel <- function (q, movementmodel = c('BVN', 'BVE', 'BVC', 'BVT','RDE', 'RDG', 'RDL'),
    move.a, move.b, truncate = Inf, lower.tail = TRUE) {
    movementmodel <- try(match.arg(movementmodel, choices = 
            .openCRstuff$kernelmodels), silent = TRUE)
    #---------------------------------------------------------------------------
    # strip "zi", record for later
    zeroinflated <- grepl('zi', movementmodel)
    if (zeroinflated) movementmodel <- gsub("zi","",movementmodel)
    #---------------------------------------------------------------------------
    if (inherits(movementmodel, 'try-error')) {
        warning("pkernel defined only for movement models ", .openCRstuff$stdmodels)
        rep(NA, length(q))
    }
    else {
        movementmodel <- stdmovement(movementmodel)   # use standard codes
        # start with upper tail area
        if (movementmodel %in% c('BVN')) {
            p <- exp(-q^2/2/move.a^2) 
        }    
        else if (movementmodel %in% c('BVE')) {
            p <- (q/move.a+1) * exp(-q/move.a)
        }
        else if (movementmodel %in% c('BVT')) {
            p <- (move.a^2/(move.a^2+q^2))^move.b
        }
        else if (movementmodel %in% c('RDE')) {
            p <- exp(-q/move.a)
        }
        else if (movementmodel %in% c('RDG')) {
            p <- pgamma(q, shape = move.b, scale = move.a, lower.tail = FALSE)
        }
        else if (movementmodel %in% c('RDL')) {
            mu <- log(move.a)
            s <- sqrt(log(1 + 1/move.b))
            p <- pnorm((log(q) - mu) / s, lower.tail = FALSE)
        }
        else if (movementmodel %in% c('UNI')) {
            p <- ifelse (q>truncate, 0, 1 - (q/truncate)^2)
        }
        else if (movementmodel %in% c('BVN2')) {
            p <- (exp(-q^2/2/move.a^2) + exp(-q^2/2/move.b^2)) / 2 
        }    
        else if (movementmodel %in% c('RDLS')) {
            mu <- log(move.a)
            s <- log(move.b)
            p <- 1 - 2/pi * atan(exp((log(q)-mu)/s))
            p[q==0] <- 0
        }
        else if (movementmodel %in% c('BVC')) {
            p <- move.a / sqrt(move.a^2 + q^2)
        }
        else {
            return(rep(NA, length(q)))
            warning('probability not available for movement kernel "', 
                movementmodel, '"')
        }
        # take complement for lower tail
        if (lower.tail) {
            p <- 1-p
            p[q==0] <- 0
        }
        if (truncate != Inf) {
            if (!lower.tail) stop ("truncation incompatible with upper tail")
            if (!(movementmodel %in% c('UNI'))) {
                plower <- pkernel(truncate, movementmodel, move.a, move.b, 
                    truncate = Inf, lower.tail = TRUE)
                p <- p / plower
            }
            p[q>truncate] <- 1.0
        }
        if (zeroinflated) {
            bz <- if (movementmodel %in% c('UNI','IND')) move.a else move.b
            if (bz>1) stop ("requested zero-inflation outside range 0-1")
            if (lower.tail)
                p <- bz + (1-bz) * p
            else
                p <- 1 - (bz + (1-bz)*(1-p))
        }
        p
    }
}
################################################################################

# quantile function for distance moved
qkernel <- function(p, movementmodel = c('BVN', 'BVE', 'BVC', 'BVT','RDE', 'RDG', 'RDL'),
    move.a, move.b, truncate = Inf, lower.tail = TRUE) {
    movementmodel <- try(match.arg(movementmodel, choices = 
            .openCRstuff$kernelmodels), silent = TRUE)
    zeroinflated <- grepl('zi', movementmodel)
    if (zeroinflated) movementmodel <- gsub("zi","",movementmodel)
    if (inherits(movementmodel, 'try-error')) {
        warning("qkernel defined only for movement models ", .openCRstuff$stdmodels)
        rep(NA, length(p))
    }
    else {
        movementmodel <- stdmovement(movementmodel)   # use standard codes
        if (zeroinflated) {
            bz <- if (movementmodel %in% c('UNI')) move.a else move.b
            p0 <- p
            if (lower.tail)
                p <- (p-bz)/(1-bz)
            else
                p <- 1 - ((1-p)-bz)/(1-bz)
        }
        else {
            bz <- 0
            p0 <- 1
        }
        
        if (truncate != Inf) {
            if (!lower.tail) stop ("cannot combine truncation and upper tail")
            ptrunc <- pkernel(truncate, movementmodel, move.a, move.b, Inf, TRUE)
            p <- p * ptrunc
        }
        if (lower.tail) {
            p <- 1-p
        }
        if (movementmodel %in% c('BVN')) {
            q <- sqrt(-log(p)*2*move.a^2)
        }    
        else if (movementmodel %in% c('BVE')) {
            # p <- (q/move.a+1) * exp(-q/move.a)
            # Doesn't work W <- VGAM::lambertW(-p/exp(1))
            # q <- (-W - 1) * move.a
            onep <- function(p) {
                onem <- function(m) {
                    # fn <- function (q, p) (q/m+1) * exp(-q/m) - p
                    fn <- function (q, p) 1 - (q/m+1) * exp(-q/m) - p
                    out <- try(uniroot(f = fn, interval = c(m*1e-3, m*1e3), p = p)$root, silent = TRUE)
                    if (inherits(out, 'try-error')) NA else out
                }
                sapply(move.a, onem)
            }
            q <- sapply(1-p, onep)   # adjusted for tail 2021-07-13
        }
        else if (movementmodel %in% c('BVT')) {
            q <- move.a * sqrt(p^(-1/move.b) - 1)
        }
        else if (movementmodel %in% c('RDE')) {
            q <- qexp(p, 1/move.a, lower.tail = FALSE)
        }
        else if (movementmodel %in% c('RDG')) {
            q <- qgamma (p, shape = move.b, scale = move.a, lower.tail = FALSE)
        }
        else if (movementmodel %in% c('RDL')) {
            mu <- log(move.a)
            s <- sqrt(log(1 + 1/move.b))
            q <- qlnorm(p, meanlog = mu, sdlog = s, lower.tail = FALSE)
        }
        else if (movementmodel %in% c('UNI')) {
            q <- truncate * sqrt(p)
        }
        else if (movementmodel %in% c('BVN2')) {
            q <- (sqrt(-log(p)*2*move.a^2) + sqrt(-log(p)*2*move.b^2))/2
        }    
        else if (movementmodel %in% c('RDLS')) {
            mu <- log(move.a)
            s <- log(move.b)
            q <- exp(2 * s / pi * log(tan(pi/2 * p)) + mu)
        }
        else if (movementmodel %in% c('BVC')) {
            q <- move.a * sqrt(p^(-2) - 1)
        }        
        else {
            q <- rep(NA, length(p))
            warning('probability not available for movement kernel "', 
                movementmodel, '"')
        }
        q[q>truncate] <- NA
        q[bz>p0] <- 1-lower.tail
        q    
    }
}

matchscale <- function(movementmodel, q = 40, expected = NULL, p = 0.5, lower = 1e-5, upper = 1e5, 
    move.b = 1, truncate = Inf) {
    if (!is.null(expected)) {
        mfn <- function(move.a) expected.d (movementmodel, move.a, move.b, truncate)-expected
    }
    else {
        mfn <- function(move.a) pkernel (q, movementmodel, move.a, move.b, truncate)-p
    }
    movementmodel <- stdmovement(movementmodel)   # use standard codes
    if (movementmodel %in% c(
        'BVN', 'BVNzi', 'BVN2',
        'BVE','BVEzi',
        'BVC', 'BVCzi',
        'BVT',
        'RDE','RDEzi',
        'RDG',
        'RDL',
        'RDLS'
        )) {
        out <- try(uniroot(mfn, c(lower,upper))$root, silent = TRUE)
    }    
    else if (movementmodel %in% c('UNIzi')) {
        upper <- min(upper, 1.0)  # do not exceed feasible limit
        if (is.finite(truncate))
            out <- try(uniroot(mfn, c(lower,upper))$root, silent = TRUE)
        else
            stop("UNIzi kernel must be truncated at finite radius")
    }
    else if (movementmodel %in% c('UNI')) {
        if (!is.null(expected)) {
            mfnt <- function(truncate) expected.d ('UNI', truncate = truncate)-expected
        }
        else {
            mfnt <- function(truncate) pkernel (q, 'UNI', truncate = truncate)-p
        }
        out <- try(uniroot(mfnt, c(lower,upper))$root, silent = TRUE)
    }
    else {
        stop(movementmodel, " not available in matchscale()")
    }
    if (inherits(out, 'try-error')) out <- NA
    out
}
