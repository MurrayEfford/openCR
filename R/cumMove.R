# cumMove.R
# 2021-06-14,16
# 2021-11-19 user may provide settlement

cumMove <- function (X, mask, kernel, edgemethod = c('truncate', 'wrap', 'none'), 
    nstep = 1, mqarray = NULL, settlecov = NULL) {
    if (ms(mask)) {
        stop ("requires single mask")
    }
    if (is.null(settlecov)) {
        settlement <- matrix(1)   # dummy for now 2021-11-21
    }
    else {
        if (!(settlecov %in% names(covariates(mask)))) {
            stop ("settlecov should be the name of a mask covariate")
        }
        if (nstep>0)
        settlement <- matrix(covariates(mask)[,settlecov], nrow=nrow(mask), ncol=nstep)
    }
    
    # cell probabilities
    kernelp <- matrix(covariates(kernel)$kernelp, nrow = nrow(kernel), 
        ncol = max(1,nstep))
    
    # integer cell positions
    # cellsize <- attr(mask,'area')^0.5 * 100   ## metres, equal mask cellsize
    cellsize <- spacing(mask)  # metres 
    kernel <- round(kernel / cellsize)
    
    # edge method will usually be 'truncate and normalise'
    edgemethod <- match.arg(edgemethod)
    edgecode <- edgemethodcode(edgemethod)  # 0 none, 1 wrap, 2 truncate
    if (attr(mask, 'type') != 'traprect' && edgemethod == 'wrap') {
        stop("edgemethod = 'wrap' requires mask of type 'traprect'")        
    }
    # generate lookup array
    if (is.null(mqarray)) {
        kernel <- linearisekernel (kernel, mask) 
        mqarray <- mqsetup (mask, kernel, cellsize, edgecode)  
    }

    # default to centre when X not specified
    if (missing(X)) X <- apply(mask, 2, mean)
    
    # initial probabilities
    if (inherits(X, 'mask') && ('pm' %in% names(covariates(X)))) {
        pm <- covariates(X)$pm
    }
    else if (inherits(X, c('SpatialPolygons', 'matrix'))) {
        pm <- pointsInPolygon(mask, X)
        pm <- pm/sum(pm)
    }
    else if (inherits(X, 'list') && all(c('x','y') %in% names(X))) {
        # polygon defined by list with x, y coordinates
        # e.g. output from buffer.contour
        X <- cbind(X$x, X$y)
        pm <- pointsInPolygon(mask, X)
        pm <- pm/sum(pm)
    }
    else {
        # isolated point
        X <- matrix(X, ncol = 2)
        index <- nearesttrap(X, as.matrix(mask))
        pm <- rep(0, nrow(mask))   
        pm[index] <- 1/nrow(X)
    }
    # iterate over steps
    pm0 <- pm
    if (nstep > 0) {
        for (j in 1:nstep) {
            pm <- convolvemqcpp(
                as.integer(j), 
                as.integer(edgecode),
                as.matrix(mqarray), 
                as.matrix(settlement),
                as.double(kernelp), 
                as.double(pm))
        }
    }
    
    # return mask object with probabilities as covariates pm0 (initial) and pm (final)
    covariates(mask) <- data.frame(pm0 = pm0, pm = pm)
    mask
}

proportionInPolygon <- function (mask, poly, cov = 'pm') {
    if (is.null(covariates(mask)) || !cov %in% names(covariates(mask))) {
        stop ("mask does not have covariate '", cov, "'")
    }
    if (inherits(poly, 'list') && !all(c('x','y') %in% names(poly))) {
        sapply(poly, proportionInPolygon, mask = mask, cov = cov)   
    }
    else {
        if (is.list(poly)) poly <- cbind(poly$x, poly$y)
        OK <- pointsInPolygon(mask, poly)
        covar <- covariates(mask)[,cov]
        sum(covar[OK], na.rm = TRUE) / sum(covar, na.rm = TRUE)
    }
}

