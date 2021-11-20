###############################################################################
## package 'openCR'
## mask.designdata.R
## Prepare design matrix for mask-level models (cf D.designdata.R in secr)
##
## 2021-11-20

###############################################################################
## NOTE does not standardize stratumcov, maskcov
###############################################################################

mask.designdata <- function (mask, maskmodel, stratumlevels, sessionlevels, 
    stratumcov = NULL, sessioncov = NULL, meanSD = NULL) 
{

    ## mask -- mask object or list of masks of the same length as stratumlevels
    ## maskmodel -- formula that may be constant ~1 or include
    ## any of the 'automatic' terms c('stratum','x','y','x2','y2','xy',
    ##     'session', 't', 'Session') or user-supplied mask-level covariates
    ## stratumlevels -- character vector of stratum names
    ## sessionlevels -- character vector of session names
    ## stratumcov -- dataframe of stratum-specific covariates
    ## sessioncov -- dataframe of session-specific covariates
    ## meanSD -- optional externally provided mean and SD (rows) for x-
    ##     and y-coordinates (columns)

    ## Output is a dataframe with one row for each combination of
    ## mask point, stratum and session. Conceptually, we use a rectangular
    ## array with enough rows to accommodate the largest mask, so some rows
    ## in the output may merely hold space to enable easy indexing.

    #--------------------------------------------------------------------------
    ## utility function

    ## standardise a numeric column from msk and pad with NA to desired length
    getcol <- function (msk, colnum) {
          mn <- attr(msk, "meanSD")[1, colnum]
          SD <- attr(msk, "meanSD")[2, colnum]
          pad1 (scale(msk[,colnum], mn, SD), nmaskrow)
    }

    #--------------------------------------------------------------------------
    ## setup
    vars  <- all.vars(maskmodel)
    nstrata <- length(stratumlevels)
    R     <- length(sessionlevels)
    if (ms(mask)) {
        maskrows <- sapply(mask, nrow)
        if (length(maskrows) != nstrata)
            stop ("number of masks does not match number of strata")
    }
    else {
        maskrows <- nrow(mask)
    }
    nmaskrow <- max(maskrows)
    dims  <- c(nmaskrow, nstrata, R)
    ## special case where new meanSD passed
    if (!is.null(meanSD)) {
        if (ms(mask))
            for (i in 1:length(mask))
                attr(mask[[i]], "meanSD") <- meanSD[[i]]
        else
            attr(mask, "meanSD") <- meanSD
    }

    #--------------------------------------------------------------------------
    ## coordinates
    ## might be condensed by operating on x,y together, but would it be clear?
    if (any (vars %in% c('x','y','x2','y2','xy'))) {
        if (ms(mask)) {
            ## stratum-specific masks
            x <- lapply(mask, getcol, 1)
            y <- lapply(mask, getcol, 2)
            dframe <- data.frame(
                x = as.vector(unlist(x)),
                y = as.vector(unlist(y))
            )
        }
        else {
            ## uniform mask across strata
            x <- getcol(mask,1)
            y <- getcol(mask,2)
            dframe <- data.frame( 
                x = rep(as.vector(unlist(x)), nstrata * R),
                y = rep(as.vector(unlist(y)), nstrata * R)
            )
        }
        #---------------------------------------------
        ## coordinates transformed for quadratic trend
        if ('x2' %in% vars) {
            dframe$x2 <- dframe$x^2
        }
        if ('y2' %in% vars) {
            dframe$y2 <- dframe$y^2
        }
        if ('xy' %in% vars) {
            dframe$xy <- dframe$x * dframe$y
        }
    }
    else {
        dframe <- data.frame(intercept = rep(1, nmaskrow * nstrata * R))
    }
    #--------------------------------------------------------------------------
    
    ## strata
    if ('stratum' %in% vars) {
        if (length(stratumlevels)<1)
            stop ("no strata specified")
        dframe$stratum <- insertdim(factor(stratumlevels), 2, dims)
    }
    #--------------------------------------------------------------------------
    
    ## sessions (synonym t)
    if ('session' %in% vars) {
        dframe$session <- insertdim(factor(sessionlevels, levels =
                sessionlevels), 3, dims)
    }
    if ('t' %in% vars) {
        dframe$t <- insertdim(factor(sessionlevels, levels =
                sessionlevels), 3, dims)
    }
    if ('Session' %in% vars) {
       dframe$Session <- insertdim(0:(R-1), 3, dims)
    }
    #--------------------------------------------------------------------------

    ## all autovars should have now been dealt with
    vars <- vars[!vars %in% c('stratum', 'x', 'y', 'x2', 'y2', 'xy',
        'session', 't', 'Session')]
    #--------------------------------------------------------------------------
    
    ## stratum covariates
    if (!is.null(stratumcov)) {
        stratumcov <- stringsAsFactors(stratumcov)   
        found <- names(stratumcov) %in% vars
        if (is.data.frame(stratumcov) & any(found)) {
            found <- names(stratumcov)[found]
            values <- as.data.frame(stratumcov[,found])
            names(values) <- found
            if (length(values)>0) {
                for (i in 1:ncol(values)) {
                    vals <- values[,i]
                    dframe[,found[i]] <- insertdim (vals, 2, dims)
                }
                vars <- vars[!(vars %in% found)]
            }
        }
    }
    #--------------------------------------------------------------------------
    
    ## session covariates
    if (!is.null(sessioncov)) {
        sessioncov <- stringsAsFactors(sessioncov)   
        found <- names(sessioncov) %in% vars
        if (is.data.frame(sessioncov) & any(found)) {
            found <- names(sessioncov)[found]
            values <- as.data.frame(sessioncov[,found])
            names(values) <- found
            if (length(values)>0) {
                for (i in 1:ncol(values)) {
                    vals <- values[,i]
                    dframe[,found[i]] <- insertdim (vals, 3, dims)
                }
                vars <- vars[!(vars %in% found)]
            }
        }
    }
    #--------------------------------------------------------------------------

    ## mask covariates
    maskcovs <- covariates(mask)
    if ((!is.null(maskcovs)) & (length(vars>0))) {
        std <- FALSE    ## no standardization
        ## requires all covariates in all masks
        expand <- function (df, n) {
            found <- names(df) %in% vars
            temp <- df[, found, drop = FALSE]
            ## pad with replicated row 1 to ensure rectangular
            if (nrow(temp) < n)
                rbind(temp, temp[rep(1,n-nrow(df)),,drop=FALSE])  ## drop=F added 2013-03-10
            else
                temp
        }
        if (ms(mask)) {
            maskcovs <- lapply(maskcovs, expand, nmaskrow)
            ncov <- sapply (maskcovs, ncol)
            if (any(ncov != ncov[1]))
                stop ("covariate missing in at least one session mask")
            maskcovs <- do.call(rbind, maskcovs)
        }
        else {
            maskcovs <- expand(maskcovs, nmaskrow)
        }
        for (i in names(maskcovs)) {
            vals <- maskcovs[,i] ## vector
            if (any(is.na(vals))) {
                warning(sum(is.na(vals)), " NA set to zero in mask covariate")
                vals[is.na(vals)] <- 0
            }
            if (!is.factor(vals) && !is.character(vals) && std)  
                vals <- scale(vals)
            if (ms(mask))
                ## vals already full length to match nrow(dframe)
                dframe[,i] <- vals
            else
                ## vals repeated across groups & sessions
                dframe[,i] <- insertdim (vals, 1, dims)
        }
        vars <- vars[!(vars %in% names(maskcovs))]
    }
    #--------------------------------------------------------------------------
    
    ## unmatched variables in model?
    if (length(vars) > 0) {
        stop (paste(vars,collapse=','), " not found")
    }
    #--------------------------------------------------------------------------
    
    ## report dimensions as an attribute
    attr(dframe, 'dimmaskdesign') <- c(nmaskrow, nstrata, R)
    attr(dframe, 'validMaskRows') <- maskrows
    dframe
}
###############################################################################


## test

# msk <- make.mask(make.grid(), nx=4)
# mask.designdata(msk, ~session, 1, 1:2)

