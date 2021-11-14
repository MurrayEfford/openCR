############################################################################################
## package 'openCR'
## openCR.make.newdata.R
## 2011 12 09
## Create (neutral) design data suitable for 'predict'
## 2015-02-06 reconciled this current version with forked 1.2.0:
## 2017-12 revamped
## 2018-04-12 allow single session
## 2018-11-22 new learned responses
## 2019-02-02 fixed bug: factor(0,1)
## 2020-10-19 agecov
## 2020-12-07 tt occasion-level time variation cf Kendall et al. 1997
## 2021-04-25 2.0.0 stratified
## 2021-05-12 fixed bug in stratified sessioncov
## 2021-07-02 fixed backward incompatibility bug details$minimumage not specified
## 2021-07-30 makeNewData method for openCR objects
############################################################################################

makeNewData.openCR <- function (object, all.levels = FALSE, ...) {
# openCR.make.newdata <- function (object, all.levels = FALSE, ...) {
        
    # 'Session', 't' are handled separately at end
    autovars <- c(.openCRstuff$learnedresponses, 'stratum', 'session','tt', 'h2','h3')
    capthist <- object$capthist
    mask <- object$mask
    
    vars <- object$vars
    
    dframe <- object$dframe
    stratanames <- factor(strata(capthist))
    nstrata <- length(stratanames)
    J <- sapply(primaryintervals(object), length)+1
    S <- if(ms(capthist)) sapply(capthist, ncol) else ncol(capthist)
    
    # fix backward compatibility bug 2021-07-02
    if (is.null(object$details$minimumage)) object$details$minimumage <- 0
    if (is.null(object$details$maximumage)) object$details$maximumage <- 1
    
    agerange <- object$details$minimumage:object$details$maximumage
    sessioncov <- stdcovlist(object$sessioncov, 'scov', nstrata, J)
    timecov    <- stdcovlist(object$timecov, 'tcov', nstrata, S)
    agecov     <- stdcovlist(object$agecov, 'acov', nstrata, diff(agerange) + 1)
    stratumcov <- stdcovlist(object$stratumcov, 'stratumcov', 1, NULL)
    
    nmix <- object$details$nmix
    if(is.null(nmix)) nmix <- 1
    mixvar <- switch(nmix, character(0),'h2','h3')
    
    #############################################################
    onestratum <- function(stratum) {
        findvars <- function (basevars, cov) {
            ## function to add covariates to a list
            ## cov should be dataframe or list of dataframes, one per stratum (R > 1),
            if (!is.data.frame(cov)) cov <- cov[[stratum]] ## assume multisession list
            if (is.null(cov) | (length(cov)==0) | (length(stratumvars)==0)) return(basevars)
            else {
                found <- ''
                for (v in stratumvars) {
                    if (v %in% names(cov)) {
                        vals <- cov[,v]
                        if (is.character(vals)) vals <- factor(vals)
                        basevars[[v]] <- if (is.factor(vals))
                            factor(levels(vals), levels = levels(vals))
                        else
                            unique(vals)
                        found <- c(found, v)
                    }
                }
                stratumvars <<- stratumvars[!(stratumvars %in% found)]
                return(basevars)
            }
        }
        if (nstrata>1) {
            capthist <- capthist[[stratum]]
            mask <- mask[[stratum]]
        }
        interv <- intervals(capthist)
        
        stratumvars <- vars
        
        # single stratum label, levels of factor apply to whole
        basevars <- list(stratum = factor(stratanames[stratum], levels=stratanames))

        # use either session or tt
        if ('tt' %in% vars) {
            basevars$tt <- factor(1:(length(interv)+1))
        }
        else {
            basevars$session <- factor(1:J[stratum])
        }
        
        mixvar <- 'h2'   ## stop gap 2018-01-22
        if (nmix>1) basevars[mixvar] <- list(as.character(1:nmix))
        for (v in stratumvars) {
            if (v=='T')  basevars$T <- 0
            for (i in .openCRstuff$learnedresponses) {
                if (v == i) basevars[[i]] <- factor(0:1)  
            }
            
            if (v=='age')  basevars$age <- factor(agerange)
            if (v=='Age')  basevars$Age <- agerange
            if (v=='Age2')  basevars$Age2 <- agerange^2
        }
        ## all autovars should now have been dealt with
        stratumvars <- stratumvars[!stratumvars %in% autovars]
        basevars <- findvars (basevars, covariates(capthist)) ## individual covariates
        
        basevars <- findvars (basevars, covariates(traps(capthist)))
        basevars <- findvars (basevars, covariates(mask))
        basevars <- findvars (basevars, timecov)
        basevars <- findvars (basevars, agecov)    ## 2020-10-19
        basevars <- findvars (basevars, covariates(traps(capthist)))
        if (!is.null(mask))
            basevars <- findvars (basevars, covariates(mask))
        if (!is.null(dframe))
            basevars <- findvars (basevars, dframe)
        
        ## revert to first level
        for (v in names(basevars)) {
            if (!all.levels & !(v %in% c('stratum', 'session', 'tt', 'h2','h3'))) {
                basevars[[v]] <- basevars[[v]][1] 
            }
        }
        basevars <- lapply(basevars, function(x) if (is.character(x)) factor(x) else x)
        out <- expand.grid(basevars)
        
        if (!'session' %in% names(out)) {
            out <- data.frame(
                stratum = out$stratum, 
                session = factor(primarysessions(interv)[as.numeric(out$tt)]),
                out[,-1, drop = FALSE])
        }
        if (!is.null(sessioncov)) {
            if (nstrata>1) {
                sessioncov <- sessioncov[[stratum]]   ## 2021-05-12
            }
            for (i in names(sessioncov)) {
                if ((i %in% vars) & !(i %in% names(out)))
                    out[,i] <- sessioncov[out$session,i]
            }
        }
        
        out
    }   # end one stratum
    
    
    newdata <- lapply(1:nstrata, onestratum)
    newdata <- do.call(rbind, newdata)
    
    if (!is.null(stratumcov)) {
        for (i in names(stratumcov)) {
            if ((i %in% vars) & !(i %in% names(newdata))) {
                cov <- stratumcov[newdata$stratum,i]
                if (is.character(cov)) cov <- factor(cov)
                newdata[,i] <- cov
            }
        }
    }
    
    if ('Session' %in% vars) {
        newdata$Session <- as.numeric(newdata$session) - 1
    }
    if ('t' %in% vars) { ## synonym 
        newdata$t <- newdata$session
    }
    
    newdata <- newdata[,names(newdata) %in% c('stratum','session',vars)]
    newdata
}
############################################################################################

