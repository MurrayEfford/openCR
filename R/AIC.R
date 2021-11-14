oneline.openCR <- function (openCR, svtol, use.rank, n = NULL) {
    
    if (is.null(n)) {
        if (ms(openCR$capthist)) {
            n <- sum(sapply (openCR$capthist, nrow))
        }
        else  {
            n  <- nrow(openCR$capthist)     # number caught
        }
    }
    
    Npar <- max(unlist(openCR$parindx))
    rank <- length(which(openCR$eigH > svtol))
    NP <- ifelse(use.rank, rank, Npar)
    ## allow for fixed beta parameters 2009 10 19
    if (!is.null(openCR$details$fixedbeta))
        NP <- NP - sum(!is.na(openCR$details$fixedbeta))
    
    AICval <- 2*(openCR$fit$value + NP)
    AICcval <- ifelse ((n - NP - 1)>0,
                       2*(openCR$fit$value + Npar) + 2 * NP * (NP+1) / (n - NP - 1), NA)
    c (
        model  = model.string(openCR$model, openCR$details$userDfn),
        npar   = Npar,
        rank   = rank,
        logLik = -openCR$fit$value,
        AIC    = round(AICval, 3),
        AICc   = round(AICcval, 3),
        fitted = openCR$fitted
    )
}
############################################################################################

logLik.openCR <- function(object, ...) {
    npar <- length(object$fit$par)
    structure (-object$fit$value, df = npar, class = 'logLik')
}

############################################################################################

AIC.openCR <- function (object, ..., sort = TRUE, k = 2, dmax = 10, use.rank = FALSE,
                        svtol = 1e-5, criterion = c('AIC','AICc'), n = NULL) {
    if (k != 2)
        stop ("'AIC.openCR' defined only for k = 2")
    
    allargs <- list(...)
    modelnames <- (c ( as.character(match.call(expand.dots=FALSE)$object),
                       as.character(match.call(expand.dots=FALSE)$...) ))
    
    allargs <- openCRlist(object, allargs)
    names(allargs) <- modelnames
    AIC(allargs, sort=sort, k=k, dmax=dmax, use.rank=use.rank, svtol=svtol, 
        criterion=criterion, n=n)
}
############################################################################################
############################################################################################

AIC.openCRlist <- function (object, ..., sort = TRUE, k = 2, dmax = 10, use.rank = FALSE,
                            svtol = 1e-5, criterion = c('AIC','AICc'), n = NULL) {
    
    if (k != 2)
        stop ("AIC.openCR defined only for k = 2")
    
    if (length(list(...)) > 0)
        warning ("... argument ignored in 'AIC.openCRlist'")

    criterion <- match.arg(criterion)
    if (any(duplicated(names(object)))) names(object) <- 1:length(object)
    modelnames <- names(object)
    allargs <- object
    if (any(sapply(allargs,class) != 'openCR'))
        stop ("components of 'object' must be 'openCR' objects")
    mat <- t(sapply(allargs, oneline.openCR, use.rank = use.rank,
                    svtol = svtol, n = n))
    output <- data.frame(mat, stringsAsFactors = FALSE)
    for (i in 3:6)
        output[,i] <- as.numeric(output[,i])

    output$delta <- output[,criterion] - min(output[,criterion])
    OK <- abs(output$delta) < abs(dmax)
    sumdelta <- sum(exp(-output$delta[OK]/2))
    output$wt <- ifelse ( OK, round(exp(-output$delta/2) / sumdelta,4), 0)
    row.names(output) <- modelnames
    if (sort) output <- output [order(output[,criterion]),]
    names(output)[7] <- paste('d',criterion,sep='')
    names(output)[8] <- paste(criterion,'wt',sep='')
    if (nrow(output)==1) { output[,8] <- NULL; output[,7] <- NULL}
    
    output
}
############################################################################################

