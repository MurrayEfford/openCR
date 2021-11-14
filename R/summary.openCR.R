# summary.openCR.R

summary.openCR <- function (object, newdata = NULL, alpha = 0.05, svtol = 1e-5, deriv = FALSE, ...) {
    
    secrmodel <- grepl("secr", object$type)
    CLtype <- grepl("CL", object$type) || grepl("PLB", object$type) || grepl("CJS", object$type) 
    ch <- object$capthist
    det <- detector(traps(ch))[1]
    out <- vector('list')
    
    # cl <- paste(names(object$call)[-1],object$call[-1], sep=' = ', collapse=', ' )
    # cl <- paste('openCR.fit(', cl, ')')
    # out$call <- strwrap(cl, getOption('width'))
    
    out$versiontime <- paste0(object$version, ', run ', object$starttime)
    if (!is.null(object$details$newdetector)) {
        out$newdetector <- object$details$newdetector
    }

    ###################
    ## Data description
    
    if (secrmodel) {
        trp <- traps(ch)
        out$traps <- data.frame (Detector = detector(trp)[1],
                                 Number = nrow(trp),
                                 Spacing = spacing(trp))
        rownames(out$traps) <- ''
        if (!is.null(usage(trp)))
            out$traps$UsagePct <- 100 * sum(usage(trp))/length(usage(trp))
        if (length(detector(trp))>1)
            out$detector <- detector(trp)
    }
    freq <- covariates(ch)$freq 
    if (is.null(freq)) freq <- rep(1,nrow(ch))
    n  <- sum(freq)                 # number caught
    ncapt <- sum(freq * apply(abs(ch),1,sum))
    nprimary <- length(object$intervals)+1
    nsecondary <- ncol(ch)
    ntrp <- nrow(traps(ch))

    ch <- unsqueeze(ch)   ## 2018-11-10    
    chsess <- suppressWarnings(split(ch, primarysessions(intervals(ch)), byoccasion = TRUE))
    bysession <- summary(chsess, terse = TRUE, moves = TRUE)
    if (is.null(ntrp)) {
        bysession <- bysession[1:3,]
        nmov <- numeric(0)
    }
    else {
        nmov <- if (nrow(bysession>4)) sum(unlist(sapply(moves(ch), function(y) y>0))) else numeric(0)
    }
    out$capthist <- cbind(bysession, Total = c(nsecondary, ncapt, n, ntrp, nmov))
    out$intervals <- object$intervals
    
    if (secrmodel) {
        out$mask <- data.frame(Cells = nrow(object$mask), Spacing = spacing(object$mask))
        if (length(maskarea(object$mask))==0)
            out$mask <- cbind(out$mask, Length = masklength(object$mask))
        else
            out$mask <- cbind(out$mask, Area = maskarea(object$mask))
    }
    out$modeldetails <- data.frame(type = object$type,
                                   fixed = fixed.string(object$fixed),
                                   distribution = if (!CLtype) object$distribution else 'none')
    if (secrmodel) {
        if (any(det %in% .openCRstuff$countdetectors)) {
            out$Countmodel <- if (object$binomN == 0) 'Poisson'
            else if (object$binomN == 1) 'Binomial, size from usage'
            else if (object$binomN > 1) paste('Binomial', object$binomN)
        }
        out$Movementmodel <- object$movementmodel
    }
    out$AICtable <- AIC(object)
    out$link <- data.frame(object$link)
    out$coef <- coef(object)
    
    if (!is.null(object$fit$hessian)) {
        out$hessian <- data.frame (rankH = length(which(object$eigH > svtol)),
                             svtol = svtol)
        out$hessian <- cbind(out$hessian, as.list(round(object$eigH, -log10(svtol))))
        neigen <- length(object$eigH)
        names(out$hessian)[3:(2+neigen)] <- paste0('Eigen', 1:neigen)
    }    
    out$predicted <- predict (object, newdata, alpha = alpha)
    
    #################################
    # Derived parameters
    #################################
    if (deriv) {
        out$derived <- derived(object, alpha=alpha, ...)
    }
    
    ## remove distracting row names
    for (i in 1:length(out)) {
        if (is.data.frame(out[[i]]))
            if (nrow(out[[i]])==1 & (!names(out)[i] %in% c('coef')))
                rownames (out[[i]]) <- ''
    }
    class (out) <- "summary.openCR"
    out
}
############################################################################################

print.summary.openCR <- function(x,...) {
    class(x) <- NULL
    attr(x,'fit') <- NULL
    attr(x,'eigH') <- NULL
    print(x)
}

############################################################################################

AIC.summary.openCR <- function (object, ..., sort = TRUE, k = 2, dmax = 10, 
                              criterion = c('AICc','AIC')) {
    ## identical to AIC.summary.secr
    criterion <- match.arg(criterion)
    allargs <- list(object, ...)
    output <- do.call(rbind, lapply(allargs, '[[', "AICtable"))
    rownames(output) <- NULL
    output$delta <- output[,criterion] - min(output[, criterion])
    ## optional sort
    if (sort) output <- output[order(output$delta),]
    ## AICwt with dmax
    OK <- abs(output$delta) < abs(dmax)
    sumdelta <- sum(exp(-output$delta[OK]/2))
    output$wt <- ifelse ( OK, round(exp(-output$delta/2) / sumdelta,4), 0)
    names(output)[7] <- paste('d',criterion,sep='')
    names(output)[8] <- paste(criterion,'wt',sep='')
    if (nrow(output)==1) { output[,8] <- NULL; output[,7] <- NULL}
    output
}
############################################################################################
