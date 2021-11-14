###############################################################################
## openCR
## 2018-05-12
## to do:
## extractFocal doesn't handle h2
## 2019-06-30 enabled sessionlabels
###############################################################################

extractFocal <- function(ocrlist, ...) {
    getcentre <- function(x, centre) 
        lapply(x, function(y) {
            cols <- c('estimate','SE.estimate','lcl','ucl')
            onerow <- as.numeric(y[centre,cols])
            names(onerow) <- cols
            onerow
        })
    pr <- predict(ocrlist, ...)
    
    tmp <- mapply(getcentre, pr, names(pr), SIMPLIFY = FALSE)
    out <- vector('list')
    for (parm in names(tmp[[1]])) {
        out[[parm]] <- do.call(rbind, lapply(tmp, '[[', parm))
    }
    out
}

moving.fit <- function(..., width = 3, centres = NULL, filestem = NULL, 
                       trace = FALSE, FUN = openCR.fit) {
    around <- function (ch, j) {
        ## subset capthist
        selectedocc <- primary %in% ((j-buff):(j+buff))
        newocc <- oldoccasion[selectedocc]
        # interv <- oldinterv[selectedocc]
        # interv <- interv[1:(length(interv)-1)]
        ch <- subset(ch, occasions = newocc)
        # intervals(ch) <- interv
        sessionlabels(ch) <- slabels[(j-buff):(j+buff)]
        ch
    }
    aroundscov <- function (scov, j) {
        ## session covariates for requested sessions
        selectedsess <- (j-buff):(j+buff)
        if (is.data.frame(scov))
            scov[selectedsess,,drop = FALSE]
        else
            scov[selectedsess]
    }
    runone <- function(j) {
        arg$capthist <- around(ch,j)
        arg$sessioncov <- aroundscov(arg$sessioncov,j)
        arg <- arg[names(arg) %in% c(allowedargs, 'capthist')]
        names(arg)[names(arg) == "capthist"] <- ""
        fit <- do.call(FUN, arg)
        if (openCRfit) fit$call <- ""
        if (!is.null(filestem)) {
            save(fit, file = paste0(filestem, j, ".RData"))
        }
        if (trace) {
            message("Completed session ", j)
        }
        fit
    }
    allowedargs <- c(names(formals(FUN)))
    openCRfit <- identical(FUN, openCR.fit)
    arg <- list(...)
    ch <- arg$capthist
    if (is.null(arg$type)) arg$type <- 'CJS'
    if (is.null(arg$movementmodel)) arg$movementmodel <- 'static'
    ## perform join if ms, etc.
    HPXpoly <- detector(traps(ch))[1] %in% c('polygon','polygonX') && (arg$detectfn == 'HPX')
    ch <- stdcapthist(ch, arg$type, arg$nclone, FALSE, HPXpoly, stratified = FALSE)   ## 2021-04-18 not yet adapted for stratified data
    oldoccasion <- 1:ncol(ch)
    # oldinterv <- intervals(ch)
    primary <- primarysessions(intervals(ch))
    slabels <- sessionlabels(ch)
    if (width>max(primary)) stop ('width exceeds number of primary sessions')
    if (width %% 2 != 1) stop ('width should be an odd integer')
    buff <- (width-1)/2
    if (is.null(centres)) {
        from <- width %/% 2 + 1
        to <- max(primary) - buff
        centres <- from : to
    }
    if (min(centres) < (1+buff))
        stop ("centres should be at least ", as.character(1+buff))
    if (max(centres) > (max(primary)-buff))
        stop ("centres should be no more than ", as.character(max(primary)-buff))
    out <- vector('list', length(unique(primary)))
    for (j in centres) out[[j]] <- runone(j)
    out <- out[!sapply(out, is.null)]
    if (openCRfit) out <- openCRlist(out)
    names(out) <- sessionlabels(ch)[centres]
    attr(out, "width") <- width
    out
}
