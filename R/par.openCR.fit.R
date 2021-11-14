## 2017-05-12
## openCR parallel fit, derived, region.N
## warning 2020-12-13

par.openCR.fit <- function (arglist, ncores = 1, seed = 123, trace = FALSE,
    logfile = NULL, prefix = "") {
    ptm  <- proc.time()
    ## 'inherits' causes R to search in enclosing frames
    if (is.character(arglist))
        arglist <- mget(arglist, inherits = TRUE)
    
    ## force 'trace' to common value across all components of arglist
    arglist <- lapply(arglist, function (x) {x$trace <- trace; x})

    ## check for capthist, mask, dframe mentioned by name
    ## objects are exported to the worker processes as required
    getnames <- function(obj = 'capthist') {
        tmpnames <- sapply(arglist, function(x) if (is.character(x[[obj]])) x[[obj]] else '')
        unique(tmpnames)
    }
    data <- c(getnames('capthist'), getnames('mask'),getnames('dframe'),getnames('details'))
    data <- data[nchar(data)>0]

    ## default details savecall to FALSE across all components of arglist
    arglist <- lapply(arglist, function (x) {
        if (is.null(x$details))
            x$details <- list(savecall = FALSE)
        else if (!('savecall' %in% names(x$details))) {
            x$details[['savecall']] <- FALSE
        }
        x
    })
    
    ## individual fits may use ncores > 1 from 1.5.0
    if (ncores > 1) {
        warning("par.secr.fit with ncores > 1 is SLOWER than par.secr.fit with ncores = 1")
        if (is.null(logfile)) {
            logfile <- tempfile("logfile", fileext = ".txt")
        }
        clust <- makeCluster(ncores, methods = FALSE, 
            useXDR = .Platform$endian=='big', outfile = logfile)
        clusterSetRNGStream(clust, seed)
        clusterExport(clust, c(data, 'openCR.fit'), environment())
        output <- parLapply(clust, arglist, do.call, what = 'openCR.fit')
        stopCluster(clust)
    }
    else {
        set.seed (seed)
        output <- lapply(arglist, do.call, what = 'openCR.fit')
    }
    
    message('Completed in ', round((proc.time() - ptm)[3]/60,3), ' minutes at ',
        format(Sys.time(), "%H:%M:%S %d %b %Y"))

    if (inherits(output[[1]], 'openCR'))
        output <- openCRlist(output)

    ## apply standard naming convention
    names(output) <- paste0(prefix, names(arglist))

    output
}