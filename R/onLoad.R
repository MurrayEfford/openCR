###############################################################################
## package 'openCR'
## onLoad.R
## 2020-11-02, 2021-03-13, 2021-04-13
## 2020-05-09 allow possibility defaultncores == 1 
###############################################################################

.onLoad <- function (libname, pkgname) {
    ## also sets environment variable RCPP_PARALLEL_NUM_THREADS
    defaultncores <- RcppParallel::defaultNumThreads()
    if (defaultncores == 1) {
        RcppParallel::setThreadOptions(1)
    }
    else {
        RcppParallel::setThreadOptions(2)
    }
    
    ## advice of Kevin Ushey 2020-03-18, 2021-04
    ## to avoid UBSAN errors from parallelFor
    ## Used in tests
    ## Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")
}

## .onLoad is preferred if actions are required for single functions 
## that may be called without attaching package
