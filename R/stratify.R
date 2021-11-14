###############################################################################
# stratify.R
## 2021-04-22, 24
###############################################################################

# Input cases
# 1. list (strata) of single-session capthist (already joined)
# 2. list (strata) of multi-session capthist; join within ms capthist
# 3. one single-session capthist, split by covariate (already joined)
# 4. one multi-session capthist, split by covariate then join

stratify <- function (..., intervals = NULL, MoreArgs = list(), 
    covariate = NULL, bytraps = FALSE) {
    splitfn <- function(x) {
        if (!bytraps)
            factor(covariates(x)[[covariate]])
        else
            factor(covariates(traps(x))[[covariate]])
    }
    input <- list(...)
    if (length(input) == 1) {
        input <- input[[1]]
        if (is.list(input)) {
            if (inherits(input, 'capthist'))
                inputcase <- 4
            else
                inputcase <- if (ms(input[[1]])) 2 else 1
        }
        else {
            if (is.null(covariate)) stop ("provide name of covariate to split capthist")
            inputcase <- if (ms(input)) 4 else 3
        }
    }
    else {
        inputcase <- if (ms(input[[1]])) 2 else 1
    }
    #-------------------------------------------------
    if (inputcase == 1) {
        stratumlist <- input
    }
    else if (inputcase == 3) {
        # assume already 'join'ed
        # with suitable intervals
        stratumlist <- split(input, splitfn(input))
    }
    else {
        # input cases 2,4
        if (inputcase == 4) input <- split(input, lapply(input, splitfn))
        if(is.null(intervals)) intervals <- lapply(input, function(x) rep(1, length(x)-1))
        stratumlist <- mapply(join, input, intervals = intervals, MoreArgs = MoreArgs)
    }
    if (!all(sapply(stratumlist, inherits, 'capthist'))) 
        stop ("stratify requires list of capthist")
    if (any(sapply(stratumlist, ms))) 
        stop ("stratify failed: a stratum cannot be multi-session")
    class(stratumlist) <- c('capthist','list')
    if (is.null(names(stratumlist))) {
        names(stratumlist) <- rep("",length(stratumlist))
        dots2 <- substitute(list(...))[-1]    
        if (length(stratumlist) == length(dots2))
            names(stratumlist) <- sapply(dots2, deparse)
    }
    if (any(duplicated(names(stratumlist))) | any(names(stratumlist)=="")) {
        warning ("session names replaced to avoid duplication or omission")
        names(stratumlist) <- 1:length(stratumlist)
    }
    MS.capthist(stratumlist)
}
