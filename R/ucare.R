################################################################################
## package 'openCR'
## ucare.cjs.R
## 2018-04-15
## 2018-12-08 get(..., asNamespace("R2ucare")) formulation 
################################################################################

ucare.cjs <- function (CH, tests = "all", by = NULL, verbose = TRUE, rounding = 3, ...) {
    
    onegroup <- function(CH) {
        ## 2-D CH
        X <- apply(CH, 1:2, sum)
        ## collapse to primary only
        X <- t(apply(X, 1, function(x) tapply(x,primarysession,max))) 
        ## ignore multiples
        X[X>1] <- 1
        freq <- rep(1,nrow(X))
        test <- function (funstr) {
            fun <- get(funstr, asNamespace("R2ucare"))
            do.call(fun, 
                    args = list(X = X, freq = freq, verbose = verbose, rounding = rounding),
                    envir = environment())
        }
        notoverall <- tests[tests != "overall_CJS"]
        out <- lapply(notoverall, test)
        names(out) <- notoverall
        if (!verbose & length(out)>0) {
            out <- lapply(out, function(x) as.data.frame(t(data.frame(x))))
            out <- rbind.fill (out)
            rownames(out) <- tests[tests != "overall_CJS"]
        }
        if ("overall_CJS" %in% tests) {
            fun <- get("overall_CJS", asNamespace("R2ucare"))
            if (verbose) {
                out <- c(out, list(fun(X, freq, rounding)))
                names(out) <- tests
            }
            else {
                out <- list(components = out, overall_CJS = fun(X, freq, rounding))
            }
        }
        if (length(out[[1]]) == 0) out <- out[-1]
        if ((length(tests) == 1) & verbose)
            out[[1]]
        else
            out
    }
    
    if (!requireNamespace("R2ucare", quietly = TRUE)) {
        stop ("Package R2ucare is required for this function; please install")
    }
    else {
        if (tolower(tests[1]) == "all") {
            tests <- c("test3sr", "test3sm", "test2ct", "test2cl", "overall_CJS")
        }
        else {
            if (!all(tests %in% c("test3sr", "test3sm", "test2ct", "test2cl", "overall_CJS")))
                stop("only these tests allowed: 'test3sr', 'test3sm', 'test2ct', 'test2cl', 'overall_CJS'")
        }
        if (ms(CH)) {
            warning("multisession capthist collapsed with 'join'")
            CH <- join(CH)
        }   
        if (is.null(intervals(CH))) intervals(CH) <- rep(1,ncol(CH)-1)
        if (sum(intervals(CH)>0) == 0) {
            stop("R2ucare is not meaningful for closed populations")
        }
        if (any(intervals(CH) == 0)) {
            warning("secondary sessions collapsed for R2ucare")
        }
        if (any(CH<0)) {
            warning("Non-release of some animals (CH<0) has been ignored")
            CH <- abs(CH)
        }
        
        ## unsqueeze
        CH <- unsqueeze(CH)
        primarysession <- primarysessions(intervals(CH)) ## map secondary to primary
        
        ## grouped
        if (!is.null(by)) {
            if (!is.character(by))
                stop ("by should be a character value naming a covariate")
            if (by %in% names(covariates(CH))) {
                CH <- split(CH, f = covariates(CH)[,by], ...)  
                lapply(CH, onegroup)
            }
            else stop ("by should name a covariate for grouping")
        }
        ## ungrouped
        else onegroup(CH)
    }
}