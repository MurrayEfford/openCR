## cyclic fixing over real parameters
## cf Schwarz and Arnason 1996 p865
## 2018-02-12, 2018-05-01
## 2018-05-02 why different eigenvalues? can judge the rank differently
## 2018-05-02 no speed gain, yet, so not published
cyclic.fit <- function (..., maxcycle = 10, tol = 1e-5, trace = FALSE) {
    ptm  <- proc.time()
    starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y")
    args <- list(...)
    if (is.null(args$method) || args$method == 'none') args$method <- 'Newton-Raphson'
    defaultdetails <- list(hessian = TRUE)
    args$details <- replacedefaults(defaultdetails, args$details)
    args$details$LLonly <- TRUE
    oldhessian <- args$details$hessian
    args$details$hessian <- FALSE
    fit0 <- do.call(openCR.fit, args)
    args$details$LLonly <- FALSE
    LL0 <- fit0[1]
    parindx <- attr(fit0, 'parindx')
    beta <- fit0[-1]
    nreal <- length(parindx)
    for (j in 1:maxcycle) {
        for (pari in 1:nreal) {
            args$details$fixedbeta <- beta
            args$details$fixedbeta[parindx[[pari]]] <- NA
            args$start <- beta
            fit1 <- do.call(openCR.fit, args)
            beta[is.na(args$details$fixedbeta)] <- fit1$fit$par
            if (trace) {
                message("Iter ", j, " ", names(parindx)[pari], " LL = ", logLik(fit1))
            }
        }
        if (abs(logLik(fit1)-LL0)<tol) break
        LL0 <- logLik(fit1)
    }
    if (abs(logLik(fit1)-LL0) >= tol) warning("reached maxcycle but tol not achieved")
    args$details$fixedbeta[] <- NA
    # args$details$hessian <- oldhessian
    args$method = 'none'
    out <- do.call(openCR.fit, args)
    out$call <- ""
    out$starttime <- starttime
    out$proctime <- proc.time()[3] - ptm[3]
    out
}

# cyclic.fit(capthist=dipperCH, model = list(p~t, phi~t), tol = 1e-5, trace = TRUE)

