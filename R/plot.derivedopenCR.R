plot.derivedopenCR <- function(x, par = 'phi', add = FALSE, xoffset = 0, ylim = NULL,
                        useintervals = TRUE, intermediate.x = TRUE,  ...) {
    if (useintervals)
        xv <- x$estimates$time
    else
        xv <- 0:nrow(x$estimates)
    sessionxv <- xv
    if (intermediate.x & (par %in% c('phi', 'f', 'lambda', 'b', 'BN','BD')))
        xv <- (xv + c(xv[-1],NA))/2
    xv <- xv + xoffset

    sessnames <- rownames(x$estimates)
    if (is.null(sessnames)) sessnames <- 1:length(xv)

    pred <- x$estimates[,par]

    xl <- range(sessionxv)
    yl <- ylim
    if (is.null(yl)) {
        yl <- c(min(pred, na.rm=TRUE)*0.8, max(pred, na.rm=TRUE)*1.05)
        if (yl[1]<0.05) yl[1] <- 0
        if (yl[2]>0.95 & yl[2]<1) yl[2] <- 1
    }
    if (!add) {
        plot (sessionxv, pred, type = 'n', xlab = 'Session', ylab = par,
              ylim = yl, xlim = xl, axes = FALSE, yaxs = 'i')
        axis(1, at = sessionxv, labels = sessnames)
        axis (2, at = pretty(yl), las=1)
        box(bty = 'l')
    }

    points (xv, pred, ...)
    invisible(xv)
}
#source ('d:/open populations/openCR/R/plot.derivedopenCR.R')
