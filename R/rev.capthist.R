## 2018-11-25

rev.capthist <- function (x) {
    if (ms(x)) {
        nsess <- length(x)
        if (is.null(intervals(x))) intervals(x) <- rep(1, nsess-1)
        out <- x
        names(out) <- names(x)[nsess:1]
        out[] <- x[nsess:1]
    }
    else {
        if (is.null(intervals(x))) {
            intervals(x) <- rep(1, ncol(x)-1)
            warning ("rev.capthist using intervals = 1")
        }
        cumss <- getcumss(x)                  ## cumulative secondary sessions per primary session
        nsess <- length(cumss)-1
        out <- x
        st <- 0 
        for (j in 1:nsess) {
            k <- nsess-j+1
            ss <- (cumss[k]+1):cumss[k+1]
            out[, st + ss-cumss[k],] <- x[,ss,]
            st <- st + length(ss)
        }
    }
    sessionlabels(out) <- rev(sessionlabels(x))
    intervals(out) <- rev(intervals(x))
    out
}
