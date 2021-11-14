###############################################################################
# JS.counts.R
## 2018-04-20 openCR 1.2.0
## 2021-04-18 stratified
###############################################################################

JS.counts <- function(object, primary.only = TRUE, stratified = FALSE) {
    if (stratified) {
        lapply(object, JS.counts)
    }
    else {
        first <- function(y) {
            w <- 1:length(y)
            as.numeric(w == min(which(y)))
        }
        last <- function(y) {
            w <- 1:length(y)
            as.numeric(w == max(which(y)))
        }
        object <- unsqueeze(object)   # 2018-02-06
        if (inherits(object, 'capthist')) {
            if (ms(object)) {
                ch <- suppressWarnings(reduce(object, by = 'all', outputdetector = 'nonspatial', verify = FALSE))
                object <- join(ch)
            }
            else {
                if (primary.only) object <- primaryonly(object)
            }
            CH <- apply(abs(object), 1:2, sum)>0   ## sum over detectors
        }
        else {
            CH <- abs(object)>0   # 0/1
        }
        nsess <- ncol(CH)
        ni <- apply(CH,2,sum)
        firsti <- as.matrix(apply(CH,1,first))
        lasti <- as.matrix(apply(CH,1,last))
        ui <- apply(firsti,1,sum)
        li <- apply(lasti,1,sum)
        mi <- ni-ui
        ri <- ni-li
        zi <- cumsum(c(0,ri[-nsess]) - mi)
        removed <- apply(object,2, function(x) sum(x<0))
        data.frame(n=ni, R=ni-removed, m=mi, r=ri, z=zi)
    }
}
############################################################################################

