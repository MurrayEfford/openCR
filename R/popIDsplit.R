popIDsplit <- function (pop) {
    if (ms(pop)) {
        nsess <- length(pop)
        ID <- unique(unlist(sapply(pop, rownames)))
        out <- array(dim = c(length(ID), nsess, 2), dimnames = list(ID, names(pop), c('x','y')))
        for (i in 1:nsess) {
            out[rownames(pop[[i]]),i,] <- unlist(pop[[i]])
        }
        out
    }
    else {
        stop("requires multisession pop object")    
    }
}