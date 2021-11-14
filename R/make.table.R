###############################################################################
## openCR
## make.table.R
## 2018-02-26 openCR 1.0.0
## 2018-03-27 openCR 1.1.0 select first nsess[i] rows from pred[i]
## 2021-04-21 stratification
###############################################################################

make.table <- function (fits, parm = 'phi', fields = 'estimate', strata = 1, 
    collapse = FALSE, ...) {
    tempname <- deparse(substitute(fits))
    if (!inherits(fits, 'openCRlist')) {
        fits <- openCRlist(fits)
    }
    nmodel <- length(fits)
    if (is.null(names(fits))) {
        names(fits) <- paste0('fit',1:length(fits))
    }
    if (nmodel == 1 && names(fits) == 'fits') {
       names(fits) <- tempname 
    }
    if (length(strata)>1) {
        all.strata <- lapply(fits, function (x) strata(x$capthist))
        if (nmodel>1) {
            for (i in 2:nmodel) {
                if (any(all.strata[[i]] != all.strata[[1]])) 
                    stop ("make.table requires strata to be the same across models")
            }
        }
        tablist <- lapply (strata, make.table, fits = fits, parm = parm, fields = fields, collapse = FALSE, ...)
        names(tablist) <- all.strata[[1]][strata] 
        
        if (collapse) {
            newtab <- do.call(rbind, tablist)
            nr <- sapply(tablist, nrow)
            stratum <- rep(names(tablist), nr)
            rownames(newtab) <- paste (rownames(newtab), stratum, sep = '.')
            names(dimnames(newtab)) <- c('model.stratum', 'session')
            as.table(newtab)
        }
        else {
            tablist
        }
    }
    else {
        # strata of length 1
        nsess <- sapply(fits, function(x) length(primaryintervals(x)[[strata]])+1) 
        pred <- predict(fits, ...)
        rown <- names(fits)
        labellist <- lapply(fits, '[[', 'sessionlabels')
        # for backward compatibility 2021-04-26
        #if (!is.list(labellist[[1]])) labellist <- list(labellist)
        # 2021-09-24 fix
        if (!is.list(labellist)) labellist <- list(labellist)
        coln <- unique(as.character(unlist(labellist)))
        tab <- matrix (nrow = length(rown), ncol = length(coln), 
            dimnames = list(model = rown, session = coln))
        tab <- as.table(tab)
        for (i in rown) { 
            if (parm %in% names(pred[[i]])) {
                preds <- pred[[i]][[parm]]
                if ('stratum' %in% names(preds)) {
                    rowcol <- cbind(rep(i, nsess[i]), unlist(labellist[[i]][[strata]]))
                    preds <- preds[preds$stratum %in% levels(preds$stratum)[strata], ]
                }
                else {
                    rowcol <- cbind(rep(i, nsess[i]), unlist(labellist[[i]]))
                }
                tab[rowcol] <- preds[1:nsess[i],fields]
            }
        }
        tab
    }
}
