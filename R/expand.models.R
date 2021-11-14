# expand.models

# all additive combinations

expand.models <- function (model, ...) {
    # model is a list of models
    # LHS, RHS defined in utility.R
    oneformula <- function(model1) {
        coll <- "+"  
        if(any(grepl('\\*', as.character(model1))))
            warning ("ignoring interactions in ", LHS(model1), RHS(model1))
        lhs <- LHS(model1)
        vars <- all.vars(RHS(model1))
        tmp <- lapply(0:length(vars), combn, x = vars)
        tmp[[which(sapply(tmp,nrow)==0)]] <- matrix('1')   ## constant
        rhs <- unlist(lapply(tmp, function(x) apply(x, 2, paste, collapse = coll)))
        lapply(paste(lhs, rhs, sep='~'), as.formula, ...)
    }
    tmp <- lapply(model, oneformula)
    names(tmp) <- sapply(model, LHS)
    apply(do.call(expand.grid, tmp), 1, as.list)
}

## expand.models(list(p~t+b+h2, phi~t))
## expand.models(list(p~1, phi~1))

## could construct input for par.openCR.fit
