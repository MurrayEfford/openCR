## 2018-02-03

cloned.fit <- function (object, nclone = 100, newdata = NULL, linkscale = FALSE) {
    
    compare <- function (pr0, pr1) {
        df <- data.frame(pr0[,c('estimate', 'SE.estimate'), drop = FALSE], 
                         pr1[,c('estimate','SE.estimate'), drop = FALSE])
        df$ratio <- round((df[,2] / df[,4]) / sqrt(nclone), 4)
        if (linkscale) est <- 'beta' else est <- 'estimate'
        names(df) <- c(paste0(c('', 'SE.'), est), 
                       paste0(paste0(c('', 'SE.'), est), '.', as.character(nclone)),
                       'SE.ratio')
        df
    }
    object$details$nclone <- nclone
    args <- object[c('capthist','type','model',
                     'mask','detectfn','binomN','movementmodel',
                     'link', 'fixed','timecov','sessioncov','agecov','details','method','trace','ncores')]
    args$start <- object$fit$par
    fit1 <- do.call(openCR.fit, args)
    if (linkscale) {
        for (i in 1:length(object$link)) object$link[[i]] <- 'identity'
        for (i in 1:length(fit1$link)) fit1$link[[i]] <- 'identity'
    }
    pred0 <- predict(object, newdata = newdata)
    pred1 <- predict(fit1, newdata = newdata)
    mapply(compare, pred0, pred1, SIMPLIFY = FALSE)
}