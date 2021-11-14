###############################################################################
# m.array.R
## 2018-04-20 openCR 1.2.0
## 2021-04-18 stratified
## 2021-07-26 bug fix
###############################################################################

m.array <- function (object, primary.only = TRUE, never.recaptured = TRUE, last.session = TRUE, stratified = FALSE) {
    if (ms(object) && stratified) {
        lapply(object, m.array, primary.only, never.recaptured, last.session, stratified = FALSE)
    }
    else {
        if (ms(object)) {
            object <- join(reduce(object, by = 'all', outputdetector = 'nonspatial', verify = FALSE))
        }
        else {
            if (primary.only) object <- primaryonly(object)
        }
        object <- unsqueeze(object)
        S <- ncol(object)
        df <- as.data.frame(object)
        df <- split(df, df$ID)
        df <- lapply(df, function(x) {
            x <- rbind(x, x[1,])
            x[nrow(x),'Occasion'] <- S+1
            cbind(x$Occasion[-nrow(x)], x$Occasion[-1])
        })
        df <- do.call(rbind, df)
        tab <- table(
            factor(df[,1], levels = 1:S),       # levels applied to fix bug 2021-07-26
            factor(df[,2], levels = 1:(S+1))
            )
        tab[lower.tri(tab, diag=TRUE)] <- NA
        if (!last.session) {
            tab <- tab[-nrow(tab),,drop = FALSE]
        }
        rows <- 1:nrow(tab)
        n <- unlist(counts(object,'n'))[rows]
        lost <- unlist(counts(object,'losses'))[rows]
        tab <- cbind(R = n - lost,
                     tab[,-c(1,S+1),drop = FALSE])
        row.names(tab) <- rows
        if (!is.null(sessionlabels(object))) {
            rownames(tab) <- sessionlabels(object)[1:nrow(tab)]
            colnames(tab)[2:ncol(tab)] <- sessionlabels(object)[2:ncol(tab)]
        }
        if (never.recaptured)
            tab <- cbind(tab, NRecap = tab[,1] - apply(tab[,-1,drop = FALSE],1,sum, na.rm = TRUE))
        as.table(tab)
    }
}
############################################################################################
