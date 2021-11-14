#############################################################################
## package 'openCR'
## read.inp.R
## 2017-12-24, 2018-01-22
#############################################################################

remove.comments <- function (x, start = '/*', end = '*/') {
    startpos <- regexpr(start, x, fixed = TRUE) 
    while (startpos>0) {
        endpos <- regexpr(end, x, fixed = TRUE) 
        x <- str_c(substring(x, 1, startpos-1), substring(x, endpos+2, nchar(x)), collapse="")
        startpos <-  regexpr(start, x, fixed = TRUE) 
    }
    return(x)
}
    
read.inp <- function(filename, ngroups = 1, grouplabel = 'group', grouplevels = NULL, covnames = NULL, skip = 0){
    if (length(grouplevels)>0) ngroups <- length(grouplevels)
    li <- readLines(filename)
    if (skip>0) li <- li[-(1:skip)] 
    all <- stringr::str_c(li, collapse ="")
    all <- remove.comments(all)
    li <- str_split(all, ';')
    li <- sapply(li, stringr::str_trim)             ## remove leading or trailing blanks
    li <- gsub('  ', ' ', li) 
    li <- li[li!=""]
    df <- do.call(rbind, stringr::str_split(li, ' '))
    ncov <- length(covnames)
    ncolch <- ncol(df)- (ngroups + ncov)
    ch <- df[, 1:ncolch, drop = FALSE]
    ch <- apply(ch, 1, paste, collapse='')
    ch <- gsub(' ','',ch)
    freq <- as.numeric(df[, ncolch+(1:ngroups)])
    gp <- rep(rep(1:ngroups, each = length(ch)),freq)
    ch <- rep(ch, ngroups)
    ch <- rep(ch,freq)

    freq <- rep(1, length(ch))
    newdf <- data.frame(ch=ch, freq=freq, stringsAsFactors = FALSE)
    if (ngroups>1) {
        if (!is.null(grouplevels)) {
            gp <- factor(grouplevels[gp], levels = grouplevels)
        }
        newdf[[grouplabel]] <- gp
    }
    if (length(covnames)>0) {
        covars <- df[,ncolch+ngroups+(1:ncov), drop = FALSE]
        names(covars) <- covnames
        newdf <- rbind(newdf, covars, )
    }
    rownames(newdf) <- 1:nrow(newdf)
    secr::unRMarkInput(newdf)
}

# alternative code 2018
# rudimentary; does not allow multiple groups
# read.inp <- function (filename, skip = 0, ...) {
#     inpstr <- readLines(filename, warn = FALSE)
#     if (skip>0) inpstr <- inpstr[-(1:skip)] 
#     strs <- strsplit(inpstr, ' ')
#     last <- sapply(strs, length)
#     freq <- sapply(strs, tail, 1)
#     freq <- as.numeric(gsub(';', '', unlist(freq)))
#     ch <- mapply('[', strs, -last, SIMPLIFY = FALSE)
#     ch <- sapply(ch, paste, collapse = "")
#     df <- data.frame(ch=ch, freq=freq, stringsAsFactors = FALSE)
#     unRMarkInput(df)
# }