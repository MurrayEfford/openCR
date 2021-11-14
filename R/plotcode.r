## plotcode.R
## 2013-05-28, 2013-10-21, 2018-03-19

chartHistories <- function (ch, seq = 1:nrow(ch), start = as.Date("1980/01/1"), finish = as.Date("2007/01/1"),
                         by = "year", IDcex = 0.6, offset = -200) {
    ## one history
    onebird <- function (i) {
        j <- seq[i]
        individual <- row.names(ch)[j]
        occ <- apply(abs(ch[j,,]),1,sum) > 0
        ntimes <- sum(occ)
        times <- cumss[occ]
        lines (times, rep(i, ntimes), type='l', col='black')
        
        # points (df2$Date[OK3], individual[OK3], pch=21, type='p', col='black', bg='white')
        # points (df2$Date[OK2], individual[OK2], pch=17, type='p')
        # text (start+off, individual+0.1, df2$Colour[1], cex=IDcex, adj=0)
        # text (start+off*1.3, individual+0.1, df2$Age[1], cex=IDcex, adj=0)
    }
    
    cumss <- c(0,cumsum(intervals(ch)))
    ## baseplot
    plot(0,0, xlim = c(0, max(cumss)), ylim = c(0,nrow(ch)),
        type = 'n', axes = F,
        xlab = '', ylab = '')
    # labs <- seq(start, finish, by = by)
    # axis(1, at=labs, lab=F)
    # labs1 <- labs[-length(labs)]  # drop last
    # axis(1, at=labs1+15, tick=F, lab=F) ## format(labs1, format='%b'))
    # mtext(side=1, line = 1, at = seq(as.Date("2009/07/1"), as.Date("2017/07/1"),by="year"),
    #     2009:2017)
    
    # par(tck=1)
    
    # dfdates <- unique(df$Date)
    # rdfdates <- range(dfdates)
    
    #axis(side=3, at=dfdates, lab=F, col='grey')
    
    # axis(side=3, at=labs, lab=F, col='grey')
    par(tck = -0.01)

    # axis(side=3, at=rdfdates, lab=format(rdfdates, format='%d/%m/%Y'))
    # mtext (side=3, line=3, at=start, text=df$Site[1], adj=0, cex=0.9)
    ## one animal at a time...
    for (i in 1:nrow(ch)) {
        onebird(i)
    }
    invisible()
}

# tmp <- join(ch8095)
# par(mfrow=c(1,1), mar=c(4,4,1,1))
# idnum <- match(as.numeric(row.names(tmp)), possums$TATTOO)
# or <- order(possums$FIRST_DATE[idnum])
# chartHistories(tmp[1:10,,], seq = or)
# axis(1)
