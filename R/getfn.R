###############################################################################
## openCR
## getfn.R
## 2018-02-26 openCR 1.0.0
## 2018-11-11 buggy if J<3; fixed for getbeta0
## 2019-04-21 getmoveargs bug fixed - moveargsi + 1
###############################################################################

#------------------------------------------------------------------------

getD <- function ( type, J, nmix, pmix, openval, PIAJ, intervals) {
    ## return superD
    if (type %in% c(7,12,13,24)) {
        # JSSAsecrf = 7
        # JSSAsecrl = 12
        # JSSAsecrb = 13
        # JSSAsecrg = 24
        return(openval[PIAJ[1,1,1,1], 4])
    }
    else {
        # JSSAsecrD = 8
        # JSSAsecrB = 14

        sumB <- 0
        for (x in 1:nmix) {
            phij <- openval[PIAJ[1,1,1:(J-1),x], 2]^intervals    ## per session
            if (type %in% 8) {
                D <- openval[PIAJ[1,1,1:J,x], 3]
                B <- D
                for (j in 1:(J-1)) {
                    B[j+1] = D[j+1] - D[j] * phij[j]
                }
                sumB <- sumB + sum(B) * pmix[x]
            }
            else { ## type == 14
                B <- openval[PIAJ[1,1,1:J,x],3] 
                sumB <- sumB + sum(B) * pmix[x]
            }
        }
        return (sumB)
    }
}
#------------------------------------------------------------------------

getN <- function ( type, ncf, J, nmix, pmix, openval, PIAJ, intervals) {
    ## return superN
    if (type %in% c(2,3,4,22,28)) {
        # JSSAb = 2
        # JSSAl = 3
        # JSSAf = 4
        # JSSAg = 22
        # JSSAk = 28
        return(openval[PIAJ[1,1,1,1], 4])
    }
    else {
        # JSSAB = 18
        # JSSAN = 19
        sumB <- 0
        for (x in 1:nmix) {
            phij <- openval[PIAJ[1,1,1:(J-1),x], 2]^intervals    ## per session
            if (type %in% 19) {
                N <- openval[PIAJ[1,1, 1:J, x], 3]
                B <- N
                for (j in 1:(J-1)) {
                    B[j+1] = N[j+1] - N[j] * phij[j]
                }
                sumB <- sumB + sum(B) * pmix[x]
            }
            else { ## type == 18
                B <- openval[PIAJ[1,1,1:J,x],3] 
                sumB <- sumB + sum(B) * pmix[x]
            }
        }
        ## if (type %in% c(2, 3, 4, 21, 22)) sumB <- sumB + ncf
        return (sumB)
    }
}
#------------------------------------------------------------------------

getp <- function (n, x, openval, PIA) {
    return(openval[PIA[1,n,,1,x], 1])    # p for each secondary session
}
#------------------------------------------------------------------------

getpj <- function (n, x, openval, PIAJ) {
    return(openval[PIAJ[1,n,,x], 1])    # p for each primary session
}
#------------------------------------------------------------------------

getphij <- function (n, x, openval, PIAJ, intervals) {
    J1 <- 1:length(intervals)
    phi <- openval[PIAJ[1,n,J1,x], 2]
    c(exp(log(phi) * intervals),0)     ## phi for each primary session, zero for last
}
#------------------------------------------------------------------------

## bug fix moveargsi + 1 2019-04-21
getmoveargs <- function (n, x, openval, PIAJ, intervals, moveargsi) {
    ## z is adjustment for intrusion of 'z' parameter between sigma and move.a
    ## in the realparameter table (openval, openval0) when detectfn has 3 parameters
    J <- length(intervals) + 1
    moveargs <- matrix(0, nrow = J, ncol = 2)
    moveargs[,1] <- openval[PIAJ[1,n,,x], moveargsi[1]+1, drop = FALSE]
    if (moveargsi[2]>0)
        moveargs[,2] <- openval[PIAJ[1,n,,x], moveargsi[2]+1, drop = FALSE]
    moveargs[nrow(moveargs),] <- 0    ## movea for each primary session, zero for last
    moveargs
}
#------------------------------------------------------------------------

getgamj <- function (n, x, openval, PIAJ, intervals) {
    J <- length(intervals)+1
    J2 <- 2:J
    g <- openval[PIAJ[1,n, J2, x],3]
    c(0, exp(log(g) * intervals))
}
#------------------------------------------------------------------------

getkapj <- function (n, x, openval, PIAJ) {
    c(1, openval[PIAJ[1,n,-1, x],3])
}
#------------------------------------------------------------------------

getgamjl <- function (n, x, openval, PIAJ, intervals) {
    J1 <- 1:length(intervals)
    phi <- openval[PIAJ[1,n, J1, x],2]
    phij <- exp(log(phi) * intervals)
    lam <- openval[PIAJ[1,n, J1, x],3]
    lamj <- exp(log(lam) * intervals)
    c(0, phij/lamj)
}
#------------------------------------------------------------------------

getfj <- function (n, x, openval, PIAJ, intervals, phij) {
    J1 <- 1:length(intervals)
    f <- openval[PIAJ[1,n, J1, x], 3]
    c(exp(log(phij[J1]+f) * intervals) - exp(log(phij[J1] * intervals)), 0)
}
#------------------------------------------------------------------------

getlj <- function (n, x, openval, PIAJ, intervals) {
    J1 <- 1:length(intervals)
    l <- openval[PIAJ[1,n, J1, x], 3]
    c(exp(log(l) * intervals), 0)
}
#------------------------------------------------------------------------

getbeta0 <- function (n, x, openval, PIAJ) {
    J <- dim(PIAJ)[3] # ncol(PIAJ)
    beta <- openval[PIAJ[1,n, , x], 3]   
    sumexp <- sum(exp(beta[-1]))
    if (J>1) {
        beta[2:J] <- exp(beta[2:J]) / (1 + sumexp)
        beta[1] <- 1 - sum(beta[2:J])
    }
    beta
}
#------------------------------------------------------------------------

# per capita recruitment cf Link & Barker 2005, Schwarz 'Gentle Intro'
getbetaf <- function (n, x, openval, PIAJ, phij, intervals) {
    fj <- getfj (n, x, openval, PIAJ, intervals, phij)
    J <- length(intervals)+1
    J1 <- 1:(J-1)
    d <- c(1, cumprod(phij+fj))[1:J]
    beta <- fj * d
    beta <- c(1, beta[J1])
    beta/sum(beta)
}
#------------------------------------------------------------------------

getbetal <- function (n, x, openval, PIAJ, phij, intervals) {
    J <- length(intervals)+1
    J1 <- 1:(J-1)
    lambdaj <- getlj (n, x, openval, PIAJ, intervals);
    fj <- ifelse(lambdaj < phij, 0, lambdaj - phij)
    d <- c(1, cumprod(phij+fj))[1:J]
    beta <- fj * d
    beta <- c(1, beta[J1])
    beta/sum(beta)
}
#------------------------------------------------------------------------

getbetag <- function (n, x, openval, PIAJ, phij, intervals) {
    J <- length(intervals)+1
    J1 <- 1:(J-1)
    gamj <- getgamj (n, x, openval, PIAJ, intervals)
    gamj1 <- gamj[-1]
    fj <- ifelse(gamj1 <= 0, 0, phij[-J] * (1/gamj1 - 1))
    d <- c(1, cumprod(phij[J1]+fj[J1]))[J1]
    beta <- fj * d
    beta <- c(1, beta)
    beta/sum(beta)
}
#------------------------------------------------------------------------

getbetak <- function (n, x, openval, PIAJ, phij, intervals) {
    J <- length(intervals)+1
    J1 <- 1:(J-1)
    f <- tau <- beta <- fprod <- numeric(J)
    kap <- getkapj (n, x, openval, PIAJ)
    p <- getpj(n, x, openval, PIAJ)
    tau[1] <- 1/p[1]
    for (j in J1) {
        f[j] <- (kap[j+1] - kap[j]/p[j] * (1 - p[j]) * phij[j] * p[j+1]) / (tau[j] * p[j+1])
        tau[j+1] <- tau[1] * prod(phij[1:j] + f[1:j])
    }
    for (j in 2:(J-1)) fprod[j] <- f[j] * prod(phij[1:(j-1)] + f[1:(j-1)])
    beta[1] <- 1/ ( 1 + f[1] + sum(fprod[2:(J-1)]))   ## beta_0
    beta[2] <- beta[1] * f[1]
    beta[3:J] <- beta[1] * fprod[2:(J-1)]
    beta
}
#------------------------------------------------------------------------

# parameterisation cf Pledger et al. 2010 p 885
getbetaB <- function (n, x, openval, PIAJ) {
    B <-  openval[PIAJ[1,n, , x],3]
    B / sum(B)
}
#------------------------------------------------------------------------

getbetaD <- function (n, x, openval, PIAJ, phij) {
    J <- length(phij)
    J1 <- 1:(J-1)
    D <- openval[PIAJ[1,n, , x],3]
    B <- D
    B[2:J] <- D[2:J] - D[J1] * phij[J1]
    B/sum(B)
}
#------------------------------------------------------------------------

getbeta <- function (type, n, x, openval, PIAJ, intervals, phij) {
    if (type %in% c(2, 17, 11, 13, 
                    30, 31, 41, 43))   ## added 2018-11-11
        getbeta0 (n, x, openval, PIAJ) 
    else if (type %in% c(4, 15, 27, 7, 9, 39))
        getbetaf (n, x, openval, PIAJ, phij, intervals)
    else if (type %in% c(3, 16, 10, 12, 20))
        getbetal (n, x, openval, PIAJ, phij, intervals)
    else if (type %in% c(14, 18))
        getbetaB (n, x, openval, PIAJ)
    else if (type %in% c(8, 19))
        getbetaD (n, x, openval, PIAJ, phij)
    else if (type %in% c(22, 23, 24, 25, 26))
        getbetag (n, x, openval, PIAJ, phij, intervals)
    else if (type %in% c(28,29))
        getbetak (n, x, openval, PIAJ, phij, intervals)
    else stop("no beta for this type")
}
#------------------------------------------------------------------------
