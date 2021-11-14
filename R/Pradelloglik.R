sumj <- function(uv, j, k) {
    if (j>k)
        return (0)
    else {
        return(sum(uv[j:k]))
    }
}

pradelloglik <- function (type, w, openval,  PIAJ, intervals) {
    s <- length(intervals) + 1
    w <- matrix(w, nrow=s)
    ni <- w[,1]       # number viewed at i
    u <- ni - w[,3]   # number viewed for first time at i
    v <- ni - w[,4]   # number viewed for last time at i
    d <- ni - w[,2]   # number removed at i
    mu <- w[,2] / ni
    mu[s] <- 1.0
    p <- getpj   (1, 1, openval, PIAJ)
    phij <- getphij (1, 1, openval, PIAJ, intervals)
    if (type==26)
        gamj <- getgamj (1, 1, openval, PIAJ, intervals)
    else
        gamj <- getgamjl (1, 1, openval, PIAJ, intervals)
    chi <- xi <- numeric(s)
    chi[s] <- 1
    for (i in (s-1):1) {
        chi[i] <- (1- phij[i]) + phij[i] *(1-p[i+1]) * chi[i+1]
    }
    xi[1] <- 1
    for (i in 2:s) {
        xi[i] <- (1-gamj[i]) + gamj[i] * (1-p[i-1]) / (1-p[i-1]*(1-mu[i-1]) ) * xi[i-1];
    }
    value <- 0
    for (i in 1:s) {
        # if ((xi[i]<=0) | (gamj[i]<=0) | (p[i]<=0) | (p[i]>=1) | (phij[i]<=0) | (mu[i]>=1) | (mu[i]<=0) | (chi[i]<=0))
        #    return(-1e10)
        if (xi[i] > 0) 
            value <- value + u[i] * log (xi[i]) 
        if (gamj[i]>0)
            value <- value + sumj(u,1,i-1) * log (gamj[i]);
        if (p[i]>0)
            value <- value + ni[i] * log (p[i])
        if (p[i]<1)
            value <- value + (sumj(u,1,i) - sumj(v,1,i-1) - ni[i]) * log (1-p[i])
        if (phij[i]>0)
            value <- value + sumj(v,i+1,s) * log (phij[i])
        if (mu[i] < 1) 
            value <- value + (ni[i] - d[i]) * log (mu[i]) + d[i] * log (1-mu[i])
        if (chi[i]>0)
            value <- value + sumj(u,i+1,s) * log (1 - p[i] * (1 - mu[i])) + (v[i]-d[i]) * log (chi[i])
    }
    # probability of detection
    BL <- 0
    for (i in 1:s) {
        prd <- xi[i]
        if (i > 1) {
            for (j in 1:(i-1)) 
                prd <- prd * phij[j] * (1-p[j] * (1-mu[j]))
        }
        if (i < s) {
            for (j in (i+1):s) prd <- prd * gamj[j]
        }
        BL <- BL + prd * p[i]
    }
    value <- c(value, - sum(u) * log (BL))
    return (value)   # return log likelihood 
    
}