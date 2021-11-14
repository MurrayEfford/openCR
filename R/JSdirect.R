JS.direct <- function (object) {
    nRmrz <- JS.counts(object)
    J <- nrow(nRmrz)
    
    est <- with(nRmrz, {
        M <- m + (R+1) * z / (r+1)  # Seber 1982 bias adjusted formula
        
        p <- m/M; p[c(1,J)] <- NA
        varp <- p^2 * (1-p)^2 * (1/r - 1/R + 1/m + 1/z)
        
        N <- ((n + 1) * M)/(m + 1); N[c(1,J)] <- NA
        pm <- M / N
        varN  <- N * (N - n) * (((M - m + R) * (1/r - 1/R)) / M + (1 - pm)/ m)
        
        phi <-  M[-1] / (M[-J] - m[-J] + R[-J]); phi[J-1] <- NA
        varphi <- phi[-J]^2 * ((M[-1]- m[-1]) * (M[-1]- m[-1] + R[-1]) *
                                     (1/r[-1] - 1/R[-1]) / M[-1]^2 +
                                     (M[-J] - m[-J]) * (1/r[-J] - 1/R[-J]) /
                                     (M[-J] - m[-J] + R[-J]) + (1 - phi[-J]) / M[-1])
        
        varphi <- c(varphi, NA)
        phi <- c(phi, NA)
        covphi <- ( - phi[-J] * phi[-1] * (N[-1] - m[-1]) * (1/r[-1] - 1/R[-1]) ) / M[-1]
        covphi <- c(covphi, NA)
        
        pm <- M/N
        B <- N[-1] - phi[-J] * (N[-J] - n[-J] + R[-J])
        B <- c(B, NA)
        varB <- B[-J]^2 * (M[-1] - m[-1]) * (M[-1] - m[-1] + R[-1]) * (1/r[-1] - 1/R[-1]) / M[-1]^2 +
         (M[-J] - m[-J]) * (phi[-J] * R[-J] * (N[-J] - M[-J]) / M[-J])^2 * (1/r[-J] - 1/R[-J]) / (M[-J] - m[-J] + R[-J]) +
         (N[-J] - n[-J]) * (N[-1] - B[-J]) * (N[-J] - M[-J]) * (1 - phi[-J]) / N[-J] / (M[-J] - m[-J] + R[-J]) +
         N[-1] * (N[-1] - n[-1]) * (N[-1] - M[-1]) / N[-1] / m[-1] +
         phi[-J]^2 * N[-J] * (N[-J] - n[-J]) * (N[-J] - M[-J]) / N[-J] / m[-1]
        varB <- c(varB, NA) 
        
        data.frame(
                   p = p, sep = varp^0.5,
                   N = N, seN = varN^0.5, 
                   phi = phi, sephi = varphi^0.5, covphi = covphi, 
                   B = B, seB = varB^0.5
                   )
    })
    cbind(nRmrz, est)
}
