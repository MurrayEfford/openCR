## JSSA models
## probability of birth at just before b and death just after d
bd.array <- function (beta, phi) {
    J <- length(beta)
    if (length(phi) != J)
        stop ("beta and phi differ in number of sessions")
    if (!all.equal(sum(beta), 1.0, tolerance = 1e-4))
        warning ("beta values do not sum to 1.0 in bd.array")
    phi[J] <- 0
    if (any (phi<0 | phi>1))
        warning ("phi outside range 0-1 in bd.array")
    pbd <- matrix(NA, nrow = J, ncol = J, dimnames = list(b=1:J, d=1:J))
    for (b in 1:J) {
        for (d in b:J) {
            if (b==d) phibd <- 1
            else phibd <- phi[b:(d-1)]
            pbd[b,d] <- beta[b] * prod(phibd) * (1-phi[d])
        }
    }
    as.table(pbd)
}
