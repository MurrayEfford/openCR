## test-kernel.R
## started 2021-06-29

library(openCR)
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

msk <- make.mask(nx = 51, ny = 51, type = 'rectangular', spacing = 10, 
    buffer = 0)
msk[,] <- msk[,] - 255   # centre on 0,0

# define full and sparse kernels
k.BVN <- make.kernel('BVN', kernelradius = 20, spacing = 10, move.a = 50, 
    clip = TRUE, sparsekernel = FALSE, r0 = 0)
k.BVT1 <- make.kernel('BVT', kernelradius = 20, spacing = 10, move.a = 50, 
    move.b = 1, clip = TRUE, sparsekernel = FALSE, r0 = 0)
k.BVN.sp <- make.kernel('BVN', kernelradius = 20, spacing = 10, move.a = 50, 
    clip = TRUE, sparsekernel = TRUE, r0 = 0)

# use traps object as basis for closed circular polygon
poly <- make.circle(n = 500, radius = 100)
poly <- as.matrix(poly[c(1:500,1),])   # close

test_that("full BVN kernel has expected proportion in 2 x move.a", {
    expect_equal(proportionInPolygon(k.BVN, poly, 'kernelp'), 0.8566384, 
        tolerance = 1e-4, check.attributes = FALSE)
    
})

# pkernel(100, 'BVT', move.a=50, move.b=1) = 0.8, but
# pkernel(100, 'BVT', move.a=50, move.b=1, truncate = 200) = 0.85 
# and for this discretized and truncated kernel...
test_that("full BVT1 kernel has expected proportion in 2 x move.a", {
    expect_equal(proportionInPolygon(k.BVT1, poly, 'kernelp'), 0.8427389, 
        tolerance = 1e-4, check.attributes = FALSE)
    
})

test_that("proportionInPolygon matches expected value after cumMove from point", {
    
    # initial distribution defaults to central point; one step
    X <- cumMove(mask = msk, kernel = k.BVN.sp, nstep = 1)
    expect_equal(proportionInPolygon(X,poly), 0.8632878, 
        tolerance = 1e-4, check.attributes = FALSE)
    
})

test_that("proportionInPolygon matches expected value after cumMove from polygon", {
    
    # initial distribution across polygon; two steps
    X <- cumMove(poly, msk, k.BVN.sp, nstep = 2)
    expect_equal(proportionInPolygon(X, poly),  0.4761272, 
        tolerance = 1e-4, check.attributes = FALSE)
    
})

test_that("qkernel medians match expected", {
    expect_equal(qkernel(0.5, 'BVN', 30), sqrt(-2*log(0.5)* 30^2), 
        tolerance = 1e-6, check.attributes = FALSE)
    expect_equal(qkernel(0.5, 'BVE', 30),  qgamma(0.5, shape=2, scale=30), 
        tolerance = 1e-6, check.attributes = FALSE)
    expect_equal(qkernel(0.5, 'BVT', 30, 3), 30*sqrt(0.5^(-1/3)-1), 
        tolerance = 1e-6, check.attributes = FALSE)
    expect_equal(qkernel(0.5, 'RDE', 30), -30*log(0.5), 
        tolerance = 1e-6, check.attributes = FALSE)
    expect_equal(qkernel(0.5, 'RDG', 30, 3), qgamma(0.5, shape=3, scale=30),  
        tolerance = 1e-6, check.attributes = FALSE)
    expect_equal(qkernel(0.5, 'RDL', 30, 3), 30, 
        tolerance = 1e-6, check.attributes = FALSE)
    
})

test_that("gkernel matches expected", {
    r <- 30
    alpha <- 30
    beta <- 3
    mu <- log(alpha)                ## frL
    sigma <- sqrt(log(1 + 1/beta))  ## frL
    expect_equal(gkernel(r, 'BVN', alpha), exp(-r^2/2/alpha^2)/2/pi/alpha^2, 
        tolerance = 1e-6, check.attributes = FALSE)
    expect_equal(gkernel(r, 'BVE', alpha),  exp(-r/alpha)/2/pi/alpha^2, 
        tolerance = 1e-6, check.attributes = FALSE)
    expect_equal(gkernel(r, 'BVT', alpha, beta), beta/pi/alpha^2 * (1 + r^2/alpha^2)^(-beta-1), 
        tolerance = 1e-6, check.attributes = FALSE)
    expect_equal(gkernel(r, 'RDE', alpha), exp(-r/alpha)/2/pi/r/alpha, 
        tolerance = 1e-6, check.attributes = FALSE)
    expect_equal(gkernel(r, 'RDG', alpha, beta), r^(beta-2) * exp(-r/alpha)/2/pi/gamma(beta)/alpha^beta,  
        tolerance = 1e-6, check.attributes = FALSE)
    expect_equal(gkernel(r, 'RDL', alpha, beta), exp(-(log(r)-mu)^2 / 2 / sigma^2) /(2*pi)^1.5/r^2/sigma, 
        tolerance = 1e-6, check.attributes = FALSE)
})

test_that("discretized kernel matches expected", {
    alpha <- 30
    # full kernel BVE
    k <- make.kernel('BVE', kernelradius = 10, spacing = 10, move.a = alpha,
        sparsekernel = FALSE, clip = TRUE, normalize = TRUE, r0=0)
    d <- sqrt(apply(k^2,1,sum))
    p <- covariates(k)$kernelp
    ci <- which(k$x==0 & k$y==0)
    expect_equal(sum(p), 1.0, tolerance = 1e-8)
    expect_equal(sum(p*d), 47.13507891, tolerance = 1e-6)
    
    # full kernel BVC, r0 = 1/sqrt(pi)
    k <- make.kernel('BVC', kernelradius = 10, spacing = 10, move.a = alpha,
        sparsekernel = FALSE, clip = TRUE, normalize = FALSE, r0=1/sqrt(pi))
    d <- sqrt(apply(k^2,1,sum))
    p <- covariates(k)$kernelp
    ci <- which(k$x==0 & k$y==0)
    expect_equal(sum(p), 0.7258447181, tolerance = 1e-6)
    p <- p/sum(p)
    expect_equal(sum(p*d), 41.59256317, tolerance = 1e-6)
    
    # sparse kernel
    k <- make.kernel('BVE', kernelradius = 10, spacing = 10, move.a = alpha,
        sparsekernel = TRUE, clip = TRUE, normalize = TRUE, r0=0)
    d <- sqrt(apply(k^2,1,sum))
    p <- covariates(k)$kernelp
    expect_equal(sum(p), 1.0, tolerance = 1e-8)
    expect_equal(sum(p*d), 45.9770124, tolerance = 1e-6)
    
    # large full kernel
    k <- make.kernel('BVE', kernelradius = 100, spacing = 5, move.a = alpha,
        sparsekernel = FALSE, clip = TRUE, normalize = TRUE, r0=0)
    d <- sqrt(apply(k^2,1,sum))
    p <- covariates(k)$kernelp
    expect_equal(sum(p*d), 59.984398563, tolerance = 1e-6)
    
})

test_that("zero-inflated quantiles matches expected", {
        expect_equal(qkernel(0.9, 'BVNzi', 30, 0.3), 59.1830910675, 
            tolerance = 1e-6, check.attributes = FALSE)
        expect_equal(qkernel(0.9, 'BVEzi', 30, 0.3), 103.0669556, 
            tolerance = 1e-6, check.attributes = FALSE)
        expect_equal(qkernel(0.9, 'RDEzi', 30, 0.3),  58.3773044717, 
            tolerance = 1e-6, check.attributes = FALSE)
        expect_equal(qkernel(0.9, 'UNIzi', 0.3, truncate=100), 100, 
            tolerance = 1e-6, check.attributes = FALSE)
        expect_equal(qkernel(0.2, 'UNIzi', 0.3, truncate=100), 0, 
            tolerance = 1e-6, check.attributes = FALSE)
})

# problem that can't find userfn 2021-08-08    
# test_that("user-defined kernel matches canned version", {
#     userfn <- function(r,a) {
#         exp(-r^2/2/a^2)
#     }
#     ku <- make.kernel(userfn, 10, 10, 30)
#     kn <- make.kernel('BVN', 10, 10, 30)
#     r <- sqrt(apply(ku^2,1,sum))
#     ku.p <- covariates(ku)$kernelp
#     kn.p <- covariates(kn)$kernelp
#     expect_equal(sum(ku.p*r), sum(kn.p*r),
#         tolerance = 1e-6, check.attributes = FALSE)
# })

# make.kernel('UNIzi', 30,10,0.5)
