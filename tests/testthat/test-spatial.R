## Started 2020-12-13

library(openCR)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

# create small working proximity dataset
suppressWarnings(smallCH <- join(subset(ovenCHp, sessions = 1:3, traps = 1:20, dropnullocc = FALSE)))

msk <- make.mask(traps(smallCH), buffer = 200, nx = 20, type = 'trapbuffer')

argssecr <- list(capthist = smallCH, mask = msk,
    start = list(lambda0 = 0.037, sigma = 65.3, phi = 0.497, f = 0.560),
    movementmodel = "static", 
    details = list(LLonly = TRUE))

argsmove <- list(capthist = smallCH, mask = msk, type = "PLBsecrf",
    start = list(lambda0 = 0.037, sigma = 65.3, phi = 0.497, f = 0.560, move.a = 100),
    movementmodel = "BVN", edgemethod = "truncate", sparsekernel = FALSE,
    kernelradius = 10, details = list(LLonly = TRUE))

argsmoves <- list(capthist = smallCH, mask = msk, type = "PLBsecrf",
    start = list(lambda0 = 0.037, sigma = 65.3, phi = 0.497, f = 0.560, move.a = 100),
    movementmodel = "BVN", edgemethod = "truncate", sparsekernel = TRUE,
    kernelradius = 20, details = list(LLonly = TRUE))

test_that("test data OK", {
    expect_equal(sum(smallCH), 60)
    expect_equal(RPSV(smallCH, CC = TRUE), 43.85065, tolerance = 1e-5)
})

test_that("correct PLBsecr likelihood", {
    argssecr$type = 'PLBsecrf'
    expect_equal(do.call(openCR.fit, argssecr)[1], -357.296968, 
        tolerance = 1e-4, check.attributes = FALSE)
})

test_that("correct JSSAsecr likelihood", {
    argssecr$start$superD <- 2.0
    argssecr$type = 'JSSAsecrf'
    expect_equal(do.call(openCR.fit, argssecr)[1], -359.89066477, 
        tolerance = 1e-4, check.attributes = FALSE)
})

test_that("correct movement likelihood", {
    expect_equal(do.call(openCR.fit, argsmove)[1], -358.092530, 
        tolerance = 1e-4, check.attributes = FALSE)
    
    argsmove$edgemethod <- "wrap"
    expect_error(do.call(openCR.fit, argsmove)) 
    
    argsmove$mask <- make.mask(traps(smallCH), buffer = 200, nx = 20, 
        type = 'traprect')
    expect_equal(do.call(openCR.fit, argsmove)[1], -358.17621329, 
        tolerance = 1e-4, check.attributes = FALSE)
    
    argsmove$movementmodel <- "BVE"
    argsmove$details$r0 <- 0   # historical
    expect_equal(do.call(openCR.fit, argsmove)[1], -358.28315775, 
        tolerance = 1e-4, check.attributes = FALSE)
    
    argsmove$movementmodel <- "BVT"
    argsmove$start$move.b <- 5
    argsmove$details$r0 <- 0.5
    expect_equal(do.call(openCR.fit, argsmove)[1], -357.34627618, 
        tolerance = 1e-4, check.attributes = FALSE)

    # sparse kernel 2.0.0
    expect_equal(do.call(openCR.fit, argsmoves)[1], -358.326036775 , 
        tolerance = 1e-4, check.attributes = FALSE)
    
})

test_that("error if non-rectangular mask for wrapped movement", {
    expect_error(openCR.fit (smallCH, type='PLBsecrf', mask = msk, 
        movementmodel='normal', edgemethod = 'wrap'))
})

## 2021-06-13 check for openCR.esa bug of Trent McDonald

## fast 'fit' without maximization gives an openCR object to play with
smallfit <- openCR.fit(smallCH, mask = msk, type = 'PLBsecrf', 
    start = list(lambda0 = 0.037, sigma = 65.3, phi = 0.497, f = 0.560), 
    method = 'none')

test_that("correct effective sampling area", {
    expect_equal(openCR.esa(smallfit)[1], 13.47202, 
        tolerance = 1e-4, check.attributes = FALSE)
})

test_that("correct derived superpopulation density", {
    der <- derived(smallfit)
    expect_equal(der$superD, 2.078382, 
        tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(der$estimates$D[3],1.079069, 
        tolerance = 1e-4, check.attributes = FALSE)
    
})
