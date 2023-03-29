## Started 2020-12-13
## 2023-03-29 test grouped-age CJS (details$agebreaks)

library(openCR)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

# create small working proximity dataset
suppressWarnings(smallCH <- join(subset(ovenCHp, sessions = 1:3, traps = 1:20, dropnullocc = FALSE)))

args <- list(capthist = smallCH,
    start = list(p = 0.4, phi = 0.497, f = 0.560),
    details = list(LLonly = TRUE))

test_that("correct non-spatial likelihood", {
    args$type = 'PLBf'
    expect_equal(do.call(openCR.fit, args)[1], -246.589008, 
        tolerance = 1e-4, check.attributes = FALSE)
})

test_that("derived.openCR works with fixed parameters", {
    args <- list(capthist = smallCH, type = 'PLBl', fixed = list(lambda = 1.0))
    fit <- do.call(openCR.fit, args)
    est <- derived(fit)$estimates
    expect_equal(est$phi, c(0.4560217,0.4560217,NA), tolerance = 1e-5)
})

test_that("warning if invalid parameters in start list", {
    expect_warning(openCR.fit (smallCH, type = 'PLBb', start= list(f = 0.6)))
})

# openCR 2.0 doubled LL and reported half variance
test_that("Pradel model results match Williams et al. Table 18.4", {
    suppressWarnings(fitmpradelg <- openCR.fit(microtusMCH, type = "Pradelg",
        model = list(p~t, phi~t, gamma~t)))
    pred <- predict(fitmpradelg)
    expect_equal(pred$gamma[3,'estimate'], 0.7052, tolerance = 1e-3)
    expect_equal(pred$gamma[3,'SE.estimate'], 0.068, tolerance = 1e-3)
})

test_that("Correct coefficients from grouped-age CJS of poss8088F", {
  # tests openCR.fit, predict.openCR and makeNewdata.openCR
  datadir <- system.file('extdata', package = 'openCR')
  CH <- readRDS(paste0(datadir,'/poss8088F.RDS'))
  fit <- openCR.fit(CH, model = list(phi ~ age), ncores = 2, details = list(
    agebreaks = c(0,3,6,Inf), initialage = 'age', maximumage = 6))
  # select Feb1980 estimates
  expect_equal(predict(fit,all.levels=TRUE)$phi$estimate[c(1,28,55)],  
    c(0.6660459, 0.9177610, 0.8211389), tolerance = 1e-4)
  
})
