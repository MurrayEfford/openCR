## test-stratified.R
## started 2021-05-03

## stratification openCR >= 2.0.0

library(openCR)
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

# create small working proximity dataset
suppressWarnings(smallCH <- join(subset(ovenCHp, sessions = 1:3, traps = 1:20, dropnullocc = FALSE)))

suppressWarnings(stratifiedch3 <- stratify (smallCH, covariate = 'Sex'))
suppressWarnings(ch2 <- subset(ovenCHp, sessions = 1:3, traps = 1:20, dropnullocc = FALSE))
stratifiedch2 <- stratify (ch2[1:3])

test_that("stratify generates suitable data for openCR.fit", {
    ## case 3: split joined data by covariate
    expect_equal(summary(stratifiedch3, terse = TRUE), 
        matrix(c(29,21,15,20,29,39,13,20), nrow = 4),
        check.attributes = FALSE)

    ## case 2: split joined data by covariate
    expect_equal(summary(stratifiedch2, terse = TRUE), 
        matrix(c(9,18,11,20,10,21,11,20,10,21,13,20), nrow = 4),
        check.attributes = FALSE)
})

test_that("correct stratified likelihood", {
    args <- list(capthist = stratifiedch3, type = 'PLBf', 
        start = list(p = 0.4, phi = 0.497, f = 0.560),
        model = p ~ stratum, details = list(LLonly = TRUE),
        stratified = TRUE)
    expect_equal(do.call(openCR.fit, args)[1], -246.589008072, 
        tolerance = 1e-4, check.attributes = FALSE)
})

