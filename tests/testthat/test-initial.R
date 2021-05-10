## Started 2020-12-13

## 2021-05-03 sparsekernel 
## 2021-05-03 stratification 
## 2021-05-04 use tinythread not TBB to avoid ASAN/UBSAN errors on CRAN

library(openCR)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

# create small working proximity dataset
suppressWarnings(smallCH <- join(subset(ovenCHp, sessions = 1:3, traps = 1:20, dropnullocc = FALSE)))

msk <- make.mask(traps(smallCH), buffer = 200, nx = 20, type = 'trapbuffer')

args <- list(capthist = smallCH,
    start = list(p = 0.4, phi = 0.497, f = 0.560),
    details = list(LLonly = TRUE))

argssecr <- list(capthist = smallCH, mask = msk,
    start = list(lambda0 = 0.037, sigma = 65.3, phi = 0.497, f = 0.560),
    movementmodel = "static", 
    details = list(LLonly = TRUE))

argsmove <- list(capthist = smallCH, mask = msk, type = "PLBsecrf",
    start = list(lambda0 = 0.037, sigma = 65.3, phi = 0.497, f = 0.560, move.a = 100),
    movementmodel = "normal", edgemethod = "truncate",
    details = list(LLonly = TRUE))

argsmoves <- list(capthist = smallCH, mask = msk, type = "PLBsecrf",
    start = list(lambda0 = 0.037, sigma = 65.3, phi = 0.497, f = 0.560, move.a = 100),
    movementmodel = "normal", edgemethod = "truncate", sparsekernel = TRUE,
    kernelradius = 20, details = list(LLonly = TRUE))

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
    
    argsmove$movementmodel <- "exponential"
    expect_equal(do.call(openCR.fit, argsmove)[1], -358.28315775, 
        tolerance = 1e-4, check.attributes = FALSE)
    
    argsmove$movementmodel <- "t2D"
    argsmove$start$move.b <- 5
    expect_equal(do.call(openCR.fit, argsmove)[1], -357.34627618, 
        tolerance = 1e-4, check.attributes = FALSE)

    # sparse kernel 2.0.0
    expect_equal(do.call(openCR.fit, argsmoves)[1], -358.326036775 , 
        tolerance = 1e-4, check.attributes = FALSE)
    
})

test_that("warning if invalid parameters in start list", {
    expect_warning(openCR.fit (smallCH, type = 'PLBb', start= list(f = 0.6)))
})

test_that("error if non-rectangular mask for wrapped movement", {
    expect_error(openCR.fit (smallCH, type='PLBsecrf', mask = msk, 
        movementmodel='normal', edgemethod = 'wrap'))
})

## stratification openCR >= 2.0.0

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
