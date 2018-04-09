library(ranger)
library(survival)
context("genabel")

## GenABEL
if (!requireNamespace("GenABEL", quietly = TRUE)) {
  stop("Package GenABEL is required for testing ranger completely. Please install it.", call. = FALSE)
} else {
  dat.gwaa <- readRDS("../test_gwaa.Rds")
  rg.gwaa <- ranger(CHD ~ ., data = dat.gwaa, verbose = FALSE, write.forest = TRUE)
}

test_that("classification gwaa rf is of class ranger with 15 elements", {
  expect_is(rg.gwaa, "ranger")
  expect_equal(length(rg.gwaa), 15)
})

test_that("GenABEL prediction works if no covariates and formula used", {
  dat <- dat.gwaa
  dat@phdata$Age <- NULL
  rf <- ranger(CHD ~ .-Sex, data = dat, num.trees = 5)
  expect_silent(predict(rf, dat))
})

test_that("SNP ordering working for binary classification", {
  rf <- ranger(CHD ~ ., data = dat.gwaa, num.trees = 5, respect.unordered.factors = "order")
  expect_is(rf$forest$snp.order, "list")
  expect_length(rf$forest$snp.order, nsnps(dat.gwaa))
})

test_that("SNP ordering working for regression", {
  rf <- ranger(Age ~ ., data = dat.gwaa, num.trees = 5, respect.unordered.factors = "order")
  expect_is(rf$forest$snp.order, "list")
  expect_length(rf$forest$snp.order, nsnps(dat.gwaa))
})

test_that("SNP ordering not working for multiclass", {
  dat <- dat.gwaa
  dat@phdata$mc <- factor(round(runif(nids(dat), 1, 5)))
  expect_error(ranger(mc ~ ., data = dat, num.trees = 5, respect.unordered.factors = "order"), 
               "Error: Ordering of SNPs currently only implemented for regression and binary outcomes.")
})

test_that("SNP ordering not working for survival", {
  dat <- dat.gwaa
  dat@phdata$time <- runif(nids(dat), 1, 100)
  dat@phdata$status <- rbinom(nids(dat), 1, 0.5)
  expect_error(ranger(Surv(time, status) ~ ., data = dat, num.trees = 5, respect.unordered.factors = "order"), 
               "Error: Ordering of SNPs currently only implemented for regression and binary outcomes.")
})
