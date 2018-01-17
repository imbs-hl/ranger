library(ranger)
library(survival)
context("ranger")

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

test_that("Matrix interface works for Probability estimation", {
  rf <- ranger(dependent.variable.name = "Species", data = data.matrix(iris), write.forest = TRUE, probability = TRUE)
  expect_equal(rf$treetype, "Probability estimation")
  expect_equal(rf$forest$independent.variable.names, colnames(iris)[1:4])
})

test_that("Matrix interface prediction works for Probability estimation", {
  dat <- data.matrix(iris)
  rf <- ranger(dependent.variable.name = "Species", data = dat, write.forest = TRUE, probability = TRUE)
  expect_silent(predict(rf, dat))
})

test_that("no warning if data.frame has two classes", {
  dat <- iris
  class(dat) <- c("data.frame", "data.table")
  expect_silent(ranger(Species ~ ., data = dat, verbose = FALSE))
})

test_that("Error if sample fraction is 0 or >1", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = 0))
  expect_error(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = 1.1))
})

test_that("Error if sample fraction is vector for regression", {
  expect_error(ranger(Sepal.Length ~ ., iris, num.trees = 5, sample.fraction = c(0.1, 0.2)), 
               "Error: Invalid value for sample\\.fraction\\. Vector values only valid for classification forests\\.")
})

test_that("Error if sample fraction is vector of wrong size", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.1, 0.2)), 
               "Error: Invalid value for sample\\.fraction\\. Expecting 3 values, provided 2\\.")
})

test_that("Error if element of sample fraction vector is <0 or >1", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.1, 1.1, 0.3)), 
               "Error: Invalid value for sample\\.fraction. Please give a value in \\(0,1\\] or a vector of values in \\[0,1\\]\\.")
  expect_error(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(-3, 0.5, 0.3)), 
               "Error: Invalid value for sample.fraction. Please give a value in \\(0,1] or a vector of values in \\[0,1\\]\\.")
})

test_that("Error if sum of sample fraction vector is 0", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0, 0, 0)), 
               "Error: Invalid value for sample\\.fraction. Sum of values must be >0\\.")
})

test_that("Error if replace=FALSE and not enough samples", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.2, 0.3, 0.4), 
                      replace = FALSE, keep.inbag = TRUE), 
               "Error: Not enough samples in class virginica; available: 50, requested: 60.")
  expect_silent(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.2, 0.3, 0.4), 
                       replace = TRUE, keep.inbag = TRUE))
})

test_that("Error if sample.fraction and case.weights", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.2, 0.3, 0.4), 
                      case.weights = rbinom(nrow(iris), 1, 0.5)), 
               "Error: Combination of case\\.weights and class-wise sampling not supported\\.")
})

test_that("Inbag counts match sample fraction, classification", {
  ## With replacement
  rf <- ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.2, 0.3, 0.4), 
               replace = TRUE, keep.inbag = TRUE)
  inbag <- do.call(cbind, rf$inbag.counts)
  expect_equal(unique(colSums(inbag[1:50, ])), 30)
  expect_equal(unique(colSums(inbag[51:100, ])), 45)
  expect_equal(unique(colSums(inbag[101:150, ])), 60)
  
  ## Without replacement
  rf <- ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.1, 0.2, 0.3), 
               replace = FALSE, keep.inbag = TRUE)
  inbag <- do.call(cbind, rf$inbag.counts)
  expect_equal(unique(colSums(inbag[1:50, ])), 15)
  expect_equal(unique(colSums(inbag[51:100, ])), 30)
  expect_equal(unique(colSums(inbag[101:150, ])), 45)
})

test_that("Inbag counts match sample fraction, probability", {
  ## With replacement
  rf <- ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.2, 0.3, 0.4), 
               replace = TRUE, keep.inbag = TRUE, probability = TRUE)
  inbag <- do.call(cbind, rf$inbag.counts)
  expect_equal(unique(colSums(inbag[1:50, ])), 30)
  expect_equal(unique(colSums(inbag[51:100, ])), 45)
  expect_equal(unique(colSums(inbag[101:150, ])), 60)
  
  ## Without replacement
  rf <- ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.1, 0.2, 0.3), 
               replace = FALSE, keep.inbag = TRUE, probability = TRUE)
  inbag <- do.call(cbind, rf$inbag.counts)
  expect_equal(unique(colSums(inbag[1:50, ])), 15)
  expect_equal(unique(colSums(inbag[51:100, ])), 30)
  expect_equal(unique(colSums(inbag[101:150, ])), 45)
})

test_that("as.factor() in formula works", {
  n <- 20
  dt <- data.frame(x = runif(n), y = rbinom(n, 1, 0.5))
  expect_silent(ranger(as.factor(y) ~ ., data = dt, num.trees = 5, write.forest = TRUE))
})

test_that("holdout mode holding out data with 0 weight", {
  weights <- rbinom(nrow(iris), 1, 0.5)
  rf <- ranger(Species ~ ., iris, num.trees = 5, importance = "permutation",  
               case.weights = weights, replace = FALSE, sample.fraction = 0.632*mean(weights), 
               holdout = TRUE, keep.inbag = TRUE)
  inbag <- data.frame(rf$inbag.counts)
  expect_true(all(inbag[weights == 0, ] == 0))
})

test_that("holdout mode uses holdout OOB data", {
  weights <- rbinom(nrow(iris), 1, 0.5)
  rf <- ranger(Species ~ ., iris, num.trees = 5, importance = "permutation",  
               case.weights = weights, replace = FALSE, sample.fraction = 0.632*mean(weights), 
               holdout = TRUE, keep.inbag = TRUE)
  expect_false(any(is.na(rf$predictions[weights == 0])))
  expect_true(all(is.na(rf$predictions[weights == 1])))
})

test_that("holdout mode not working if no weights", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, importance = "permutation", holdout = TRUE))
})

test_that("holdout mode: no OOB prediction if no 0 weights", {
  weights <- runif(nrow(iris))
  rf <- ranger(Species ~ ., iris, num.trees = 5, importance = "permutation",  
               case.weights = weights, replace = FALSE, 
               holdout = TRUE, keep.inbag = TRUE)
  expect_true(all(is.na(rf$predictions)))
})

test_that("OOB error is correct for 1 tree, classification", {
  n <- 50
  dat <- data.frame(y = factor(rbinom(n, 1, .5)), x = rnorm(n))
  rf <- ranger(y ~ ., dat, num.trees = 1)
  expect_equal(rf$prediction.error, mean(rf$predictions != dat$y, na.rm = TRUE))
})

test_that("OOB error is correct for 1 tree, probability prediction", {
  n <- 50
  dat <- data.frame(y = factor(rbinom(n, 1, .5)), x = rnorm(n))
  rf <- ranger(y ~ ., dat, num.trees = 1, probability = TRUE)
  prob <- c(rf$predictions[dat$y == "0", 1], rf$predictions[dat$y == "1", 2])
  expect_equal(rf$prediction.error, mean((1 - prob)^2, na.rm = TRUE))
})

test_that("OOB error is correct for 1 tree, regression", {
  n <- 50
  dat <- data.frame(y = rbinom(n, 1, .5), x = rnorm(n))
  rf <- ranger(y ~ ., dat, num.trees = 1)
  expect_equal(rf$prediction.error, mean((dat$y - rf$predictions)^2, na.rm = TRUE))
})

test_that("Missing value columns detected in training", {
  dat <- iris
  dat[4, 5] <- NA
  dat[25, 1] <- NA
  expect_error(ranger(Species ~ ., dat, num.trees = 5), "Missing data in columns: Species, Sepal.Length")
})

test_that("No error if missing value in irrelevant column, training", {
  dat <- iris
  dat[1, "Sepal.Width"] <- NA
  expect_silent(ranger(Species ~ Sepal.Length, dat, num.trees = 5))
})

test_that("No error if missing value in irrelevant column, prediction", {
  rf <- ranger(Species ~ Sepal.Length, iris, num.trees = 5)
  dat <- iris
  dat[1, "Sepal.Width"] <- NA
  expect_silent(predict(rf, dat))
})

test_that("Split points are at (A+B)/2 for numeric features, regression variance splitting", {
  dat <- data.frame(y = rbinom(100, 1, .5), x = rbinom(100, 1, .5))
  rf <- ranger(y ~ x, dat, num.trees = 10)
  split_points <- mapply(function(varID, value) {
    value[varID > 0]
    }, 
    rf$forest$split.varIDs, 
    rf$forest$split.values
  )
  expect_equal(split_points, rep(0.5, rf$num.trees))
})

test_that("Split points are at (A+B)/2 for numeric features, regression maxstat splitting", {
  dat <- data.frame(y = rbinom(100, 1, .5), x = rbinom(100, 1, .5))
  rf <- ranger(y ~ x, dat, num.trees = 10, splitrule = "maxstat", alpha = 1)
  split_points <- mapply(function(varID, value) {
    value[varID > 0]
    }, 
    rf$forest$split.varIDs, 
    rf$forest$split.values
  )
  expect_equal(split_points, rep(0.5, rf$num.trees))
})

test_that("Split points are at (A+B)/2 for numeric features, classification", {
  dat <- data.frame(y = factor(rbinom(100, 1, .5)), x = rbinom(100, 1, .5))
  rf <- ranger(y ~ x, dat, num.trees = 10)
  split_points <- mapply(function(varID, value) {
    value[varID > 0]
  }, 
  rf$forest$split.varIDs, 
  rf$forest$split.values
  )
  expect_equal(split_points, rep(0.5, rf$num.trees))
})

test_that("Split points are at (A+B)/2 for numeric features, probability", {
  dat <- data.frame(y = factor(rbinom(100, 1, .5)), x = rbinom(100, 1, .5))
  rf <- ranger(y ~ x, dat, num.trees = 10, probability = TRUE)
  split_points <- mapply(function(varID, value) {
    value[varID > 0]
  }, 
  rf$forest$split.varIDs, 
  rf$forest$split.values
  )
  expect_equal(split_points, rep(0.5, rf$num.trees))
})

test_that("Split points are at (A+B)/2 for numeric features, survival logrank splitting", {
  dat <- data.frame(time = runif(100, 1, 10), status = rbinom(100, 1, .5), x = rbinom(100, 1, .5))
  rf <- ranger(Surv(time, status) ~ x, dat, num.trees = 10, splitrule = "logrank")
  split_points <- mapply(function(varID, value) {
    value[varID > 0]
  }, 
  rf$forest$split.varIDs, 
  rf$forest$split.values
  )
  expect_equal(split_points, rep(0.5, rf$num.trees))
})

test_that("Split points are at (A+B)/2 for numeric features, survival C-index splitting", {
  dat <- data.frame(time = runif(100, 1, 10), status = rbinom(100, 1, .5), x = rbinom(100, 1, .5))
  rf <- ranger(Surv(time, status) ~ x, dat, num.trees = 10, splitrule = "C")
  split_points <- mapply(function(varID, value) {
    value[varID > 0]
  }, 
  rf$forest$split.varIDs, 
  rf$forest$split.values
  )
  expect_equal(split_points, rep(0.5, rf$num.trees))
})

test_that("Split points are at (A+B)/2 for numeric features, survival maxstat splitting", {
  dat <- data.frame(time = runif(100, 1, 10), status = rbinom(100, 1, .5), x = rbinom(100, 1, .5))
  rf <- ranger(Surv(time, status) ~ x, dat, num.trees = 10, splitrule = "maxstat", alpha = 1)
  split_points <- mapply(function(varID, value) {
    value[varID > 0]
  }, 
  rf$forest$split.varIDs, 
  rf$forest$split.values
  )
  expect_equal(split_points, rep(0.5, rf$num.trees))
})

test_that("No error if variable named forest", {
  dat <- iris
  dat$forest <- rnorm(150)
  rf <- ranger(Species ~ ., dat, num.trees = 5)
  expect_silent(predict(rf, dat))
})

test_that("GenABEL prediction works if no covariates and formula used", {
  dat <- dat.gwaa
  dat@phdata$Age <- NULL
  rf <- ranger(CHD ~ .-Sex, data = dat, num.trees = 5)
  expect_silent(predict(rf, dat))
})
