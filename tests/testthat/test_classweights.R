## Tests for class weights

library(ranger)
context("ranger_classweights")

test_that("No error if class weights used", {
  expect_silent(ranger(Species ~ ., iris, num.trees = 5, class.weights = c(0.5, 1, 0.1)))
})

test_that("Prediction accuracy for minority class increases with higher weight", {
  n <- 100
  x <- rnorm(n)
  beta0 <- -3
  beta <- 1
  y <- factor(rbinom(n, 1, plogis(beta0 + beta * x)))
  dat <- data.frame(y = y, x)

  rf <- ranger(y ~ ., dat, num.trees = 5, min.node.size = 50, class.weights = c(1, 1))
  acc_major <- mean((rf$predictions == dat$y)[dat$y == "0"], na.rm = TRUE)
  acc_minor <- mean((rf$predictions == dat$y)[dat$y == "1"], na.rm = TRUE)

  rf <- ranger(y ~ ., dat, num.trees = 5, min.node.size = 50, class.weights = c(0.01, 0.99))
  acc_major_weighted <- mean((rf$predictions == dat$y)[dat$y == "0"], na.rm = TRUE)
  acc_minor_weighted <- mean((rf$predictions == dat$y)[dat$y == "1"], na.rm = TRUE)

  expect_gt(acc_minor_weighted, acc_minor)
})

# test_that("Prediction error worse if class weights 0", {
#   rf <- ranger(Species ~ ., iris, num.trees = 5)
#   rf_weighted <- ranger(Species ~ ., iris, num.trees = 5, class.weights = c(1, 0, 0))
#   expect_lt(rf$prediction.error, rf_weighted$prediction.error)
# })

test_that("Error if class weights of wrong size", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, class.weights = c(0.5, 1)),
               "Error: Number of class weights not equal to number of classes.")
})

test_that("Error if class weights NA", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, class.weights = c(0.5, 1, NA)),
               "missing value where TRUE/FALSE needed")
})

test_that("Error if class weights not numeric", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, class.weights = c(0.5, 1, "a")),
               "Error: Invalid value for class.weights. Please give a vector of non-negative values.")
})


