library(ranger)
library(survival)
context("ranger_deforest")


test_that("deforest works as expected for probability estimation", {
  rfo <- ranger(Species ~ ., data = iris, num.trees = 10, probability = TRUE)
  dfo <- deforest(rfo, which.trees = c(1, 3, 5), warn = FALSE)
  pred.rfo <- predict(rfo, data = iris, predict.all = TRUE)$predictions
  pred.dfo <- predict(dfo, data = iris, predict.all = TRUE)$predictions
  expect_identical(pred.rfo[, , -c(1, 3, 5)], pred.dfo)
})

test_that("deforest works as expected for classification", {
  rfo <- ranger(Species ~ ., data = iris, num.trees = 10)
  dfo <- deforest(rfo, which.trees = c(1, 3, 5), warn = FALSE)
  pred.rfo <- predict(rfo, data = iris, predict.all = TRUE)$predictions
  pred.dfo <- predict(dfo, data = iris, predict.all = TRUE)$predictions
  expect_identical(pred.rfo[, -c(1, 3, 5)], pred.dfo)
})

test_that("deforest works as expected for regression", {
  n <- 50
  x <- runif(n, min = 0, max = 2*pi)
  dat <- data.frame(x = x, y = sin(x) + rnorm(n, sd = 0.1))
  rfo <- ranger(y ~ ., data = dat, num.trees = 10)
  dfo <- deforest(rfo, which.trees = c(1, 3, 5), warn = FALSE)
  pred.rfo <- predict(rfo, data = dat, predict.all = TRUE)$predictions
  pred.dfo <- predict(dfo, data = dat, predict.all = TRUE)$predictions
  expect_identical(pred.rfo[, -c(1, 3, 5)], pred.dfo)
})

test_that("deforest works as expected for censored outcomes", {
  dat <- data.frame(time = runif(100, 1, 10), status = rbinom(100, 1, .5), 
                    x = rbinom(100, 1, .5))
  rfo <- ranger(Surv(time, status) ~ x, data = dat, num.trees = 10, 
                splitrule = "logrank")
  dfo <- deforest(rfo, which.trees = c(1, 3, 5), warn = FALSE)
  pred.rfo <- predict(rfo, data = dat, predict.all = TRUE)
  pred.dfo <- predict(dfo, data = dat, predict.all = TRUE)
  expect_identical(pred.rfo$chf[, , -c(1, 3, 5)], pred.dfo$chf)
  expect_identical(pred.rfo$survival[, , -c(1, 3, 5)], pred.dfo$survival)
})
