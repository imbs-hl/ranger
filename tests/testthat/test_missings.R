library(ranger)
library(survival)
context("ranger_unordered")

test_that("Third child for missings only there if missings in data", {
  rf1 <- ranger(Species ~ ., iris, num.trees = 5)
  expect_length(rf1$forest$child.nodeIDs[[1]], 2)
  
  dat <- iris
  dat[1, 1] <- NA
  rf2 <- ranger(Species ~ ., dat, num.trees = 5)
  expect_length(rf2$forest$child.nodeIDs[[1]], 3)
})

test_that("Training works with missing values in x but not in y", {
  dat <- iris
  dat[25, 1] <- NA
  expect_silent(ranger(Species ~ ., dat, num.trees = 5))
  expect_silent(ranger(Petal.Width ~ ., dat, num.trees = 5))
  expect_error(ranger(Sepal.Length ~ ., dat, num.trees = 5), "Missing data in dependent variable.")
  
  dat <- iris
  dat[4, 5] <- NA
  expect_error(ranger(Species ~ ., dat, num.trees = 5), "Missing data in dependent variable.")
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

test_that("Prediction works with missing values, classification", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  
  dat <- iris
  dat[4, 4] <- NA
  dat[25, 1] <- NA
  expect_silent(predict(rf, dat))
})

test_that("Prediction works with missing values, regression", {
  rf <- ranger(Sepal.Width ~ ., iris, num.trees = 5, write.forest = TRUE)
  
  dat <- iris
  dat[4, 4] <- NA
  dat[25, 1] <- NA
  expect_silent(predict(rf, dat))
})

test_that("Order splitting working with missing values for classification", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D", NA), n, replace = TRUE), 
                   y = factor(rbinom(n, 1, 0.5)),
                   stringsAsFactors = FALSE)
  
  rf <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, respect.unordered.factors = 'order')
  expect_true(all(rf$forest$is.ordered))
})

test_that("Order splitting working with missing values for multiclass classification", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D", NA), n, replace = TRUE), 
                   y = factor(sample(c("A", "B", "C", "D"), n, replace = TRUE)),
                   stringsAsFactors = FALSE)
  
  rf <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, respect.unordered.factors = 'order')
  expect_true(all(rf$forest$is.ordered))
})

test_that("Missing values for survival not yet working", {
  dat <- veteran
  dat[1, 1] <- NA
  
  expect_error(ranger(Surv(time, status) ~ ., dat, num.trees = 5), "Error: Missing value handling not yet implemented for survival forests\\.")
})
