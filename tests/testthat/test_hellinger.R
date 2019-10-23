library(ranger)
library(survival)
context("ranger_hellinger")

test_that("Hellinger splitting works for classification", {
  rf <- ranger(Species ~ ., droplevels(iris[1:100, ]), splitrule = "hellinger", num.trees = 5)
  expect_is(rf, "ranger")
  expect_lt(rf$prediction.error, 0.3)
})

test_that("Hellinger splitting works for classification with non-factor response", {
  dat <- iris[1:100, ]
  dat$Species <- as.numeric(dat$Species)
  rf <- ranger(Species ~ ., dat, num.trees = 5, classification = TRUE, splitrule = "hellinger")
  expect_is(rf, "ranger")
  expect_lt(rf$prediction.error, 0.3)
})

test_that("Hellinger splitting works for probability estimation", {
  rf <- ranger(Species ~ ., droplevels(iris[1:100, ]), splitrule = "hellinger", 
               num.trees = 5, probability = TRUE)
  expect_is(rf, "ranger")
  expect_lt(rf$prediction.error, 0.3)
})

test_that("Hellinger splitting not working for multiclass", {
  expect_error(ranger(Species ~ ., iris, splitrule = "hellinger"), 
               "Error: Hellinger splitrule only implemented for binary classification\\.")
})

test_that("Hellinger splitting not working for regression", {
  expect_error(ranger(Sepal.Length ~ ., iris, splitrule = "hellinger"), 
               "Error: Hellinger splitrule only implemented for binary classification\\.")
})
