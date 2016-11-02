library(ranger)
library(survival)
context("ranger_extratrees")

test_that("extratrees splitting works for classification", {
  rf <- ranger(Species ~ ., iris, splitrule = "extratrees")
  expect_is(rf, "ranger")
  expect_lt(rf$prediction.error, 0.2)
})

test_that("extratrees unordered splitting works for classification", {
  rf <- ranger(Species ~ ., iris, splitrule = "extratrees", respect.unordered.factors = "partition")
  expect_is(rf, "ranger")
  expect_lt(rf$prediction.error, 0.2)
})

test_that("extratrees splitting works for probability estimation", {
  rf <- ranger(Species ~ ., iris, probability = TRUE, splitrule = "extratrees")
  expect_is(rf, "ranger")
  expect_lt(rf$prediction.error, 0.2)
})

test_that("extratrees unordered splitting works for probability estimation", {
  rf <- ranger(Species ~ ., iris, probability = TRUE, splitrule = "extratrees", respect.unordered.factors = "partition")
  expect_is(rf, "ranger")
  expect_lt(rf$prediction.error, 0.2)
})

test_that("extratrees splitting works for regression", {
  rf <- ranger(Sepal.Length ~ ., iris, splitrule = "extratrees")
  expect_is(rf, "ranger")
  expect_gt(rf$r.squared, 0.5)
})

test_that("extratrees unordered splitting works for regression", {
  rf <- ranger(Sepal.Length ~ ., iris, splitrule = "extratrees", respect.unordered.factors = "partition")
  expect_is(rf, "ranger")
  expect_gt(rf$r.squared, 0.5)
})

test_that("extratrees splitting works for survival", {
  rf <- ranger(Surv(time, status) ~ ., veteran, splitrule = "extratrees")
  expect_is(rf, "ranger")
  expect_lt(rf$prediction.error, 0.4)
})

test_that("extratrees unordered splitting works for survival", {
  rf <- ranger(Surv(time, status) ~ ., veteran, splitrule = "extratrees", respect.unordered.factors = "partition")
  expect_is(rf, "ranger")
  expect_lt(rf$prediction.error, 0.4)
})
