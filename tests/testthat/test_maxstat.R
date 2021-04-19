library(rangerts)
library(survival)
context("rangerts_maxstat")

test_that("maxstat splitting works for survival", {
  rf <- rangerts(Surv(time, status) ~ ., veteran, splitrule = "maxstat")
  expect_is(rf, "rangerts")
  expect_lt(rf$prediction.error, 0.4)
})

test_that("maxstat splitting works for regression", {
  rf <- rangerts(Sepal.Length ~ ., iris, splitrule = "maxstat")
  expect_is(rf, "rangerts")
  expect_gt(rf$r.squared, 0.5)
})

test_that("maxstat splitting, alpha or minprop out of range throws error", {
  expect_error(rangerts(Surv(time, status) ~ ., veteran, splitrule = "maxstat", alpha = -1))
  expect_error(rangerts(Surv(time, status) ~ ., veteran, splitrule = "maxstat", alpha = 2))
  expect_error(rangerts(Surv(time, status) ~ ., veteran, splitrule = "maxstat", minprop = -1))
  expect_error(rangerts(Surv(time, status) ~ ., veteran, splitrule = "maxstat", minprop = 1))
})

test_that("maxstat splitting not working for classification", {
  expect_error(rangerts(Species ~ ., iris, splitrule = "maxstat"))
})

test_that("maxstat impurity importance is positive", {
  rf <- rangerts(Surv(time, status) ~ ., veteran, num.trees = 5,
               splitrule = "maxstat", importance = "impurity")
  expect_gt(mean(rf$variable.importance), 0)

  rf <- rangerts(Sepal.Length ~ ., iris, num.trees = 5,
               splitrule = "maxstat", importance = "impurity")
  expect_gt(mean(rf$variable.importance), 0)
})

test_that("maxstat corrected impurity importance is positive (on average)", {
  rf <- rangerts(Surv(time, status) ~ ., veteran, num.trees = 50,
               splitrule = "maxstat", importance = "impurity_corrected")
  expect_gt(mean(rf$variable.importance), 0)

  rf <- rangerts(Sepal.Length ~ ., iris, num.trees = 5,
               splitrule = "maxstat", importance = "impurity_corrected")
  expect_gt(mean(rf$variable.importance), 0)
})
