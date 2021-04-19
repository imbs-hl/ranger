library(rangerts)
library(survival)
context("rangerts_extratrees")

test_that("extratrees splitting works for classification", {
  rf <- rangerts(Species ~ ., iris, splitrule = "extratrees")
  expect_is(rf, "rangerts")
  expect_lt(rf$prediction.error, 0.2)
})

test_that("extratrees unordered splitting works for classification", {
  n <- 20
  dat <- data.frame(x = sample(c("A", "B", "C", "D"), n, replace = TRUE),
                    y = factor(rbinom(n, 1, 0.5)),
                    stringsAsFactors = FALSE)
  rf <- rangerts(y ~ ., data = dat, num.trees = 5, min.node.size = n/2,
               splitrule = "extratrees", respect.unordered.factors = 'partition')

  expect_is(rf, "rangerts")
})

test_that("extratrees splitting works for probability estimation", {
  rf <- rangerts(Species ~ ., iris, probability = TRUE, splitrule = "extratrees")
  expect_is(rf, "rangerts")
  expect_lt(rf$prediction.error, 0.2)
})

test_that("extratrees unordered splitting works for probability estimation", {
  n <- 20
  dat <- data.frame(x = sample(c("A", "B", "C", "D"), n, replace = TRUE),
                    y = factor(rbinom(n, 1, 0.5)),
                    stringsAsFactors = FALSE)
  rf <- rangerts(y ~ ., data = dat, num.trees = 5, min.node.size = n/2,
               splitrule = "extratrees", respect.unordered.factors = 'partition',
               probability = TRUE)
  expect_is(rf, "rangerts")
})

test_that("extratrees splitting works for regression", {
  rf <- rangerts(Sepal.Length ~ ., iris, splitrule = "extratrees")
  expect_is(rf, "rangerts")
  expect_gt(rf$r.squared, 0.5)
})

test_that("extratrees unordered splitting works for regression", {
  rf <- rangerts(Sepal.Length ~ ., iris, splitrule = "extratrees", respect.unordered.factors = "partition")
  expect_is(rf, "rangerts")
  expect_gt(rf$r.squared, 0.5)
})

test_that("extratrees splitting works for survival", {
  rf <- rangerts(Surv(time, status) ~ ., veteran, splitrule = "extratrees")
  expect_is(rf, "rangerts")
  expect_lt(rf$prediction.error, 0.4)
})

test_that("extratrees unordered splitting works for survival", {
  rf <- rangerts(Surv(time, status) ~ ., veteran, splitrule = "extratrees", respect.unordered.factors = "partition")
  expect_is(rf, "rangerts")
  expect_lt(rf$prediction.error, 0.4)
})

test_that("extratrees splitting works for large number of random splits", {
  expect_silent(rangerts(Species ~ ., iris, splitrule = "extratrees", num.random.splits = 100))
})
