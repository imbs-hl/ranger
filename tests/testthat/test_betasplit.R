library(ranger)
context("ranger_betasplit")

# Generate data with 0..1 outcome
n <- 100
p <- 4
beta <- c(0, 1, 2, 3)
x <- replicate(p, runif(n))
y <- as.vector(x %*% beta)
y <- (y-min(y))/(max(y)-min(y))
dat <- data.frame(y = y, x)

test_that("beta splitting works for regression", {
  rf <- ranger(y ~ ., dat, splitrule = "beta", num.trees = 50)
  expect_is(rf, "ranger")
  expect_lt(rf$prediction.error, 0.2)
})

test_that("beta splitting not working for non 0..1 outcome", {
  expect_error(ranger(Sepal.Length ~ ., iris, splitrule = "beta", num.trees = 50))
})