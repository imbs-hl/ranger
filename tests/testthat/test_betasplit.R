library(rangerts)
context("rangerts_betasplit")

# Generate data with 0..1 outcome
n <- 100
p <- 4
beta <- c(0, 1, 2, 3)
x <- replicate(p, runif(n))
y <- as.vector(x %*% beta)
y <- (y-min(y))/(max(y)-min(y))
dat <- data.frame(y = y, x)

test_that("beta splitting works for regression", {
  rf <- rangerts(y ~ ., dat, splitrule = "beta", num.trees = 50)
  expect_is(rf, "rangerts")
  expect_lt(rf$prediction.error, 0.2)
})

test_that("beta splitting not working for non 0..1 outcome", {
  expect_error(rangerts(Sepal.Length ~ ., iris, splitrule = "beta", num.trees = 50))
})
