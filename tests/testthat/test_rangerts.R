library(rangerts)
context("rangerts for time series")

test_that("Fonctional when time series parameters are active", {
  # create a periodical dataset
  dat <- iris
  dat$Sepal.Length <- dat$Sepal.Length + rep(1:10, 15) * 2
  dat$Time <- rep(1:10, 15)

  # ---- every time series bootstrap mode goes well ----

  # when bootstrap.ts is absent
  # but other time series parameters are given
  # the function takes by default iid setting
  expect_silent(rangerts(Sepal.Length ~ ., data = dat,
                       block.size = 10))

  # non overlapping boostrap
  expect_silent(rangerts(Sepal.Length ~ ., data = dat,
                       bootstrap.ts = "nonoverlapping",
                       block.size = 10))

  # moving block boostrap
  expect_silent(rangerts(Sepal.Length ~ ., data = dat,
                       bootstrap.ts = "moving",
                       block.size = 10))

  # circular block boostrap
  expect_silent(rangerts(Sepal.Length ~ ., data = dat,
                       bootstrap.ts = "circular",
                       block.size = 10))

  # seasonal block boostrap
  expect_silent(rangerts(Sepal.Length ~ ., data = dat,
                       bootstrap.ts = "seasonal",
                       block.size = 10))

  # ---- other details ----
  # block bootstrapped forest is different from iid one
  rf_iid <- rangerts(Sepal.Length ~ ., data = dat)
  rf_no  <- rangerts(Sepal.Length ~ ., data = dat,
                   bootstrap.ts = "nonoverlapping",
                   block.size = 10)
  expect_true(any(rf_iid$predictions != rf_no$predictions))

  # if bootstrap takes temporal order of observations
  # or takes the inversed order
  rf_temp <- rangerts(Sepal.Length ~ ., data = dat,
                   bootstrap.ts = "nonoverlapping",
                   block.size = 10, by.end = F)
  rf_inv  <- rangerts(Sepal.Length ~ ., data = dat,
                   bootstrap.ts = "nonoverlapping",
                   block.size = 10)
  expect_true(any(rf_inv$predictions != rf_temp$predictions))

})
