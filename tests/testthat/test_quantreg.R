library(ranger)
context("ranger_quantreg")

rf.quant <- ranger(mpg ~ ., mtcars[1:26, ], quantreg = TRUE, 
                   keep.inbag = TRUE, num.trees = 50)
pred.quant <- predict(rf.quant, mtcars[27:32, ], type = "quantiles")

test_that("Quantile prediction is of correct size", {
  expect_equal(dim(pred.quant$predictions), c(pred.quant$num.samples, 3))
})

test_that("Quantile OOB prediction is of correct size", {
  pred <- predict(rf.quant, data = NULL, type = "quantiles")
  expect_equal(dim(pred$predictions), c(rf.quant$num.samples, 3))
})

test_that("Lower quantile smaller central smaller upper", {
  expect_true(all(pred.quant$predictions[, 1] < pred.quant$predictions[, 2]))
  expect_true(all(pred.quant$predictions[, 2] < pred.quant$predictions[, 3]))
})

test_that("Working for single quantile", {
  expect_silent(pred <- predict(rf.quant, data = mtcars[27:32, ], type = "quantiles", quantiles = .5))
  expect_silent(pred <- predict(rf.quant, data = NULL, type = "quantiles", quantiles = .5))
})

test_that("Working for single new observation", {
  expect_silent(pred <- predict(rf.quant, mtcars[27, ], type = "quantiles"))
})

test_that("Working for constant variables", {
  dat <- data.frame(x1 = 1, x2 = seq(1,10), y = seq(1,10))
  rf <- ranger(y ~ ., dat, quantreg = TRUE)
  expect_silent(predict(rf, dat, type = "quantiles"))
})

test_that("Error message if no regression forest", {
  rf <- ranger(Species ~ ., iris, num.trees = 5)
  expect_error(predict(rf, iris, type = "quantiles"), 
               "Error\\: Quantile prediction implemented only for regression outcomes\\.")
})

test_that("Error message if not grown with quantreg=TRUE", {
  rf <- ranger(mpg ~ ., mtcars[1:26, ], quantreg = FALSE, num.trees = 50)
  expect_error(predict(rf, mtcars[27:32, ], type = "quantiles"), 
               "Error\\: Set quantreg\\=TRUE in ranger\\(\\.\\.\\.\\) for quantile prediction\\.")
})

test_that("User specified function works as expected", {
  pred_sample <- predict(rf.quant, mtcars[27:32, ], type = "quantiles", what = function(x) sample(x, 10, replace = TRUE))
  expect_equal(dim(pred_sample$predictions), c(pred_sample$num.samples, 10))
})

test_that("Working for factor variables", {
  expect_silent(rf <- ranger(Sepal.Length ~ ., iris, quantreg = TRUE))
  expect_silent(predict(rf, iris, type = "quantiles"))
})

test_that("Working for unordered factor variables", {
  expect_silent(rf <- ranger(Sepal.Length ~ ., iris, quantreg = TRUE, respect.unordered.factors = "order"))
  expect_silent(predict(rf, iris, type = "quantiles"))
})
