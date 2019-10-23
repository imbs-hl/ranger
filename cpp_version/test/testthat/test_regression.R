
context("ranger_cpp_regression")

test_that("Prediction is equal to R version", {
  # C++ version
  rf <- ranger_cpp(data = iris, depvarname = "Sepal.Length", treetype = 3, ntree = 5, write = "", seed = 10)
  pred <- ranger_cpp(data = iris, treetype = 3, predict = "ranger_out.forest", seed = 20)
  preds_cpp <- as.data.frame(fread("ranger_out.prediction"))[, 1]
  
  # R version
  rf <- ranger(Sepal.Length ~ ., iris, num.trees = 5, seed = 10)
  preds_r <- as.numeric(predict(rf, iris, seed = 20)$predictions)
  
  expect_equal(round(preds_cpp, 4), round(preds_r, 4))
})

test_that("Predictions are in range of original data", {
  rf <- ranger_cpp(data = iris, depvarname = "Sepal.Width", treetype = 3, ntree = 5, write = "", seed = 10)
  pred <- ranger_cpp(data = iris, predict = "ranger_out.forest", treetype = 3)
  preds_cpp <- as.data.frame(fread("ranger_out.prediction"))[, 1]
  expect_is(preds_cpp, "numeric")
  expect_true(all(preds_cpp >= min(iris$Sepal.Width)))
  expect_true(all(preds_cpp <= max(iris$Sepal.Width)))
})

test_that("Same result with default splitting", {
  # C++ version
  rf <- ranger_cpp(data = iris, depvarname = "Petal.Length", treetype = 3, ntree = 5, seed = 10)
  err_cpp <- grep("Overall OOB prediction error", rf, value = TRUE)
  err_cpp <- as.numeric(gsub("[^0-9.]", "", err_cpp))
  
  # R version
  rf <- ranger(Petal.Length ~ ., iris, num.trees = 5, seed = 10)
  err_r <- rf$prediction.error
  
  expect_equal(round(err_cpp, 4), round(err_r, 4))
})

test_that("Same result with extratrees splitting", {
  # C++ version
  rf <- ranger_cpp(data = iris, depvarname = "Sepal.Length", treetype = 3, ntree = 5, splitrule = 5, catvars = "Species", seed = 10)
  err_cpp <- grep("Overall OOB prediction error", rf, value = TRUE)
  err_cpp <- as.numeric(gsub("[^0-9.]", "", err_cpp))
  
  # R version
  rf <- ranger(Sepal.Length ~ ., iris, num.trees = 5, splitrule = "extratrees", seed = 10)
  err_r <- rf$prediction.error
  
  expect_equal(round(err_cpp, 4), round(err_r, 4))
})

test_that("Same result with maxstat splitting", {
  # C++ version
  rf <- ranger_cpp(data = iris, depvarname = "Sepal.Length", treetype = 3, ntree = 5, splitrule = 4, seed = 10)
  err_cpp <- grep("Overall OOB prediction error", rf, value = TRUE)
  err_cpp <- as.numeric(gsub("[^0-9.]", "", err_cpp))
  
  # R version
  rf <- ranger(Sepal.Length ~ ., iris, num.trees = 5, splitrule = "maxstat", seed = 10)
  err_r <- rf$prediction.error
  
  expect_equal(round(err_cpp, 4), round(err_r, 4))
})

test_that("Same result with beta splitting", {
  # Generate data with 0..1 outcome
  n <- 100
  p <- 4
  beta <- c(0, 1, 2, 3)
  x <- round(replicate(p, runif(n)), 3)
  y <- as.vector(x %*% beta)
  y <- round((y-min(y))/(max(y)-min(y)), 3)
  dat <- data.frame(y = y, x)
  
  # C++ version
  rf <- ranger_cpp(data = dat, depvarname = "y", treetype = 3, ntree = 5, splitrule = 6, seed = 10)
  err_cpp <- grep("Overall OOB prediction error", rf, value = TRUE)
  err_cpp <- as.numeric(gsub("[^0-9.]", "", err_cpp))
  
  # R version
  rf <- ranger(y ~ ., dat, num.trees = 5, splitrule = "beta", seed = 10)
  err_r <- rf$prediction.error
  
  expect_equal(round(err_cpp, 4), round(err_r, 4))
})

test_that("Same result with unordered splitting", {
  # C++ version
  rf <- ranger_cpp(data = iris, depvarname = "Sepal.Length", treetype = 3, ntree = 5, catvars = "Species", seed = 10)
  err_cpp <- grep("Overall OOB prediction error", rf, value = TRUE)
  err_cpp <- as.numeric(gsub("[^0-9.]", "", err_cpp))
  
  # R version
  rf <- ranger(Sepal.Length ~ ., iris, num.trees = 5, respect.unordered.factors = "partition", seed = 10)
  err_r <- rf$prediction.error
  
  expect_equal(round(err_cpp, 4), round(err_r, 4))
})
