
context("ranger_cpp_classification")

test_that("Prediction is equal to R version", {
  # C++ version
  rf <- ranger_cpp(data = iris, depvarname = "Species", ntree = 5, write = "", seed = 10)
  pred <- ranger_cpp(data = iris, predict = "ranger_out.forest", seed = 20)
  preds_cpp <- as.data.frame(fread("ranger_out.prediction"))[, 1]
  
  # R version
  rf <- ranger(Species ~ ., iris, num.trees = 5, seed = 10)
  preds_r <- as.numeric(predict(rf, iris, seed = 20)$predictions)
  
  expect_equal(preds_cpp, preds_r)
})

test_that("Predictions are positive numbers", {
  rf <- ranger_cpp(data = iris, depvarname = "Species", ntree = 5, write = "", seed = 10)
  pred <- ranger_cpp(data = iris, predict = "ranger_out.forest")
  preds_cpp <- as.data.frame(fread("ranger_out.prediction"))[, 1]
  expect_is(preds_cpp, "integer")
  expect_true(all(preds_cpp > 0))
})

test_that("Same result with default splitting", {
  # C++ version
  rf <- ranger_cpp(data = iris, depvarname = "Species", ntree = 5, seed = 10)
  err_cpp <- grep("Overall OOB prediction error", rf, value = TRUE)
  err_cpp <- as.numeric(gsub("[^0-9.]", "", err_cpp))
  
  # R version
  rf <- ranger(Species ~ ., iris, num.trees = 5, seed = 10)
  err_r <- rf$prediction.error
  
  expect_equal(round(err_cpp, 4), round(err_r, 4))
})

test_that("Same result with extratrees splitting", {
  # C++ version
  rf <- ranger_cpp(data = iris, depvarname = "Species", ntree = 5, splitrule = 5, seed = 10)
  err_cpp <- grep("Overall OOB prediction error", rf, value = TRUE)
  err_cpp <- as.numeric(gsub("[^0-9.]", "", err_cpp))
  
  # R version
  rf <- ranger(Species ~ ., iris, num.trees = 5, splitrule = "extratrees", seed = 10)
  err_r <- rf$prediction.error
  
  expect_equal(round(err_cpp, 4), round(err_r, 4))
})
