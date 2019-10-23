
context("ranger_cpp_probability")

test_that("Prediction is equal to R version", {
  # C++ version
  rf <- ranger_cpp(data = iris, depvarname = "Species", probability = "", ntree = 5, write = "", seed = 10)
  pred <- ranger_cpp(data = iris, predict = "ranger_out.forest", probability = "", seed = 20)
  preds_cpp <- as.matrix(fread("ranger_out.prediction"))
  colnames(preds_cpp) <- NULL
  
  # R version
  rf <- ranger(Species ~ ., iris, probability = TRUE, num.trees = 5, seed = 10)
  preds_r <- predict(rf, iris, seed = 20)$predictions
  colnames(preds_r) <- NULL
  
  expect_equal(round(preds_cpp, 4), round(preds_r, 4))
})

test_that("Predictions are probabilites", {
  rf <- ranger_cpp(data = iris, depvarname = "Species", probability = "", ntree = 5, write = "", seed = 10)
  pred <- ranger_cpp(data = iris, predict = "ranger_out.forest", probability = "")
  preds_cpp <- as.matrix(fread("ranger_out.prediction"))
  expect_is(preds_cpp, "matrix")
  expect_equal(dim(preds_cpp), c(150, 3))
  expect_true(all(preds_cpp >= 0))
  expect_true(all(preds_cpp <= 1))
})

test_that("Same result with default splitting", {
  # C++ version
  rf <- ranger_cpp(data = iris, depvarname = "Species", probability = "", ntree = 5, seed = 10)
  err_cpp <- grep("Overall OOB prediction error", rf, value = TRUE)
  err_cpp <- as.numeric(gsub("[^0-9.]", "", err_cpp))
  
  # R version
  rf <- ranger(Species ~ ., iris, probability = TRUE, num.trees = 5, seed = 10)
  err_r <- rf$prediction.error
  
  expect_equal(round(err_cpp, 4), round(err_r, 4))
})

test_that("Same result with extratrees splitting", {
  # C++ version
  rf <- ranger_cpp(data = iris, depvarname = "Species", probability = "", ntree = 5, splitrule = 5, seed = 10)
  err_cpp <- grep("Overall OOB prediction error", rf, value = TRUE)
  err_cpp <- as.numeric(gsub("[^0-9.]", "", err_cpp))
  
  # R version
  rf <- ranger(Species ~ ., iris, probability = TRUE, num.trees = 5, splitrule = "extratrees", seed = 10)
  err_r <- rf$prediction.error
  
  expect_equal(round(err_cpp, 4), round(err_r, 4))
})
