
context("ranger_cpp_regression")

test_that("Prediction is equal to R version", {
  # C++ version
  rf <- ranger_cpp(data = veteran, depvarname = "time", statusvarname = "status", treetype = 5, ntree = 5, write = "", seed = 10)
  pred <- ranger_cpp(data = veteran, treetype = 5, predict = "ranger_out.forest", seed = 20)
  preds_cpp <- as.matrix(fread("ranger_out.prediction", skip = 4))
  dimnames(preds_cpp) <- NULL
  
  # R version
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, seed = 10)
  preds_r <- predict(rf, veteran, seed = 20)$chf
  
  expect_equal(round(preds_cpp, 2), round(preds_r, 2))
})

test_that("Predictions are positive numerics", {
  rf <- ranger_cpp(data = veteran, depvarname = "time", statusvarname = "status", treetype = 5, ntree = 5, write = "", seed = 10)
  pred <- ranger_cpp(data = veteran, predict = "ranger_out.forest", treetype = 5)
  preds_cpp <- as.matrix(fread("ranger_out.prediction", skip = 4))
  expect_is(preds_cpp, "matrix")
  expect_true(all(preds_cpp >= 0))
})

test_that("Same result with default splitting", {
  # C++ version
  rf <- ranger_cpp(data = veteran, depvarname = "time", statusvarname = "status", treetype = 5, ntree = 5, seed = 10)
  err_cpp <- grep("Overall OOB prediction error", rf, value = TRUE)
  err_cpp <- as.numeric(gsub("[^0-9.]", "", err_cpp))
  
  # R version
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, seed = 10)
  err_r <- rf$prediction.error
  
  expect_equal(round(err_cpp, 4), round(err_r, 4))
})

test_that("Same result with extratrees splitting", {
  # C++ version
  rf <- ranger_cpp(data = veteran, depvarname = "time", statusvarname = "status", treetype = 5, ntree = 5, splitrule = 5,  catvars = "celltype", seed = 10)
  err_cpp <- grep("Overall OOB prediction error", rf, value = TRUE)
  err_cpp <- as.numeric(gsub("[^0-9.]", "", err_cpp))
  
  # R version
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, splitrule = "extratrees", seed = 10)
  err_r <- rf$prediction.error
  
  expect_equal(round(err_cpp, 4), round(err_r, 4))
})

test_that("Same result with maxstat splitting", {
  # C++ version
  rf <- ranger_cpp(data = veteran, depvarname = "time", statusvarname = "status", treetype = 5, ntree = 5, splitrule = 4, seed = 10)
  err_cpp <- grep("Overall OOB prediction error", rf, value = TRUE)
  err_cpp <- as.numeric(gsub("[^0-9.]", "", err_cpp))
  
  # R version
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, splitrule = "maxstat", seed = 10)
  err_r <- rf$prediction.error
  
  expect_equal(round(err_cpp, 4), round(err_r, 4))
})

test_that("Same result with unordered splitting", {
  # C++ version
  rf <- ranger_cpp(data = veteran, depvarname = "time", statusvarname = "status", treetype = 5, ntree = 5, catvars = "celltype", seed = 10)
  err_cpp <- grep("Overall OOB prediction error", rf, value = TRUE)
  err_cpp <- as.numeric(gsub("[^0-9.]", "", err_cpp))
  
  # R version
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, respect.unordered.factors = "partition", seed = 10)
  err_r <- rf$prediction.error
  
  expect_equal(round(err_cpp, 4), round(err_r, 4))
})
