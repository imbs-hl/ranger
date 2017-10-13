## Tests for the (infitesimal) jackknife for standard error prediction

library(ranger)
library(survival)
context("ranger_jackknife")

test_that("jackknife standard error prediction working for regression", {
  idx <- sample(nrow(iris), 10)
  test <- iris[idx, ]
  train <- iris[-idx, ]
  
  rf <- ranger(Petal.Length ~ ., train, num.trees = 5, keep.inbag = TRUE)
  pred <- predict(rf, test, type = "se", se.method = "jack")
  
  expect_equal(length(pred$predictions), nrow(test))
})

test_that("jackknife standard error prediction not working for other tree types", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE)
  expect_error(predict(rf, iris, type = "se", se.method = "jack"), 
               "Error: Jackknife standard error prediction currently only available for regression.")
  
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, probability = TRUE)
  expect_error(predict(rf, iris, type = "se", se.method = "jack"), 
               "Error: Jackknife standard error prediction currently only available for regression.")
  
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, keep.inbag = TRUE)
  expect_error(predict(rf, veteran, type = "se", se.method = "jack"), 
               "Error: Jackknife standard error prediction currently only available for regression.")
})

test_that("IJ standard error prediction working for regression", {
  idx <- sample(nrow(iris), 10)
  test <- iris[idx, ]
  train <- iris[-idx, ]
  
  rf <- ranger(Petal.Length ~ ., train, num.trees = 5, keep.inbag = TRUE)
  pred <- predict(rf, test, type = "se", se.method = "infjack")
  
  expect_equal(length(pred$predictions), nrow(test))
})

test_that("IJ standard error prediction working for probability", {
  idx <- sample(nrow(iris), 25)
  test <- iris[idx, ]
  train <- iris[-idx, ]
  
  rf <- ranger(Species ~ ., train, num.trees = 5, keep.inbag = TRUE, probability = TRUE)
  pred <- predict(rf, test, type = "se", se.method = "infjack")
  
  expect_equal(nrow(pred$predictions), nrow(test))
})

test_that("IJ standard error prediction not working for other tree types", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE)
  expect_error(predict(rf, iris, type = "se", se.method = "infjack"), 
               "Error: Not a probability forest. Set probability=TRUE to use the infinitesimal jackknife standard error prediction for classification.")
  
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, keep.inbag = TRUE)
  expect_error(predict(rf, veteran, type = "se", se.method = "infjack"), 
               "Error: Infinitesimal jackknife standard error prediction not yet available for survival.")
})

test_that("standard error prediction not working if keep.inbag = FALSE", {
  rf <- ranger(Petal.Length ~ ., iris, num.trees = 5)
  expect_error(predict(rf, iris, type = "se"), 
               "Error: No saved inbag counts in ranger object. Please set keep.inbag=TRUE when calling ranger.")
})

test_that("standard error prediction not working if no OOB observations", {
  test <- iris[-1, ]
  train <- iris[1, ]
  rf <- ranger(Petal.Length ~ ., train, num.trees = 5, keep.inbag = TRUE)
  expect_error(predict(rf, iris, type = "se"), 
               "Error: No OOB observations found, consider increasing num.trees or reducing sample.fraction.")
})

test_that("standard error prediction working for single testing observation", {
  test <- iris[1, ]
  train <- iris[-1, ]
  
  rf <- ranger(Petal.Length ~ ., train, num.trees = 5, keep.inbag = TRUE)
  
  pred_jack <- predict(rf, test, type = "se", se.method = "jack")
  expect_equal(length(pred_jack$predictions), nrow(test))
  
  pred_ij <- predict(rf, test, type = "se", se.method = "infjack")
  expect_equal(length(pred_ij$predictions), nrow(test))
})

test_that("standard error response prediction is the same as response prediction", {
  idx <- sample(nrow(iris), 25)
  test <- iris[idx, ]
  train <- iris[-idx, ]
  
  set.seed(100)
  rf_se <- ranger(Petal.Length ~ ., train, num.trees = 5, keep.inbag = TRUE)
  pred_se <- predict(rf_se, test, type = "se")
  
  set.seed(100)
  rf_resp <- ranger(Petal.Length ~ ., train, num.trees = 5)
  pred_resp <- predict(rf_resp, test, type = "response")
  
  expect_equal(pred_se$predictions, pred_resp$predictions)
})

test_that("standard error is larger for fewer trees, regression", {
  idx <- sample(nrow(iris), 25)
  test <- iris[idx, ]
  train <- iris[-idx, ]
  
  rf5 <- ranger(Petal.Length ~ ., train, num.trees = 5, keep.inbag = TRUE)
  pred5_jack <- predict(rf5, test, type = "se", se.method = "jack")
  pred5_ij <- predict(rf5, test, type = "se", se.method = "infjack")
  
  rf50 <- ranger(Petal.Length ~ ., train, num.trees = 50, keep.inbag = TRUE)
  pred50_jack <- predict(rf50, test, type = "se", se.method = "jack")
  pred50_ij <- predict(rf50, test, type = "se", se.method = "infjack")
  
  expect_lt(mean(pred50_jack$se), mean(pred5_jack$se))
  expect_lt(mean(pred50_ij$se), mean(pred5_ij$se))
})

test_that("standard error is larger for fewer trees, probability", {
  idx <- sample(nrow(iris), 25)
  test <- iris[idx, ]
  train <- iris[-idx, ]
  
  rf5 <- ranger(Species ~ ., train, num.trees = 5, keep.inbag = TRUE, probability = TRUE)
  pred5 <- predict(rf5, test, type = "se", se.method = "infjack")
  
  rf50 <- ranger(Petal.Length ~ ., train, num.trees = 50, keep.inbag = TRUE)
  pred50 <- predict(rf50, test, type = "se", se.method = "infjack")
  
  expect_lt(mean(pred50$se), mean(pred5$se))
})

test_that("Warning for few observations with IJ", {
  idx <- sample(nrow(iris), 10)
  test <- iris[idx, ]
  train <- iris[-idx, ]
  
  rf5 <- ranger(Species ~ ., train, num.trees = 5, keep.inbag = TRUE, probability = TRUE)
  expect_warning(predict(rf5, test, type = "se", se.method = "infjack"), 
                 "Sample size <=20, no calibration performed.")
})
