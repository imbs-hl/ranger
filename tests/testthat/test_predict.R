## Tests for predictions

library(ranger)
library(survival)
context("ranger_pred")

## Tests
test_that("predict returns good prediction", {
  rf <- ranger(Species ~ ., iris, write.forest = TRUE)
  pred <- predict(rf, iris)
  expect_gt(mean(iris$Species == predictions(pred)), 0.9)
})

test_that("case weights work", {
  expect_silent(ranger(Species ~ ., iris, num.trees = 5, case.weights = rep(1, nrow(iris))))
  ## Should only predict setosa now
  weights <- c(rep(1, 50), rep(0, 100))
  rf <- ranger(Species ~ ., iris, num.trees = 5, case.weights = weights, write.forest = TRUE)
  pred <- predict(rf, iris)$predictions
  expect_true(all(pred == "setosa"))
})

test_that("Prediction works correctly if dependent variable is not first or last", {
  dat <- iris[, c(1:2, 5, 3:4)]
  rf <- ranger(Species ~ ., dat, num.trees = 5, write.forest = TRUE)
  expect_gte(mean(predictions(predict(rf, dat)) == dat$Species), 0.9)
  
  ## No response column
  expect_gte(mean(predictions(predict(rf, dat[, -3])) == dat$Species), 0.9)
})

test_that("Prediction works correctly if dependent variable is not first or last, alternative interface", {
  dat <- iris[, c(1:2, 5, 3:4)]
  rf <- ranger(dependent.variable.name = "Species", data = dat, num.trees = 5, write.forest = TRUE)
  expect_gte(mean(predictions(predict(rf, dat)) == dat$Species), 0.9)
  
  ## No response column
  expect_gte(mean(predictions(predict(rf, dat[, -3])) == dat$Species), 0.9)
})

test_that("Missing value columns detected in predict", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  
  dat <- iris
  dat[4, 4] <- NA
  dat[25, 1] <- NA
  expect_error(predict(rf, dat), "Missing data in columns: Sepal.Length, Petal.Width.")
})

test_that("If num.trees set, these number is used for predictions", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred <- predict(rf, iris, predict.all = TRUE, num.trees = 3)
  expect_equal(pred$num.trees, 3)
  expect_equal(dim(pred$predictions), c(nrow(iris), 3))
})

test_that("If num.trees not set, all trees are used for prediction", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred <- predict(rf, iris, predict.all = TRUE)
  expect_equal(pred$num.trees, 5)
  expect_equal(dim(pred$predictions), c(nrow(iris), 5))
})

test_that("Error if unknown value for type", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  expect_error(predict(rf, iris, type = "class"))
})

test_that("Terminal nodes returned by predict are node ids, classification", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred <- predict(rf, iris, type = "terminalNodes")
  
  expect_equal(dim(pred$predictions), c(nrow(iris), rf$num.trees))
  expect_true(all(pred$predictions > 0))
  expect_true(all(pred$predictions < max(sapply(rf$forest$split.varIDs, length))))
})


test_that("Terminal nodes returned by predict are node ids, probability", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE, probability = TRUE)
  pred <- predict(rf, iris, type = "terminalNodes")
  
  expect_equal(dim(pred$predictions), c(nrow(iris), rf$num.trees))
  expect_true(all(pred$predictions > 0))
  expect_true(all(pred$predictions < max(sapply(rf$forest$split.varIDs, length))))
})

test_that("Terminal nodes returned by predict are node ids, regression", {
  rf <- ranger(Sepal.Length ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred <- predict(rf, iris, type = "terminalNodes")
  
  expect_equal(dim(pred$predictions), c(nrow(iris), rf$num.trees))
  expect_true(all(pred$predictions > 0))
  expect_true(all(pred$predictions < max(sapply(rf$forest$split.varIDs, length))))
})

test_that("Terminal nodes returned by predict are node ids, survival", {
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, write.forest = TRUE)
  pred <- predict(rf, veteran, type = "terminalNodes")
  
  expect_equal(dim(pred$predictions), c(nrow(veteran), rf$num.trees))
  expect_true(all(pred$predictions > 0))
  expect_true(all(pred$predictions < max(sapply(rf$forest$split.varIDs, length))))
})

test_that("Same result with warning if getTerminalNodeIDs() used", {
  rf <- ranger(Species ~ ., iris, num.trees = 5)
  pred <- predict(rf, iris, type = "terminalNodes")
  expect_warning(expect_equal(getTerminalNodeIDs(rf, iris), pred$predictions))
})

test_that("predict.all works for single observation", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred <- predict(rf, iris[1, ], predict.all = TRUE)
  
  expect_equal(dim(pred$predictions), c(1, rf$num.trees))
})

test_that("Warning if predicting with corrected impurity importance", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, importance = "impurity_corrected")
  expect_warning(predict(rf, iris))
})

test_that("standard error prediction working for regression", {
  idx <- sample(nrow(iris), 10)
  test <- iris[idx, ]
  train <- iris[-idx, ]
  
  rf <- ranger(Petal.Length ~ ., train, num.trees = 5, keep.inbag = TRUE)
  pred <- predict(rf, test, type = "se")
  
  expect_equal(length(pred$predictions), nrow(test))
})

test_that("standard error prediction not working for other tree types", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE)
  expect_error(predict(rf, iris, type = "se"), 
               "Error: Standard error prediction currently only available for regression.")
  
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, probability = TRUE)
  expect_error(predict(rf, iris, type = "se"), 
               "Error: Standard error prediction currently only available for regression.")
  
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, keep.inbag = TRUE)
  expect_error(predict(rf, veteran, type = "se"), 
               "Error: Standard error prediction currently only available for regression.")
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
  pred <- predict(rf, test, type = "se")
  
  expect_equal(length(pred$predictions), nrow(test))
})

test_that("standard error response prediction is the same as response prediction", {
  idx <- sample(nrow(iris), 10)
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

test_that("standard error is larger for fewer trees", {
  idx <- sample(nrow(iris), 10)
  test <- iris[idx, ]
  train <- iris[-idx, ]
  
  rf5 <- ranger(Petal.Length ~ ., train, num.trees = 5, keep.inbag = TRUE)
  pred5 <- predict(rf5, test, type = "se")
  
  rf50 <- ranger(Petal.Length ~ ., train, num.trees = 50, keep.inbag = TRUE)
  pred50 <- predict(rf50, test, type = "se")
  
  expect_lt(mean(pred50$se), mean(pred5$se))
})
