## Tests for predictions

library(ranger)
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

test_that("Terminal nodes returned by predict are node ids", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred <- predict(rf, iris, type = "terminalNodeIds")
  
  expect_equal(dim(pred), c(nrow(iris), rf$num.trees))
})

test_that("Terminal nodes returned by predict are the same as by getTerminalNodeIds", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred <- predict(rf, iris, type = "terminalNodeIds")
  
  nodeIds <- getTerminalNodeIDs(rf, iris)
  
  expect_equal(pred$predictions, nodeIds)
})
