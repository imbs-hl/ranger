## Tests for hierarchical shrinkage

library(ranger)
context("ranger_hshrink")

## Tests
test_that("horizontal shrinkage gives an error when node.stats=FALSE", {
  rf <- ranger(Sepal.Length ~ ., iris, num.trees = 1, node.stats = FALSE)
  expect_error(hshrink(rf, lambda = 5))
})

test_that("horizontal shrinkage does not work for hard classification", {
  rf <- ranger(Species ~ ., iris, num.trees = 1, node.stats = TRUE, probability = FALSE)
  expect_error(hshrink(rf, lambda = 5))
})

test_that("horizontal shrinkage with lambda=0 doesn't change leafs and prediction, regression", {
  rf <- ranger(Sepal.Length ~ ., iris, num.trees = 1, node.stats = TRUE)
  split_values_before <- rf$forest$split.values[[1]]
  pred_before <- predict(rf, iris)$predictions
  hshrink(rf, lambda = 0)
  split_values_after <- rf$forest$split.values[[1]]
  pred_after <- predict(rf, iris)$predictions
  expect_equal(split_values_before, split_values_after)
  expect_equal(pred_before, pred_after)
})

test_that("horizontal shrinkage with lambda=0 doesn't change leafs and prediction, probability", {
  rf <- ranger(Species ~ ., iris, num.trees = 1, node.stats = TRUE, probability = TRUE)
  class_freq_before <- simplify2array(rf$forest$terminal.class.counts[[1]])
  pred_before <- predict(rf, iris)$predictions
  hshrink(rf, lambda = 0)
  class_freq_after <- simplify2array(rf$forest$terminal.class.counts[[1]])
  pred_after <- predict(rf, iris)$predictions
  expect_equal(class_freq_before, class_freq_after)
  expect_equal(pred_before, pred_after)
})

test_that("horizontal shrinkage with lambda>0 does change leafs and prediction, regression", {
  rf <- ranger(Sepal.Length ~ ., iris, num.trees = 1, replace = FALSE, sample.fraction = 1, node.stats = TRUE)
  split_values_before <- rf$forest$split.values[[1]]
  pred_before <- predict(rf, iris)$predictions
  split_values_before[1] <- 0 # Modify to create deep copy
  hshrink(rf, lambda = 100)
  split_values_after <- rf$forest$split.values[[1]]
  split_values_after[1] <- 0 # Also modify here
  pred_after <- predict(rf, iris)$predictions
  expect_false(all(split_values_before == split_values_after))
  
  # Shrinkage reduces variance
  expect_lt(var(pred_after), var(pred_before))
  
})

test_that("horizontal shrinkage with lambda>0 does change leafs and prediction, probability", {
  rf <- ranger(Species ~ ., iris, num.trees = 1, node.stats = TRUE, probability = TRUE)
  class_freq_before <- simplify2array(rf$forest$terminal.class.counts[[1]])
  pred_before <- predict(rf, iris)$predictions
  hshrink(rf, lambda = 100)
  class_freq_after <- simplify2array(rf$forest$terminal.class.counts[[1]])
  pred_after <- predict(rf, iris)$predictions
  expect_false(all(class_freq_before == class_freq_after))
  
  # Shrinkage reduces variance
  expect_lt(var(pred_after[, 1]), var(pred_before[, 1]))
  expect_lt(var(pred_after[, 2]), var(pred_before[, 2]))
  expect_lt(var(pred_after[, 3]), var(pred_before[, 3]))
})
