## Tests for importance measures

library(ranger)
library(survival)
context("ranger_imp")

## Initialize the random forests
rg.imp <- ranger(Species ~ ., data = iris, num.trees = 10,
                 importance = "impurity")
rg.unbiased <- ranger(Species ~ ., data = iris, num.trees = 10,
                      importance = "impurity_unbiased")
rg.perm <- ranger(Species ~ ., data = iris, num.trees = 10,
                 importance = "permutation")
rg.scale.perm <- ranger(Species ~ ., data = iris, num.trees = 50,
                 importance = "permutation", scale.permutation.importance = TRUE)

## Tests
test_that("importance measures work", {
  expect_is(rg.imp$variable.importance, "numeric")
  expect_is(rg.perm$variable.importance, "numeric")
  expect_is(rg.scale.perm$variable.importance, "numeric")
  expect_is(rg.unbiased$variable.importance, "numeric")
})

test_that("gini importance is larger than 1", {
  expect_gt(rg.imp$variable.importance[1], 1)
})

test_that("Unbiased gini importance is larger than 1", {
  expect_gt(rg.unbiased$variable.importance[1], 1)
})

test_that("unscaled importance is smaller than 1", {
  expect_lt(rg.perm$variable.importance[1], 1)
})

test_that("scaled importance is larger than 1", {
  expect_gt(rg.scale.perm$variable.importance[1], 1)
})

test_that("error thrown if no importance in object", {
  rf <- ranger(Species ~ ., data = iris, num.trees = 5)
  expect_error(importance(rf), "No variable importance found. Please use 'importance' option when growing the forest.")
})

test_that("Error thrown if Unbiased gini importance used with split.select.weights", {
  expect_error(ranger(Species ~ ., data = iris, num.trees = 5, 
                      split.select.weights = rep(.5, 4), importance = "impurity_unbiased"), 
               "Unbiased impurity importance not supported in combination with split.select.weights.")
})

test_that("Survival permutation importance is smaller than 1", {
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, importance = "permutation")
  expect_lt(rf$variable.importance[1], 1)
})

test_that("Survival impurity importance is larger than 1", {
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, importance = "impurity")
  expect_gt(rf$variable.importance[1], 1)
})

test_that("Survival corrected impurity importance is smaller than 1", {
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, importance = "impurity_unbiased")
  expect_lt(rf$variable.importance[1], 1)
})

