## Tests for importance measures

library(ranger)
context("ranger_imp")

## Classification
rg.imp.class <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                       num.trees = 5, importance = "impurity")
rg.perm.class <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                        num.trees = 5, importance = "permutation")
rg.scale.perm.class <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                              num.trees = 5, importance = "permutation", scale.permutation.importance = TRUE)

## Probability estimation
rg.imp.prob <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                      num.trees = 5, importance = "impurity", probability = TRUE)
rg.perm.prob <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                      num.trees = 5, importance = "permutation", probability = TRUE)
rg.scale.perm.prob <- ranger(Species ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                             num.trees = 5, importance = "permutation", scale.permutation.importance = TRUE, probability = TRUE)

## Regression
rg.imp.regr <- ranger(Sepal.Length ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                      num.trees = 5, importance = "impurity")
rg.perm.regr <- ranger(Sepal.Length ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                       num.trees = 5, importance = "permutation")
rg.scale.perm.regr <- ranger(Sepal.Length ~ ., data = iris, verbose = FALSE, write.forest = TRUE,
                             num.trees = 5, importance = "permutation", scale.permutation.importance = TRUE)

## Survival
rg.perm.surv <- ranger(Surv(time, status) ~ ., data = veteran, verbose = FALSE, write.forest = TRUE,
                       num.trees = 5, importance = "permutation")
rg.scale.perm.surv <- ranger(Surv(time, status) ~ ., data = veteran, verbose = FALSE, write.forest = TRUE,
                             num.trees = 5, importance = "permutation", scale.permutation.importance = TRUE)


## Tests
test_that("importance measures work, classification", {
  expect_is(rg.imp.class$variable.importance, "numeric")
  expect_is(rg.perm.class$variable.importance, "numeric")
  expect_is(rg.scale.perm.class$variable.importance, "numeric")
})

test_that("importance measures work, probability", {
  expect_is(rg.imp.prob$variable.importance, "numeric")
  expect_is(rg.perm.prob$variable.importance, "numeric")
  expect_is(rg.scale.perm.prob$variable.importance, "numeric")
})

test_that("importance measures work, regression", {
  expect_is(rg.imp.regr$variable.importance, "numeric")
  expect_is(rg.perm.regr$variable.importance, "numeric")
  expect_is(rg.scale.perm.regr$variable.importance, "numeric")
})

test_that("importance measures work, survival", {
  expect_is(rg.perm.surv$variable.importance, "numeric")
  expect_is(rg.scale.perm.surv$variable.importance, "numeric")
})

test_that("impurity importance is larger than 1", {
  expect_gt(rg.imp.class$variable.importance[1], 1)
  expect_gt(rg.imp.prob$variable.importance[1], 1)
  expect_gt(rg.imp.regr$variable.importance[1], 1)
})

test_that("unscaled importance is smaller than 1", {
  expect_lt(rg.perm.class$variable.importance[1], 1)
  expect_lt(rg.perm.prob$variable.importance[1], 1)
  expect_lt(rg.perm.regr$variable.importance[1], 1)
  expect_lt(rg.perm.surv$variable.importance[3], 1)
})

test_that("scaled importance is larger than 1", {
  expect_gt(rg.scale.perm.class$variable.importance[1], 1)
  expect_gt(rg.scale.perm.prob$variable.importance[1], 1)
  expect_gt(rg.scale.perm.regr$variable.importance[1], 1)
  expect_gt(rg.scale.perm.surv$variable.importance[3], 1)
})

test_that("error thrown if no importance in object", {
  rf <- ranger(Species ~ ., data = iris, num.trees = 5)
  expect_error(importance(rf), "No variable importance found. Please use 'importance' option when growing the forest.")
})

