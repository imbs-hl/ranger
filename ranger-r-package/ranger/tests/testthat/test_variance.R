## Tests for variance estimation

library(ranger)
context("infinitesimal jackknife")

## Set dataset
train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
iris.train <- iris[train.idx, ]
iris.test <- iris[-train.idx, ]

## Train random forest
rf  <- ranger(Species ~ ., data = iris.train, num.trees = 5, write.forest = TRUE, 
             keep.inbag = TRUE)
rf2 <- ranger(Species ~ ., data = iris.train, num.trees = 5, write.forest = TRUE, 
             keep.inbag = FALSE)

test_that("var_IJ_U returns variance of prediction", {
  pred <- predict(rf, iris.test, predict.all = TRUE, var.IJ = TRUE)
  expect_is(pred$var.ij.u, "numeric")
  expect_equal(length(pred$var.ij.u), nrow(iris.test))
})

test_that("Error if keep.inbag equals false", {
  expect_error(predict(rf2, iris.test, predict.all = TRUE, var.IJ = TRUE))
})

