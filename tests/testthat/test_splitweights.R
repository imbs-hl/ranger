## Tests for split select weights

library(ranger)
context("ranger_splitweights")

## Tests
test_that("split select weights work", {
  expect_silent(ranger(Species ~ ., iris, num.trees = 5, split.select.weights = c(0.1, 0.2, 0.3, 0.4)))
  expect_error(ranger(Species ~ ., iris, num.trees = 5, split.select.weights = c(0.1, 0.2, 0.3)))
})

test_that("split select weights work with 0s and 1s", {
  num.trees <- 5
  weights <- replicate(num.trees, sample(c(0, 0, 1, 1)), simplify = FALSE)
  rf <- ranger(Species ~ ., iris, num.trees = num.trees, split.select.weights = weights)
  selected_correctly <- sapply(1:rf$num.trees, function(i) {
    all(treeInfo(rf, i)[,"splitvarID"] %in% c(which(weights[[i]] > 0) - 1, NA))
  })
  expect_true(all(selected_correctly))
})

test_that("Tree-wise split select weights work", {
  num.trees <- 5
  weights <- replicate(num.trees, runif(ncol(iris)-1), simplify = FALSE)
  expect_silent(ranger(Species ~ ., iris, num.trees = num.trees, split.select.weights = weights))
  
  weights <- replicate(num.trees+1, runif(ncol(iris)-1), simplify = FALSE)
  expect_error(ranger(Species ~ ., iris, num.trees = num.trees, split.select.weights = weights))
})

test_that("always split variables work", {
  expect_silent(ranger(Species ~ ., iris, num.trees = 10, 
                       always.split.variables = c("Petal.Length", "Petal.Width"), mtry = 2))
  expect_silent(ranger(Species ~ ., iris, num.trees = 10, 
                       always.split.variables = c("Petal.Width", "Petal.Length"), mtry = 2))
  expect_silent(ranger(dependent.variable.name = "Species", data = iris, num.trees = 10, 
                       always.split.variables = c("Petal.Length", "Petal.Width"), mtry = 2))
})

test_that("Tree-wise split select weights work with 0s", {
  num.trees <- 5
  weights <- replicate(num.trees, sample(c(0, 0, 0.5, 0.5)), simplify = FALSE)
  rf <- ranger(Species ~ ., iris, mtry = 2, num.trees = num.trees, 
               split.select.weights = weights)
  selected_correctly <- sapply(1:num.trees, function(i) {
    all(treeInfo(rf, i)[,"splitvarID"] %in% c(which(weights[[i]] > 0) - 1, NA))
  })
  expect_true(all(selected_correctly))
})

test_that("always split variables respect split select weights", {
    iris_vars <- setdiff(names(iris), 'Species')
    n_vars <- length(iris_vars)
    last_var <- iris_vars[n_vars]
    with_last_zero <- c(rep(1, n_vars-1), 0)
    expect_silent(
        ranger(Species ~ ., iris, num.trees=5,
               always.split.variables=last_var, mtry=n_vars-1,
               split.select.weights=with_last_zero)
    )
})

