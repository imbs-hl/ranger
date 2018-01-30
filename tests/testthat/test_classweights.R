## Tests for class weights

library(ranger)
context("ranger_classweights")

test_that("No error if class weights used", {
  expect_silent(ranger(Species ~ ., iris, num.trees = 5, class.weights = c(0.5, 1, 0.1)))
})

test_that("Prediction error worse if class weights 0", {
  rf <- ranger(Species ~ ., iris, num.trees = 5)
  rf_weighted <- ranger(Species ~ ., iris, num.trees = 5, class.weights = c(1, 0, 0))
  expect_lt(rf$prediction.error, rf_weighted$prediction.error)
})

test_that("Error if class weights of wrong size", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, class.weights = c(0.5, 1)), 
               "Error: Number of class weights not equal to number of classes.")
})

test_that("Error if class weights NA", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, class.weights = c(0.5, 1, NA)), 
               "missing value where TRUE/FALSE needed")
})

test_that("Error if class weights not numeric", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, class.weights = c(0.5, 1, "a")), 
               "Error: Invalid value for class.weights. Please give a vector of non-negative values.")
})


