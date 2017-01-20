library(ranger)
library(survival)
context("ranger_interface")

## Formula interface
test_that("All variables included if . in formula", {
  rf <- ranger(Species ~ ., iris, num.trees = 5)
  expect_equal(sort(rf$forest$independent.variable.names), 
               sort(colnames(iris)[1:4]))
})

test_that("Variable excluded if - in formula", {
  rf <- ranger(Species ~ . -Petal.Length, iris, num.trees = 5)
  expect_equal(sort(rf$forest$independent.variable.names), 
               sort(c("Sepal.Length", "Sepal.Width", "Petal.Width")))
})

test_that("Interaction included if : in formula", {
  rf <- ranger(Species ~ Petal.Length + Sepal.Length:Sepal.Width, iris, num.trees = 5)
  expect_equal(sort(rf$forest$independent.variable.names), 
               sort(c("Petal.Length", "Sepal.Length:Sepal.Width")))
})

test_that("Interaction included if * in formula", {
  rf <- ranger(Species ~ Petal.Length + Sepal.Length*Sepal.Width, iris, num.trees = 5)
  expect_equal(sort(rf$forest$independent.variable.names), 
               sort(c("Petal.Length", "Sepal.Length", "Sepal.Width", "Sepal.Length:Sepal.Width")))
})
