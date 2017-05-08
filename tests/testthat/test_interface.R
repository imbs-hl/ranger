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

## Formula interface, survival
test_that("All variables included if . in formula", {
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5)
  expect_equal(sort(rf$forest$independent.variable.names), 
               sort(colnames(veteran)[c(1:2, 5:8)]))
})

test_that("Variable excluded if - in formula", {
  rf <- ranger(Surv(time, status) ~ . - celltype - age, veteran, num.trees = 5)
  expect_equal(sort(rf$forest$independent.variable.names), 
               sort(c("trt", "karno", "diagtime", "prior")))
})

test_that("Interaction included if : in formula", {
  rf <- ranger(Surv(time, status) ~ celltype + age:prior, veteran, num.trees = 5)
  expect_equal(sort(rf$forest$independent.variable.names), 
               sort(c("celltype", "age:prior")))
})

test_that("Interaction included if * in formula", {
  rf <- ranger(Surv(time, status) ~ celltype + age*prior, veteran, num.trees = 5)
  expect_equal(sort(rf$forest$independent.variable.names), 
               sort(c("celltype", "age", "prior", "age:prior")))
})

test_that("Error if interaction of factor variable included", {
  expect_error(ranger(Surv(time, status) ~ celltype*prior, veteran, num.trees = 5), 
               "Error: Only numeric columns allowed in interaction terms.")
})
