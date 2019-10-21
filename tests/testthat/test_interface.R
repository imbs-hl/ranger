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

test_that("Working if dependent variable has attributes other than names", {
  iris2 <- iris
  attr(iris2$Sepal.Width, "aaa") <- "bbb"
  expect_silent(ranger(data = iris2, dependent.variable = "Sepal.Width"))
})

test_that("Working if dependent variable is matrix with one column", {
  iris2 <- iris
  iris2$Sepal.Width = scale(iris$Sepal.Width)
  expect_silent(ranger(data = iris2, dependent.variable = "Sepal.Width"))
})

test_that("Same result with x/y interface, classification", {
  set.seed(300)
  rf_formula <- ranger(Species ~ ., iris, num.trees = 5)
  
  set.seed(300)
  rf_xy <- ranger(y = iris[, 5], x = iris[, -5], num.trees = 5)
  
  expect_equal(rf_formula$prediction.error, rf_xy$prediction.error)
  expect_equal(rf_formula$predictions, rf_xy$predictions)
})

test_that("Same result with x/y interface, regression", {
  set.seed(300)
  rf_formula <- ranger(Sepal.Length ~ ., iris, num.trees = 5)
  
  set.seed(300)
  rf_xy <- ranger(y = iris[, 1], x = iris[, -1], num.trees = 5)
  
  expect_equal(rf_formula$prediction.error, rf_xy$prediction.error)
  expect_equal(rf_formula$predictions, rf_xy$predictions)
})

test_that("Same result with x/y interface, survival", {
  set.seed(300)
  rf_formula <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5)
  
  set.seed(300)
  rf_xy <- ranger(y = veteran[, c(3, 4)], x = veteran[, c(-3, -4)], num.trees = 5)
  
  expect_equal(rf_formula$prediction.error, rf_xy$prediction.error)
  expect_equal(rf_formula$predictions, rf_xy$predictions)
})

test_that("Column order does not change prediction", {
  dat <- iris[, c(sample(1:4), 5)]
  rf <- ranger(dependent.variable.name = "Species", data = iris)
  
  set.seed(42)
  pred1 <- predict(rf, iris)$predictions
  
  set.seed(42)
  pred2 <- predict(rf, dat)$predictions
  
  expect_equal(pred1, pred2)
})

# Tibbles
# This is failing on Rdevel. Possible without suggesting tibble package?
# if (requireNamespace("tibble", quietly = TRUE)) {
#   tb <- tibble::as_tibble(iris)
# }
# test_that("Training works with tibbles, formula interface", {
#   skip_if_not_installed("tibble")
#   set.seed(1000)
#   rf1 <- ranger(Species ~ ., tb, num.trees = 5)
#   
#   set.seed(1000)
#   rf2 <- ranger(Species ~ ., iris, num.trees = 5)
#   
#   expect_equal(rf1$prediction.error, rf2$prediction.error)
#   
#   pred1 <- levels(iris$Species)[rf1$predictions[!is.na(rf1$predictions)]]
#   pred2 <- as.character(rf2$predictions[!is.na(rf2$predictions)])
#   expect_equal(pred1, pred2)
# })
# 
# test_that("Training works with tibbles, alternative interface", {
#   skip_if_not_installed("tibble")
#   set.seed(1000)
#   rf1 <- ranger(dependent.variable.name = "Species", data = tb, num.trees = 5)
#   
#   set.seed(1000)
#   rf2 <- ranger(dependent.variable.name = "Species", data = iris, num.trees = 5)
#   
#   expect_equal(rf1$prediction.error, rf2$prediction.error)
#   
#   pred1 <- levels(iris$Species)[rf1$predictions[!is.na(rf1$predictions)]]
#   pred2 <- as.character(rf2$predictions[!is.na(rf2$predictions)])
#   expect_equal(pred1, pred2)
# })
# 
# test_that("Prediction works with tibbles, formula interface", {
#   skip_if_not_installed("tibble")
#   set.seed(1000)
#   rf1 <- ranger(Species ~ ., tb, num.trees = 5)
#   
#   set.seed(1000)
#   rf2 <- ranger(Species ~ ., iris, num.trees = 5)
#   
#   set.seed(1000)
#   pred1 <- predict(rf1, tb)
#   set.seed(1000)
#   pred2 <- predict(rf1, iris)
#   set.seed(1000)
#   pred3 <- predict(rf2, tb)
#   set.seed(1000)
#   pred4 <- predict(rf2, iris)
#   
#   expect_equal(pred1$predictions, pred2$predictions)
#   expect_equal(pred2$predictions, pred3$predictions)
#   expect_equal(pred3$predictions, pred4$predictions)
# })
# 
# test_that("Prediction works with tibbles, alternative interface", {
#   skip_if_not_installed("tibble")
#   set.seed(1000)
#   rf1 <- ranger(dependent.variable.name = "Species", data = tb, num.trees = 5)
#   
#   set.seed(1000)
#   rf2 <- ranger(dependent.variable.name = "Species", data = iris, num.trees = 5)
#   
#   set.seed(1000)
#   pred1 <- predict(rf1, tb)
#   set.seed(1000)
#   pred2 <- predict(rf1, iris)
#   set.seed(1000)
#   pred3 <- predict(rf2, tb)
#   set.seed(1000)
#   pred4 <- predict(rf2, iris)
#   
#   expect_equal(pred1$predictions, pred2$predictions)
#   expect_equal(pred2$predictions, pred3$predictions)
#   expect_equal(pred3$predictions, pred4$predictions)
# })
