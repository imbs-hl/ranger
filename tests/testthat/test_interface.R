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

# Tibbles
if (requireNamespace("tibble", quietly = TRUE)) {
  attachNamespace("tibble")
  tb <- as_tibble(iris)
}
test_that("Training works with tibbles, formula interface", {
  skip_if_not_installed("tibble")
  set.seed(1000)
  rf1 <- ranger(Species ~ ., tb, num.trees = 5)
  
  set.seed(1000)
  rf2 <- ranger(Species ~ ., iris, num.trees = 5)
  
  expect_equal(rf1$prediction.error, rf2$prediction.error)
  
  pred1 <- levels(iris$Species)[rf1$predictions[!is.na(rf1$predictions)]]
  pred2 <- as.character(rf2$predictions[!is.na(rf2$predictions)])
  expect_equal(pred1, pred2)
})

test_that("Training works with tibbles, alternative interface", {
  skip_if_not_installed("tibble")
  set.seed(1000)
  rf1 <- ranger(dependent.variable.name = "Species", data = tb, num.trees = 5)
  
  set.seed(1000)
  rf2 <- ranger(dependent.variable.name = "Species", data = iris, num.trees = 5)
  
  expect_equal(rf1$prediction.error, rf2$prediction.error)
  
  pred1 <- levels(iris$Species)[rf1$predictions[!is.na(rf1$predictions)]]
  pred2 <- as.character(rf2$predictions[!is.na(rf2$predictions)])
  expect_equal(pred1, pred2)
})

test_that("Prediction works with tibbles, formula interface", {
  skip_if_not_installed("tibble")
  set.seed(1000)
  rf1 <- ranger(Species ~ ., tb, num.trees = 5)
  
  set.seed(1000)
  rf2 <- ranger(Species ~ ., iris, num.trees = 5)
  
  set.seed(1000)
  pred1 <- predict(rf1, tb)
  set.seed(1000)
  pred2 <- predict(rf1, iris)
  set.seed(1000)
  pred3 <- predict(rf2, tb)
  set.seed(1000)
  pred4 <- predict(rf2, iris)
  
  expect_equal(pred1$predictions, pred2$predictions)
  expect_equal(pred2$predictions, pred3$predictions)
  expect_equal(pred3$predictions, pred4$predictions)
})

test_that("Prediction works with tibbles, alternative interface", {
  skip_if_not_installed("tibble")
  set.seed(1000)
  rf1 <- ranger(dependent.variable.name = "Species", data = tb, num.trees = 5)
  
  set.seed(1000)
  rf2 <- ranger(dependent.variable.name = "Species", data = iris, num.trees = 5)
  
  set.seed(1000)
  pred1 <- predict(rf1, tb)
  set.seed(1000)
  pred2 <- predict(rf1, iris)
  set.seed(1000)
  pred3 <- predict(rf2, tb)
  set.seed(1000)
  pred4 <- predict(rf2, iris)
  
  expect_equal(pred1$predictions, pred2$predictions)
  expect_equal(pred2$predictions, pred3$predictions)
  expect_equal(pred3$predictions, pred4$predictions)
})
