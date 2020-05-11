library(ranger)
library(survival)

test_that("LHS formula", {
  callRanger <- function() {
    myTransformation <- function(x) { x }
    ranger(myTransformation(Species) ~ ., data = iris)
  }
  
  expect_error(callRanger(), NA)
})

test_that("Formula works with matrix data", {
  set.seed(10)
  rf1 <- ranger(Species ~ ., iris, num.trees = 5)
  pred1 <- rf1$predictions
  
  set.seed(10)
  rf2 <- ranger(Species ~ ., data.matrix(iris), num.trees = 5, classification = TRUE)
  pred2 <- factor(levels(iris$Species)[rf2$predictions], levels = levels(iris$Species))
  
  expect_equal(pred1, pred2)
})

test_that("Formula works with matrix data, survival", {
  set.seed(10)
  rf1 <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5)
  
  set.seed(10)
  rf2 <- ranger(Surv(time, status) ~ ., data.matrix(veteran), num.trees = 5)
  
  expect_equal(rf1$chf, rf2$chf)
})

