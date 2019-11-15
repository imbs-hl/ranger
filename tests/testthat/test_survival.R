## Tests for random forests for survival analysis

library(ranger)
library(survival)
context("ranger_surv")

## Initialize the random forest for survival analysis
rg.surv <- ranger(Surv(time, status) ~ ., data = veteran, num.trees = 10)

## Basic tests (for all random forests equal)
test_that("survival result is of class ranger with 15 elements", {
  expect_is(rg.surv, "ranger")
  expect_equal(length(rg.surv), 15)
})

test_that("results have right number of trees", {
  expect_equal(rg.surv$num.trees, 10)
})

test_that("results have right number of independent variables", {
  expect_equal(rg.surv$num.independent.variables, ncol(veteran) - 2)
})

test_that("Alternative interface works for survival", {
  rf <- ranger(dependent.variable.name = "time", status.variable.name = "status", data = veteran, num.trees = 10)
  expect_equal(rf$treetype, "Survival")
})

test_that("Alternative interface prediction works for survival", {
  rf <- ranger(dependent.variable.name = "time", status.variable.name = "status", data = veteran, num.trees = 10)
  expect_equal(predict(rf, veteran)$num.independent.variables, ncol(veteran) - 2) 
  expect_equal(predict(rf, veteran[, setdiff(names(veteran), c("time", "status"))])$num.independent.variables, ncol(veteran) - 2) 
})

test_that("Matrix interface works for survival", {
  rf <- ranger(dependent.variable.name = "time", status.variable.name = "status", data = data.matrix(veteran), write.forest = TRUE, num.trees = 10)
  expect_equal(rf$treetype, "Survival")
  expect_equal(rf$forest$independent.variable.names, colnames(veteran)[c(1:2, 5:8)])
})

test_that("Matrix interface prediction works for survival", {
  dat <- data.matrix(veteran)
  rf <- ranger(dependent.variable.name = "time", status.variable.name = "status", data = dat, write.forest = TRUE, num.trees = 10)
  expect_silent(predict(rf, dat))
})

test_that("growing works for single observations, survival", {
  rf <- ranger(Surv(time, status) ~ ., veteran[1, ], write.forest = TRUE, num.trees = 10)
  expect_is(rf$survival, "matrix")
})

test_that("predict works for single observations, survival", {
  rf <- ranger(Surv(time, status) ~ ., veteran, write.forest = TRUE, num.trees = 10)
  pred <- predict(rf, head(veteran, 1))
  expect_equal(length(pred$survival), length(rf$unique.death.times))
})

## Special tests for random forests for survival analysis
test_that("unique death times in survival result is right", {
  expect_equal(rg.surv$unique.death.times, sort(unique(veteran$time)))
})

test_that("C-index splitting works", {
  rf <- ranger(Surv(time, status) ~ ., data = veteran, 
               splitrule = "C", num.trees = 10)
  expect_equal(rf$treetype, "Survival")
})

test_that("C-index splitting not working on classification data", {
  expect_error(ranger(Species ~ ., iris, splitrule = "C", num.trees = 10))
})

test_that("Logrank splitting not working on classification data", {
  expect_error(ranger(Species ~ ., iris, splitrule = "logrank", num.trees = 10))
})

test_that("No error if survival tree without OOB observations", {
  dat <- data.frame(time = c(1,2), status = c(0,1), x = c(1,2))
  expect_silent(ranger(Surv(time, status) ~ ., dat, num.trees = 1, num.threads = 1))
})

test_that("predict.all for survival returns 3d array of size samples x times x trees", {
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5)
  pred <- predict(rf, veteran, predict.all = TRUE)
  
  expect_is(pred$survival, "array")
  expect_equal(dim(pred$survival), 
               c(nrow(veteran), length(pred$unique.death.times), rf$num.trees))
  expect_is(pred$chf, "array")
  expect_equal(dim(pred$chf), 
               c(nrow(veteran), length(pred$unique.death.times), rf$num.trees))
})

test_that("Mean of predict.all for survival is equal to forest prediction", {
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5)
  pred_forest <- predict(rf, veteran, predict.all = FALSE)
  pred_trees <- predict(rf, veteran, predict.all = TRUE)
  expect_equal(apply(pred_trees$chf, 1:2, mean), pred_forest$chf)
})

test_that("timepoints() function returns timepoints", {
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5)
  expect_equal(timepoints(rf), rf$unique.death.times)
  
  pred <- predict(rf, veteran)
  expect_equal(timepoints(pred), rf$unique.death.times)
})

test_that("timepoints() working on survival forest only", {
  rf <- ranger(Species ~ ., iris, num.trees = 5)
  expect_error(timepoints(rf), "No timepoints found. Object is no Survival forest.")
  
  pred <- predict(rf, iris)
  expect_error(timepoints(pred), "No timepoints found. Object is no Survival prediction object.")
})

test_that("Survival error without covariates", {
  expect_error(ranger(Surv(time, status) ~ ., veteran[, c("time", "status")], num.trees = 5), 
               "Error: No covariates found.")
})

test_that("Survival error for competing risk data", {
  sobj <- Surv(veteran$time, factor(sample(1:3, nrow(veteran), replace = TRUE)))
  expect_error(ranger(y = sobj, x = veteran[, 1:2], num.trees = 5), 
               "Error: Competing risks not supported yet\\. Use status=1 for events and status=0 for censoring\\.")
})
