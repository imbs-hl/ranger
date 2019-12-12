library(ranger)
library(survival)

context("importance_pvalues")

## GenABEL data
if (requireNamespace("GenABEL", quietly = TRUE)) {
  dat_gwaa <- readRDS("../test_gwaa.rds")
}

## 0 noise variables
rf_p0 <- ranger(Species ~., iris, num.trees = 5, 
                importance = "permutation", write.forest = TRUE)
holdout_p0 <- holdoutRF(Species ~., iris, num.trees = 5)

## 100 noise variables
n <- nrow(iris)
p <- 100
noise <- replicate(p, rnorm(n))
colnames(noise) <- paste0("noise", 1:p)
dat_n100 <- cbind(iris, noise)

rf_p100 <- ranger(Species ~., dat_n100, num.trees = 5,
                  importance = "permutation", write.forest = TRUE)
holdout_p100 <- holdoutRF(Species ~., dat_n100, num.trees = 5)

## General
test_that("Importance p-values Janitza: Error if impurity importance", {
  rf <- ranger(Species ~., iris, num.trees = 5, importance = "impurity")
  expect_error(importance_pvalues(rf, method = "janitza"))
})

## Janitza
test_that("Importance p-values Janitza: warning if few negative importance values", {
  expect_warning(importance_pvalues(rf_p100, method = "janitza"))
})

test_that("Importance p-values Janitza: returns correct dimensions", {
  expect_warning(vimp <- importance_pvalues(rf_p100, method = "janitza"))
  expect_is(vimp, "matrix")
  expect_equal(dim(vimp), c(104, 2))
})

test_that("Importance p-values Janitza: error if no importance", {
  rf_none <- ranger(Species ~., iris, num.trees = 5, importance = "none", write.forest = TRUE)
  expect_error(importance_pvalues(rf_none, method = "janitza"))
})
  
test_that("Importance p-values Janitza: error if Gini importance", {
  rf_imp <- ranger(Species ~., iris, num.trees = 5, importance = "impurity", write.forest = TRUE)
  expect_error(importance_pvalues(rf_imp, method = "janitza"))
})

test_that("Importance p-values Janitza: error if no unimportant variables", {
  # Error only when all importance values positive
  skip_if(any(rf_p0$variable.importance <= 0))

  expect_warning(expect_error(importance_pvalues(rf_p0, method = "janitza")))
})

test_that("Importance p-values Janitza: warning for regression", {
  rf <- ranger(Sepal.Length ~., dat_n100, num.trees = 5, importance = "permutation", write.forest = TRUE)
  expect_warning(importance_pvalues(rf, method = "janitza"))
})

test_that("Importance p-values Janitza-Holdout: returns correct dimensions", {
  expect_warning(vimp <- importance_pvalues(holdout_p100, method = "janitza"))
  expect_is(vimp, "matrix")
  expect_equal(dim(vimp), c(104, 2))
})

## Altmann
test_that("Importance p-values Altmann: returns correct dimensions", {
  vimp <- importance_pvalues(rf_p0, method = "altmann", formula = Species ~ ., data = iris)
  expect_is(vimp, "matrix")
  expect_equal(dim(vimp), c(4, 2))
})

test_that("Importance p-values Altmann: error if no importance", {
  rf_none <- ranger(Species ~., iris, num.trees = 5, importance = "none", write.forest = TRUE)
  expect_error(importance_pvalues(rf_none, method = "altmann", formula = Species ~ ., data = iris))
})

test_that("Importance p-values Altmann: not working for holdoutRF", {
  expect_error(importance_pvalues(holdout_p0, method = "altmann", formula = Species ~ ., data = iris))
})

test_that("Importance p-values Altmann: No zero p-values", {
  vimp <- importance_pvalues(rf_p0, method = "altmann", formula = Species ~ ., data = iris)
  expect_false(any(vimp[, "pvalue"] == 0))
})

test_that("Importance p-values Altmann: working with character formula", {
  vimp <- importance_pvalues(rf_p0, method = "altmann", formula = "Species ~ .", data = iris)
  expect_is(vimp, "matrix")
  expect_equal(dim(vimp), c(4, 2))
})

## Hold-out RF
test_that("HoldoutRF working", {
  expect_is(holdout_p0, "holdoutRF")
})

test_that("HoldoutRF working with GenABEL data", {
  skip_if_not_installed("GenABEL")
  holdout_gwaa <- holdoutRF(CHD ~., dat_gwaa, num.trees = 5)
  expect_is(holdout_p0, "holdoutRF")
})

test_that("HoldoutRF ... argument working", {
  rf <- holdoutRF(Species ~., iris, num.trees = 5)
  expect_equal(rf$rf1$num.trees, 5)
})

test_that("HoldoutRF working with formula", {
  rf <- holdoutRF(Species ~., iris, num.trees = 5)
  expect_equal(rf$rf1$treetype, "Classification")
  
  rf <- holdoutRF(Species ~., data = iris, num.trees = 5)
  expect_equal(rf$rf1$treetype, "Classification")
  
  rf <- holdoutRF(formula = Species ~., iris, num.trees = 5)
  expect_equal(rf$rf1$treetype, "Classification")
  
  rf <- holdoutRF(data = iris, formula = Species ~., num.trees = 5)
  expect_equal(rf$rf1$treetype, "Classification")
})

test_that("HoldoutRF working with dependent.variable.name", {
  rf <- holdoutRF(dependent.variable.name = "Species", data = iris, num.trees = 5)
  expect_equal(rf$rf1$treetype, "Classification")
  
  rf <- holdoutRF(data = iris, dependent.variable.name = "Species", num.trees = 5)
  expect_equal(rf$rf1$treetype, "Classification")
})

test_that("HoldoutRF not working if importance argument used", {
  expect_error(holdoutRF(Species ~., iris, num.trees = 5, importance = "impurity"), 
               "Error: Argument 'importance' not supported in holdoutRF.")
})

test_that("HoldoutRF not working if replace argument used", {
  expect_error(holdoutRF(Species ~., iris, num.trees = 5, replace = TRUE), 
               "Error: Argument 'replace' not supported in holdoutRF.")
})

## Survival, 0 noise variables
rf_p0_surv <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, 
                     importance = "permutation", write.forest = TRUE)
#holdout_p0_surv <- holdoutRF(Surv(time, status) ~ ., veteran, num.trees = 5)

## Survival, 100 noise variables
p <- 100
noise <- replicate(p, rnorm(nrow(veteran)))
colnames(noise) <- paste0("noise", 1:p)
dat_n100_surv <- cbind(veteran, noise)

rf_p100_surv <- ranger(Surv(time, status) ~., dat_n100_surv, num.trees = 5,
                       importance = "permutation", write.forest = TRUE)
#holdout_p100_surv <- holdoutRF(Surv(time, status) ~., dat_n100_surv, num.trees = 5)

test_that("Survival importance p-values Janitza: returns correct dimensions", {
  expect_warning(vimp <- importance_pvalues(rf_p100_surv, method = "janitza"))
  expect_is(vimp, "matrix")
  expect_equal(dim(vimp), c(106, 2))
})

test_that("Survival importance p-values Altmann: returns correct dimensions", {
  vimp <- importance_pvalues(rf_p0_surv, method = "altmann", formula = Surv(time, status) ~ ., data = veteran)
  expect_is(vimp, "matrix")
  expect_equal(dim(vimp), c(6, 2))
})

test_that("Survival importance p-values Altmann working with corrected impurity importance", {
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, 
               importance = "impurity_corrected")
  
  vimp <- importance_pvalues(rf, method = "altmann", formula = Surv(time, status) ~ ., data = veteran)
  expect_is(vimp, "matrix")
  expect_equal(dim(vimp), c(6, 2))
})

test_that("Survival importance p-values Janitza working with corrected impurity importance", {
  rf <- ranger(Surv(time, status) ~ ., dat_n100_surv, num.trees = 5, 
               importance = "impurity_corrected")
  
  expect_warning(vimp <- importance_pvalues(rf, method = "janitza"))
  expect_is(vimp, "matrix")
  expect_equal(dim(vimp), c(106, 2))
})

