library(ranger)
library(survival)
context("ranger_unordered")

test_that("Third child for missings only there if missings in data", {
  rf1 <- ranger(Species ~ ., iris, num.trees = 5)
  expect_length(rf1$forest$child.nodeIDs[[1]], 2)
  
  dat <- iris
  dat[1, 1] <- NA
  rf2 <- ranger(Species ~ ., dat, num.trees = 5)
  expect_length(rf2$forest$child.nodeIDs[[1]], 3)
})

test_that("Training works with missing values in x but not in y", {
  dat <- iris
  dat[25, 1] <- NA
  expect_silent(ranger(Species ~ ., dat, num.trees = 5))
  expect_silent(ranger(Petal.Width ~ ., dat, num.trees = 5))
  expect_error(ranger(Sepal.Length ~ ., dat, num.trees = 5), "Missing data in dependent variable.")
  
  dat <- iris
  dat[4, 5] <- NA
  expect_error(ranger(Species ~ ., dat, num.trees = 5), "Missing data in dependent variable.")
})

test_that("No error if missing value in irrelevant column, training", {
  dat <- iris
  dat[1, "Sepal.Width"] <- NA
  expect_silent(ranger(Species ~ Sepal.Length, dat, num.trees = 5))
})

test_that("No error if missing value in irrelevant column, prediction", {
  rf <- ranger(Species ~ Sepal.Length, iris, num.trees = 5)
  dat <- iris
  dat[1, "Sepal.Width"] <- NA
  expect_silent(predict(rf, dat))
})

test_that("Prediction works with missing values, classification", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  
  dat <- iris
  dat[4, 4] <- NA
  dat[25, 1] <- NA
  expect_silent(predict(rf, dat))
})

test_that("Prediction works with missing values, regression", {
  rf <- ranger(Sepal.Width ~ ., iris, num.trees = 5, write.forest = TRUE)
  
  dat <- iris
  dat[4, 4] <- NA
  dat[25, 1] <- NA
  expect_silent(predict(rf, dat))
})

test_that("Order splitting working with missing values for classification", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D", NA), n, replace = TRUE), 
                   y = factor(rbinom(n, 1, 0.5)),
                   stringsAsFactors = FALSE)
  
  rf <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, respect.unordered.factors = 'order')
  expect_true(all(rf$forest$is.ordered))
})

test_that("Order splitting working with missing values for multiclass classification", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D", NA), n, replace = TRUE), 
                   y = factor(sample(c("A", "B", "C", "D"), n, replace = TRUE)),
                   stringsAsFactors = FALSE)
  
  rf <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, respect.unordered.factors = 'order')
  expect_true(all(rf$forest$is.ordered))
})

test_that("Missing values for survival not yet working", {
  dat <- veteran
  dat[1, 1] <- NA
  
  expect_error(ranger(Surv(time, status) ~ ., dat, num.trees = 5), "Error: Missing value handling not yet implemented for survival forests\\.")
})

test_that("na.omit leads to same result as manual removal, classification", {
  dat <- iris
  dat[1, 1] <- NA
  rf1 <- ranger(Species ~ ., dat, num.trees = 5, seed = 10, na.action = "na.omit")
  
  dat2 <- na.omit(dat)
  rf2 <- ranger(Species ~ ., dat2, num.trees = 5, seed = 10)
  
  expect_equal(rf1$predictions, rf2$predictions)
})

test_that("na.omit leads to same result as manual removal, probability", {
  dat <- iris
  dat[1, 1] <- NA
  rf1 <- ranger(Species ~ ., dat, num.trees = 5, probability = TRUE, seed = 10, na.action = "na.omit")
  
  dat2 <- na.omit(dat)
  rf2 <- ranger(Species ~ ., dat2, num.trees = 5, probability = TRUE, seed = 10)
  
  expect_equal(rf1$predictions, rf2$predictions)
})

test_that("na.omit leads to same result as manual removal, regression", {
  dat <- iris
  dat[1, 1] <- NA
  rf1 <- ranger(Sepal.Width ~ ., dat, num.trees = 5, seed = 10, na.action = "na.omit")
  
  dat2 <- na.omit(dat)
  rf2 <- ranger(Sepal.Width ~ ., dat2, num.trees = 5, seed = 10)
  
  expect_equal(rf1$predictions, rf2$predictions)
})

test_that("na.omit leads to same result as manual removal, survival", {
  dat <- veteran
  dat[1, 1] <- NA
  rf1 <- ranger(Surv(time, status) ~ ., dat, num.trees = 5, seed = 10, na.action = "na.omit")
  
  dat2 <- na.omit(dat)
  rf2 <- ranger(Surv(time, status) ~ ., dat2, num.trees = 5, seed = 10)
  
  expect_equal(rf1$chf, rf2$chf)
})

test_that("na.omit not working if no observations left", {
  dat <- iris
  dat[1:150, 1] <- NA
  expect_error(ranger(Species ~ ., dat, num.trees = 5, na.action = "na.omit"), "Error: No observations left after removing missing values\\.")
})

test_that("NA vs no-NA split is done with MIA but not without, classification", {
  x <- c(rep(NA, 10), rep(1, 10))
  y <- factor(c(rep(0, 10), rep(1, 10)))
  dat <- data.frame(x = x, y = y)
  
  # Should not split at all
  rf <- ranger(y~., dat, num.trees = 1, sample.fraction = 1, replace = FALSE, na.action = "na.learn", mia = FALSE)
  expect_length(rf$forest$split.values[[1]], 1)
  
  # Should split NA vs non-NA (inf. split value)
  rf <- ranger(y~., dat, num.trees = 1, sample.fraction = 1, replace = FALSE, na.action = "na.learn", mia = TRUE)
  expect_length(rf$forest$split.values[[1]], 3)
  expect_equal(rf$forest$split.values[[1]][[1]], Inf)
})

