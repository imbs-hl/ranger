## Tests for character data

library(ranger)
context("ranger_char")

## Initialize random forests
dat <- iris
dat$Test <- paste0("AA",as.character(1:nrow(dat)))

## Tests
test_that("no warning or error if character vector in data", {
  expect_silent(rf <- ranger(Species ~ ., dat, 
                             num.trees = 5, write.forest = TRUE))
  expect_silent(predict(rf, dat))
})

test_that("no warning or error if character vector in data, alternative interface", {
  expect_silent(rf <- ranger(dependent.variable.name = "Species", data = dat, 
                             num.trees = 5, write.forest = TRUE))
  expect_silent(predict(rf, dat))
})

test_that("same result for single row", {
  rf <- ranger(Species ~ ., dat, num.trees = 5, probability = TRUE)
  pred1 <- predict(rf, dat)$predictions
  pred2 <- t(sapply(1:nrow(dat), function(i) {
    predict(rf, dat[i, ])$predictions
  }))
  colnames(pred2) <- colnames(pred1)
  expect_equal(pred1, pred2)
})
