## Tests for cross-platform reproducibility

library(ranger)
context("ranger_repro")

test_that("returns same predictions as before, R seed", {
  set.seed(2023)
  rf <- ranger(Species ~ ., iris, probability = TRUE)
  pred <- predict(rf, iris)
  
  expect_equal(pred$predictions[148:150, 3], c(1.0, 0.9871778, 0.9511746), tolerance = 1e-5)
})

test_that("returns same predictions as before, cpp seed", {
  set.seed(2023)
  rf <- ranger(Species ~ ., iris, probability = TRUE, seed = 2023)
  pred <- predict(rf, iris, seed = 2023)
  
  expect_equal(pred$predictions[148:150, 3], c(0.9997500, 0.9889944, 0.9401063), tolerance = 1e-5)
})



