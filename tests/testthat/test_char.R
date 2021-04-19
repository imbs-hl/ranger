## Tests for character data

library(rangerts)
context("rangerts_char")

## Initialize random forests
dat <- iris
dat$Test <- paste0("AA",as.character(1:nrow(dat)))

## Tests
test_that("no warning or error if character vector in data", {
  expect_silent(rf <- rangerts(Species ~ ., dat,
                             num.trees = 5, write.forest = TRUE))
  expect_silent(predict(rf, dat))
})

test_that("no warning or error if character vector in data, alternative interface", {
  expect_silent(rf <- rangerts(dependent.variable.name = "Species", data = dat,
                             num.trees = 5, write.forest = TRUE))
  expect_silent(predict(rf, dat))
})
