library(ranger)

test_that("LHS formula", {
  callRanger <- function() {
    myTransformation <- function(x) { x }
    ranger(myTransformation(Species) ~ ., data = iris)
  }
  
  expect_error(callRanger(), NA)
})
