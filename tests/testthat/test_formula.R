library(rangerts)

test_that("LHS formula", {
  callrangerts <- function() {
    myTransformation <- function(x) { x }
    rangerts(myTransformation(Species) ~ ., data = iris)
  }

  expect_error(callrangerts(), NA)
})
