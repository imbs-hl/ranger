## Tests for print function

library(rangerts)
context("rangerts_print")

## Initialize the random forest
rf <- rangerts(Species ~ ., iris, num.trees = 5, write.forest = TRUE)

## Test print rangerts function
expect_that(print(rf), prints_text("rangerts result"))

## Test print forest function
expect_that(print(rf$forest), prints_text("rangerts forest object"))

## Test print prediction function
expect_that(print(predict(rf, iris)), prints_text("rangerts prediction"))

## Test str rangerts function
expect_that(str(rf), prints_text("List of 14"))

## Test str forest function
expect_that(str(rf$forest), prints_text("List of 9"))
