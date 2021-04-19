## Tests for using seeds

library(rangerts)
context("rangerts_seed")

## Initialize the random forests
ind = 1:150 %in% sample(150, 100)

set.seed(2)
mod1 = rangerts(Species ~ ., data = iris[ind, ], write.forest = TRUE, num.trees = 50)
pred1 = predict(mod1, data = iris[!ind, ])

set.seed(2)
mod2 = rangerts(Species ~ ., data = iris[ind, ], write.forest = TRUE, num.trees = 50)
pred2 = predict(mod2, data = iris[!ind, ])

set.seed(2)
mod3 = rangerts(dependent.variable.name = "Species", data = iris[ind, ], write.forest = TRUE, num.trees = 50)
pred3 = predict(mod3, data = iris[!ind, ])

## Tests
test_that("same result with same seed", {
  expect_equal(pred1$predictions, pred2$predictions)
})

test_that("same result with same seed, different interface", {
  expect_equal(pred1$predictions, pred3$predictions)
})
