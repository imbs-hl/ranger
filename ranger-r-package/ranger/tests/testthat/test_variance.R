library(ranger)

context("infinitesimal jackknife")

test_that("var_IJ_U returns variance of prediction", {
  train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
  iris.train <- iris[train.idx, ]
  iris.test <- iris[-train.idx, ]
  
  rf <- ranger(Species ~ ., data = iris.train, num.trees = 5, write.forest = TRUE, 
               keep.inbag = TRUE)
  
  pred <- predict(rf, iris.test, predict.all = TRUE, var_IJ = TRUE)
  
  expect_that(pred$var_ij_u, is_a("numeric"))
  expect_that(length(pred$var_ij_u), equals(nrow(iris.test)))
})