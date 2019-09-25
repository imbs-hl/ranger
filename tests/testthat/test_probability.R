## Tests for random forests for probability estimation

library(ranger)
context("ranger_prob")

## Initialize random forest
train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
iris.train <- iris[train.idx, ]
iris.test <- iris[-train.idx, ]

rg.prob <- ranger(Species ~ ., data = iris.train, write.forest = TRUE, probability = TRUE)
prob <- predict(rg.prob, iris.test)

## Tests
test_that("probability estimations are a matrix with correct size", {
  expect_is(prob$predictions, "matrix")
  expect_equal(nrow(prob$predictions), nrow(iris.test))
  expect_equal(ncol(prob$predictions), length(rg.prob$forest$levels))
})

test_that("growing works for single observations, probability prediction", {
  expect_warning(rf <- ranger(Species ~ ., iris[1, ], write.forest = TRUE, probability = TRUE), 
                 "Dropped unused factor level\\(s\\) in dependent variable\\: versicolor\\, virginica\\.")
  expect_is(rf$predictions, "matrix")
})

test_that("probability estimations are between 0 and 1 and sum to 1", {
  expect_true(all(prob$predictions > -1e-5 & prob$predictions <= 1 + 1e-5))
  expect_equal(rowSums(prob$predictions), rep(1, nrow(prob$predictions)))
})

test_that("save.memory option works for probability", {
  rf <- ranger(Species ~ ., data = iris, probability = TRUE, save.memory = TRUE)
  expect_equal(rf$treetype, "Probability estimation")
})

test_that("predict works for single observations, probability prediction", {
  rf <- ranger(Species ~ ., iris, write.forest = TRUE, probability = TRUE)
  pred <- predict(rf, head(iris, 1))
  expect_is(pred$predictions, "matrix")
  expect_equal(names(which.max(pred$predictions[1, ])), as.character(iris[1,"Species"]))
  
  dat <- iris
  dat$Species <- as.numeric(dat$Species)
  rf <- ranger(Species ~ ., dat, write.forest = TRUE, probability = TRUE)
  pred <- predict(rf, head(dat, 1))
  expect_is(pred$predictions, "matrix")
  expect_equal(which.max(pred$predictions[1, ]), as.numeric(iris[1,"Species"]))
})

test_that("Probability estimation works correctly if labels are reversed", {
  ## Simulate data
  n <- 50
  a1 <- c(rnorm(n, 3, sd = 2), rnorm(n, 8, sd = 2))
  a2 <- c(rnorm(n, 8, sd = 2), rnorm(n, 3, sd = 2))
  
  ## create labels for data
  labels <- as.factor(c(rep("0", n), rep("1", n)))
  dat <- data.frame(label = labels, a1, a2)
  
  labels.rev <- as.factor(c(rep("1", n), rep("0", n))) 
  dat.rev <- data.frame(label = labels.rev, a1, a2)
  
  ## Train
  rf <- ranger(dependent.variable.name = "label", data = dat, probability = TRUE, 
               write.forest = TRUE, num.trees = 5)
  rf.rev <- ranger(dependent.variable.name = "label", data = dat.rev, probability = TRUE, 
                   write.forest = TRUE, num.trees = 5)
  
  ## Check OOB predictions
  expect_gte(mean(rf$predictions[1:n, "0"], na.rm = TRUE), 0.5)
  expect_gte(mean(rf$predictions[(n+1):(2*n), "1"], na.rm = TRUE), 0.5)
  
  expect_gte(mean(rf.rev$predictions[1:n, "1"], na.rm = TRUE), 0.5)
  expect_gte(mean(rf.rev$predictions[(n+1):(2*n), "0"], na.rm = TRUE), 0.5)
  
  ## Check predict() predictions
  pred <- predict(rf, dat)
  expect_gte(mean(pred$predictions[1:n, "0"], na.rm = TRUE), 0.5)
  expect_gte(mean(pred$predictions[(n+1):(2*n), "1"], na.rm = TRUE), 0.5)
  
  pred.rev <- predict(rf.rev, dat.rev)
  expect_gte(mean(pred.rev$predictions[1:n, "1"], na.rm = TRUE), 0.5)
  expect_gte(mean(pred.rev$predictions[(n+1):(2*n), "0"], na.rm = TRUE), 0.5)
})

test_that("Probability estimation works correctly if first or second factor level empty", {
  expect_warning(rf <- ranger(Species ~ ., iris[51:150, ], probability = TRUE), 
                 "^Dropped unused factor level\\(s\\) in dependent variable\\: setosa\\.")
  expect_silent(pred <- predict(rf, iris[101:150, ]))
  expect_gte(mean(pred$predictions[1:50, "virginica"], na.rm = TRUE), 0.9)
  
  expect_warning(rf <- ranger(Species ~ ., iris[c(101:150, 51:100), ], probability = TRUE), 
                 "^Dropped unused factor level\\(s\\) in dependent variable\\: setosa\\.")
  expect_silent(pred <- predict(rf, iris[c(101:150, 51:100), ]))
  expect_gte(mean(pred$predictions[1:50, "virginica"], na.rm = TRUE), 0.9)
  expect_gte(mean(pred$predictions[51:100, "versicolor"], na.rm = TRUE), 0.9)
})

test_that("No error if unused factor levels in outcome", {
  expect_warning(rf <- ranger(Species ~ ., iris[1:100, ], num.trees = 5, probability = TRUE),
                 "^Dropped unused factor level\\(s\\) in dependent variable\\: virginica\\.")
  pred <- predict(rf, iris)
  expect_equal(ncol(pred$predictions), 2)
})

test_that("predict.all for probability returns 3d array of size samples x classes x trees", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE, probability = TRUE)
  pred <- predict(rf, iris, predict.all = TRUE)
  expect_is(pred$predictions, "array")
  expect_equal(dim(pred$predictions), 
               c(nrow(iris), nlevels(iris$Species), rf$num.trees))
})

test_that("Mean of predict.all for probability is equal to forest prediction", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE, probability = TRUE)
  pred_forest <- predict(rf, iris, predict.all = FALSE)
  pred_trees <- predict(rf, iris, predict.all = TRUE)
  expect_equivalent(apply(pred_trees$predictions, 1:2, mean), pred_forest$predictions)
})
