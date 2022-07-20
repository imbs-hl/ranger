library(ranger)
library(survival)
context("ranger_unordered")

test_that("Old parameters still work", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D"), n, replace = TRUE), 
                   y = rbinom(n, 1, 0.5), 
                   stringsAsFactors = FALSE)
  
  rf.false <- ranger(y ~ ., data = dt, num.trees = 5, write.forest = TRUE, 
                     respect.unordered.factors = FALSE)
  rf.true <- ranger(y ~ ., data = dt, num.trees = 5, write.forest = TRUE, 
                    respect.unordered.factors = TRUE)
  
  expect_equal(rf.false$forest$covariate.levels$x, levels(factor(dt$x)))
  expect_equal(length(rf.true$forest$covariate.levels), 1)
})

test_that("If respect.unordered.factors='partition', regard characters as unordered", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D"), n, replace = TRUE), 
                   y = rbinom(n, 1, 0.5), 
                   stringsAsFactors = FALSE)
  
  set.seed(2)
  rf.char <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, respect.unordered.factors = 'partition')
  
  dt$x <- factor(dt$x, ordered = FALSE)
  set.seed(2)
  rf.fac <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, respect.unordered.factors = 'partition')
  
  expect_equal(rf.char$prediction.error, rf.fac$prediction.error)
})

test_that("If respect.unordered.factors='ignore', regard characters as ordered", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D"), n, replace = TRUE), 
                   y = rbinom(n, 1, 0.5), 
                   stringsAsFactors = FALSE)
  
  set.seed(2)
  rf.char <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, respect.unordered.factors = 'ignore')
  
  dt$x <- factor(dt$x, ordered = FALSE)
  set.seed(2)
  rf.fac <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, respect.unordered.factors = 'ignore')
  
  expect_equal(rf.char$prediction.error, rf.fac$prediction.error)
})

test_that("Error if other value for respect.unordered.factors", {
  expect_error(ranger(y ~ ., iris, num.trees = 5, respect.unordered.factors = NULL))
})

test_that("Same results if no unordered factors", {
  set.seed(100)
  rf1 <- ranger(Species ~ ., iris, num.trees = 5, respect.unordered.factors = 'ignore')
  set.seed(100)
  rf2 <- ranger(Species ~ ., iris, num.trees = 5, respect.unordered.factors = 'order')
  set.seed(100)
  rf3 <- ranger(Species ~ ., iris, num.trees = 5, respect.unordered.factors = 'partition')
  
  expect_equal(rf1$confusion.matrix, 
               rf2$confusion.matrix)
  expect_equal(rf1$confusion.matrix, 
               rf3$confusion.matrix)
})

test_that("Unordered splitting working for classification", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D"), n, replace = TRUE), 
                   y = factor(rbinom(n, 1, 0.5)),
                   stringsAsFactors = FALSE)
  
  rf <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, respect.unordered.factors = 'partition')
  expect_true(any(!rf$forest$is.ordered))
})

test_that("Unordered splitting working for probability", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D"), n, replace = TRUE), 
                   y = factor(rbinom(n, 1, 0.5)),
                   stringsAsFactors = FALSE)
  
  rf <- ranger(y ~ ., data = dt, num.trees = 5, min.node.size = n/2, 
               respect.unordered.factors = 'partition', probability = TRUE)
  expect_true(any(!rf$forest$is.ordered))
})

test_that("Order splitting working for multiclass classification", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D"), n, replace = TRUE), 
                   y = factor(sample(c("A", "B", "C", "D"), n, replace = TRUE)),
                   stringsAsFactors = FALSE)
  
  rf <- ranger(y ~ ., data = dt, num.trees = 5, 
               respect.unordered.factors = 'order', probability = TRUE)
  expect_true(all(rf$forest$is.ordered))
})

test_that("Order splitting working for multiclass probability", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D"), n, replace = TRUE), 
                   y = factor(sample(c("A", "B", "C", "D"), n, replace = TRUE)),
                   stringsAsFactors = FALSE)
  
  rf <- ranger(y ~ ., data = dt, num.trees = 5, 
               respect.unordered.factors = 'order', probability = TRUE)
  expect_true(all(rf$forest$is.ordered))
})

test_that("Order splitting working with alternative interface", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C", "D"), n, replace = TRUE), 
                   y = factor(sample(c("A", "B", "C", "D"), n, replace = TRUE)),
                   stringsAsFactors = FALSE)
  
  rf <- ranger(dependent.variable.name = "y", data = dt, num.trees = 5, respect.unordered.factors = 'order')
  expect_true(all(rf$forest$is.ordered))
  
  rf <- ranger(dependent.variable.name = "y", data = dt, num.trees = 5, respect.unordered.factors = 'order', probability = TRUE)
  expect_true(all(rf$forest$is.ordered))
})

test_that("Order splitting working with single level factor", {
  n <- 20
  
  # Binary classification
  dt_class <- data.frame(x = sample(c("A"), n, replace = TRUE), 
                         y = factor(sample(c("A", "B"), n, replace = TRUE)),
                         stringsAsFactors = FALSE)
  expect_silent(ranger(y ~ ., data = dt_class, num.trees = 5, 
                       respect.unordered.factors = 'order', probability = FALSE))
  expect_silent(ranger(y ~ ., data = dt_class, num.trees = 5, 
                       respect.unordered.factors = 'order', probability = TRUE))
  
  # Multiclass classification
  dt_mult <- data.frame(x = sample(c("A"), n, replace = TRUE), 
                        y = factor(sample(c("A", "B", "C", "D"), n, replace = TRUE)),
                        stringsAsFactors = FALSE)
  expect_silent(ranger(y ~ ., data = dt_class, num.trees = 5, 
                       respect.unordered.factors = 'order', probability = FALSE))
  expect_silent(ranger(y ~ ., data = dt_class, num.trees = 5, 
                       respect.unordered.factors = 'order', probability = TRUE))
  
  # Regression
  dt_cont <- data.frame(x = sample(c("A"), n, replace = TRUE), 
                        y = rnorm(n),
                        stringsAsFactors = FALSE)
  expect_silent(ranger(y ~ ., data = dt_cont, num.trees = 5, 
                       respect.unordered.factors = 'order'))
  
  # Survival
  dt_surv <- data.frame(x = sample(c("A"), n, replace = TRUE), 
                        time = rnorm(n),
                        status = rbinom(n, 1, .5),
                        stringsAsFactors = FALSE)
  expect_silent(ranger(Surv(time, status) ~ ., data = dt_surv, num.trees = 5, 
                       respect.unordered.factors = 'order'))
})

test_that("Unordered splitting working for survival", {
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, min.node.size = 50, respect.unordered.factors = 'partition')
  expect_true(any(!rf$forest$is.ordered))
})

test_that("Order splitting working for survival", {
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, min.node.size = 50, respect.unordered.factors = 'order')
  expect_true(all(rf$forest$is.ordered))
})

test_that("Order splitting working for survival with alternative interface", {
  rf <- ranger(dependent.variable.name = "time", status.variable.name = "status", data = veteran, num.trees = 5, min.node.size = 50, respect.unordered.factors = 'order')
  expect_true(all(rf$forest$is.ordered))
})

test_that("Error if too many factors in 'partition' mode", {
  n <- 100
  dt <- data.frame(x = factor(1:100, ordered = FALSE),  
                   y = rbinom(n, 1, 0.5))
  
  expect_error(ranger(y ~ ., data = dt, num.trees = 5, respect.unordered.factors = 'partition'))
})

test_that("Survival forest with 'order' mode works", {
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, 
                              write.forest = TRUE, respect.unordered.factors = 'order')
  expect_equal(sort(rf$forest$covariate.levels$celltype), 
               sort(levels(veteran$celltype)))
    
  pred <- predict(rf, veteran)
  expect_is(pred, "ranger.prediction")
})

test_that("maxstat splitting not working with unordered factors", {
  expect_error(ranger(Sepal.Length ~ ., iris, splitrule = "maxstat", respect.unordered.factors = "partition"))
  expect_error(ranger(Surv(time, status) ~ ., veteran, splitrule = "maxstat", respect.unordered.factors = "partition"))
})

test_that("C splitting not working with unordered factors", {
  expect_error(ranger(Surv(time, status) ~ ., veteran, splitrule = "C", respect.unordered.factors = "partition"))
})

test_that("Warning for maxstat with 'order' mode", {
  expect_warning(ranger(Sepal.Length ~ ., iris, num.trees = 5, splitrule = "maxstat",
                        respect.unordered.factors = 'order'))
})

test_that("No warning for multiclass classification/probability and survival with 'order' mode", {
  expect_silent(ranger(Species ~ ., iris, num.trees = 5, 
                        respect.unordered.factors = 'order'))
  expect_silent(ranger(Species ~ ., iris, num.trees = 5, probability = TRUE,
                        respect.unordered.factors = 'order'))
  expect_silent(ranger(Surv(time, status) ~ ., veteran, num.trees = 5, 
                        respect.unordered.factors = 'order'))
})

test_that("No error if new levels in predict, 1 column", {
  set.seed(1)
  n <- 20
  train <- data.frame(x = sample(c("A", "B", "C"), n, replace = TRUE), 
                   y = rbinom(n, 1, 0.5), 
                   stringsAsFactors = FALSE)
  
  test <- data.frame(x = sample(c("A", "B", "C", "D"), n, replace = TRUE), 
                      y = rbinom(n, 1, 0.5), 
                      stringsAsFactors = FALSE)
  
  ## ignore
  rf.ignore <- ranger(y ~ ., data = train, num.trees = 5, respect.unordered.factors = 'ignore')
  expect_silent(predict(rf.ignore, test))

  ## partition  
  rf.partition <- ranger(y ~ ., data = train, num.trees = 5, respect.unordered.factors = 'partition')
  expect_silent(predict(rf.partition, test))
  
  ## order
  rf.order <- ranger(y ~ ., data = train, num.trees = 5, respect.unordered.factors = 'order')
  expect_silent(predict(rf.order, test))
})

test_that("No error if new levels in predict, 2 columns", {
  set.seed(1)
  n <- 20
  train <- data.frame(x1 = sample(c("A", "B", "C"), n, replace = TRUE), 
                      x2 = sample(c("A", "B", "C"), n, replace = TRUE),
                      y = rbinom(n, 1, 0.5), 
                      stringsAsFactors = FALSE)
  
  test <- data.frame(x1 = sample(c("A", "B", "C", "D"), n, replace = TRUE), 
                     x2 = sample(c("A", "B", "C", "D"), n, replace = TRUE), 
                     y = rbinom(n, 1, 0.5), 
                     stringsAsFactors = FALSE)
  
  ## ignore
  rf.ignore <- ranger(y ~ ., data = train, num.trees = 5, respect.unordered.factors = 'ignore')
  expect_silent(predict(rf.ignore, test))
  
  ## partition  
  rf.partition <- ranger(y ~ ., data = train, num.trees = 5, respect.unordered.factors = 'partition')
  expect_silent(predict(rf.partition, test))
  
  ## order
  rf.order <- ranger(y ~ ., data = train, num.trees = 5, respect.unordered.factors = 'order')
  expect_silent(predict(rf.order, test))
})

test_that("No error if NA factor levels and order", {
  df <- data.frame(x = addNA(factor(c("a", "a", NA, NA, "b", "b"))),
                   y = c(1, 2, 3, 4, 5, 6))
  expect_silent(ranger(dependent.variable.name = "y", data = df, 
                respect.unordered.factors = "order"))
})

test_that("Order splitting working when numerics in data", {
  n <- 20
  
  # Binary classification
  dt_class <- data.frame(x1 = sample(c("A", "B", "C"), n, replace = TRUE), 
                         x2 = sample(1:3, n, replace = TRUE),
                         y = factor(sample(c("A", "B"), n, replace = TRUE)),
                         stringsAsFactors = FALSE)
  rf <- expect_silent(ranger(y ~ ., data = dt_class, num.trees = 5, 
                       respect.unordered.factors = 'order', probability = FALSE))
  expect_silent(predict(rf, dt_class))
  rf <- expect_silent(ranger(y ~ ., data = dt_class, num.trees = 5, 
                       respect.unordered.factors = 'order', probability = TRUE))
  expect_silent(predict(rf, dt_class))
  
  # Multiclass classification
  dt_mult <- data.frame(x1 = sample(c("A", "B", "C"), n, replace = TRUE), 
                        x2 = sample(1:3, n, replace = TRUE),
                        y = factor(sample(c("A", "B", "C", "D"), n, replace = TRUE)),
                        stringsAsFactors = FALSE)
  rf <- expect_silent(ranger(y ~ ., data = dt_class, num.trees = 5, 
                       respect.unordered.factors = 'order', probability = FALSE))
  expect_silent(predict(rf, dt_mult))
  rf <- expect_silent(ranger(y ~ ., data = dt_class, num.trees = 5, 
                       respect.unordered.factors = 'order', probability = TRUE))
  expect_silent(predict(rf, dt_mult))
  
  # Regression
  dt_cont <- data.frame(x1 = sample(c("A", "B", "C"), n, replace = TRUE), 
                        x2 = sample(1:3, n, replace = TRUE),
                        y = rnorm(n),
                        stringsAsFactors = FALSE)
  rf <- expect_silent(ranger(y ~ ., data = dt_cont, num.trees = 5, 
                       respect.unordered.factors = 'order'))
  expect_silent(predict(rf, dt_cont))
  
  # Survival
  dt_surv <- data.frame(x1 = sample(c("A", "B", "C"), n, replace = TRUE), 
                        x2 = as.numeric(sample(1:3, n, replace = TRUE)),
                        time = rnorm(n),
                        status = rbinom(n, 1, .5),
                        stringsAsFactors = FALSE)
  rf <- expect_silent(ranger(Surv(time, status) ~ ., data = dt_surv, num.trees = 5, 
                       respect.unordered.factors = 'order'))
  expect_silent(predict(rf, dt_surv))
  
  # Survival with Surv() in data
  dt_surv <- data.frame(x1 = sample(c("A", "B", "C"), n, replace = TRUE), 
                        x2 = as.numeric(sample(1:3, n, replace = TRUE)),
                        y = Surv(rnorm(n), rbinom(n, 1, .5)),
                        stringsAsFactors = FALSE)
  rf <- expect_silent(ranger(y ~ ., data = dt_surv, num.trees = 5, 
                             respect.unordered.factors = 'order'))
  expect_silent(predict(rf, dt_surv))
})

test_that("Partition splitting working for large number of levels", {
  n <- 43
  dt <- data.frame(x = factor(1:n, ordered = FALSE),  
                   y = rbinom(n, 1, 0.5))
  
  rf <- ranger(y ~ ., data = dt, num.trees = 10, splitrule = "extratrees")
  
  max_split <- max(sapply(1:rf$num.trees, function(i) {
    max(log2(rf$forest$split.values[[i]]), na.rm = TRUE)
  }))

  expect_lte(max_split, n)
})

test_that("Order splitting working for quantreg forests", {
  n <- 20
  dt <- data.frame(x = sample(c("A", "B", "C"), n, replace = TRUE),
                   y = rbinom(n, 1, 0.5), 
                   stringsAsFactors = TRUE)
  expect_silent(ranger(y ~ ., dt, num.trees = 5,
                       quantreg = TRUE, respect.unordered.factors = 'order'))
})
