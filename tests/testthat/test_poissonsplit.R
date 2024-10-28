library(ranger)
context("ranger_poisson")

# Generate poisson distributed outcome
set.seed(42)
n <- 1000
p <- 4
beta <- c(0, 0.1, -0.2, 0.3)
x <- replicate(p, runif(n))
# Use exp(..) to keep it positive (and adds interactions).
# Add -1 to make it small as Poisson should be better for small frequencies.
lambda <- exp(-1 + as.vector(x %*% beta))
y <- rpois(n, lambda)
df <- data.frame(y = y, x)

# And a simple dataset with zero outcomes
df2 <- data.frame(y = c(0, 0, 0, 0, 0, 1, 2, 3, 4, 5),
                  x1 = c("a", "a", "a", "a", "a", "b", "b", "b", "b", "b"),
                  x2 = c(0, 0, 0, 0, 0, 1, 1, 1, 2, 2))

poisson_deviance <- function(y_true, y_pred) {
  if (any(y_true == 0 & y_pred == 0)) {
    stop("Error: Poisson deviance does not exist for y_pred == y_true == 0.")
  }
  pos <- y_true > 0
  dev <- y_pred
  dev[pos] <- y_true[pos] * log(y_true[pos] / y_pred[pos]) - y_true[pos] + y_pred[pos]
  return(2 * mean(dev))
}


test_that("poisson splitting works on poisson distributed data", {
  n_train = 1:(4*n %/% 5)
  n_test = (max(n_train)+1):n
  df_train = df[n_train, ]
  df_test = df[n_test, ]
  rf_poi <- ranger(y ~ ., df_train, splitrule = "poisson", num.trees = 50, min.node.size = 50, poisson.tau = 1, seed = 123)
  rf_mse <- ranger(y ~ ., df_train, splitrule = "variance", num.trees = 50, min.node.size = 50, seed = 123)
  
  expect_is(rf_poi, "ranger")
  # deviance on test set
  expect_lt(poisson_deviance(df_test$y, predict(rf_poi, df_test)$predictions),
            poisson_deviance(df_test$y, predict(rf_mse, df_test)$predictions))
})

test_that("poisson splitting not working for negative outcome", {
  expect_error(ranger(y ~ ., data.frame(y = c(-1.5, 2), x = c(1, 2)), splitrule = "poisson"))
  expect_error(ranger(y ~ ., data.frame(y = c(0, 0), x = c(1, 2)), splitrule = "poisson"))
})

test_that("poisson.tau <= 0 throws error", {
  expect_error(ranger(y ~ ., df2, poisson.tau = 0))
})

test_that("poisson splitting predicts positive even on nodes with all values equal 0", {
  rf <- ranger(y ~ ., df2, splitrule = "poisson", poisson.tau = 0.1, mtry = 2, num.trees = 2,
               min.node.size = 1, seed = 123)
  expect_true(all(c(predict(rf, df2, predict.all = TRUE)$predictions) > 0))
})

test_that("poisson splitting gives larger predictions for larger values of poisson.tau on pure nodes with y = 0", {
  rf1 <- ranger(y ~ ., df2, splitrule = "poisson", poisson.tau = 0.1, mtry = 2, num.trees = 2,
               min.node.size = 1, seed = 123)
  rf2 <- ranger(y ~ ., df2, splitrule = "poisson", poisson.tau = 10, mtry = 2, num.trees = 2,
               min.node.size = 1, seed = 123)
  expect_true(all(predict(rf2, df2)$predictions[df2$y == 0] > predict(rf1, df2)$predictions[df2$y == 0]))
})
