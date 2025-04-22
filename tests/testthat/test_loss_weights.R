
library(ranger)
context("ranger_loss_weights")


# Regression: tiny dataset with known outcome
test_that("loss.weights yields correct regression predictions", {
  
  train <- data.frame(
    x = c(0, 1, 1),
    y = c(0, 10, 20)
  )
  test <- data.frame(x = c(0, 1))
  
  # 1) uniform weights: mean(10,20) = 15
  m1 <- ranger(
    formula = y ~ x, data = train, replace = FALSE,
    num.trees = 1, mtry = 1, min.node.size = 1, sample.fraction = 1, 
    loss.weights = c(1,1,1))
  
  p1 <- predict(m1, test)$predictions
  expect_equal(p1, c(0, 15))
  
  # 2) up‐weight the y=20 point by factor 4
  m2 <- ranger(
    formula = y ~ x, data = train, replace = FALSE,
    num.trees = 1, mtry = 1, min.node.size = 1, sample.fraction = 1, 
    loss.weights = c(1,1,4))
  
  p2 <- predict(m2, test)$predictions
  # now leaf at x=1 is (10*1 + 20*4)/(1+4) = 18
  expect_equal(p2, c(0, 18))
})


# Classification: tiny dataset with known classes
test_that("loss.weights affects classification majority vote", {
  
  #  x=0: one A (leaf always A)
  #  x=1: three samples: two A, one B
  train <- data.frame(
    x = c(0, 1, 1, 1),
    y = factor(c("A", "A", "A", "B"), levels = c("A", "B"))
  )
  test <- data.frame(x = c(0, 1))
  
  # 1) Uniform weights: in leaf x=1: 2*A vs 1*B -> predict A
  m1 <- ranger(
    formula = y ~ x, data = train, replace = FALSE,
    num.trees = 1, mtry = 1, min.node.size = 1, sample.fraction = 1,
    loss.weights = c(1,1,1,1), probability = TRUE)
  
  p1 <- predict(m1, test)$predictions
  # Uniform weights: both leaves predict A
  # A         B
  # 1.0000000 0.0000000
  # 0.6666667 0.3333333
  expect_equal(c(p1), c(1, 2/3, 0, 1/3))
  
  # 2) Up‑weight the lone B in the x=1 leaf -> vote flips to B
  m2 <- ranger(
    formula = y ~ x, data = train, replace = FALSE,
    num.trees = 1, mtry = 1, min.node.size = 1, sample.fraction = 1,
    loss.weights = c(1,1,1,4), probability = TRUE)
  
  p2 <- predict(m2, test)$predictions
  # Weighted towards B: second leaf predicts B
  # A         B
  # 1.0000000 0.0000000
  # 0.3333333 0.6666667
  expect_equal(c(p2), c(1, 1/3, 0, 2/3))
})
