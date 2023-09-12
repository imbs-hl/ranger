## Tests for node statistics

library(ranger)
library(survival)
context("ranger_nodestats")

test_that("if node.stats FALSE, no nodestats saved, classification", {
  rf <- ranger(Species ~ ., iris, num.trees = 5)
  expect_null(rf$forest$num.samples.nodes)
  expect_null(rf$forest$node.predictions)
})

test_that("if node.stats FALSE, no nodestats saved, probability", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, probability = TRUE)
  expect_null(rf$forest$num.samples.nodes)
  expect_null(rf$forest$node.predictions)
  expect_length(rf$forest$terminal.class.counts[[1]][[1]], 0)
})

test_that("if node.stats FALSE, no nodestats saved, regression", {
  rf <- ranger(Sepal.Length ~ ., iris, num.trees = 5)
  expect_null(rf$forest$num.samples.nodes)
  expect_null(rf$forest$node.predictions)
})

test_that("if node.stats FALSE, no nodestats saved, survival", {
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5)
  expect_null(rf$forest$num.samples.nodes)
  expect_null(rf$forest$node.predictions)
  expect_length(rf$forest$chf[[1]][[1]], 0)
})

test_that("if node.stats TRUE, nodestats saved, classification", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, node.stats = TRUE)
  
  expect_is(rf$forest$num.samples.nodes, "list")
  expect_length(rf$forest$num.samples.nodes, rf$num.trees)
  expect_equal(rf$forest$num.samples.nodes[[1]][1], nrow(iris))
  
  expect_is(rf$forest$node.predictions, "list")
  expect_length(rf$forest$node.predictions, rf$num.trees)
  expect_is(rf$forest$node.predictions[[1]], "numeric")
})

test_that("if node.stats TRUE, nodestats saved, probability", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, probability = TRUE, node.stats = TRUE)
  
  expect_is(rf$forest$num.samples.nodes, "list")
  expect_length(rf$forest$num.samples.nodes, rf$num.trees)
  expect_equal(rf$forest$num.samples.nodes[[1]][1], nrow(iris))
  
  expect_null(rf$forest$node.predictions)
  
  expect_is(rf$forest$terminal.class.counts, "list")
  expect_length(rf$forest$terminal.class.counts, rf$num.trees)
  expect_length(rf$forest$terminal.class.counts[[1]][[1]], nlevels(iris$Species))
})

test_that("if node.stats TRUE, nodestats saved, regression", {
  rf <- ranger(Sepal.Length ~ ., iris, num.trees = 5, node.stats = TRUE)
  
  expect_is(rf$forest$num.samples.nodes, "list")
  expect_length(rf$forest$num.samples.nodes, rf$num.trees)
  expect_equal(rf$forest$num.samples.nodes[[1]][1], nrow(iris))
  
  expect_is(rf$forest$node.predictions, "list")
  expect_length(rf$forest$node.predictions, rf$num.trees)
  expect_is(rf$forest$node.predictions[[1]], "numeric")
})

test_that("if node.stats TRUE, nodestats saved, survival", {
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, node.stats = TRUE)
  
  expect_is(rf$forest$num.samples.nodes, "list")
  expect_length(rf$forest$num.samples.nodes, rf$num.trees)
  expect_equal(rf$forest$num.samples.nodes[[1]][1], nrow(veteran))
  
  expect_null(rf$forest$node.predictions)
  
  expect_is(rf$forest$chf, "list")
  expect_length(rf$forest$chf, rf$num.trees)
  expect_is(rf$forest$chf[[1]], "list")
  expect_is(rf$forest$chf[[1]][[1]], "numeric")
  expect_length(rf$forest$chf[[1]][[1]], length(rf$unique.death.times))
})



rf <- ranger(Species ~ ., iris, num.trees = 10, probability = TRUE, node.stats = TRUE)
rf$forest$num.samples.nodes
rf$forest$node.predictions
rf$forest$terminal.class.counts


rf <- ranger(Sepal.Length ~ ., iris, num.trees = 10, node.stats = TRUE)
rf$forest$num.samples.nodes
rf$forest$node.predictions

# Survival

rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 10, node.stats = TRUE)
rf$forest$num.samples.nodes
rf$forest$node.predictions
rf$forest$chf



