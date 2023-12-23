library(ranger)
library(survival)
context("ranger_treeInfo")

## Classification
rf.class.formula <- ranger(Species ~ ., iris, num.trees = 5)
rf.class.first <- ranger(dependent.variable.name = "Species", data = iris[, c(5, 1:4)], num.trees = 5)
rf.class.mid <- ranger(dependent.variable.name = "Species", data = iris[, c(1:2, 5, 3:4)], num.trees = 5)
rf.class.last <- ranger(dependent.variable.name = "Species", data = iris, num.trees = 5)

ti.class.formula <- treeInfo(rf.class.formula)
ti.class.first <- treeInfo(rf.class.first)
ti.class.mid <- treeInfo(rf.class.mid)
ti.class.last <- treeInfo(rf.class.last)

test_that("Terminal nodes have only prediction, non-terminal nodes all others, classification formula", {
  expect_true(all(is.na(ti.class.formula[ti.class.formula$terminal, 2:6])))
  expect_true(all(!is.na(ti.class.formula[ti.class.formula$terminal, c(1, 7:8)])))
  expect_true(all(!is.na(ti.class.formula[!ti.class.formula$terminal, -8])))
  expect_true(all(is.na(ti.class.formula[!ti.class.formula$terminal, 8])))
})
test_that("Terminal nodes have only prediction, non-terminal nodes all others, classification depvarname first", {
  expect_true(all(is.na(ti.class.first[ti.class.first$terminal, 2:6])))
  expect_true(all(!is.na(ti.class.first[ti.class.first$terminal, c(1, 7:8)])))
  expect_true(all(!is.na(ti.class.first[!ti.class.first$terminal, -8])))
  expect_true(all(is.na(ti.class.first[!ti.class.first$terminal, 8])))
})
test_that("Terminal nodes have only prediction, non-terminal nodes all others, classification depvarname mid", { 
  expect_true(all(is.na(ti.class.mid[ti.class.mid$terminal, 2:6])))
  expect_true(all(!is.na(ti.class.mid[ti.class.mid$terminal, c(1, 7:8)])))
  expect_true(all(!is.na(ti.class.mid[!ti.class.mid$terminal, -8])))
  expect_true(all(is.na(ti.class.mid[!ti.class.mid$terminal, 8])))
})
test_that("Terminal nodes have only prediction, non-terminal nodes all others, classification depvarname last", { 
  expect_true(all(is.na(ti.class.last[ti.class.last$terminal, 2:6])))
  expect_true(all(!is.na(ti.class.last[ti.class.last$terminal, c(1, 7:8)])))
  expect_true(all(!is.na(ti.class.last[!ti.class.last$terminal, -8])))
  expect_true(all(is.na(ti.class.last[!ti.class.last$terminal, 8])))
})

test_that("Names in treeInfo match, classification", {
  varnames <- colnames(iris)[1:4]
  expect_true(all(is.na(ti.class.formula$splitvarName) | ti.class.formula$splitvarName %in% varnames))
  expect_true(all(is.na(ti.class.first$splitvarName) | ti.class.first$splitvarName %in% varnames))
  expect_true(all(is.na(ti.class.mid$splitvarName) | ti.class.mid$splitvarName %in% varnames))
  expect_true(all(is.na(ti.class.last$splitvarName) | ti.class.last$splitvarName %in% varnames))
})

test_that("Prediction for classification is factor with correct levels", {
  expect_is(ti.class.formula$prediction, "factor")
  expect_equal(levels(ti.class.formula$prediction), levels(iris$Species))
})

test_that("Prediction for classification is same as class prediction", {
  dat <- iris[sample(nrow(iris)), ]
  rf <- ranger(dependent.variable.name = "Species", data = dat, num.trees = 1, 
               replace = FALSE, sample.fraction = 1)
  pred_class <- predict(rf, dat)$predictions
  nodes <- predict(rf, dat, type = "terminalNodes")$predictions[, 1]
  ti <- treeInfo(rf, 1)
  pred_ti <- sapply(nodes, function(x) {
    ti[ti$nodeID == x, "prediction"]
  })
  expect_equal(pred_ti, pred_class)
})

test_that("Prediction for classification is same as class prediction, new factor", {
  dat <- iris[sample(nrow(iris)), ]
  dat$Species <- factor(dat$Species, levels = sample(levels(dat$Species)))
  rf <- ranger(dependent.variable.name = "Species", data = dat, num.trees = 1, 
               replace = FALSE, sample.fraction = 1)
  pred_class <- predict(rf, dat)$predictions
  nodes <- predict(rf, dat, type = "terminalNodes")$predictions[, 1]
  ti <- treeInfo(rf, 1)
  pred_ti <- sapply(nodes, function(x) {
    ti[ti$nodeID == x, "prediction"]
  })
  expect_equal(pred_ti, pred_class)
})

test_that("Prediction for classification is same as class prediction, unused factor levels", {
  dat <- iris[c(101:150, 51:100), ]
  expect_warning(rf <- ranger(dependent.variable.name = "Species", data = dat, num.trees = 1, 
    replace = FALSE, sample.fraction = 1))
  pred_class <- predict(rf, dat)$predictions
  nodes <- predict(rf, dat, type = "terminalNodes")$predictions[, 1]
  ti <- treeInfo(rf, 1)
  pred_ti <- sapply(nodes, function(x) {
    ti[ti$nodeID == x, "prediction"]
  })
  expect_equal(pred_ti, pred_class)
})

test_that("Prediction for probability is same as probability prediction", {
  dat <- iris[sample(nrow(iris)), ]
  rf <- ranger(dependent.variable.name = "Species", data = dat, num.trees = 1, 
               sample.fraction = 1, replace = FALSE, probability = TRUE)
  ti <- treeInfo(rf)
  pred_prob <- predict(rf, dat)$predictions
  nodes <- predict(rf, dat, type = "terminalNodes")$predictions[, 1]
  pred_ti <- t(sapply(nodes, function(x) {
    as.matrix(ti[ti$nodeID == x, 8:10])
  }))
  colnames(pred_ti) <- gsub("pred\\.", "", colnames(ti)[8:10])
  expect_equal(pred_prob, pred_ti)
})

test_that("Prediction for probability is same as probability prediction, new factor", {
  dat <- iris[sample(nrow(iris)), ]
  dat$Species <- factor(dat$Species, levels = sample(levels(dat$Species)))
  rf <- ranger(dependent.variable.name = "Species", data = dat, num.trees = 1, 
               sample.fraction = 1, replace = FALSE, probability = TRUE)
  ti <- treeInfo(rf)
  pred_prob <- predict(rf, dat)$predictions
  nodes <- predict(rf, dat, type = "terminalNodes")$predictions[, 1]
  pred_ti <- t(sapply(nodes, function(x) {
    as.matrix(ti[ti$nodeID == x, 8:10])
  }))
  colnames(pred_ti) <- gsub("pred\\.", "", colnames(ti)[8:10])
  expect_equal(pred_prob, pred_ti)
})

test_that("Prediction for probability is same as probability prediction, unused factor levels", {
  dat <- iris[c(101:150, 51:100), ]
  dat$Species <- factor(dat$Species, levels = sample(levels(dat$Species)))
  expect_warning(rf <- ranger(dependent.variable.name = "Species", data = dat, num.trees = 1, 
    sample.fraction = 1, replace = FALSE, probability = TRUE))
  ti <- treeInfo(rf)
  pred_prob <- predict(rf, dat)$predictions
  nodes <- predict(rf, dat, type = "terminalNodes")$predictions[, 1]
  pred_ti <- t(sapply(nodes, function(x) {
    as.matrix(ti[ti$nodeID == x, 8:9])
  }))
  colnames(pred_ti) <- gsub("pred\\.", "", colnames(ti)[8:9])
  expect_equal(pred_prob, pred_ti)
})

test_that("Prediction for matrix classification is integer with correct values", {
  rf <- ranger(dependent.variable.name = "Species", data = data.matrix(iris), 
               num.trees = 5, classification = TRUE)
  ti <- treeInfo(rf, 1)
  expect_is(ti$prediction, "numeric")
  expect_equal(sort(unique(ti$prediction)), 1:3)
})

## Regression
n <- 20
dat <- data.frame(y = rnorm(n), 
                  replicate(2, runif(n)), 
                  replicate(2, rbinom(n, size = 1, prob = .5)))

rf.regr.formula <- ranger(y ~ ., dat, num.trees = 5)
rf.regr.first <- ranger(dependent.variable.name = "y", data = dat, num.trees = 5)
rf.regr.mid <- ranger(dependent.variable.name = "y", data = dat[, c(2:3, 1, 4:5)], num.trees = 5)
rf.regr.last <- ranger(dependent.variable.name = "y", data = dat[, c(2:5, 1)], num.trees = 5)

ti.regr.formula <- treeInfo(rf.regr.formula)
ti.regr.first <- treeInfo(rf.regr.first)
ti.regr.mid <- treeInfo(rf.regr.mid)
ti.regr.last <- treeInfo(rf.regr.last)

test_that("Terminal nodes have only prediction, non-terminal nodes all others, regression formula", {
  expect_true(all(is.na(ti.regr.formula[ti.regr.formula$terminal, 2:6])))
  expect_true(all(!is.na(ti.regr.formula[ti.regr.formula$terminal, c(1, 7:8)])))
  expect_true(all(!is.na(ti.regr.formula[!ti.regr.formula$terminal, -8])))
  expect_true(all(is.na(ti.regr.formula[!ti.regr.formula$terminal, 8])))
})
test_that("Terminal nodes have only prediction, non-terminal nodes all others, regression depvarname first", {
  expect_true(all(is.na(ti.regr.first[ti.regr.first$terminal, 2:6])))
  expect_true(all(!is.na(ti.regr.first[ti.regr.first$terminal, c(1, 7:8)])))
  expect_true(all(!is.na(ti.regr.first[!ti.regr.first$terminal, -8])))
  expect_true(all(is.na(ti.regr.first[!ti.regr.first$terminal, 8])))
})
test_that("Terminal nodes have only prediction, non-terminal nodes all others, regression depvarname mid", { 
  expect_true(all(is.na(ti.regr.mid[ti.regr.mid$terminal, 2:6])))
  expect_true(all(!is.na(ti.regr.mid[ti.regr.mid$terminal, c(1, 7:8)])))
  expect_true(all(!is.na(ti.regr.mid[!ti.regr.mid$terminal, -8])))
  expect_true(all(is.na(ti.regr.mid[!ti.regr.mid$terminal, 8])))
})
test_that("Terminal nodes have only prediction, non-terminal nodes all others, regression depvarname last", { 
  expect_true(all(is.na(ti.regr.last[ti.regr.last$terminal, 2:6])))
  expect_true(all(!is.na(ti.regr.last[ti.regr.last$terminal, c(1, 7:8)])))
  expect_true(all(!is.na(ti.regr.last[!ti.regr.last$terminal, -8])))
  expect_true(all(is.na(ti.regr.last[!ti.regr.last$terminal, 8])))
})

test_that("Names in treeInfo match, regression", {
  varnames <- c("X1", "X2", "X1.1", "X2.1")
  expect_true(all(is.na(ti.regr.formula$splitvarName) | ti.regr.formula$splitvarName %in% varnames))
  expect_true(all(is.na(ti.regr.first$splitvarName) | ti.regr.first$splitvarName %in% varnames))
  expect_true(all(is.na(ti.regr.mid$splitvarName) | ti.regr.mid$splitvarName %in% varnames))
  expect_true(all(is.na(ti.regr.last$splitvarName) | ti.regr.last$splitvarName %in% varnames))
})

test_that("Prediction for regression is numeric in correct range", {
  expect_is(ti.regr.formula$prediction, "numeric")
  expect_true(all(is.na(ti.regr.formula$prediction) | ti.regr.formula$prediction >= min(dat$y)))
  expect_true(all(is.na(ti.regr.formula$prediction) | ti.regr.formula$prediction <= max(dat$y)))
})

## Probability estimation
rf.prob.formula <- ranger(Species ~ ., iris, num.trees = 5, probability = TRUE)
rf.prob.first <- ranger(dependent.variable.name = "Species", data = iris[, c(5, 1:4)], num.trees = 5, probability = TRUE)
rf.prob.mid <- ranger(dependent.variable.name = "Species", data = iris[, c(1:2, 5, 3:4)], num.trees = 5, probability = TRUE)
rf.prob.last <- ranger(dependent.variable.name = "Species", data = iris, num.trees = 5, probability = TRUE)

ti.prob.formula <- treeInfo(rf.prob.formula)
ti.prob.first <- treeInfo(rf.prob.first)
ti.prob.mid <- treeInfo(rf.prob.mid)
ti.prob.last <- treeInfo(rf.prob.last)

test_that("Terminal nodes have only prediction, non-terminal nodes all others, probability formula", {
  expect_true(all(is.na(ti.prob.formula[ti.prob.formula$terminal, 2:6])))
  expect_true(all(!is.na(ti.prob.formula[ti.prob.formula$terminal, c(1, 7:10)])))
  expect_true(all(!is.na(ti.prob.formula[!ti.prob.formula$terminal, c(-8, -9, -10)])))
  expect_true(all(is.na(ti.prob.formula[!ti.prob.formula$terminal, 8:10])))
})
test_that("Terminal nodes have only prediction, non-terminal nodes all others, probability depvarname first", {
  expect_true(all(is.na(ti.prob.first[ti.prob.first$terminal, 2:6])))
  expect_true(all(!is.na(ti.prob.first[ti.prob.first$terminal, c(1, 7:8)])))
  expect_true(all(!is.na(ti.prob.first[!ti.prob.first$terminal, c(-8, -9, -10)])))
  expect_true(all(is.na(ti.prob.first[!ti.prob.first$terminal, 8:10])))
})
test_that("Terminal nodes have only prediction, non-terminal nodes all others, probability depvarname mid", { 
  expect_true(all(is.na(ti.prob.mid[ti.prob.mid$terminal, 2:6])))
  expect_true(all(!is.na(ti.prob.mid[ti.prob.mid$terminal, c(1, 7:8)])))
  expect_true(all(!is.na(ti.prob.mid[!ti.prob.mid$terminal, c(-8, -9, -10)])))
  expect_true(all(is.na(ti.prob.mid[!ti.prob.mid$terminal, 8:10])))
})
test_that("Terminal nodes have only prediction, non-terminal nodes all others, probability depvarname last", { 
  expect_true(all(is.na(ti.prob.last[ti.prob.last$terminal, 2:6])))
  expect_true(all(!is.na(ti.prob.last[ti.prob.last$terminal, c(1, 7:8)])))
  expect_true(all(!is.na(ti.prob.last[!ti.prob.last$terminal, c(-8, -9, -10)])))
  expect_true(all(is.na(ti.prob.last[!ti.prob.last$terminal, 8:10])))
})

test_that("Names in treeInfo match, probability", {
  varnames <- colnames(iris)[1:4]
  expect_true(all(is.na(ti.prob.formula$splitvarName) | ti.prob.formula$splitvarName %in% varnames))
  expect_true(all(is.na(ti.prob.first$splitvarName) | ti.prob.first$splitvarName %in% varnames))
  expect_true(all(is.na(ti.prob.mid$splitvarName) | ti.prob.mid$splitvarName %in% varnames))
  expect_true(all(is.na(ti.prob.last$splitvarName) | ti.prob.last$splitvarName %in% varnames))
})

test_that("Prediction for probability is one probability per class, sum to 1", {
  expect_equal(ncol(ti.prob.formula), 10)
  expect_is(ti.prob.formula$pred.setosa, "numeric")
  expect_true(all(!ti.prob.formula$terminal | rowSums(ti.prob.formula[, 8:10]) == 1))
})

test_that("Prediction for probability has correct factor levels", {
  dat <- iris[c(101:150, 1:100), ]
  rf <- ranger(dependent.variable.name = "Species", data = dat, num.trees = 5, probability = TRUE)
  
  # Predict
  pred_rf <- predict(rf, dat, num.trees = 1)$predictions
  
  # Predict with treeInfo
  ti <- treeInfo(rf)
  terminal_nodes <- predict(rf, dat, type = "terminalNodes")$predictions[, 1]
  pred_ti <- as.matrix(ti[terminal_nodes + 1, grep("pred", colnames(ti))])
  colnames(pred_ti) <- gsub("pred\\.", "", colnames(pred_ti))
  rownames(pred_ti) <- NULL
  
  expect_equal(pred_rf, pred_ti)
})

## Survival
rf.surv.formula <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5)
rf.surv.first <- ranger(dependent.variable.name = "time", status.variable.name = "status", data = veteran[, c(3:4, 1:2, 5:8)], num.trees = 5)
rf.surv.mid <- ranger(dependent.variable.name = "time", status.variable.name = "status", data = veteran, num.trees = 5)
rf.surv.last <- ranger(dependent.variable.name = "time", status.variable.name = "status", data = veteran[, c(2, 1, 5:8, 3:4)], num.trees = 5)

ti.surv.formula <- treeInfo(rf.surv.formula)
ti.surv.first <- treeInfo(rf.surv.first)
ti.surv.mid <- treeInfo(rf.surv.mid)
ti.surv.last <- treeInfo(rf.surv.last)

test_that("Terminal nodes have only nodeID, non-terminal nodes all, survival formula", {
  expect_true(all(is.na(ti.surv.formula[ti.surv.formula$terminal, 2:6])))
  expect_true(all(!is.na(ti.surv.formula[ti.surv.formula$terminal, c(1, 7)])))
  expect_true(all(!is.na(ti.surv.formula[!ti.surv.formula$terminal, ])))
})
test_that("Terminal nodes have only prediction, non-terminal nodes all others, survival depvarname first", {
  expect_true(all(is.na(ti.surv.first[ti.surv.first$terminal, 2:6])))
  expect_true(all(!is.na(ti.surv.first[ti.surv.first$terminal, c(1, 7)])))
  expect_true(all(!is.na(ti.surv.first[!ti.surv.first$terminal, ])))
})
test_that("Terminal nodes have only prediction, non-terminal nodes all others, survival depvarname mid", { 
  expect_true(all(is.na(ti.surv.mid[ti.surv.mid$terminal, 2:6])))
  expect_true(all(!is.na(ti.surv.mid[ti.surv.mid$terminal, c(1, 7)])))
  expect_true(all(!is.na(ti.surv.mid[!ti.surv.mid$terminal, ])))
})
test_that("Terminal nodes have only prediction, non-terminal nodes all others, survival depvarname last", { 
  expect_true(all(is.na(ti.surv.last[ti.surv.last$terminal, 2:6])))
  expect_true(all(!is.na(ti.surv.last[ti.surv.last$terminal, c(1, 7)])))
  expect_true(all(!is.na(ti.surv.last[!ti.surv.last$terminal, ])))
})

test_that("Names in treeInfo match, survival", {
  varnames <- colnames(veteran)[c(1:2, 5:8)]
  expect_true(all(is.na(ti.surv.formula$splitvarName) | ti.surv.formula$splitvarName %in% varnames))
  expect_true(all(is.na(ti.surv.first$splitvarName) | ti.surv.first$splitvarName %in% varnames))
  expect_true(all(is.na(ti.surv.mid$splitvarName) | ti.surv.mid$splitvarName %in% varnames))
  expect_true(all(is.na(ti.surv.last$splitvarName) | ti.surv.last$splitvarName %in% varnames))
})

test_that("No prediction for Survival", {
  expect_equal(ncol(ti.surv.formula), 7)
})

## General
test_that("Error if no saved forest", {
  expect_error(treeInfo(ranger(Species ~ ., iris, write.forest = FALSE)), 
               "Error\\: No saved forest in ranger object\\. Please set write.forest to TRUE when calling ranger\\.")
})

## Unordered splitting
test_that("Spitting value is comma separated list for partition splitting", {
  n <- 50
  dat <- data.frame(x = sample(c("A", "B", "C", "D", "E"), n, replace = TRUE), 
                    y = rbinom(n, 1, 0.5), 
                    stringsAsFactors = FALSE)
  rf.partition <- ranger(y ~ ., dat, num.trees = 5, respect.unordered.factors = "partition")
  ti.partition <- treeInfo(rf.partition)
  
  expect_is(ti.partition$splitval, "character")
  expect_true(all(is.na(ti.partition$splitval) | grepl("^\\d+(?:,\\d+)*$", ti.partition$splitval)))
})

test_that("Spitting value is numeric for order splitting", {
  set.seed(100)
  rf.order <- ranger(Sepal.Length ~ ., iris, num.trees = 5, respect.unordered.factors = "order")
  ti.order <- treeInfo(rf.order)
  expect_is(ti.order$splitval[!ti.order$terminal & ti.order$splitvarName == "Species"], "numeric")
})

test_that("treeInfo works for 31 unordered factor levels but not for 32", {
  n <- 31
  dt <- data.frame(x = factor(1:n, ordered = FALSE),  
                   y = rbinom(n, 1, 0.5))
  rf <- ranger(y ~ ., data = dt, num.trees = 10, splitrule = "extratrees")
  expect_silent(treeInfo(rf))
  
  n <- 32
  dt <- data.frame(x = factor(1:n, ordered = FALSE),  
                   y = rbinom(n, 1, 0.5))
  rf <- ranger(y ~ ., data = dt, num.trees = 10, splitrule = "extratrees")
  expect_warning(treeInfo(rf), "Unordered splitting levels can only be shown for up to 31 levels.")
})



