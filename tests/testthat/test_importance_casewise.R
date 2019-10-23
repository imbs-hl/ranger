context("test_casewise_importances")

test_that("casewise importance works, classification", {
  n <- 1000
  data <- 
    data.frame(
      x = round(runif(n), 1),
      y = round(rnorm(n, mean = 1), 1),
      z = round(rnorm(n, mean = 2), 1)
    )
  rownames(data) <- paste0("case_", seq_len(nrow(data)))
  data$a <- factor(ifelse(ifelse(data$x < .5, data$y, data$z) > 1.5, "left", "right"))
  
  rf <- ranger(
    data = data,
    dependent.variable.name = "a",
    importance = "permutation",
    local.importance = TRUE,
    num.trees = 5
  )
  vic <- rf$variable.importance.local
  
  # should see clear pattern here:
  # pheatmap::pheatmap(vic[order(data$x),], cluster_cols = FALSE, cluster_rows = FALSE)
  
  expect_equal(rownames(vic), rownames(data))
  expect_equal(colnames(vic), colnames(data)[1:3])
  
  expect_lte(wilcox.test(vic[data$x < .5, 2], vic[data$x >= .5, 2], "greater")$p.value, .01)
  expect_gte(wilcox.test(vic[data$x < .5, 2], vic[data$x >= .5, 2], "less")$p.value, .99)
  expect_lte(wilcox.test(vic[data$x < .5, 3], vic[data$x >= .5, 3], "less")$p.value, .01)
  expect_gte(wilcox.test(vic[data$x < .5, 3], vic[data$x >= .5, 3], "greater")$p.value, .99)
})

test_that("casewise importance works, regression", {
  n <- 1000
  data <- 
    data.frame(
      x = round(runif(n), 1),
      y = round(rnorm(n, mean = 1), 1),
      z = round(rnorm(n, mean = 2), 1)
    )
  rownames(data) <- paste0("case_", seq_len(nrow(data)))
  data$a <- ifelse(data$x < .5, data$y, data$z)
  
  rf <- ranger(
    data = data,
    dependent.variable.name = "a",
    importance = "permutation",
    local.importance = TRUE,
    num.trees = 5
  )
  vic <- rf$variable.importance.local
  
  # should see clear pattern here:
  # pheatmap::pheatmap(vic[order(data$x),], cluster_cols = FALSE, cluster_rows = FALSE)
  
  expect_equal(rownames(vic), rownames(data))
  expect_equal(colnames(vic), colnames(data)[1:3])
  
  expect_lte(wilcox.test(vic[data$x < .5, 2], vic[data$x >= .5, 2], "greater")$p.value, .01)
  expect_gte(wilcox.test(vic[data$x < .5, 2], vic[data$x >= .5, 2], "less")$p.value, .99)
  expect_lte(wilcox.test(vic[data$x < .5, 3], vic[data$x >= .5, 3], "less")$p.value, .01)
  expect_gte(wilcox.test(vic[data$x < .5, 3], vic[data$x >= .5, 3], "greater")$p.value, .99)
})

test_that("casewise importance works, probability", {
  n <- 1000
  data <- 
    data.frame(
      x = round(runif(n), 1),
      y = round(rnorm(n, mean = 1), 1),
      z = round(rnorm(n, mean = 2), 1)
    )
  rownames(data) <- paste0("case_", seq_len(nrow(data)))
  # data$a <- ifelse(data$x < .5, data$y, data$z)
  data$a <- factor(ifelse(ifelse(data$x < .5, data$y, data$z) > 1.5, "left", "right"))
  
  rf <- ranger(
    data = data,
    dependent.variable.name = "a",
    importance = "permutation",
    probability = TRUE,
    local.importance = TRUE,
    num.trees = 5
  )
  vic <- rf$variable.importance.local
  
  # should see clear pattern here:
  # pheatmap::pheatmap(vic[order(data$x),], cluster_cols = FALSE, cluster_rows = FALSE)
  
  expect_equal(rownames(vic), rownames(data))
  expect_equal(colnames(vic), colnames(data)[1:3])
  
  expect_lte(wilcox.test(vic[data$x < .5, 2], vic[data$x >= .5, 2], "greater")$p.value, .01)
  expect_gte(wilcox.test(vic[data$x < .5, 2], vic[data$x >= .5, 2], "less")$p.value, .99)
  expect_lte(wilcox.test(vic[data$x < .5, 3], vic[data$x >= .5, 3], "less")$p.value, .01)
  expect_gte(wilcox.test(vic[data$x < .5, 3], vic[data$x >= .5, 3], "greater")$p.value, .99)
})

test_that("casewise importance works, survival", {
  n <- 1000
  data <- 
    data.frame(
      x = round(runif(n), 1),
      y = round(rnorm(n, mean = 1), 1),
      z = round(rnorm(n, mean = 2), 1),
      surv = rbinom(n, 1, .8)
    )
  rownames(data) <- paste0("case_", seq_len(nrow(data)))
  data$a <- (ifelse(data$x < .5, data$y, data$z))
  
  rf <- ranger(
    data = data,
    dependent.variable.name = "a",
    status.variable.name = "surv",
    importance = "permutation",
    local.importance = TRUE,
    num.trees = 5
  )
  vic <- rf$variable.importance.local
  
  # should see clear pattern here:
  # pheatmap::pheatmap(vic[order(data$x),], cluster_cols = FALSE, cluster_rows = FALSE)
  
  expect_equal(rownames(vic), rownames(data))
  expect_equal(colnames(vic), colnames(data)[1:3])
  
  expect_lte(wilcox.test(vic[data$x < .5, 2], vic[data$x >= .5, 2], "greater")$p.value, .1)
  expect_gte(wilcox.test(vic[data$x < .5, 2], vic[data$x >= .5, 2], "less")$p.value, .9)
  expect_lte(wilcox.test(vic[data$x < .5, 3], vic[data$x >= .5, 3], "less")$p.value, .1)
  expect_gte(wilcox.test(vic[data$x < .5, 3], vic[data$x >= .5, 3], "greater")$p.value, .9)
})
