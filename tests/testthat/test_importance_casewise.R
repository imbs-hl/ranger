context("test_casewise_importances")

test_that("casewise importance works, classification", {
  data <- 
    data.frame(
      x = runif(1000),
      y = rnorm(1000, mean = 1),
      z = rnorm(1000, mean = 2)
    )
  rownames(data) <- paste0("case_", seq_len(nrow(data)))
  data$a <- factor(ifelse(ifelse(data$x < .5, data$y, data$z) > 1.5, "left", "right"))
  
  rf <- ranger(
    data = data,
    dependent.variable.name = "a",
    importance = "permutation",
    local.importance = TRUE,
    num.threads = 1
  )
  vic <- rf$variable.importance.casewise
  
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
  data <- 
    data.frame(
      x = runif(1000),
      y = rnorm(1000, mean = 1),
      z = rnorm(1000, mean = 2)
    )
  rownames(data) <- paste0("case_", seq_len(nrow(data)))
  data$a <- ifelse(data$x < .5, data$y, data$z)
  
  rf <- ranger(
    data = data,
    dependent.variable.name = "a",
    importance = "permutation",
    local.importance = TRUE,
    num.threads = 1
  )
  vic <- rf$variable.importance.casewise
  
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
  data <- 
    data.frame(
      x = runif(1000),
      y = rnorm(1000, mean = 1),
      z = rnorm(1000, mean = 2)
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
    num.threads = 1
  )
  vic <- rf$variable.importance.casewise
  
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
  data <- 
    data.frame(
      x = runif(1000),
      y = rnorm(1000, mean = 1),
      z = rnorm(1000, mean = 2),
      surv = sample(0:1, 1000, replace = TRUE)
    )
  rownames(data) <- paste0("case_", seq_len(nrow(data)))
  data$a <- ifelse(data$x < .5, data$y, data$z) * ifelse(data$surv, runif(1000, min = .8, max = 1), 1)
  
  rf <- ranger(
    data = data,
    dependent.variable.name = "a",
    status.variable.name = "surv",
    importance = "permutation",
    local.importance = TRUE,
    num.threads = 1
  )
  vic <- rf$variable.importance.casewise
  
  # should see clear pattern here:
  # pheatmap::pheatmap(vic[order(data$x),], cluster_cols = FALSE, cluster_rows = FALSE)
  
  expect_equal(rownames(vic), rownames(data))
  expect_equal(colnames(vic), colnames(data)[1:3])
  
  expect_lte(wilcox.test(vic[data$x < .5, 2], vic[data$x >= .5, 2], "greater")$p.value, .01)
  expect_gte(wilcox.test(vic[data$x < .5, 2], vic[data$x >= .5, 2], "less")$p.value, .99)
  expect_lte(wilcox.test(vic[data$x < .5, 3], vic[data$x >= .5, 3], "less")$p.value, .01)
  expect_gte(wilcox.test(vic[data$x < .5, 3], vic[data$x >= .5, 3], "greater")$p.value, .99)
})