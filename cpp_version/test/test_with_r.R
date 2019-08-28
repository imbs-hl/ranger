
library(data.table)
library(testthat)

# Function to call C++ version from R
ranger_cpp <- function(data, ...) {
  if (is.data.frame(data) && any(sapply(data, is.numeric))) {
    idx_numeric <- sapply(data, is.numeric)
    data[, !idx_numeric] <- lapply(data[, !idx_numeric, drop = FALSE], as.numeric)
  }
  fwrite(data, "temp_data.csv")
  ret <- system2("../build/ranger", 
                 args = c("--verbose", "--file temp_data.csv", paste0("--", names(list(...)), " ", list(...))), 
                 stdout = TRUE, stderr = TRUE)
  if (length(ret) == 1 &&  nchar(ret) >= 5 && substr(ret, 1, 5) == "Error") {
    stop(ret)
  }
  #unlink("temp_data.csv")
  ret
}

test_that("Prediction is equal to R version", {
  # C++ version
  rf <- ranger_cpp(data = iris, depvarname = "Species", ntree = 5, write = "", seed = 10)
  pred <- ranger_cpp(data = iris, depvarname = "Species", predict = "ranger_out.forest", seed = 20)
  preds_cpp <- as.data.frame(fread("ranger_out.prediction"))[, 1]
  
  # R version
  rf <- ranger(Species ~ ., iris, num.trees = 5, seed = 10)
  preds_r <- as.numeric(predict(rf, iris, seed = 20)$predictions)
  
  expect_equal(preds_cpp, preds_r)
})

test_that("Error if sample fraction is 0 or >1", {
  expect_error(ranger_cpp(data = iris, depvarname = "Species", ntree = 5, fraction = 0), 
               "Error: Illegal argument for option 'fraction'\\. Please give a value in \\(0,1\\]\\. See '--help' for details\\. Ranger will EXIT now\\.")
  expect_error(ranger_cpp(data = iris, depvarname = "Species", ntree = 5, fraction = 1.1), 
               "Error: Illegal argument for option 'fraction'\\. Please give a value in \\(0,1\\]\\. See '--help' for details\\. Ranger will EXIT now\\.")
})

test_that("Same result with beta splitting", {
  # Generate data with 0..1 outcome
  n <- 100
  p <- 4
  beta <- c(0, 1, 2, 3)
  x <- round(replicate(p, runif(n)), 3)
  y <- as.vector(x %*% beta)
  y <- round((y-min(y))/(max(y)-min(y)), 3)
  dat <- data.frame(y = y, x)
  
  # C++ version
  rf <- ranger_cpp(data = dat, depvarname = "y", treetype = 3, ntree = 5, splitrule = 6, seed = 10)
  err_cpp <- grep("Overall OOB prediction error", rf, value = TRUE)
  err_cpp <- as.numeric(gsub("[^0-9.]", "", err_cpp))
  
  # R version
  rf <- ranger(y ~ ., dat, num.trees = 5, splitrule = "beta", seed = 10)
  err_r <- rf$prediction.error

  expect_equal(round(err_cpp, 4), round(err_r, 4))
})
