
context("ranger_cpp_arguments")

test_that("Error if sample fraction is 0 or >1", {
  expect_warning(
    expect_error(ranger_cpp(data = iris, depvarname = "Species", ntree = 5, fraction = 0), 
      "Error: Illegal argument for option 'fraction'\\. Please give a value in \\(0,1\\]\\. See '--help' for details\\. Ranger will EXIT now\\."))
  expect_warning(
    expect_error(ranger_cpp(data = iris, depvarname = "Species", ntree = 5, fraction = 1.1), 
      "Error: Illegal argument for option 'fraction'\\. Please give a value in \\(0,1\\]\\. See '--help' for details\\. Ranger will EXIT now\\."))
})
