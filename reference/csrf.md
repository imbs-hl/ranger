# Case-specific random forests.

In case-specific random forests (CSRF), random forests are built
specific to the cases of interest. Instead of using equal probabilities,
the cases are weighted according to their difference to the case of
interest.

## Usage

``` r
csrf(
  formula,
  training_data,
  test_data,
  params1 = list(),
  params2 = list(),
  verbose = FALSE
)
```

## Arguments

- formula:

  Object of class `formula` or `character` describing the model to fit.

- training_data:

  Training data of class `data.frame`.

- test_data:

  Test data of class `data.frame`.

- params1:

  Parameters for the proximity random forest grown in the first step.

- params2:

  Parameters for the prediction random forests grown in the second step.

- verbose:

  Logical indicating whether or not to print computation progress.

## Value

Predictions for the test dataset.

## Details

The algorithm consists of 3 steps:

1.  Grow a random forest on the training data

2.  For each observation of interest (test data), the weights of all
    training observations are computed by counting the number of trees
    in which both observations are in the same terminal node.

3.  For each test observation, grow a weighted random forest on the
    training data, using the weights obtained in step 2. Predict the
    outcome of the test observation as usual.

In total, n+1 random forests are grown, where n is the number
observations in the test dataset. For details, see Xu et al. (2014).

## References

Xu, R., Nettleton, D. & Nordman, D.J. (2014). Case-specific random
forests. J Comp Graph Stat 25:49-65.
[doi:10.1080/10618600.2014.983641](https://doi.org/10.1080/10618600.2014.983641)
.

## Author

Marvin N. Wright

## Examples

``` r
## Split in training and test data
train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
iris.train <- iris[train.idx, ]
iris.test <- iris[-train.idx, ]

## Run case-specific RF
csrf(Species ~ ., training_data = iris.train, test_data = iris.test, 
     params1 = list(num.trees = 50, mtry = 4), 
     params2 = list(num.trees = 5))
#>  [1] setosa     setosa     setosa     setosa     setosa     setosa    
#>  [7] setosa     setosa     setosa     setosa     setosa     versicolor
#> [13] versicolor versicolor versicolor versicolor versicolor virginica 
#> [19] versicolor versicolor versicolor versicolor versicolor versicolor
#> [25] versicolor versicolor versicolor versicolor versicolor versicolor
#> [31] versicolor virginica  virginica  virginica  virginica  virginica 
#> [37] virginica  virginica  versicolor virginica  virginica  virginica 
#> [43] versicolor virginica  virginica  virginica  virginica  virginica 
#> [49] virginica  virginica 
#> Levels: setosa versicolor virginica
```
