# Ranger prediction

Prediction with new data and a saved forest from Ranger.

## Usage

``` r
# S3 method for class 'ranger'
predict(
  object,
  data = NULL,
  predict.all = FALSE,
  num.trees = object$num.trees,
  type = "response",
  se.method = "infjack",
  quantiles = c(0.1, 0.5, 0.9),
  what = NULL,
  seed = NULL,
  num.threads = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  Ranger `ranger` object.

- data:

  New test data of class `data.frame` or `gwaa.data` (GenABEL).

- predict.all:

  Return individual predictions for each tree instead of aggregated
  predictions for all trees. Return a matrix (sample x tree) for
  classification and regression, a 3d array for probability estimation
  (sample x class x tree) and survival (sample x time x tree).

- num.trees:

  Number of trees used for prediction. The first `num.trees` in the
  forest are used.

- type:

  Type of prediction. One of 'response', 'se', 'terminalNodes',
  'quantiles' with default 'response'. See below for details.

- se.method:

  Method to compute standard errors. One of 'jack', 'infjack' with
  default 'infjack'. Only applicable if type = 'se'. See below for
  details.

- quantiles:

  Vector of quantiles for quantile prediction. Set `type = 'quantiles'`
  to use.

- what:

  User specified function for quantile prediction used instead of
  `quantile`. Must return numeric vector, see examples.

- seed:

  Random seed. Default is `NULL`, which generates the seed from `R`. Set
  to `0` to ignore the `R` seed. The seed is used in case of ties in
  classification mode.

- num.threads:

  Number of threads. Use 0 for all available cores. Default is 2 if not
  set by options/environment variables (see below).

- verbose:

  Verbose output on or off.

- ...:

  further arguments passed to or from other methods.

## Value

Object of class `ranger.prediction` with elements

|                             |                                                                           |
|-----------------------------|---------------------------------------------------------------------------|
| `predictions`               | Predicted classes/values (only for classification and regression)         |
| `unique.death.times`        | Unique death times (only for survival).                                   |
| `chf`                       | Estimated cumulative hazard function for each sample (only for survival). |
| `survival`                  | Estimated survival function for each sample (only for survival).          |
| `num.trees`                 | Number of trees.                                                          |
| `num.independent.variables` | Number of independent variables.                                          |
| `treetype`                  | Type of forest/tree. Classification, regression or survival.              |
| `num.samples`               | Number of samples.                                                        |

## Details

For `type = 'response'` (the default), the predicted classes
(classification), predicted numeric values (regression), predicted
probabilities (probability estimation) or survival probabilities
(survival) are returned. For `type = 'se'`, the standard error of the
predictions are returned (regression only). The
jackknife-after-bootstrap or infinitesimal jackknife for bagging is used
to estimate the standard errors based on out-of-bag predictions. See
Wager et al. (2014) for details. For `type = 'terminalNodes'`, the IDs
of the terminal node in each tree for each observation in the given
dataset are returned. For `type = 'quantiles'`, the selected quantiles
for each observation are estimated. See Meinshausen (2006) for details.

If `type = 'se'` is selected, the method to estimate the variances can
be chosen with `se.method`. Set `se.method = 'jack'` for
jackknife-after-bootstrap and `se.method = 'infjack'` for the
infinitesimal jackknife for bagging.

For classification and `predict.all = TRUE`, a factor levels are
returned as numerics. To retrieve the corresponding factor levels, use
`rf$forest$levels`, if `rf` is the ranger object.

By default, ranger uses 2 threads. The default can be changed with: (1)
`num.threads` in ranger/predict call, (2) environment variable
R_RANGER_NUM_THREADS, (3) `options(ranger.num.threads = N)`, (4)
`options(Ncpus = N)`, with precedence in that order.

## References

- Wright, M. N. & Ziegler, A. (2017). ranger: A Fast Implementation of
  Random Forests for High Dimensional Data in C++ and R. J Stat Softw
  77:1-17.
  [doi:10.18637/jss.v077.i01](https://doi.org/10.18637/jss.v077.i01) .

- Wager, S., Hastie T., & Efron, B. (2014). Confidence Intervals for
  Random Forests: The Jackknife and the Infinitesimal Jackknife. J Mach
  Learn Res 15:1625-1651. <https://jmlr.org/papers/v15/wager14a.html>.

- Meinshausen (2006). Quantile Regression Forests. J Mach Learn Res
  7:983-999. <https://www.jmlr.org/papers/v7/meinshausen06a.html>.

## See also

[`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md)

## Author

Marvin N. Wright

## Examples

``` r
## Classification forest
ranger(Species ~ ., data = iris)
#> Ranger result
#> 
#> Call:
#>  ranger(Species ~ ., data = iris) 
#> 
#> Type:                             Classification 
#> Number of trees:                  500 
#> Sample size:                      150 
#> Number of independent variables:  4 
#> Mtry:                             2 
#> Target node size:                 1 
#> Variable importance mode:         none 
#> Splitrule:                        gini 
#> OOB prediction error:             4.00 % 
train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
iris.train <- iris[train.idx, ]
iris.test <- iris[-train.idx, ]
rg.iris <- ranger(Species ~ ., data = iris.train)
pred.iris <- predict(rg.iris, data = iris.test)
table(iris.test$Species, pred.iris$predictions)
#>             
#>              setosa versicolor virginica
#>   setosa         15          0         0
#>   versicolor      0         16         2
#>   virginica       0          2        15

## Quantile regression forest
rf <- ranger(mpg ~ ., mtcars[1:26, ], quantreg = TRUE)
pred <- predict(rf, mtcars[27:32, ], type = "quantiles", quantiles = c(0.1, 0.5, 0.9))
pred$predictions
#>      quantile= 0.1 quantile= 0.5 quantile= 0.9
#> [1,]          21.0          24.4          32.4
#> [2,]          21.0          22.8          32.4
#> [3,]          13.3          17.8          30.4
#> [4,]          15.5          21.0          22.8
#> [5,]          13.3          14.3          19.2
#> [6,]          21.0          22.8          32.4

## Quantile regression forest with user-specified function
rf <- ranger(mpg ~ ., mtcars[1:26, ], quantreg = TRUE)
pred <- predict(rf, mtcars[27:32, ], type = "quantiles", 
                what = function(x) sample(x, 10, replace = TRUE))
pred$predictions
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,] 21.0 33.9 27.3 27.3 32.4 33.9 22.8 19.2 14.3  27.3
#> [2,] 33.9 21.5 21.5 22.8 21.5 22.8 21.0 22.8 21.0  27.3
#> [3,] 32.4 15.2 21.0 30.4 32.4 13.3 21.0 15.2 10.4  14.3
#> [4,] 13.3 22.8 21.0 17.8 19.2 21.0 21.0 19.2 22.8  21.0
#> [5,] 10.4 13.3 10.4 13.3 15.2 14.3 13.3 13.3 13.3  19.2
#> [6,] 19.2 21.5 21.4 22.8 27.3 22.8 22.8 24.4 17.8  27.3
```
