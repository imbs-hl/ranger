# Ranger prediction

Prediction with new data and a saved forest from Ranger.

## Usage

``` r
# S3 method for class 'ranger.forest'
predict(
  object,
  data,
  predict.all = FALSE,
  num.trees = object$num.trees,
  type = "response",
  se.method = "infjack",
  seed = NULL,
  num.threads = NULL,
  verbose = TRUE,
  inbag.counts = NULL,
  ...
)
```

## Arguments

- object:

  Ranger `ranger.forest` object.

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

- seed:

  Random seed. Default is `NULL`, which generates the seed from `R`. Set
  to `0` to ignore the `R` seed. The seed is used in case of ties in
  classification mode.

- num.threads:

  Number of threads. Use 0 for all available cores. Default is 2 if not
  set by options/environment variables (see below).

- verbose:

  Verbose output on or off.

- inbag.counts:

  Number of times the observations are in-bag in the trees.

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
dataset are returned.

If `type = 'se'` is selected, the method to estimate the variances can
be chosen with `se.method`. Set `se.method = 'jack'` for jackknife after
bootstrap and `se.method = 'infjack'` for the infinitesimal jackknife
for bagging.

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

## See also

[`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md)

## Author

Marvin N. Wright
