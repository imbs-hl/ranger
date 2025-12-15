# Ranger

Ranger is a fast implementation of random forests (Breiman 2001) or
recursive partitioning, particularly suited for high dimensional data.
Classification, regression, and survival forests are supported.
Classification and regression forests are implemented as in the original
Random Forest (Breiman 2001), survival forests as in Random Survival
Forests (Ishwaran et al. 2008). Includes implementations of extremely
randomized trees (Geurts et al. 2006) and quantile regression forests
(Meinshausen 2006).

## Usage

``` r
ranger(
  formula = NULL,
  data = NULL,
  num.trees = 500,
  mtry = NULL,
  importance = "none",
  write.forest = TRUE,
  probability = FALSE,
  min.node.size = NULL,
  min.bucket = NULL,
  max.depth = NULL,
  replace = TRUE,
  sample.fraction = ifelse(replace, 1, 0.632),
  case.weights = NULL,
  class.weights = NULL,
  splitrule = NULL,
  num.random.splits = 1,
  alpha = 0.5,
  minprop = 0.1,
  poisson.tau = 1,
  split.select.weights = NULL,
  always.split.variables = NULL,
  respect.unordered.factors = NULL,
  scale.permutation.importance = FALSE,
  local.importance = FALSE,
  regularization.factor = 1,
  regularization.usedepth = FALSE,
  keep.inbag = FALSE,
  inbag = NULL,
  holdout = FALSE,
  quantreg = FALSE,
  time.interest = NULL,
  oob.error = TRUE,
  num.threads = NULL,
  save.memory = FALSE,
  verbose = TRUE,
  node.stats = FALSE,
  seed = NULL,
  na.action = "na.learn",
  dependent.variable.name = NULL,
  status.variable.name = NULL,
  classification = NULL,
  x = NULL,
  y = NULL,
  ...
)
```

## Arguments

- formula:

  Object of class `formula` or `character` describing the model to fit.
  Interaction terms supported only for numerical variables.

- data:

  Training data of class `data.frame`, `matrix`, `dgCMatrix` (Matrix) or
  `gwaa.data` (GenABEL).

- num.trees:

  Number of trees.

- mtry:

  Number of variables to possibly split at in each node. Default is the
  (rounded down) square root of the number variables. Alternatively, a
  single argument function returning an integer, given the number of
  independent variables.

- importance:

  Variable importance mode, one of 'none', 'impurity',
  'impurity_corrected', 'permutation'. The 'impurity' measure is the
  Gini index for classification, the variance of the responses for
  regression and the sum of test statistics (see `splitrule`) for
  survival.

- write.forest:

  Save `ranger.forest` object, required for prediction. Set to `FALSE`
  to reduce memory usage if no prediction intended.

- probability:

  Grow a probability forest as in Malley et al. (2012).

- min.node.size:

  Minimal node size to split at. Default 1 for classification, 5 for
  regression, 3 for survival, and 10 for probability. For
  classification, this can be a vector of class-specific values.

- min.bucket:

  Minimal terminal node size. No nodes smaller than this value can
  occur. Default 3 for survival and 1 for all other tree types. For
  classification, this can be a vector of class-specific values.

- max.depth:

  Maximal tree depth. A value of NULL or 0 (the default) corresponds to
  unlimited depth, 1 to tree stumps (1 split per tree).

- replace:

  Sample with replacement.

- sample.fraction:

  Fraction of observations to sample. Default is 1 for sampling with
  replacement and 0.632 for sampling without replacement. For
  classification, this can be a vector of class-specific values.

- case.weights:

  Weights for sampling of training observations. Observations with
  larger weights will be selected with higher probability in the
  bootstrap (or subsampled) samples for the trees.

- class.weights:

  Weights for the outcome classes (in order of the factor levels) in the
  splitting rule (cost sensitive learning). Classification and
  probability prediction only. For classification the weights are also
  applied in the majority vote in terminal nodes.

- splitrule:

  Splitting rule. For classification and probability estimation "gini",
  "extratrees" or "hellinger" with default "gini". For regression
  "variance", "extratrees", "maxstat", "beta" or "poisson" with default
  "variance". For survival "logrank", "extratrees", "C" or "maxstat"
  with default "logrank".

- num.random.splits:

  For "extratrees" splitrule.: Number of random splits to consider for
  each candidate splitting variable.

- alpha:

  For "maxstat" splitrule: Significance threshold to allow splitting.

- minprop:

  For "maxstat" splitrule: Lower quantile of covariate distribution to
  be considered for splitting.

- poisson.tau:

  For "poisson" splitrule: The coefficient of variation of the
  (expected) frequency is \\1/\tau\\. If a terminal node has only 0
  responses, the estimate is set to \\\alpha 0 + (1-\alpha)
  mean(parent)\\ with \\\alpha = samples(child) mean(parent) / (\tau +
  samples(child) mean(parent))\\.

- split.select.weights:

  Numeric vector with weights between 0 and 1, used to calculate the
  probability to select variables for splitting. Alternatively, a list
  of size num.trees, containing split select weight vectors for each
  tree can be used.

- always.split.variables:

  Character vector with variable names to be always selected in addition
  to the `mtry` variables tried for splitting.

- respect.unordered.factors:

  Handling of unordered factor covariates. One of 'ignore', 'order' and
  'partition'. For the "extratrees" splitrule the default is "partition"
  for all other splitrules 'ignore'. Alternatively TRUE (='order') or
  FALSE (='ignore') can be used. See below for details.

- scale.permutation.importance:

  Scale permutation importance by standard error as in (Breiman 2001).
  Only applicable if permutation variable importance mode selected.

- local.importance:

  Calculate and return local importance values as in (Breiman 2001).
  Only applicable if `importance` is set to 'permutation'.

- regularization.factor:

  Regularization factor (gain penalization), either a vector of length p
  or one value for all variables.

- regularization.usedepth:

  Consider the depth in regularization.

- keep.inbag:

  Save how often observations are in-bag in each tree.

- inbag:

  Manually set observations per tree. List of size num.trees, containing
  inbag counts for each observation. Can be used for stratified
  sampling.

- holdout:

  Hold-out mode. Hold-out all samples with case weight 0 and use these
  for variable importance and prediction error.

- quantreg:

  Prepare quantile prediction as in quantile regression forests
  (Meinshausen 2006). Regression only. Set `keep.inbag = TRUE` to
  prepare out-of-bag quantile prediction.

- time.interest:

  Time points of interest (survival only). Can be `NULL` (default, use
  all observed time points), a vector of time points or a single number
  to use as many time points (grid over observed time points).

- oob.error:

  Compute OOB prediction error. Set to `FALSE` to save computation time,
  e.g. for large survival forests.

- num.threads:

  Number of threads. Use 0 for all available cores. Default is 2 if not
  set by options/environment variables (see below).

- save.memory:

  Use memory saving (but slower) splitting mode. No effect for survival
  and GWAS data. Warning: This option slows down the tree growing, use
  only if you encounter memory problems.

- verbose:

  Show computation status and estimated runtime.

- node.stats:

  Save node statistics. Set to `TRUE` to save prediction, number of
  observations and split statistics for each node.

- seed:

  Random seed. Default is `NULL`, which generates the seed from `R`. Set
  to `0` to ignore the `R` seed.

- na.action:

  Handling of missing values. Set to "na.learn" to internally handle
  missing values (default, see below), to "na.omit" to omit observations
  with missing values and to "na.fail" to stop if missing values are
  found.

- dependent.variable.name:

  Name of dependent variable, needed if no formula given. For survival
  forests this is the time variable.

- status.variable.name:

  Name of status variable, only applicable to survival data and needed
  if no formula given. Use 1 for event and 0 for censoring.

- classification:

  Set to `TRUE` to grow a classification forest. Only needed if the data
  is a matrix or the response numeric.

- x:

  Predictor data (independent variables), alternative interface to data
  with formula or dependent.variable.name.

- y:

  Response vector (dependent variable), alternative interface to data
  with formula or dependent.variable.name. For survival use a
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) object or a
  matrix with time and status.

- ...:

  Further arguments passed to or from other methods (currently ignored).

## Value

Object of class `ranger` with elements

- `forest`:

  Saved forest (If write.forest set to TRUE). Note that the variable IDs
  in the `split.varIDs` object do not necessarily represent the column
  number in R.

- `predictions`:

  Predicted classes/values, based on out-of-bag samples (classification
  and regression only).

- `variable.importance`:

  Variable importance for each independent variable.

- `variable.importance.local`:

  Variable importance for each independent variable and each sample, if
  `local.importance` is set to TRUE and `importance` is set to
  'permutation'.

- `prediction.error`:

  Overall out-of-bag prediction error. For classification this is
  accuracy (proportion of misclassified observations), for probability
  estimation the Brier score, for regression the mean squared error and
  for survival one minus Harrell's C-index.

- `r.squared`:

  R squared. Also called explained variance or coefficient of
  determination (regression only). Computed on out-of-bag data.

- `confusion.matrix`:

  Contingency table for classes and predictions based on out-of-bag
  samples (classification only).

- `unique.death.times`:

  Unique death times (survival only).

- `chf`:

  Estimated cumulative hazard function for each sample (survival only).

- `survival`:

  Estimated survival function for each sample (survival only).

- `call`:

  Function call.

- `num.trees`:

  Number of trees.

- `num.independent.variables`:

  Number of independent variables.

- `mtry`:

  Value of mtry used.

- `min.node.size`:

  Value of minimal node size used.

- `treetype`:

  Type of forest/tree. classification, regression or survival.

- `importance.mode`:

  Importance mode used.

- `num.samples`:

  Number of samples.

- `inbag.counts`:

  Number of times the observations are in-bag in the trees.

- `dependent.variable.name`:

  Name of the dependent variable. This is NULL when x/y interface is
  used.

- `status.variable.name`:

  Name of the status variable (survival only). This is NULL when x/y
  interface is used.

## Details

The tree type is determined by the type of the dependent variable. For
factors classification trees are grown, for numeric values regression
trees and for survival objects survival trees. The Gini index is used as
default splitting rule for classification. For regression, the estimated
response variances or maximally selected rank statistics (Wright et al.
2016) can be used. For Survival the log-rank test, a C-index based
splitting rule (Schmid et al. 2015) and maximally selected rank
statistics (Wright et al. 2016) are available. For all tree types,
forests of extremely randomized trees (Geurts et al. 2006) can be grown.

With the `probability` option and factor dependent variable a
probability forest is grown. Here, the node impurity is used for
splitting, as in classification forests. Predictions are class
probabilities for each sample. In contrast to other implementations,
each tree returns a probability estimate and these estimates are
averaged for the forest probability estimate. For details see Malley et
al. (2012).

Note that nodes with size smaller than `min.node.size` can occur because
`min.node.size` is the minimal node size *to split at*, as in original
Random Forests. To restrict the size of terminal nodes, set
`min.bucket`. Variables selected with `always.split.variables` are tried
additionally to the mtry variables randomly selected. In
`split.select.weights`, weights do not need to sum up to 1, they will be
normalized later. The weights are assigned to the variables in the order
they appear in the formula or in the data if no formula is used. Names
of the `split.select.weights` vector are ignored. Weights assigned by
`split.select.weights` to variables in `always.split.variables` are
ignored. The usage of `split.select.weights` can increase the
computation times for large forests.

Unordered factor covariates can be handled in 3 different ways by using
`respect.unordered.factors`: For 'ignore' all factors are regarded
ordered, for 'partition' all possible 2-partitions are considered for
splitting. For 'order' and 2-class classification the factor levels are
ordered by their proportion falling in the second class, for regression
by their mean response, as described in Hastie et al. (2009), chapter
9.2.4. For multiclass classification the factor levels are ordered by
the first principal component of the weighted covariance matrix of the
contingency table (Coppersmith et al. 1999), for survival by the median
survival (or the largest available quantile if the median is not
available). The use of 'order' is recommended, as it computationally
fast and can handle an unlimited number of factor levels. Note that the
factors are only reordered once and not again in each split.

The 'impurity_corrected' importance measure is unbiased in terms of the
number of categories and category frequencies and is almost as fast as
the standard impurity importance. It is a modified version of the method
by Sandri & Zuccolotto (2008), which is faster and more memory
efficient. See Nembrini et al. (2018) for details. This importance
measure can be combined with the methods to estimate p-values in
[`importance_pvalues`](http://imbs-hl.github.io/ranger/reference/importance_pvalues.md).
We recommend not to use the 'impurity_corrected' importance when making
predictions since the feature permutation step might reduce predictive
performance (a warning is raised when predicting on new data).

Note that ranger has different default values than other packages. For
example, our default for `mtry` is the square root of the number of
variables for all tree types, whereas other packages use different
values for regression. Also, changing one hyperparameter does not change
other hyperparameters (where possible). For example,
`splitrule="extratrees"` uses randomized splitting but does not disable
bagging as in Geurts et al. (2006). To disable bagging, use
`replace = FALSE, sample.fraction = 1`. This can also be used to grow a
single decision tree without bagging and feature subsetting:
`ranger(..., num.trees = 1, mtry = p, replace = FALSE, sample.fraction = 1)`,
where p is the number of independent variables.

While random forests are known for their robustness, default
hyperparameters not always work well. For example, for high dimensional
data, increasing the `mtry` value and the number of trees `num.trees` is
recommended. For more details and recommendations, see Probst et al.
(2019). To find the best hyperparameters, consider hyperparameter tuning
with the `tuneRanger` or `mlr3` packages.

Out-of-bag prediction error is calculated as accuracy (proportion of
misclassified observations) for classification, as Brier score for
probability estimation, as mean squared error (MSE) for regression and
as one minus Harrell's C-index for survival. Harrell's C-index is
calculated based on the sum of the cumulative hazard function (CHF) over
all timepoints, i.e., `rowSums(chf)`, where `chf` is the the out-of-bag
CHF; for details, see Ishwaran et al. (2008). Calculation of the
out-of-bag prediction error can be turned off with `oob.error = FALSE`.

Regularization works by penalizing new variables by multiplying the
splitting criterion by a factor, see Deng & Runger (2012) for details.
If `regularization.usedepth=TRUE`, \\f^d\\ is used, where *f* is the
regularization factor and *d* the depth of the node. If regularization
is used, multithreading is deactivated because all trees need access to
the list of variables that are already included in the model.

Missing values can be internally handled by setting
`na.action = "na.learn"` (default), by omitting observations with
missing values with `na.action = "na.omit"` or by stopping if missing
values are found with `na.action = "na.fail"`. With
`na.action = "na.learn"`, in each node either all missings go left or
all missings go right. The direction is chosen based on the split
criterion value (i.e., decrease of impurity). For prediction, this
direction is saved as the "default" direction. If a missing occurs in
prediction at a node where there is no default direction, it goes left.

For a large number of variables and data frames as input data the
formula interface can be slow or impossible to use. Alternatively
`dependent.variable.name` (and `status.variable.name` for survival) or
`x` and `y` can be used. Use `x` and `y` with a matrix for `x` to avoid
conversions and save memory. Consider setting `save.memory = TRUE` if
you encounter memory problems for very large datasets, but be aware that
this option slows down the tree growing.

For GWAS data consider combining `ranger` with the `GenABEL` package.
See the Examples section below for a demonstration using `Plink` data.
All SNPs in the `GenABEL` object will be used for splitting. To use only
the SNPs without sex or other covariates from the phenotype file, use
`0` on the right hand side of the formula. Note that missing values are
treated as an extra category while splitting.

By default, ranger uses 2 threads. The default can be changed with: (1)
`num.threads` in ranger/predict call, (2) environment variable
R_RANGER_NUM_THREADS, (3) `options(ranger.num.threads = N)`, (4)
`options(Ncpus = N)`, with precedence in that order.

See <https://github.com/imbs-hl/ranger> for the development version.

## References

- Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of
  random forests for high dimensional data in C++ and R. J Stat Softw
  77:1-17.
  [doi:10.18637/jss.v077.i01](https://doi.org/10.18637/jss.v077.i01) .

- Schmid, M., Wright, M. N. & Ziegler, A. (2016). On the use of
  Harrell's C for clinical risk prediction via random survival forests.
  Expert Syst Appl 63:450-459.
  [doi:10.1016/j.eswa.2016.07.018](https://doi.org/10.1016/j.eswa.2016.07.018)
  .

- Wright, M. N., Dankowski, T. & Ziegler, A. (2017). Unbiased split
  variable selection for random survival forests using maximally
  selected rank statistics. Stat Med 36:1272-1284.
  [doi:10.1002/sim.7212](https://doi.org/10.1002/sim.7212) .

- Nembrini, S., Koenig, I. R. & Wright, M. N. (2018). The revival of the
  Gini Importance? Bioinformatics.
  [doi:10.1093/bioinformatics/bty373](https://doi.org/10.1093/bioinformatics/bty373)
  .

- Breiman, L. (2001). Random forests. Mach Learn, 45:5-32.
  [doi:10.1023/A:1010933404324](https://doi.org/10.1023/A%3A1010933404324)
  .

- Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S.
  (2008). Random survival forests. Ann Appl Stat 2:841-860.
  [doi:10.1097/JTO.0b013e318233d835](https://doi.org/10.1097/JTO.0b013e318233d835)
  .

- Malley, J. D., Kruppa, J., Dasgupta, A., Malley, K. G., & Ziegler, A.
  (2012). Probability machines: consistent probability estimation using
  nonparametric learning machines. Methods Inf Med 51:74-81.
  [doi:10.3414/ME00-01-0052](https://doi.org/10.3414/ME00-01-0052) .

- Hastie, T., Tibshirani, R., Friedman, J. (2009). The Elements of
  Statistical Learning. Springer, New York. 2nd edition.

- Geurts, P., Ernst, D., Wehenkel, L. (2006). Extremely randomized
  trees. Mach Learn 63:3-42.
  [doi:10.1007/s10994-006-6226-1](https://doi.org/10.1007/s10994-006-6226-1)
  .

- Meinshausen (2006). Quantile Regression Forests. J Mach Learn Res
  7:983-999. <https://www.jmlr.org/papers/v7/meinshausen06a.html>.

- Sandri, M. & Zuccolotto, P. (2008). A bias correction algorithm for
  the Gini variable importance measure in classification trees. J Comput
  Graph Stat, 17:611-628.
  [doi:10.1198/106186008X344522](https://doi.org/10.1198/106186008X344522)
  .

- Coppersmith D., Hong S. J., Hosking J. R. (1999). Partitioning nominal
  attributes in decision trees. Data Min Knowl Discov 3:197-217.
  [doi:10.1023/A:1009869804967](https://doi.org/10.1023/A%3A1009869804967)
  .

- Deng & Runger (2012). Feature selection via regularized trees. The
  2012 International Joint Conference on Neural Networks (IJCNN),
  Brisbane, Australia.
  [doi:10.1109/IJCNN.2012.6252640](https://doi.org/10.1109/IJCNN.2012.6252640)
  .

- Probst, P., Wright, M. N. & Boulesteix, A-L. (2019). Hyperparameters
  and tuning strategies for random forest. WIREs Data Mining Knowl
  Discov
  9:e1301.[doi:10.1002/widm.1301](https://doi.org/10.1002/widm.1301) .

## See also

[`predict.ranger`](http://imbs-hl.github.io/ranger/reference/predict.ranger.md)

## Author

Marvin N. Wright

## Examples

``` r
## Classification forest with default settings
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

## Prediction
train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
iris.train <- iris[train.idx, ]
iris.test <- iris[-train.idx, ]
rg.iris <- ranger(Species ~ ., data = iris.train)
pred.iris <- predict(rg.iris, data = iris.test)
table(iris.test$Species, pred.iris$predictions)
#>             
#>              setosa versicolor virginica
#>   setosa         11          0         0
#>   versicolor      0         17         2
#>   virginica       0          3        17

## Quantile regression forest
rf <- ranger(mpg ~ ., mtcars[1:26, ], quantreg = TRUE)
pred <- predict(rf, mtcars[27:32, ], type = "quantiles")
pred$predictions
#>      quantile= 0.1 quantile= 0.5 quantile= 0.9
#> [1,]          21.0          24.4          33.9
#> [2,]          21.0          22.8          32.4
#> [3,]          13.3          15.5          30.4
#> [4,]          15.2          21.0          22.8
#> [5,]          13.3          14.3          19.2
#> [6,]          21.0          22.8          32.4

## Variable importance
rg.iris <- ranger(Species ~ ., data = iris, importance = "impurity")
rg.iris$variable.importance
#> Sepal.Length  Sepal.Width Petal.Length  Petal.Width 
#>     9.244839     2.440269    45.990477    41.568461 

## Survival forest
require(survival)
#> Loading required package: survival
rg.veteran <- ranger(Surv(time, status) ~ ., data = veteran)
plot(rg.veteran$unique.death.times, rg.veteran$survival[1,])


## Alternative interfaces (same results)
ranger(dependent.variable.name = "Species", data = iris)
#> Ranger result
#> 
#> Call:
#>  ranger(dependent.variable.name = "Species", data = iris) 
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
ranger(y = iris[, 5], x = iris[, -5])
#> Ranger result
#> 
#> Call:
#>  ranger(y = iris[, 5], x = iris[, -5]) 
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

if (FALSE) { # \dontrun{
## Use GenABEL interface to read Plink data into R and grow a classification forest
## The ped and map files are not included
library(GenABEL)
convert.snp.ped("data.ped", "data.map", "data.raw")
dat.gwaa <- load.gwaa.data("data.pheno", "data.raw")
phdata(dat.gwaa)$trait <- factor(phdata(dat.gwaa)$trait)
ranger(trait ~ ., data = dat.gwaa)
} # }
```
