# Deforesting a random forest

The main purpose of this function is to allow for post-processing of
ensembles via L2 regularized regression (i.e., the LASSO), as described
in Friedman and Popescu (2003). The basic idea is to use the LASSO to
post-process the predictions from the individual base learners in an
ensemble (i.e., decision trees) in the hopes of producing a much smaller
model without sacrificing much in the way of accuracy, and in some
cases, improving it. Friedman and Popescu (2003) describe conditions
under which tree-based ensembles, like random forest, can potentially
benefit from such post-processing (e.g., using shallower trees trained
on much smaller samples of the training data without replacement).
However, the computational benefits of such post-processing can only be
realized if the base learners "zeroed out" by the LASSO can actually be
removed from the original ensemble, hence the purpose of this function.
A complete example using
[`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md) can be
found at <https://github.com/imbs-hl/ranger/issues/568>.

## Usage

``` r
deforest(object, which.trees = NULL, ...)

# S3 method for class 'ranger'
deforest(object, which.trees = NULL, warn = TRUE, ...)
```

## Arguments

- object:

  A fitted random forest (e.g., a
  [`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md)
  object).

- which.trees:

  Vector giving the indices of the trees to remove.

- ...:

  Additional (optional) arguments. (Currently ignored.)

- warn:

  Logical indicating whether or not to warn users that some of the
  standard output of a typical
  [`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md) object
  or no longer available after deforestation. Default is `TRUE`.

## Value

An object of class `"deforest.ranger"`; essentially, a
[`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md) object
with certain components replaced with `NA`s (e.g., out-of-bag (OOB)
predictions, variable importance scores (if requested), and OOB-based
error metrics).

## Note

This function is a generic and can be extended by other packages.

## References

Friedman, J. and Popescu, B. (2003). Importance sampled learning
ensembles, Technical report, Stanford University, Department of
Statistics. <https://jerryfriedman.su.domains/ftp/isle.pdf>.

## Author

Brandon M. Greenwell

## Examples

``` r
## Example of deforesting a random forest
rfo <- ranger(Species ~ ., data = iris, probability = TRUE, num.trees = 100)
dfo <- deforest(rfo, which.trees = c(1, 3, 5))
#> Warning: Many of the components of a typical "ranger" object are not available after deforestation and are instead replaced with `NA` (e.g., out-of-bag (OOB) predictions, variable importance scores (if requested), and OOB-based error metrics).
dfo  # same as `rfo` but with trees 1, 3, and 5 removed
#> Ranger (deforested) result
#> 
#> Note that many of the components of a typical "ranger" object are not available after deforestation and are instead replaced with `NA` (e.g., out-of-bag (OOB) predictions, variable importance scores (if requested), and OOB-based error metrics) 
#> 
#> Type:                             Probability estimation 
#> Number of trees:                  97 
#> Sample size:                      150 
#> Number of independent variables:  4 
#> Mtry:                             2 
#> Target node size:                 10 
#> Variable importance mode:         none 
#> Splitrule:                        gini 
#> OOB prediction error (Brier s.):  NA 

## Sanity check
preds.rfo <- predict(rfo, data = iris, predict.all = TRUE)$predictions
preds.dfo <- predict(dfo, data = iris, predict.all = TRUE)$predictions
identical(preds.rfo[, , -c(1, 3, 5)], y = preds.dfo)
#> [1] TRUE
```
