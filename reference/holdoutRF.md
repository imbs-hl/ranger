# Hold-out random forests

Grow two random forests on two cross-validation folds. Instead of
out-of-bag data, the other fold is used to compute permutation
importance. Related to the novel permutation variable importance by
Janitza et al. (2015).

## Usage

``` r
holdoutRF(...)
```

## Arguments

- ...:

  All arguments are passed to
  [`ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md)
  (except `importance`, `case.weights`, `replace` and `holdout`.).

## Value

Hold-out random forests with variable importance.

## References

Janitza, S., Celik, E. & Boulesteix, A.-L., (2015). A computationally
fast variable importance test for random forests for high-dimensional
data. Adv Data Anal Classif
[doi:10.1007/s11634-016-0276-4](https://doi.org/10.1007/s11634-016-0276-4)
.  

## See also

[`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md)

## Author

Marvin N. Wright
