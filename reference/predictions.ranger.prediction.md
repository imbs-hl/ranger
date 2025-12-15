# Ranger predictions

Extract predictions of Ranger prediction object.

## Usage

``` r
# S3 method for class 'ranger.prediction'
predictions(x, ...)
```

## Arguments

- x:

  Ranger prediction object.

- ...:

  Further arguments passed to or from other methods.

## Value

Predictions: Classes for Classification forests, Numerical values for
Regressions forests and the estimated survival functions for all
individuals for Survival forests.

## See also

[`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md)

## Author

Marvin N. Wright
