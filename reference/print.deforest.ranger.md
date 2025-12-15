# Print deforested ranger summary

Print basic information about a deforested
[`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md) object.

## Usage

``` r
# S3 method for class 'deforest.ranger'
print(x, ...)
```

## Arguments

- x:

  A [`deforest`](http://imbs-hl.github.io/ranger/reference/deforest.md)
  object (i.e., an object that inherits from class `"deforest.ranger"`).

- ...:

  Further arguments passed to or from other methods.

## Note

Many of the components of a typical
[`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md) object
are not available after deforestation and are instead replaced with `NA`
(e.g., out-of-bag (OOB) predictions, variable importance scores (if
requested), and OOB-based error metrics).

## See also

[`deforest`](http://imbs-hl.github.io/ranger/reference/deforest.md).

## Author

Brandon M. Greenwell
