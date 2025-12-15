# Hierarchical shrinkage

Apply hierarchical shrinkage to a ranger object. Hierarchical shrinkage
is a regularization technique that recursively shrinks node predictions
towards parent node predictions. For details see Agarwal et al. (2022).

## Usage

``` r
hshrink(rf, lambda)
```

## Arguments

- rf:

  ranger object, created with `node.stats = TRUE`.

- lambda:

  Non-negative shrinkage parameter.

## Value

The ranger object is modified in-place.

## References

- Agarwal, A., Tan, Y.S., Ronen, O., Singh, C. & Yu, B. (2022).
  Hierarchical Shrinkage: Improving the accuracy and interpretability of
  tree-based models. Proceedings of the 39th International Conference on
  Machine Learning, PMLR 162:111-135.

## Author

Marvin N. Wright

## Examples

``` r
## Hierarchical shrinkage for a probablity forest
rf <- ranger(Species ~ ., iris, node.stats = TRUE, probability = TRUE)
hshrink(rf, lambda = 5)
```
