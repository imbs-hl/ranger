# Tree information in human readable format

Extract tree information of a `ranger` object.

## Usage

``` r
treeInfo(object, tree = 1)
```

## Arguments

- object:

  `ranger` object.

- tree:

  Number of the tree of interest.

## Value

A data.frame with the columns

|                |                                                                                                                                                                                                                                                                             |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `nodeID`       | The nodeID, 0-indexed.                                                                                                                                                                                                                                                      |
| `leftChild`    | ID of the left child node, 0-indexed.                                                                                                                                                                                                                                       |
| `rightChild`   | ID of the right child node, 0-indexed.                                                                                                                                                                                                                                      |
| `splitvarID`   | ID of the splitting variable, 0-indexed. Caution, the variable order changes if the formula interface is used.                                                                                                                                                              |
| `splitvarName` | Name of the splitting variable.                                                                                                                                                                                                                                             |
| `splitval`     | The splitting value. For numeric or ordinal variables, all values smaller or equal go to the left, larger values to the right. For unordered factor variables see above.                                                                                                    |
| `terminal`     | Logical, TRUE for terminal nodes.                                                                                                                                                                                                                                           |
| `prediction`   | One column with the predicted class (factor) for classification and the predicted numerical value for regression. One probability per class for probability estimation in several columns. Nothing for survival, refer to `object$forest$chf` for the CHF node predictions. |
| `numSamples`   | Number of samples in the node (only if ranger called with `node.stats = TRUE`).                                                                                                                                                                                             |
| `splitStat`    | Split statistics, i.e., value of the splitting criterion (only if ranger called with `node.stats = TRUE`).                                                                                                                                                                  |

## Details

Node and variable ID's are 0-indexed, i.e., node 0 is the root node. If
the formula interface is used in the `ranger` call, the variable ID's
are usually different to the original data used to grow the tree. Refer
to the variable name instead to be sure.

Splitting at unordered factors (nominal variables) depends on the option
`respect.unordered.factors` in the `ranger` call. For the "ignore" and
"order" approaches, all values smaller or equal the `splitval` value go
to the left and all values larger go to the right, as usual. However,
with "order" the values correspond to the order in
`object$forest$covariate.levels` instead of the original order (usually
alphabetical). In the "partition" mode, the `splitval` values for
unordered factor are comma separated lists of values, representing the
factor levels (in the original order) going to the right.

## See also

[`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md)

## Author

Marvin N. Wright

## Examples

``` r
rf <- ranger(Species ~ ., data = iris)
treeInfo(rf, 1)
#>    nodeID leftChild rightChild splitvarID splitvarName splitval terminal
#> 1       0         1          2          2 Petal.Length     2.45    FALSE
#> 2       1        NA         NA         NA         <NA>       NA     TRUE
#> 3       2         3          4          0 Sepal.Length     6.15    FALSE
#> 4       3         5          6          2 Petal.Length     4.80    FALSE
#> 5       4         7          8          1  Sepal.Width     2.40    FALSE
#> 6       5        NA         NA         NA         <NA>       NA     TRUE
#> 7       6         9         10          0 Sepal.Length     5.95    FALSE
#> 8       7        NA         NA         NA         <NA>       NA     TRUE
#> 9       8        11         12          0 Sepal.Length     6.85    FALSE
#> 10      9        NA         NA         NA         <NA>       NA     TRUE
#> 11     10        13         14          2 Petal.Length     5.05    FALSE
#> 12     11        15         16          3  Petal.Width     1.75    FALSE
#> 13     12        NA         NA         NA         <NA>       NA     TRUE
#> 14     13        NA         NA         NA         <NA>       NA     TRUE
#> 15     14        NA         NA         NA         <NA>       NA     TRUE
#> 16     15        NA         NA         NA         <NA>       NA     TRUE
#> 17     16        NA         NA         NA         <NA>       NA     TRUE
#>    prediction
#> 1        <NA>
#> 2      setosa
#> 3        <NA>
#> 4        <NA>
#> 5        <NA>
#> 6  versicolor
#> 7        <NA>
#> 8  versicolor
#> 9        <NA>
#> 10  virginica
#> 11       <NA>
#> 12       <NA>
#> 13  virginica
#> 14  virginica
#> 15 versicolor
#> 16 versicolor
#> 17  virginica
```
