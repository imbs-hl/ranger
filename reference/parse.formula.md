# Parse formula

Parse formula and return dataset containing selected columns.
Interactions are supported for numerical columns only. An interaction
column is the product of all interacting columns.

## Usage

``` r
parse.formula(formula, data, env = parent.frame())
```

## Arguments

- formula:

  Object of class `formula` or `character` describing the model to fit.

- data:

  Training data of class `data.frame`.

- env:

  The environment in which the left hand side of `formula` is evaluated.

## Value

Dataset including selected columns and interactions.
