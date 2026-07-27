# Coerce the `Rceattle` plotting argument to a list of fits

Accepts a single `Rceattle` fit, a list of fits (multi-model overlay),
or an MSE object, and returns a plain list of `Rceattle` fits.

## Usage

``` r
.as_model_list(Rceattle, mse = FALSE, OM = TRUE)
```

## Arguments

- Rceattle:

  A fit, list of fits, or (when `mse = TRUE`) MSE object.

- mse, OM:

  When `mse = TRUE`, pull the operating model (`OM = TRUE`) or the
  terminal estimation model from each MSE element.

## Value

A list of `Rceattle` fits.
