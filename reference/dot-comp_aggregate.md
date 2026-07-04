# Aggregate composition proportions across hindcast years

Sums observed and fitted proportions over years (for one fleet x type
panel) and renormalizes to the mean composition; joint-sex groups keep
their shared normalization (females + males sum to 1).

## Usage

``` r
.comp_aggregate(d)
```

## Arguments

- d:

  One panel's rows from
  [`.comp_resid_long()`](https://grantdadams.github.io/Rceattle/reference/dot-comp_resid_long.md).

## Value

A data frame with `bin`, `sex_grp`, `obs`, `hat`, `y_obs`, `y_hat`.
