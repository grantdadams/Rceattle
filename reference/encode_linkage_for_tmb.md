# Encode a linkage table into TMB-friendly parallel vectors

Encode a linkage table into TMB-friendly parallel vectors

## Usage

``` r
encode_linkage_for_tmb(table, X)
```

## Arguments

- table:

  an `Rceattle_linkage_table` (typically the output of
  [`pool_linkages()`](https://grantdadams.github.io/Rceattle/reference/pool_linkages.md)).

- X:

  the global design matrix from
  [`pool_linkages()`](https://grantdadams.github.io/Rceattle/reference/pool_linkages.md)
  (passed through unchanged so callers can stash it alongside).

## Value

A named list with components:

- n_linkage:

  `integer(1)`. Number of coefficients (`nrow(table)`).

- linkage_process, linkage_param, linkage_species, linkage_sex,
  linkage_age_bin, linkage_X_col, linkage_link, linkage_prior_family:

  Parallel `integer(n_linkage)` vectors.

- linkage_prior_p1, linkage_prior_p2:

  Parallel `numeric(n_linkage)` vectors.

- linkage_X:

  The design matrix as passed in.

- linkage_design_names:

  `character` colnames of `linkage_X`.
