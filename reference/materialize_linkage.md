# Materialize a linkage spec into linkage-table rows

Expands an
[`linkage_spec()`](https://grantdadams.github.io/Rceattle/reference/linkage_spec.md)
into rows of the canonical linkage table. One row is produced for every
combination of:

## Usage

``` r
materialize_linkage(spec, process, env_data, strata = list())
```

## Arguments

- spec:

  an `Rceattle_linkage_spec`.

- process:

  one of
  [LINKAGE_PROCESSES](https://grantdadams.github.io/Rceattle/reference/LINKAGE_PROCESSES.md).

- env_data:

  data frame of environmental covariates (one row per model year). Must
  contain every variable named on the RHS of `spec$formula`.

- strata:

  named list giving the discrete levels for each stratifying factor
  named in `spec$by`. For example, for `by = ~species` the user must
  supply `strata = list(species = 1:3)`. Each element should be a
  1-based integer vector of stratum ids. Allowed names are `"species"`,
  `"sex"`, and `"age_bin"`.

## Value

An `Rceattle_linkage_table` with one row per coefficient.

## Details

- design-matrix column from `model.matrix(spec$formula, env_data)`,

- stratum implied by `spec$by` (e.g. one per species when
  `by = ~species`).

The `X_col` column is initially set to a *local* column index into the
per-spec design matrix; it is the caller's responsibility (in
[`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)'s
pooling step) to remap these to global column indices once all specs
have been combined into a single shared `X` matrix.
