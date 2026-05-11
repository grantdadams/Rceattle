# Build a single linkage row

Convenience constructor that returns a one-row linkage table with
default values for optional columns. Useful in tests and for incremental
table assembly.

## Usage

``` r
linkage_row(
  process,
  param,
  X_col,
  species = NA_integer_,
  sex = NA_integer_,
  age_bin = NA_integer_,
  design_col = NA_character_,
  link = "identity",
  init = 0,
  init_supplied = FALSE,
  lower = -Inf,
  upper = Inf,
  prior_family = "none",
  prior_p1 = NA_real_,
  prior_p2 = NA_real_,
  re_group = NA_character_,
  est_phase = 1L
)
```

## Arguments

- process, param, X_col:

  required identifying fields.

- species, sex, age_bin:

  stratum ids; `NA` = shared across the dimension.

- design_col:

  name of the design matrix column.

- link:

  link function; one of
  [LINKAGE_LINKS](https://grantdadams.github.io/Rceattle/reference/LINKAGE_PROCESSES.md).

- init:

  initial value (default `0`).

- lower, upper:

  bounds (default `-Inf`, `Inf`).

- prior_family:

  one of
  [PRIOR_FAMILIES](https://grantdadams.github.io/Rceattle/reference/PRIOR_FAMILIES.md).
  `"none"` = no prior.

- prior_p1, prior_p2:

  family-specific prior parameters; ignored when
  `prior_family == "none"`.

- re_group:

  random-effect grouping label; `NA` = fixed.

- est_phase:

  estimation phase ordinal; `0` = fix at `init`.

## Value

A one-row `Rceattle_linkage_table`.
