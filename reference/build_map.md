# Main function to construct the TMB map argument for CEATTLE

Orchestrates the building of the TMB map object by calling specialized
helper functions for each parameter block (Recruitment, M1, Predation,
Selectivity, Catchability, etc.).

## Usage

``` r
build_map(
  data_list,
  params,
  debug = FALSE,
  random_rec = FALSE,
  random_sel = FALSE
)
```

## Arguments

- data_list:

  an Rceattle data_list

- params:

  A parameter list created from
  [`build_params`](https://grantdadams.github.io/Rceattle/reference/build_params.md).

- debug:

  Logical. If TRUE, sets all map values to NA except the dummy
  parameter, running the model without parameter estimation.

- random_rec:

  Logical. If TRUE, treats recruitment deviations as random effects,
  meaning the variance parameter (`R_log_sd`) is estimated.

- random_sel:

  Logical. If TRUE, treats selectivity deviations as random effects,
  meaning the variance parameter (`sel_dev_log_sd`) is estimated.

## Value

A list containing the factorized TMB map (`mapFactor`) and the original
map matrix/array list (`mapList`).
