# Helper to set map for Recruitment parameters

Maps the recruitment deviations (`rec_dev`, `init_dev`),
stock-recruitment parameters (`rec_pars`), their variances (`R_log_sd`),
and environmental linkages (`beta_rec_pars`).

see
[`build_srr()`](https://grantdadams.github.io/Rceattle/reference/build_srr.md)
for options,

## Usage

``` r
build_map_recruitment(map_list, data_list, nyrs_hind, nyrs_proj, random_rec)
```

## Arguments

- map_list:

  The current TMB map list.

- data_list:

  The data list containing model settings.

- nyrs_hind:

  Number of historical years.

- nyrs_proj:

  Total number of years (historical + projected).

- random_rec:

  Logical indicating if recruitment deviations are random effects.

## Value

Updated `map_list`.
