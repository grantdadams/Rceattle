# Helper to set map for Natural Mortality (M1) parameters

Maps the fixed parameters (`log_M1`) and the random effects parameters
(`log_M1_dev`, `M1_dev_log_sd`, `M1_rho`) based on `M1_model` and
`M1_re` settings.

## Usage

``` r
build_map_m1(map_list, data_list, nyrs_hind)
```

## Arguments

- map_list:

  The current TMB map list.

- data_list:

  The data list containing model settings.

- nyrs_hind:

  Number of historical years.

## Value

Updated `map_list`.
