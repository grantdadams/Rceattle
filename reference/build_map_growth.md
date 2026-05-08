# Helper to set map for growth parameters

Maps the fixed parameters (`log_growth_pars`) and the random effects
parameters (`log_growth_par_devs`, `growth_dev_log_sd`, `growth_rho`)
based on `growth_model` and `growth_re` settings.

## Usage

``` r
build_map_growth(map_list, data_list, nyrs_hind)
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
