# Helper to set map for Fishing Mortality and Data Weights

Maps fishing mortality parameters (`log_F`) and related targets. Also
maps data weighting parameters (`catch_log_sd`, `comp_weights`,
`caal_weights`).

## Usage

``` r
build_map_f_and_data_weights(map_list, data_list, nyrs_hind)
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
