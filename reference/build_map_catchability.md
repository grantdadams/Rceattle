# Helper to set map for Catchability parameters

Maps survey catchability base parameters (`index_log_q`), time-varying
deviations (`index_q_dev`), and environmental linkages (`index_q_beta`,
`index_q_rho`).

## Usage

``` r
build_map_catchability(map_list, data_list, nyrs_hind)
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
