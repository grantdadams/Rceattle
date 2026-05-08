# Helper to set map for Predation Mortality (M2) parameters

Maps predation suitability parameters (`log_gam_a`, `log_gam_b`,
`log_phi`) and diet weight parameters based on `msmMode` and `suitMode`.

## Usage

``` r
build_map_predation(map_list, data_list)
```

## Arguments

- map_list:

  The current TMB map list.

- data_list:

  The data list containing model settings.

## Value

Updated `map_list`.
