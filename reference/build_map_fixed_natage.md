# Helper to set map for Fixed N-at-Age models

Turns off (sets to NA) most population and fleet parameters for species
where the dynamics are fixed (`estDynamics > 0`).

## Usage

``` r
build_map_fixed_natage(map_list, data_list)
```

## Arguments

- map_list:

  The current TMB map list.

- data_list:

  The data list containing model settings.

## Value

Updated `map_list`.
