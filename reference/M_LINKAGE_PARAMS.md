# Allowed M-parameter names for `linkages` in [`build_M1()`](https://grantdadams.github.io/Rceattle/reference/build_M1.md)

Natural-scale names of the underlying natural-mortality parameters that
the linkage system can address. Currently just `M1` – with the default
log link the offset is added on the log scale to `log_M1` (broadcast
across age unless the linkage row pins a specific `age_bin`).

## Usage

``` r
M_LINKAGE_PARAMS
```

## Format

An object of class `character` of length 1.

## Details

Note: predation mortality `M2` is a derived quantity in CEATTLE (a
function of predator abundance, suitability, and ration), not a
parameter. There is no `M2` linkage target; environmental effects on
predation are mediated upstream via recruitment, growth, suitability, or
ration inputs.
