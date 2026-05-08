# Sentinel meaning "applies to every level of this stratum"

Stored where the user passed `NA` for `species`, `sex`, or `age_bin`.
The TMB template expands these by iterating over `1:nspp` / `1:nsex` /
`1:nages` rather than indexing a single cell.

## Usage

``` r
LINKAGE_STRATUM_ALL
```

## Format

An object of class `integer` of length 1.
