# Simulate Rceattle data

Simulates data used in Rceattle from the expected values estimated from
an existing Rceattle model. The variances and uncertainty are consistent
with those used in the operating model. The function simulates: survey
biomass (log-normal), catch-at-age/length composition (multinomial or
dirichlet-multinomial), conditional-age-at-length (CAAL; multinomial or
dirichlet-multinomial), and total catch (log-normal).

## Usage

``` r
sim_mod(Rceattle, simulate = FALSE)
```

## Arguments

- Rceattle:

  A CEATTLE model object exported from `Rceattle`.

- simulate:

  Logical. If `TRUE`, simulates data from distributions. If `FALSE`,
  returns the expected values (hats).

## Value

A `data_list` object containing the simulated or expected data values,
formatted for use in `Rceattle`.
