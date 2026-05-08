# Evaluate simulation performance

Function to evaluate the simulation performance with regard to bias
using the median relative error (MRE) and precision using the
coefficient of variation.

## Usage

``` r
compare_sim(operating_mod, simulation_mods, object = "quantities")
```

## Arguments

- operating_mod:

  CEATTLE model object exported from `Rceattle` to be used as the
  operating model

- simulation_mods:

  List of CEATTLE model objects exported from `Rceattle` fit to
  simulated data

- object:

  character string specifying which part of the model to compare
  (default = "quantities")

## Value

A data frame summarising simulation performance metrics
