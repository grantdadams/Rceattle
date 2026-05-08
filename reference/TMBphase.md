# Run TMB using phases

This function runs TMB with ADMB-like phasing of parameter estimation.
Function with normal inputs, passed via "...", plus two additional
arguments, "phase" Optimizer by default is nlminb phase is a tagged list
where missing elements are populated with a vector of 1s, and
non-missing elements are integers, and where the optimizer loops through
values of phase while progressively changing map to turn on parameters

## Usage

``` r
TMBphase(
  data,
  parameters,
  map,
  random,
  phases,
  model_name,
  silent,
  use_gradient = TRUE,
  control = list(eval.max = 1e+09, iter.max = 1e+09, trace = 0)
)
```

## Arguments

- data:

  A list to be passed to TMB

- parameters:

  A list of parameters of the model

- map:

  a list of map object from the model

- random:

  A character vector of names of parameters that are random effects

- phases:

  A list of the phases for the parameters of the model (same structure
  as your parameter list)

- model_name:

  A string describing the model name. Must be the name of your .cpp file

- silent:

  logical. If TRUE, suppresses output from TMB (default = TRUE).

- use_gradient:

  logical. If TRUE, uses gradient in optimization (default = TRUE).

- control:

  A list of control parameters. For details see
  [`?nlminb`](https://rdrr.io/r/stats/nlminb.html)

## Value

A list of parameter estimates and their standard errors

## Author

Gavin Fay
https://github.com/kaskr/TMB_contrib_R/blob/master/TMBphase/R/TMBphase.R
