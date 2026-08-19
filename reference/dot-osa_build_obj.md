# Rebuild a fitted Rceattle TMB object in OSA mode at the fitted parameters

Returns a TMB ADFun object equivalent to `fit$obj` but with
`osa_mode = 1`, built at the fitted parameter values and the same map /
random-effect structure, ready for
[`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html).
In OSA mode the composition likelihoods read their counts from `obsvec`
and use unweighted proper densities; the aggregate catch/index
likelihood is identical in both modes, so this single object serves
every observation type.

## Usage

``` r
.osa_build_obj(fit, osa_dat = NULL, force = FALSE)
```

## Arguments

- fit:

  A fitted `Rceattle` object.

- osa_dat:

  Optional pre-built OSA observation data (the list returned by
  [`build_osa_data()`](https://grantdadams.github.io/Rceattle/reference/build_osa_data.md)
  with `build_osa = TRUE`) to reuse instead of rebuilding it. `NULL`
  (the default) rebuilds it from `fit$obj$env$data`.

- force:

  Build a new object even when `fit` is already in OSA mode, where this
  otherwise returns `fit$obj` itself. The retry after a failed parallel
  one-step-ahead loop needs a genuinely new one.

## Value

A TMB ADFun object with `osa_mode = 1`.
