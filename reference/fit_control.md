# Bundle the optimizer / sdreport / phasing controls for `fit_mod()`

[`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
carries roughly a dozen optimizer- and reporting-related arguments
(`bias.correct`, `getsd`, `loopnum`, `newtonsteps`, ...). That is a lot
of surface area when the user mostly cares about "what model am I
fitting" rather than "how is it being fit." `fit_control()` collects
those knobs into a single object so calls to
[`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
can stay focused on the model spec:

## Usage

``` r
fit_control(
  bias.correct = FALSE,
  getsd = TRUE,
  getJointPrecision = TRUE,
  getReportCovariance = FALSE,
  use_gradient = TRUE,
  rel_tol = 1,
  loopnum = 5,
  newtonsteps = 0,
  phase = FALSE,
  TMBfilename = NULL,
  verbose = 1,
  nlminb_control = list(eval.max = 1e+09, iter.max = 1e+09, trace = 0)
)
```

## Arguments

- bias.correct:

  logical. If `TRUE`, applies bias correction via
  [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html).
  Default `FALSE`.

- getsd:

  logical. If `TRUE`, run
  [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html) after
  optimization. Default `TRUE`.

- getJointPrecision:

  logical. Return the full Hessian of fixed and random effects. Default
  `TRUE` (matches
  [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
  default).

- getReportCovariance:

  logical. Return the variance-covariance of `ADREPORT` variables.
  Default `FALSE`.

- use_gradient:

  logical. Use the analytic gradient during phasing. Default `TRUE`.

- rel_tol:

  Numeric tolerance used to flag discontinuous likelihood warnings
  (compares the TMB and `nlminb` objectives). Default `1`.

- loopnum:

  Integer. Number of times to re-start optimization (`loopnum = 3`
  sometimes achieves a lower final gradient than `loopnum = 1`). Default
  `5`.

- newtonsteps:

  Integer. Number of extra Newton steps to take after optimization
  (alternative to `loopnum`). Default `0`.

- phase:

  `TRUE`/`FALSE` or a list. If `FALSE`, the model is not phased. If
  `TRUE`, default phasing is used. Can also accept a list of parameter
  object names with corresponding phase. Default `FALSE`.

- TMBfilename:

  Optional character. Path (without `.cpp`) to an alternate TMB template
  for development. Default `NULL` (use the bundled `ceattle_v01_11`).

- verbose:

  `0` = silent, `1` = print updates of model fit, `2` = print updates of
  model fit and TMB estimation progress. Default `1`.

- nlminb_control:

  A list of control parameters passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html). See
  [`?nlminb`](https://rdrr.io/r/stats/nlminb.html). Default
  `list(eval.max = 1e9, iter.max = 1e9, trace = 0)`.

## Value

A list of class `"Rceattle_fit_control"`.

## Details

    fit <- fit_mod(
      data_list   = BS2017SS,
      msmMode     = 0,
      fit_control = fit_control(loopnum = 1, getsd = FALSE)
    )

Pass the result via the `fit_control` argument to
[`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md).
When supplied, the values in the `fit_control` object override the
corresponding individual arguments to
[`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md).
Individual arguments are kept for backward compatibility.

## Examples

``` r
# Quick-and-dirty fit: skip sdreport, single optimizer pass
ctl <- fit_control(getsd = FALSE, loopnum = 1)

# Production fit with bias correction and joint precision
ctl <- fit_control(bias.correct = TRUE, getJointPrecision = TRUE)
```
