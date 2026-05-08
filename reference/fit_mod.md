# This function runs CEATTLE

This function estimates population parameters of CEATTLE using maximum
likelihood in TMB.

## Usage

``` r
fit_mod(
  data_list = NULL,
  inits = NULL,
  map = NULL,
  bounds = NULL,
  file = NULL,
  estimateMode = 0,
  projection_uncertainty = FALSE,
  random_rec = FALSE,
  random_q = FALSE,
  random_sel = FALSE,
  HCR = build_hcr(),
  niter = 3,
  recFun = build_srr(),
  M1Fun = build_M1(),
  growthFun = build_growth(),
  msmMode = 0,
  avgnMode = 0,
  initMode = "NonEquilibrium",
  suitMode = 0,
  suit_styr = NULL,
  suit_endyr = NULL,
  fit_control = NULL,
  ...
)
```

## Arguments

- data_list:

  A data list read in via
  [`read_data`](https://grantdadams.github.io/Rceattle/reference/read_data.md)
  or built directly in R; see
  [`vignette("data-without-excel", package = "Rceattle")`](https://grantdadams.github.io/Rceattle/articles/data-without-excel.md).

- inits:

  (Optional) A named list of initial parameter values, as returned by
  [`build_params`](https://grantdadams.github.io/Rceattle/reference/build_params.md)
  or extracted from a previous fit (`model$estimated_params`). If
  `NULL`, parameters are initialized from scratch via
  [`build_params`](https://grantdadams.github.io/Rceattle/reference/build_params.md).

- map:

  (Optional) A map object from
  [`build_map`](https://grantdadams.github.io/Rceattle/reference/build_map.md).

- bounds:

  (Optional) A bounds object from
  [`build_bounds`](https://grantdadams.github.io/Rceattle/reference/build_bounds.md).

- file:

  (Optional) Filename where files will be saved. If NULL, no file is
  saved.

- estimateMode:

  0 = Fit the hindcast model and projection with HCR specified via
  `HCR`. 1 = Fit the hindcast model only (no projection). 2 = Run the
  projection only with HCR specified via `HCR` given the initial
  parameters in `inits`. 3 = debug mode 1: runs the model through
  MakeADFun, but not nlminb, 4 = runs the model through MakeADFun and
  nlminb (will all parameters mapped out).

- projection_uncertainty:

  logical. If TRUE, accounts for hindcast parameter uncertainty in
  projections when using an HCR. Default is FALSE for speed.

- random_rec:

  logical. If TRUE, treats recruitment deviations as random effects
  using the laplace approximation.The default is FALSE.

- random_q:

  logical. If TRUE, treats annual catchability deviations as random
  effects using the laplace approximation.The default is FALSE.

- random_sel:

  logical. If TRUE, treats annual selectivity deviations as random
  effects using the laplace approximation.The default is FALSE.

- HCR:

  HCR list object from
  [`build_hcr`](https://grantdadams.github.io/Rceattle/reference/build_hcr.md)

- niter:

  Number of iterations for multispecies model

- recFun:

  The stock recruit-relationship parameterization from
  [`build_srr`](https://grantdadams.github.io/Rceattle/reference/build_srr.md).

- M1Fun:

  M1 parameterizations and priors. Use `build_M1`.

- growthFun:

  The weight-at-age parameterization from
  [`build_growth`](https://grantdadams.github.io/Rceattle/reference/build_growth.md).

- msmMode:

  The predation mortality functions to used. Defaults to no predation
  mortality used.

- avgnMode:

  The average abundance-at-age approximation to be used for predation
  mortality equations. 0 (default) is the \\N/Z ( 1 - exp(-Z) )\\, 1 is
  \\N exp(-Z/2)\\, 2 is \\N\\.

- initMode:

  how the population is initialized. 0 = initial age-structure estimated
  as free parameters; 1 = equilibrium age-structure estimated out from
  R0 + mortality (M1); 2 = non-equilibrium age-structure estimated out
  from R0, mortality (M1), and initial population deviates; 3 =
  non-equilibrium age-structure estimated out from initial fishing
  mortality (Finit), R0, mortality (M1), and initial population
  deviates; 4 = non-equilibrium age-structure version 2 where initial
  fishing mortality (Finit) scales R0.

- suitMode:

  Switch for suitability derivation for each predator (single value or
  vector). 0 = empirical based on diet data (Holsman et al. 2015), 1 =
  length-based gamma suitability, 2 = weight-based gamma suitability, 3
  = length-based lognormal suitability, 4 = weight-based lognormal
  suitability, 5 = length-based normal suitability, 6 = weight-based
  normal suitability.

- suit_styr:

  Integer. The first year used to calculate mean suitability. Defaults
  to \$styr\$ in \$data_list\$. Used when diet data were sampled from a
  subset of years.

- suit_endyr:

  Integer. The last year used to calculate mean suitability. Defaults to
  \$endyr\$ in \$data_list\$. Used when diet data were sampled from a
  subset of years.

- fit_control:

  A list returned by
  [`fit_control()`](https://grantdadams.github.io/Rceattle/reference/fit_control.md)
  that bundles the optimizer / sdreport / phasing knobs (`phase`,
  `bias.correct`, `getsd`, `getJointPrecision`, `getReportCovariance`,
  `use_gradient`, `rel_tol`, `loopnum`, `newtonsteps`, `TMBfilename`,
  `verbose`, `nlminb_control`). Defaults to
  [`fit_control()`](https://grantdadams.github.io/Rceattle/reference/fit_control.md).
  See
  [`fit_control()`](https://grantdadams.github.io/Rceattle/reference/fit_control.md)
  for the meaning and defaults of each field.

- ...:

  Deprecated optimizer / sdreport / phasing arguments (e.g. `phase`,
  `getsd`, `bias.correct`, `use_gradient`, `rel_tol`, `control`,
  `getJointPrecision`, `getReportCovariance`, `loopnum`, `newtonsteps`,
  `verbose`, `TMBfilename`). These are forwarded into `fit_control` with
  a deprecation warning; pass them via
  [`fit_control()`](https://grantdadams.github.io/Rceattle/reference/fit_control.md)
  instead.

## Value

A list of class "Rceattle" including:

- data_list: List of data inputs

- initial_params: List of starting parameters

- bounds: Parameter bounds used for estimation

- map: List of map used in TMB

- obj: TMB model object

- opt: Optimized model object from `nlminb`

- sdrep: Object of class `sdreport` exported by TMB including the
  standard errors of estimated parameters

- estimated_params: List of estimated parameters

- quantities: Derived quantities from CEATTLE

- run_time: Model run time

## Details

CEATTLE is an age-structured population dynamics model that can be fit
with or without predation mortality. The default is to exclude predation
mortality by setting `msmMode` to 0. Predation mortality can be included
by setting `msmMode` with the following options:

- 0\. Single species mode

- 1\. Holsman et al. 2015 predation based on multi-species virtual
  population analysis (MSVPA) based predation formation.

- 2\. MSVPA Holling Type III

Values 3 through 9 (Kinzey & Punt 2009 functional responses – Holling
Type I/II/III, predator interference, predator preemption,
Hassell-Varley, Ecosim) are blocked at runtime by
[`data_check()`](https://grantdadams.github.io/Rceattle/reference/data_check.md)
because the implementations have not been validated against the current
parameter set. See `src/TMB/predation.hpp`.

## Examples

``` r
# \donttest{
data(BS2017SS)
ss_run <- fit_mod(
  data_list    = BS2017SS,
  estimateMode = 0,
  msmMode      = 0,
  fit_control  = fit_control(phase = FALSE, verbose = 0)
)
#> 'Diet_loglike' are not included in data, assuming 'Multinomial'
#> 'Selectivity_dimension' not specified in 'fleet_control', assuming 'Age'
#> 'CAAL_weights' not specified in 'fleet_control', assuming 1
#> `age_trans_matrix` data does not span range of age for species 1 will fill with 0s
# }
```
