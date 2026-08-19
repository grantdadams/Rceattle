# Fit the CEATTLE assessment model

Estimate CEATTLE population parameters by maximum likelihood, and
optionally project the stock and apply a harvest control rule.

## Usage

``` r
fit_mod(
  data_list = NULL,
  inits = NULL,
  map = NULL,
  bounds = NULL,
  file = NULL,
  estimateMode = 0,
  random_rec = FALSE,
  random_q = FALSE,
  random_sel = FALSE,
  HCR = build_hcr(),
  niter = 3,
  recFun = build_srr(),
  M1Fun = build_M1(),
  growthFun = build_growth(),
  qFun = build_catchability(),
  selFun = build_selectivity(),
  compFun = build_composition(),
  msmMode = 0,
  avgnMode = 0,
  initMode = "NonEquilibrium",
  suitMode = 0,
  suit_styr = NULL,
  suit_endyr = NULL,
  fit_control = NULL,
  config = NULL,
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

  What to fit, given as a string alias or the integer code:
  `"Estimate"` (0) = fit the hindcast model and the HCR projection
  (`HCR`); `"Hindcast"` (1) = fit the hindcast only (no fitting
  BRPs/HCR/projection); `"Projection"` (2) = fit the BRPs/HCR/projection
  only, from the initial parameters in `inits`; `"DebugBuild"` (3) =
  build through `MakeADFun` but not `nlminb` – the returned `obj`
  carries the real objective and gradient, so `obj$fn()` / `obj$gr()`
  are usable for diagnosing a model before committing to a fit;
  `"DebugOptimize"` (4) = optimize with all parameters mapped out, so
  the objective is a placeholder (`dummy^2`), not a likelihood. Defaults
  to `"Estimate"`.

- random_rec:

  logical. If TRUE, treats recruitment deviations as random effects
  using the Laplace approximation. The default is FALSE.

- random_q:

  logical. If TRUE, treats annual catchability deviations as random
  effects using the Laplace approximation. The default is FALSE.

- random_sel:

  logical. If TRUE, treats annual selectivity deviations as random
  effects using the Laplace approximation. The default is FALSE.

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

- qFun:

  Catchability specification from
  [`build_catchability`](https://grantdadams.github.io/Rceattle/reference/build_catchability.md),
  carrying any environmental linkages on q.

- selFun:

  Selectivity specification from
  [`build_selectivity`](https://grantdadams.github.io/Rceattle/reference/build_selectivity.md),
  carrying any environmental linkages on selectivity parameters.

- compFun:

  Composition-weighting specification from
  [`build_composition`](https://grantdadams.github.io/Rceattle/reference/build_composition.md),
  carrying any priors on the Dirichlet-multinomial weights.

- msmMode:

  The predation-mortality mode, as a string alias or integer code:
  `"SingleSpecies"` (0, the default, no predation), `"MSVPA"` (1, the
  Type-II MSVPA predation of Holsman et al. 2015) or `"TypeIIIMSVPA"`
  (2). Higher integer codes (Kinzey-Punt, Holling forms) are not yet
  implemented.

- avgnMode:

  The average abundance-at-age approximation used in the
  predation-mortality equations. Only mode 0, \\N/Z ( 1 - exp(-Z) )\\
  (the MSVPA form), is currently active; the model always uses it. The
  alternatives \\N exp(-Z/2)\\ (1) and \\N\\ (2) are not currently
  implemented and setting them has no effect.

- initMode:

  how the population is initialized, as a string alias or integer code.
  `"FreeParams"` (0) = initial age-structure estimated as free
  parameters; `"Equilibrium"` (1) = unfished (Finit = 0) equilibrium
  age-structure estimated out from R0 + mortality (M1);
  `"NonEquilibrium"` (2, the default) = non-equilibrium age-structure
  estimated out from R0, mortality (M1), and initial population
  deviates; `"FishedNonEquilibrium"` (3) = non-equilibrium age-structure
  estimated out from initial fishing mortality (Finit), R0, mortality
  (M1), and initial population deviates;
  `"FishedNonEquilibriumScaled"` (4) = non-equilibrium age-structure
  version 2 where initial fishing mortality (Finit) scales R0;
  `"OffsetEquilibrium"` (5) = unfished (Finit = 0) equilibrium
  age-structure seeded by the first-year recruitment
  (`R_init * exp(rec_dev[year 1])`) decayed by residual natural
  mortality M1, closed with the usual geometric plus group at the
  maximum age, with initial deviates turned off and no init-dev penalty
  (the Cole Monnahan / AFSC GOA pollock convention). Mode 5 differs from
  `"Equilibrium"` by exactly one term: both start from the initial
  equilibrium recruitment `R_init`, but 5 displaces it by the year-1
  recruitment deviation. Under a random-about-mean stock-recruit
  relationship `R_init` equals `R0`; under Beverton-Holt or Ricker it is
  the equilibrium recruitment implied by the curve. Note the decay uses
  M1 only, so in multispecies mode predation mortality (M2) does not
  enter the initial age structure.

- suitMode:

  Switch for suitability derivation for each predator (single value or
  vector). 0 = empirical based on diet data (Holsman et al. 2015), 1 =
  length-based gamma suitability, 2 = weight-based gamma suitability, 3
  = length-based lognormal suitability, 4 = weight-based lognormal
  suitability, 5 = length-based normal suitability, 6 = weight-based
  normal suitability. The length-based modes (1, 3, 5) are not yet
  implemented and are rejected by
  [`data_check()`](https://grantdadams.github.io/Rceattle/reference/data_check.md);
  use a weight-based mode (2, 4, 6) or empirical suitability (0).

- suit_styr:

  The first year used to calculate mean suitability. A single integer is
  applied to every predator, or a vector of length `nspp` sets a
  distinct start year per predator. Defaults to `styr` in `data_list`.
  Used when diet data were sampled from a subset of years.

- suit_endyr:

  The last year used to calculate mean suitability. A single integer is
  applied to every predator, or a vector of length `nspp` sets a
  distinct end year per predator. Defaults to `endyr` in `data_list`.
  Used when diet data were sampled from a subset of years.

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

- config:

  (Optional) An `Rceattle_run_config` from
  [`load_config()`](https://grantdadams.github.io/Rceattle/reference/load_config.md)
  (or
  [`run_config()`](https://grantdadams.github.io/Rceattle/reference/run_config.md)).
  Its stored `model_config` structure and estimation controls
  (`estimateMode`, `random_rec`/`random_q`/`random_sel`, `suit_styr`/
  `suit_endyr`, `fit_control`) overlay only the arguments the caller did
  *not* pass – an explicit argument always wins. `NULL` (default)
  applies no configuration. Example:
  `fit_mod(data_list, config = load_config("run.yaml"))`.

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
#> `age_trans_matrix` data does not span range of age for species 1 will fill with 0s
# }
```
