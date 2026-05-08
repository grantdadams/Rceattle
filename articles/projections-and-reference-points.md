# 4. Projections and reference points

Once a model has been estimated, Rceattle can project the population
forward under alternative fishing mortality scenarios and recruitment
assumptions. This vignette covers:

- Fixed-F projections (including mean historical F)
- Manual and stochastic recruitment in the projection period
- Random recruitment effects
- Visualising projection uncertainty across scenarios

## Setup

``` r

library(Rceattle)

data(BS2017SS)
data(BS2017MS)

# Set the terminal projection year
BS2017SS$projyr <- 2060
BS2017MS$projyr <- 2060

# Plot input data
plot_data(BS2017SS)

# Fit baseline single- and multi-species models
ss_run <- fit_mod(
  data_list    = BS2017SS,
  inits        = NULL,
  file         = NULL,
  estimateMode = 0,
  random_rec   = FALSE,
  msmMode      = 0,
  fit_control = fit_control(
    phase        = TRUE,
    verbose      = 1))

ms_run <- fit_mod(
  data_list    = BS2017MS,
  inits        = ss_run$estimated_params,
  file         = NULL,
  estimateMode = 0,
  niter        = 5,
  random_rec   = FALSE,
  msmMode      = 1,
  suitMode     = 0,
  fit_control = fit_control(
    verbose      = 1))
```

## Fixed-F projections

`estimateMode = 2` runs only the projection period (holding all
estimated parameters fixed) under a specified harvest control rule.
`build_hcr(HCR = 2, Ftarget = ...)` sets a constant F for each species.

The `proj_F_prop` column in `fleet_control` controls how projected F is
split across fleets when there are multiple fisheries targeting the same
species (values must sum to 1 per species).

``` r

BS2017MS$fleet_control$proj_F_prop <- 1  # One fishery fleet per species

ms_run_proj <- fit_mod(
  data_list    = BS2017MS,
  inits        = ms_run$estimated_params,
  file         = NULL,
  estimateMode = 2,
  HCR          = build_hcr(HCR = "ConstantF",
                            Ftarget = c(0.2342936, 0.513, 0.0774777)),
  niter        = 5,
  random_rec   = FALSE,
  msmMode      = 1,
  suitMode     = 0,
  fit_control = fit_control(
    verbose      = 1))

plot_catch(ms_run_proj, incl_proj = TRUE)
```

## Recruitment in the projection period

### Manual recruitment deviations

Projection recruitment deviations are stored in the parameter list and
can be replaced directly before a forward-projection run.
`estimateMode = 3` evaluates the model without re-estimating any
parameters, so only the updated recruitment deviations affect the
projection.

``` r

nyrs      <- BS2017MS$endyr - BS2017MS$styr + 1
nyrs_proj <- BS2017MS$projyr - BS2017MS$styr + 1
yrs_proj  <- (nyrs + 1):nyrs_proj

# Replace future rec_devs with draws from N(0, 0.707)
ms_run$estimated_params$rec_dev[, yrs_proj] <- stats::rnorm(
  n    = length(ms_run$estimated_params$rec_dev[, yrs_proj]),
  mean = 0,
  sd   = 0.707
)

ms_run_proj2 <- fit_mod(
  data_list    = BS2017MS,
  inits        = ms_run$estimated_params,
  file         = NULL,
  estimateMode = 3,   # Evaluate only — do not re-estimate
  niter        = 5,
  random_rec   = FALSE,
  msmMode      = 1,
  suitMode     = 0,
  fit_control = fit_control(
    verbose      = 1))
```

### Stochastic recruitment via `sample_rec()`

[`sample_rec()`](https://grantdadams.github.io/Rceattle/reference/sample_rec.md)
bootstraps historical recruitment deviations into the projection period.
The optional `rec_trend` argument adds a linear trend:
`rec_trend = -0.5` imposes a 50% decline in mean recruitment by the end
of the projection.

``` r

ms_run_proj3 <- sample_rec(
  Rceattle     = ms_run,
  sample_rec   = TRUE,
  update_model = TRUE,
  rec_trend    = -0.5   # 50% decline in mean recruitment
)
```

## Comparing projection scenarios

``` r

mod_list  <- list(ms_run, ms_run_proj, ms_run_proj2, ms_run_proj3)
mod_names <- c(
  "No fishing",
  "Mean historical F",
  "Mean F + stochastic recruitment",
  "Mean F + declining recruitment"
)

plot_biomass(Rceattle = mod_list, model_names = mod_names, incl_proj = TRUE)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, incl_proj = TRUE)
plot_catch(Rceattle = mod_list, model_names = mod_names, incl_proj = TRUE)
```

## Random recruitment effects

Setting `random_rec = TRUE` treats log-recruitment deviations as random
effects marginalised via the Laplace approximation rather than as
penalised fixed effects. This propagates process uncertainty into
projections automatically.

``` r

ss_re <- fit_mod(
  data_list    = BS2017SS,
  inits        = NULL,
  file         = NULL,
  estimateMode = 0,
  random_rec   = TRUE,
  msmMode      = 0,
  fit_control = fit_control(
    phase        = FALSE,
    verbose      = 1))

# Compare fixed vs. random effects recruitment uncertainty
plot_recruitment(
  Rceattle    = list(ss_run, ss_re),
  model_names = c("Fixed effects", "Random effects"),
  add_ci      = TRUE
)
```

## Available harvest control rules

The `HCR` argument to
[`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
accepts the output of
[`build_hcr()`](https://grantdadams.github.io/Rceattle/reference/build_hcr.md).
See
[`?build_hcr`](https://grantdadams.github.io/Rceattle/reference/build_hcr.md)
for full parameter details.

| Value | String | Description | Multi-species? |
|----|----|----|----|
| 0 | `"NoFishing"` | No fishing — estimate hindcast only | Yes |
| 1 | `"CMSY"` | CMSY: maximise catch across all species (can constrain depletion ≥ `Plimit`) | Yes |
| 2 | `"ConstantF"` | Constant F set at `Ftarget` for each species | Yes |
| 3 | `"ConstantFSSB"` | F that achieves `Ftarget`% of SSB₀ at the end of the projection | Yes |
| 4 | `"ConstantFSPR"` | Constant F_SPR set at `Ftarget`; can be scaled by `Fmult` (NEFSC convention) | No |
| 5 | `"NPFMC"` | NPFMC Tier 3 SPR-based HCR (`Ftarget`, `Flimit`, `Ptarget`, `Alpha`) | No |
| 6 | `"PFMC"` | PFMC Category 1 40-10 ACL rule with uncertainty buffer (`Pstar`, `Sigma`) | Yes |
| 7 | `"SESSF"` | SESSF Tier 1 SPR-based HCR | No |

``` r

# Example: NPFMC Tier 3 control rule (single-species only)
ss_run_tier3 <- fit_mod(
  data_list    = BS2017SS,
  inits        = ss_run$estimated_params,
  file         = NULL,
  estimateMode = 2,
  HCR          = build_hcr(HCR = "NPFMC",
                            Ftarget = 0.40,   # F40%
                            Flimit  = 0.35,   # F35%
                            Ptarget = 0.40,
                            Plimit  = 0.20,
                            Alpha   = 0.05),
  random_rec   = FALSE,
  msmMode      = 0,
  fit_control = fit_control(
    verbose      = 1))
```

For closed-loop testing of control rules across many simulated
trajectories, see
[`vignette("hcrs-and-mses")`](https://grantdadams.github.io/Rceattle/articles/hcrs-and-mses.md).
