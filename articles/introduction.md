# 1. Rceattle: An Introduction

## Overview

In this vignette we walk through fitting three models using the
`Rceattle` package to run single- and multi-species stock assessment
models. Rceattle is based on Holsman et al. (2016) and Adams et
al. (2022) to fit single- and multi-species age structured models in
TMB. The code base has been inspired by VAST, WHAM, and other TMB based
packages and generalized for flexibly fitting models with a variable
number of species, surveys, and fishing fleets. In this example we apply
Rceattle to walleye pollock, arrowtooth flounder, and Pacific cod in the
Gulf of Alaska and Eastern Bering Sea. All functions have associated
helper text using the `?`. Here, we demonstrate the basic workflow:

1.  Install `Rceattle`

2.  Load `Rceattle` and associated excel data

3.  Fit assessment models:

- model 1: a single-species assessment model estimated with recruitment
  estimated as penalized effects;

- model 2: three single-species assessment models jointly estimated with
  recruitment estimated as penalized effects;

- model 3: as model 2, but with natural mortality estimated for each
  species;

- model 4: as model 3, but with recruitment estimated as random effects;

- model 5: full multi-species assessment models with time- and age-
  varying natural mortality due to predation.

4.  Compare models by AIC and Mohn’s rho (retrospective analysis)

5.  Review plots of input data, diagnostics, and results.

## 1. Installation

`Rceattle` and associated dependencies and can be installed as follows:

``` r

# Rceattle (pulls CRAN dependencies automatically; e.g. TMB, Matrix, dplyr)
install.packages("remotes")
remotes::install_github("grantdadams/Rceattle")

# Optional: TMBhelper provides richer optimization diagnostics.
# Rceattle falls back to plain nlminb + sdreport if it's not installed.
# devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
```

## 2. Data

To run `Rceattle` a data file must first be loaded. The data files from
2018 Gulf of Alaska and 2017 Bering Sea single- and multi-species
assessments is already loaded within the `Rceattle` package:

``` r

library(Rceattle)
data("GOApollock") # Data for GOA-pollock (1-species) ?GOApollock for more information on the data 
data("BS2017SS") # Combined data for 3 species. ?BS2017SS for more information on the data 
data("BS2017MS") # Multi-species data. Note: the only difference is the residual mortality (M1_base) is lower
```

Data can be visualized using the `plot_data` function:

``` r

plot_data(GOApollock)
```

The data files can be written to excel, modified within excel, and read
back into `R` as follows:

``` r

Rceattle::write_data(data_list = BS2017SS, file = "BS2017SS.xlsx")
mydata <- Rceattle::read_data( file = "BS2017SS.xlsx")
```

A description of the inputs and file structure can be found on the first
sheet (`meta_data`) and the associated help function
[`?BS2017SS`](https://grantdadams.github.io/Rceattle/reference/BS2017SS.md).

## 3. Fit models

Models can be fit using the `fit_mod` function. For the first model we
use fit a single-species model with no stock-recruit to data on Gulf of
Alaska pollock. The model has 8 fleets: 1 fishing fleet, 6 surveys, and
1 fleet turned off and all data are excluded from the likelihood (see
`model_1$fleet_control`).

``` r

model_1 <- Rceattle::fit_mod(data_list = GOApollock,
                            inits = NULL, # Initial parameters = 0
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            fit_control = fit_control(
                              phase = TRUE,
                              verbose = 1))
summary(model_1)
```

For the second model we use fit three single-species model jointly with
no stock-recruit curve and natural mortality is fixed as the input
values (`mydata$M1_base`). The model has 7 fleets: 3 fishing fleets that
target each species individually and 4 surveys (see
`mydata$fleet_control`). The single-species model can be fit by the
default setting `msmMode = 0` using the `fit_mod` function:

``` r

model_2 <- Rceattle::fit_mod(data_list = mydata,
                            inits = NULL, # Initial parameters = 0
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Penalized likelihood
                            fit_control = fit_control(
                              phase = TRUE,
                              verbose = 1))
summary(model_2)
```

Alternatively, we can estimate natural mortality using the `build_M1`
function, where setting `M1_model = 1` estimates a single M for each
species. Here, we also begin estimation at the maximum likelihood
estimates from `model_1` using the `inits` argument:

``` r

model_3 <- Rceattle::fit_mod(data_list = mydata,
                              estimateMode = 0, # Estimate
                              M1Fun = build_M1(M1_model = "sex_age_invariant"), # Estimate sex-invariant M)
                              random_rec = FALSE, # Penalized likelihood
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(
                                phase = TRUE,
                                verbose = 1))
summary(model_3)
```

To estimate `model_4` where recruitment is treated as random effects we
can use the `random_rec` argument. For this model we initialized at the
previous model’s MLEs and don’t phase estimation to decrease run-time:

``` r

model_4 <- Rceattle::fit_mod(data_list = mydata,
                              inits = model_3$estimated_params, # Start estimation at model 3's MLEs
                              estimateMode = 0, # Estimate
                              M1Fun = build_M1(M1_model = "sex_age_invariant"),
                              random_rec = TRUE, # Random recruitment
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(
                                phase = FALSE,
                                verbose = 1))
summary(model_4)
```

To estimate time- and age-varying natural mortality due to predation we
can estimate the model in multi-species mode by setting `msmMode = 1`.
This follows the MSVPA parameterization described in [Magnusson
(1995)](https://link.springer.com/article/10.1007/BF00179756).

``` r

model_5 <- Rceattle::fit_mod(data_list = BS2017MS,
                            inits = model_3$estimated_params, # Initial parameters from single species estimates
                            M1Fun = build_M1(M1_model = "sex_age_invariant"),
                            estimateMode = 0, # Estimate
                            niter = 3, # 3 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # Empirical suitability
                            fit_control = fit_control(
                              verbose = 1))
summary(model_5)
```

## 4. Compare models

We can compare models by evaluating the marginal or joint negative
log-likelihood, AIC, and retrospective analysis. NOTE, we can not
compare models estimated using penalized likelihood vs marginal
likelihood using AIC. The likelihood based components all can be found
in the model objects build above:

``` r

# Evaluate 1 model
model_1$opt$AIC # AIC
model_1$opt$objective # Negative log-likelihood
model_1$quantities$jnll_comp # Negative log-likelihood components

# Compare AIC across 3-species penalized likelihood models
sapply(list(model_2, model_3, model_5), function(x) x$opt$AIC)
```

We can also look at retrospective bias in each assessment using the
`retrospective` function, which creates a list of model objects and a
data.frame with mohn’s rho calculations.

``` r

model_1_retro <- retrospective(model_1, peels = 7)

plot_biomass(model_1_retro$Rceattle_list, model_names = paste("Pollock; mohns =", round(model_1_retro$mohns[1,2], 3)))
```

## 5. Model plots and diagnostics

Various plots and diagnostics are built into `Rceattle`. Specifically,
time-series hindcast and forecast biomass, SSB, recruitment, biomass
consumed as prey, mortality-at-age, and depletion with uncertainty. The
time-series plots can accept a list of models (see above).

``` r

model_list <- list(model_2, model_3, model_4, model_5)
model_names <- list("single-spp", "Single-spp M", "Single-spp M - RE", "Multi-spp")

# Biomass
plot_biomass(model_list, model_names = model_names, incl_proj = TRUE, add_ci = TRUE)
```

``` r

# Recruitment
plot_recruitment(model_list, model_names = model_names, incl_proj = TRUE, add_ci = TRUE)

# SSB
plot_ssb(model_list, model_names = model_names, incl_proj = TRUE, add_ci = TRUE)

# Biomass depletion
plot_depletion(model_list, model_names = model_names, incl_proj = FALSE, add_ci = FALSE)

# SSB depletion
plot_depletion(model_list, model_names = model_names, incl_proj = FALSE, add_ci = FALSE)
```

Diagnostics plots include plotting selectivity, fits to composition
data, index fit, and catch fit. Selectivity and composition plot
functions can only input one model. Index and catch fit plots can take
multiple model inputs.

``` r

# Selectivity
plot_selectivity(model_2)

# Plot composition fit
plot_comp(model_2)

# Plot index
plot_index(model_list) # - All models

# Plot catch data
plot_catch(model_list)
```

Outputs of the model can queried using the `$` sign:

- `model_1$initial_params` - list of initial parameter values prior to
  estimation

- `model_1$bounds` - list of parameter bounds

- `model_1$map` - list of parameter map in factor form
  (`model_1$map$mapFactor`) or list form (`model_1$map$mapList`)

- `model_1$estimated_params` - list of estimated parameters post
  estimation

- `model_1$quantities` - list of derived quantities (SSB, biomass,
  selectivity, etc)

- `model_1$data_list` - the input data

- `model_1$run_time` - model run time

- `model_1$obj` - the TMB makeadfun object.

- `model_1$opt` - the fit_tmb optimization object.

- `model_1$sdrep` - the TMB sdreport object.
