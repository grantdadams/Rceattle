# 2. HCRs and MSEs: an introduction

### Fitting OMs

An example MSE using the northern rockfish stock assessment.

``` r

library(Rceattle)
data("NorthernRockfish2022")
nrdata <- NorthernRockfish2022
nrdata$fleet_control$proj_F_prop <- 1

om <- Rceattle::fit_mod(
  data_list = nrdata,
  estimateMode = 0,  # Estimate
  initMode = 1,      # Assume unfished equilibrium
  M1Fun = build_M1(
    M1_model = "sex_age_invariant",
    linkages = list(
      M1 = linkage_spec(
        formula = ~ 1,
        priors = list(
          "(Intercept)" = lognormal(log(0.06), 0.05)
        )
      )
    )
  ),
  fit_control = fit_control(phase = TRUE)
)
```

------------------------------------------------------------------------

### Fitting EMs (Tier 3 HCR)

We can use the `HCR = build_hcr()` argument in `fit_mod` to specify the
HCR and BRPs of our management strategy. A variety of HCRs and BRPs are
included and can be found by
[`?build_hcr`](https://grantdadams.github.io/Rceattle/reference/build_hcr.md).

``` r

em <- Rceattle::fit_mod(
  data_list = nrdata,
  estimateMode = 0, # Estimate
  initMode = 1,      # Assume unfished equilibrium
  M1Fun = build_M1(
    M1_model = "sex_age_invariant",
    linkages = list(
      M1 = linkage_spec(
        formula = ~ 1,
        priors = list(
          "(Intercept)" = lognormal(log(0.06), 0.05)
        )
      )
    )
  ),
  HCR = build_hcr(HCR = "NPFMC", # Tier3 HCR
                  Ftarget = 0.4, # F40%
                  Flimit = 0.35, # F35%
                  Plimit = 0.2,  # No fishing when SB<SB20
                  Alpha = 0.05),
  fit_control = fit_control(verbose = 1, phase = TRUE)
)
```

------------------------------------------------------------------------

### Running MSEs

``` r

mse1 <- run_mse(
  om = om, em = em,
  nsim = 10,
  assessment_period = 2,
  # Fishery = fleet 1 (annual), BTS = fleet 2 (biannual)
  sampling_period = c(1, 2),
  seed = 666 # default; per-sim seed = seed + sim
)
```

#### Reproducibility

[`run_mse()`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)
accepts two seed arguments, both with default `seed = 666`:

- `seed` controls the random draws *inside* each simulated trajectory
  (recruitment, observation error, composition resampling). Each
  simulation worker calls `set.seed(seed + sim)`, so simulation `1` uses
  seed `667`, simulation `2` uses `668`, and so on. Re-running the same
  [`run_mse()`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)
  call with the same `seed` reproduces the trajectory bit-for-bit.
- `regenerate_seed` (defaults to `seed`) controls the optional
  `regenerate_past = TRUE` step that re-fits the EM to OM-simulated
  historical data. Set this separately if you want to vary past-data
  realizations while holding the projection draws fixed.

The MSE uses a PSOCK cluster
([`parallel::parLapply`](https://rdrr.io/r/parallel/clusterApply.html))
by default, so worker processes are clean and only see the explicit
objects shipped by `clusterExport()` — there is no leakage from the
parent’s `.GlobalEnv` RNG state. As a result the seed plumbing above is
the *only* thing that controls reproducibility: there is no need to call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before
[`run_mse()`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)
for the simulation results to be deterministic. (You may still want
[`set.seed()`](https://rdrr.io/r/base/Random.html) for any plotting code
that draws random colours or jitter.)

For management-facing runs, record `packageVersion("Rceattle")` and the
`seed` alongside the saved `.rds` so a later reviewer can rerun the
exact same trajectories from the same package tag.

------------------------------------------------------------------------

### Evaluating MSEs via plots

``` r

plot_depletionSSB(lapply(mse1, function(x) x$OM))
```

``` r

plot_depletionSSB(
  c(mse1$Sim_1$EM,           # EMs
    list(mse1$Sim_1$OM,      # OM
         mse1$Sim_1$OM_no_F  # OM w/ no fishing
    )),    
  line_col = c(
    paste0("grey", round(seq(80, 50, length.out = 15))),
    1, 3
  )
)
```

------------------------------------------------------------------------

### Summary statistics

``` r

mse_summ <- mse_summary(mse1)
knitr::kable(mse_summ, digits = 2)
```

------------------------------------------------------------------------

## Alternative scenarios

------------------------------------------------------------------------

### Different BRPs

``` r

em2 <- Rceattle::fit_mod(
  data_list = nrdata,
  estimateMode = 0, # Estimate
  initMode = 1,      # Assume unfished equilibrium
  M1Fun = build_M1(
    M1_model = "sex_age_invariant",
    linkages = list(
      M1 = linkage_spec(
        formula = ~ 1,
        priors = list(
          "(Intercept)" = lognormal(log(0.06), 0.05)
        )
      )
    )
  ),
  HCR = build_hcr(HCR = "NPFMC", # Tier3 HCR
                  Ftarget = 0.30,# F40%
                  Flimit = 0.10, # F35%
                  Plimit = 0,    # No limit
                  Alpha = 0.1),
  fit_control = fit_control(verbose = 1, phase = TRUE)
)

mse2 <- run_mse(
  om = om, em = em2, 
  nsim = 10, 
  assessment_period = 2,
  sampling_period = c(1, 2) # Fishery = fleet 1, BTS = fleet 2
)
```

------------------------------------------------------------------------

### Triannual BTS

``` r

mse3 <- run_mse(
  om = om, em = em, 
  nsim = 10, 
  assessment_period = 2,
  sampling_period = c(1, 3) # Fishery = fleet 1, BTS = fleet 2
)
```

------------------------------------------------------------------------

### 4-year cycle and 50% sampling

``` r

mse4 <- run_mse(
  om = om, em = em, 
  nsim = 10, 
  fut_sample = 0.5,
  assessment_period = 4,
  sampling_period = c(1, 2) # Fishery = fleet 1, BTS = fleet 2
)
```

------------------------------------------------------------------------

### Climate MSE

#### Climate linked OM

``` r


# Generate random env data
projyrs <- nrdata$styr:nrdata$projyr
projnyrs <- length(projyrs)

# Env data
nrdata$env_data <- merge(nrdata$env_data, data.frame(Year = projyrs, EnvIndex = seq(0,1, length.out = projnyrs)))

# Fit OM
om2 <- Rceattle::fit_mod(
  data_list = nrdata,
  estimateMode = 0,  # Estimate
  initMode = 1,      # Assume unfished equilibrium
  M1Fun = build_M1(
    M1_model = "sex_age_invariant",
    linkages = list(
      M1 = linkage_spec(
        formula = ~ 1,
        priors = list(
          "(Intercept)" = lognormal(log(0.06), 0.05)
        )
      )
    )
  ),
  recFun = build_srr(
    srr_fun = "mean",    # Mean R with covariates
    proj_mean_rec = FALSE,
    linkages = list(
      R0 = linkage_spec(formula = ~ EnvIndex)
    )),
  fit_control = fit_control(phase = TRUE)
)

plot_biomass(om2, incl_proj = TRUE)
```

#### Run MSE

``` r

mse5 <- run_mse(
  om = om2, em = em, 
  nsim = 10, 
  fut_sample = 0.5,
  assessment_period = 4,
  sampling_period = c(1, 2) # Fishery = fleet 1, BTS = fleet 2
)
```

------------------------------------------------------------------------

## Overall summary

``` r

mse_summ$MSE = 1
mse_summ2 <- mse_summary(mse2); mse_summ2$MSE = 2

knitr::kable(do.call("rbind", list(mse_summ, mse_summ2)) |>
               dplyr::relocate(MSE), digits = 2)
```
