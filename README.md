# Rceattle

Rceattle: an R package for fitting and testing climate-linked, single- and multi-species age-structured stock assessment models via diagnostics, simulation, and management strategy evaluation.

`Rceattle` implements models in R using Template Model Builder (`TMB`; Kristensen et al., 2015). Data can be supplied via the bundled Excel template or constructed programmatically as an R list (see *Building a data object without Excel* vignette). Capabilities include:

- **Single-species** (`msmMode = 0`) and **multispecies** (`msmMode > 0`) configurations, with one- or two-sex population dynamics
- **One or multiple stocks** can be fit jointly.
- **Multiple fisheries and surveys** with flexible **catchability** and **selectivity** parameterizations (logistic, double-logistic, non-parametric, time-varying)
- **Stock–recruitment** options (Beverton–Holt, Ricker, mean-recruitment, environmentally-driven)
- **Estimable growth** (von Bertalanffy / empirical weight-at-age)
- **Environmental linkages and priors** on recruitment, natural mortality, and growth
- **Bioenergetics-based predation** with temperature-specific consumption (multispecies mode)
- **Forward projections** under alternative **harvest control rules** and climate scenarios
- **Closed-loop MSE**, **retrospective**, **jitter**, and **simulation testing**
- **Tidy outputs** via S3 methods (`as.data.frame`, `coef`, `logLik`, `vcov`, `residuals`, `plot`)

See `browseVignettes("Rceattle")` or the [package website](https://grantdadams.github.io/Rceattle/) for full documentation.


**Installation**

```r
# Rceattle (pulls CRAN dependencies automatically)
install.packages("remotes")
remotes::install_github("grantdadams/Rceattle")

# Optional: TMBhelper provides richer optimization diagnostics.
# Rceattle falls back to plain nlminb + sdreport if it's not installed.
# devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
```

For operational / management use, pin a specific tagged release rather than
tracking `main`, e.g.:

```r
devtools::install_github("grantdadams/Rceattle@v4.3.0")
```

The maintainer email in `DESCRIPTION` (`grant.adams@noaa.gov`) is the
preferred contact for questions about the package; please open an
[issue](https://github.com/grantdadams/Rceattle/issues) for bug reports
and feature requests so the discussion is publicly searchable.

**Getting started**

Fit a single-species model to the bundled Bering Sea data, inspect the
fit, and plot it:

```r
library(Rceattle)

# Bundled Eastern Bering Sea single-species data (1979-2017)
data(BS2017SS)

# Fit the hindcast (single-species mode)
fit <- fit_mod(
  data_list    = BS2017SS,
  msmMode      = 0,                              # 0 = single-species
  estimateMode = 0,                              # 0 = fit hindcast + projection
  fit_control  = fit_control(phase = TRUE)       # optimizer / sdreport / phasing knobs
)

# S3 methods: the fit behaves like an R model object
fit                       # compact print
summary(fit)              # same compact summary
coef(fit)                 # estimated fixed-effect parameters
logLik(fit)               # logLik with df attribute (AIC works)
AIC(fit)
vcov(fit)                 # fixed-effect covariance from sdreport
as.data.frame(fit)        # return data.frame of derived quantities

# Plot dispatcher: pick a panel with `what`
plot(fit, what = "biomass")
plot(fit, what = "ssb")
plot(fit, what = "recruitment")
plot(fit, what = "index")

# Residuals across data sources (long-format data frame)
ix_resid   <- residuals(fit, type = "index")    # log-scale by default
comp_resid <- residuals(fit, type = "comp")     # Pearson, with Age0_Length1
all_resid  <- residuals(fit, type = "all")
```

For a multispecies fit, swap in `BS2017MS` and set `msmMode = 1`. See the
package vignettes (`browseVignettes("Rceattle")`) for projections,
alternative harvest control rules, MSE testing, model diagnostics, and
non-Excel data construction. For deeper context, see the
[onboarding document](https://github.com/grantdadams/Rceattle/wiki/Onboarding)
and Wiki. The model can be updated following instructions
[here](https://github.com/grantdadams/Rceattle/wiki/Workflow-for-updating-the-Rceattle).

**Examples**
Additional code and function examples using data from the Bering Sea and Gulf of Alaska groundfish application can be found in the [examples](https://github.com/grantdadams/Rceattle/tree/master/examples) folder and include:
* [Fitting single-species models](https://github.com/grantdadams/Rceattle/blob/master/examples/Fit_2018_GOA_single-species_models.R)
* [Fitting multi-species models](https://github.com/grantdadams/Rceattle/blob/master/examples/Fit_2018_GOA_multi-species_model.R)
* [Alternative HCRs and MSE testing](https://github.com/grantdadams/Rceattle/blob/master/examples/HCRs_and_MSE_testing.R)
* [Simulation](https://github.com/grantdadams/Rceattle/blob/master/examples/Simulation_testing.R)
* [Model diagnostics](https://github.com/grantdadams/Rceattle/blob/master/examples/Model_diagnostics.R)

**References**

Adams, G.D., Holsman, K.K., Barbeaux, S.J., Dorn, M.W., Ianelli, J.N., Spies, I., Stewart, I.J., Punt, A.E., 2022. An ensemble approach to understand predation mortality for groundfish in the Gulf of Alaska. Fish. Res. 251, 106303. doi:10.1016/j.fishres.2022.106303

Holsman, K. K., Ianelli, J., Aydin, K., Punt, A. E., & Moffitt, E. A. (2016). A comparison of fisheries biological reference points estimated from temperature-specific multi-species and single-species climate-enhanced stock assessment models. Deep Sea Research Part II: Topical Studies in Oceanography, 134, 360-378.

Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., & Bell, B. (2015). TMB: automatic differentiation and Laplace approximation. arXiv preprint arXiv:1509.00660.

Wassermann, S. N., Adams, G. D., Haltuch, M. A., Kaplan, I. C., Marshall, K. N., & Punt, A. E. (2025). Even low levels of cannibalism can bias population estimates for Pacific hake. ICES Journal of Marine Science, 82(1), fsae064.
