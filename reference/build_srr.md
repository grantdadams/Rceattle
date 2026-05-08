# Specify the stock-recruit relationship (SRR) for Rceattle

**Stock recruitment relationships currently implemented in Rceattle:**

- `srr_fun = 0` or `"mean"`: No stock recruit relationship. Recruitment
  is a function of R0 and annual deviates (i.e. steepness = 0.99).
  \$\$R_y = exp(R0 + R\_{dev,y})\$\$

- `srr_fun = 2` or `"BevertonHolt"`: Beverton-holt stock-recruitment
  relationship \$\$R_y = \frac{\alpha\_{srr} \*
  SB\_{y-minage}}{1+\beta\_{srr} \* SB\_{y-minage}}\$\$

- `srr_fun = 4` or `"Ricker"`: Ricker stock-recruitment relationship
  \$\$R_y = \alpha\_{srr} \* SB\_{y-minage} \* exp(-\beta\_{srr} \*
  SB\_{y-minage})\$\$

When `srr_pred_fun > 0` and `srr_fun = 0` recruitment in the hindcast is
estimated as in `srr_fun = 0` \$\$R_y = exp(R0 + R\_{dev,y})\$\$, but an
additional stock recruitment relationship defined by `srr_pred_fun` is
estimated between `srr_hat_styr` and `srr_hat_endyr` and treated as an
additional penalty. The stock recruitment relationship defined by
`srr_pred_fun` is then used in the projection.

## Usage

``` r
build_srr(
  srr_fun = 0,
  srr_pred_fun = srr_fun,
  proj_mean_rec = TRUE,
  srr_mse_switchyr = NULL,
  srr_hat_styr = NULL,
  srr_hat_endyr = NULL,
  srr_est_mode = 1,
  srr_prior = 4,
  srr_prior_sd = 1,
  srr_indices = NA,
  Bmsy_lim = NA,
  linkages = NULL
)
```

## Arguments

- srr_fun:

  Stock recruit function to be used for hindcast estimation of Rceattle
  (see @description below). Default = 0

- srr_pred_fun:

  stock recruit function for projection, reference points, and penalties
  to be used for Rceattle (see below). When `srr_fun == 0`, it treats
  the stock-recruit curve as an additional penalty onto the annualy
  estimated recruitment from the hindcast (sensu AMAK and Jim Ianelli's
  pollock model). If `srr_fun > 0` then `srr_pred_fun = srr_fun` and no
  additional penalty is included.

- proj_mean_rec:

  Project the model using: 0 = mean recruitment (average R of hindcast)
  or 1 = SRR(omega, srr_devs)

- srr_mse_switchyr:

  is used for MSEs to deal with AMAK and Jim Ianelli's estimation where
  a stock recruit function is estimated as an additional penalty
  (srr_fun = 0 and srr_pred_fun \> 0). It tells the model in what year
  to switch to the stock recruit function.

- srr_hat_styr:

  Integer. The first year used for estimating the recruitment function
  as an additional penalty. It will add additional penalties sensu AMAK
  and Jim Ianelli's pollock model when `srr_pred_fun > 0` and
  `srr_fun = 0`, starting at `styr` + 1. Defaults to \$styr + 1\$ in
  \$data_list\$. Useful if environmental data used to condition
  stock-recruit relationships is not available until end-year, but
  projections are desired.

- srr_hat_endyr:

  Integer. The last year used for estimating the recruitment function as
  an additional penalty. It will add an additional penalties sensu AMAK
  and Jim Ianelli's pollock model when `srr_pred_fun > 0` and
  `srr_fun = 0`. Recruitment Defaults to \$endyr\$ in \$data_list\$.
  Useful if environmental data used to condition stock-recruit
  relationships is not available for the full time-series, but
  projections are desired.

- srr_est_mode:

  Switch to determine estimation mode. 0 = fix alpha to prior mean, 1 =
  freely estimate R0, alpha, and/or beta (default), 2 = use lognormally
  distributed prior for alpha (Ricker) or steepness (Beverton), 3 = use
  beta distributed prior for steepness (Beverton) given mean and sd.

- srr_prior:

  mean for normally distributed prior for stock-recruit parameter

- srr_prior_sd:

  Prior standard deviation for stock-recruit parameter

- srr_indices:

  Soft-deprecated. Use the `linkages` argument instead. See
  `vignette("environmental-linkages")`.

- Bmsy_lim:

  Upper limit for Ricker based SSB-MSY (e.g 1/Beta). Will add a
  likelihood penalty if beta is estimated above this limit. Default `NA`
  is not used.

- linkages:

  Optional named list of
  [`linkage_spec()`](https://grantdadams.github.io/Rceattle/reference/linkage_spec.md)
  objects keyed by recruitment parameter name (must be one of
  `"log_R0"`, `"log_alpha"`, `"log_beta"`). Each spec describes how that
  parameter depends on environmental covariates and on stratifying
  factors (species, sex). The offset enters additively (on the log
  scale) inside the recruitment compute. See
  `vignette("environmental-linkages")` for details.

## Value

A `list` containing the stock recruitment relationship settings
