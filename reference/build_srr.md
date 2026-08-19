# Specify the stock-recruit relationship (SRR) for Rceattle

**Stock recruitment relationships currently implemented in Rceattle:**

- `srr_fun = 0` or `"mean"`: No stock recruit relationship. Recruitment
  is a function of \\R0\\ (on the log scale) and annual deviates (i.e.
  steepness = 0.99). \$\$R_y = exp(R0 + R\_{dev,y})\$\$

- `srr_fun = 2` or `"BevertonHolt"`: Beverton-holt stock-recruitment
  relationship \$\$R_y = \frac{\alpha\_{srr} \*
  SB\_{y-minage}}{1+\beta\_{srr} \* SB\_{y-minage}}\$\$

- `srr_fun = 4` or `"Ricker"`: Ricker stock-recruitment relationship
  \$\$R_y = \alpha\_{srr} \* SB\_{y-minage} \* exp(-\beta\_{srr} \*
  SB\_{y-minage})\$\$

The Beverton-Holt and Ricker curves above are the deterministic mean;
realized recruitment applies the annual log deviation, \\R_y \cdot
exp(R\_{dev,y})\\, as in the mean form. For numerical stability the
Ricker \\\beta\_{srr}\\ is estimated on a scale divided by 1,000,000, so
the fitted `beta` is 1e6 times the density-dependence coefficient in the
equation above; `Bmsy_lim` (\\\approx 1/\beta\_{srr}\\) carries the same
scaling.

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
  srr_alpha_init = NULL,
  srr_beta_init = NULL,
  srr_indices = NA,
  Bmsy_lim = NA,
  linkages = NULL
)
```

## Arguments

- srr_fun:

  Stock-recruit function used in the hindcast estimation (see the list
  below). Default = 0

- srr_pred_fun:

  Stock-recruit function used for projection, reference points, and
  penalties (see below). When `srr_fun == 0`, the stock-recruit curve is
  added as a penalty on the annually estimated hindcast recruitment
  (following AMAK and Jim Ianelli's pollock model). If `srr_fun > 0`,
  then `srr_pred_fun = srr_fun` and no extra penalty is added.

- proj_mean_rec:

  Project the model using: 0 = mean recruitment (average R of hindcast)
  or 1 = SRR(omega, srr_devs)

- srr_mse_switchyr:

  Year at which an MSE switches from the annual recruitment-penalty
  estimate to the stock-recruit function (the `srr_fun = 0`,
  `srr_pred_fun > 0` case).

- srr_hat_styr:

  Integer. First year used to estimate the recruitment-penalty function
  (the AMAK/Ianelli penalty, active when `srr_pred_fun > 0` and
  `srr_fun = 0`), starting at `styr + 1`. Defaults to `styr + 1` in
  `data_list`. Useful when the environmental data conditioning the
  stock-recruit relationship is not available until the terminal year
  but projections are still wanted.

- srr_hat_endyr:

  Integer. Last year used to estimate the recruitment-penalty function
  (the AMAK/Ianelli penalty, active when `srr_pred_fun > 0` and
  `srr_fun = 0`). Defaults to `endyr` in `data_list`. Useful when the
  environmental data conditioning the stock-recruit relationship does
  not span the full time series but projections are still wanted.

- srr_est_mode:

  Switch to determine estimation mode. Accepts integer codes or the
  equivalent readable strings: 0 / "Fixed" = fix alpha to prior mean; 1
  / "Estimated" = freely estimate R0, alpha, and/or beta (default); 2 /
  "LognormalPrior" = lognormally distributed prior for alpha (Ricker) or
  steepness (Beverton); 3 / "BetaPrior" = beta distributed prior for
  steepness (Beverton) given mean and sd.

- srr_prior:

  mean for normally distributed prior for stock-recruit parameter. Note
  this is not the same quantity for every curve: for Ricker
  (`srr_pred_fun` 4 / 5) it is a prior on \\\alpha\\ itself, while for
  Beverton-Holt (2 / 3) it is a prior on **steepness** and so must lie
  in (0, 1). See `srr_est_mode`.

- srr_prior_sd:

  Prior standard deviation for stock-recruit parameter

- srr_alpha_init, srr_beta_init:

  Optional starting values for the stock-recruit \\\alpha\\ and
  \\\beta\\ parameters, on the natural (not log) scale, one value per
  species. Only used when the curve actually estimates them (`srr_fun`
  or `srr_pred_fun` above 1); ignored for mean recruitment, where both
  are mapped out.

  The package defaults (\\\alpha = e^3\\, \\\beta = 3\\) are
  placeholders with no knowledge of the stock's scale. \\\beta\\ in
  particular sets the density dependence in \\R = \alpha S / (1 + \beta
  S)\\, so it must be on the order of \\(\alpha - 1/\phi_0) / R_0\\ –
  typically \\10^{-3}\\ or smaller for a stock measured in tonnes.
  Starting three orders of magnitude away drives predicted recruitment
  to near zero and the optimizer returns `NA/NaN gradient evaluation`.
  For a Beverton-Holt seeded from a steepness \\h\\ and unfished
  spawning biomass per recruit \\\phi_0\\: \$\$\alpha = \frac{4h}{\phi_0
  (1 - h)}, \qquad \beta = \frac{\alpha - 1/\phi_0}{R_0}.\$\$

- srr_indices:

  Soft-deprecated. Use the `linkages` argument instead. See
  [`vignette("environmental-linkages-and-priors")`](https://grantdadams.github.io/Rceattle/articles/environmental-linkages-and-priors.md).

- Bmsy_lim:

  Upper limit for Ricker based SSB-MSY (e.g 1/Beta). Will add a
  likelihood penalty if beta is estimated above this limit. Default `NA`
  is not used.

- linkages:

  Optional named list of
  [`linkage_spec()`](https://grantdadams.github.io/Rceattle/reference/linkage_spec.md)
  objects keyed by recruitment parameter name (must be one of `"R0"`,
  `"alpha"`, `"beta"`). Each spec describes how that parameter depends
  on environmental covariates and on stratifying factors (species, sex).
  The offset enters additively (on the log scale) inside the recruitment
  compute. See
  [`vignette("environmental-linkages-and-priors")`](https://grantdadams.github.io/Rceattle/articles/environmental-linkages-and-priors.md)
  for details.

## Value

A `list` containing the stock recruitment relationship settings
