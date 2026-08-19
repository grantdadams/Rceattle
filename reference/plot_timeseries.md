# Plot time-series

Function that plots the time-series (SSB/B/R/Depletion) 95% CI trends as
estimated from Rceattle

## Usage

``` r
plot_timeseries(
  Rceattle,
  output = "biomass",
  ylab = NULL,
  file = NULL,
  model_names = NULL,
  line_col = NULL,
  species = NULL,
  spnames = NULL,
  add_ci = FALSE,
  lwd = 3,
  save = FALSE,
  legend.pos = "topright",
  right_adj = 0,
  width = 7,
  height = 6.5,
  minyr = NULL,
  maxyr = NULL,
  incl_proj = FALSE,
  mod_cex = 1,
  lty = rep(1, length(Rceattle)),
  alpha = 0.4,
  mse = FALSE,
  OM = TRUE,
  reference = NULL,
  zero_y = FALSE,
  mod_avg = rep(FALSE, length(Rceattle))
)
```

## Arguments

- Rceattle:

  Single or list of Rceattle model objects exported from `Rceattle`

- output:

  derived quantity of interest: recruitment, biomass, ssb, depletion, or
  ssb_depletion. Uses same name as ".cpp" file.

- ylab:

  Y-axis label. `NULL` (default) derives one from `output` and the
  model's `minage`.

- file:

  name of a file to identified the files exported by the function.

- model_names:

  Names of models to be used in legend

- line_col:

  Colors of models to be used for line color

- species:

  What species to include 1:nspp

- spnames:

  Species names for legend

- add_ci:

  If the confidence interval is to be added (see Details for how it is
  constructed)

- lwd:

  Line width as specified by user

- save:

  Save derived quantity?

- legend.pos:

  Position of the legend as used by
  [`legend`](https://rdrr.io/r/graphics/legend.html) (default =
  "topright").

- right_adj:

  Multiplier for to add to the right side of the figure for fitting the
  legend.

- width:

  Figure width in inches

- height:

  Figure height in inches

- minyr:

  First year to plot

- maxyr:

  max year to plot

- incl_proj:

  TRUE/FALSE, include projection years

- mod_cex:

  Cex of text for model name legend

- lty:

  line type

- alpha:

  shading for confidence intervals

- mse:

  Is an MSE object from
  [`load_mse`](https://grantdadams.github.io/Rceattle/reference/load_mse.md)
  or
  [`run_mse`](https://grantdadams.github.io/Rceattle/reference/run_mse.md)

- OM:

  if mse == TRUE, use the OM (TRUE) or EM (FALSE) for plotting?

- reference:

  Reference model

- zero_y:

  Anchor the y-axis at zero. TRUE for the absolute series (biomass, SSB,
  recruitment), where a truncated axis exaggerates change; FALSE for the
  depletions, which are already on a relative scale.

- mod_avg:

  TRUE/FALSE

## Value

Returns and saves a figure with the population trajectory.

## Units

The model carries numbers-at-age in **thousands** and weight-at-age in
**kg**, so every biomass series (`biomass`, `ssb`,
`exploitable_biomass`) comes out of the model in **mt** and recruitment
comes out in **thousands of fish**. For display these are divided by 1e6
(million mt) and 1e3 (millions of recruits) respectively; depletion is a
ratio and is not rescaled. Supply the model's inputs on that convention
– catch and index in mt, weight-at-age in kg – or the axis labels will
not describe what is plotted.

## Confidence intervals

`add_ci = TRUE` draws a 95% interval. Every strictly positive series
takes it on the log scale as `exp(log(x) +/- qnorm(0.975) * sd_log)`.
These quantities are built multiplicatively, so `log(x)` is close to
linear in the estimated parameters where `x` is exponential in them: the
interval is both a better linearization and right-skewed, and it cannot
cross zero the way a symmetric natural-scale interval does for weak year
classes and depleted stocks.

`sd_log` comes from the model's own `log_biomass` / `log_ssb` / `log_R`
where those are reported, and is otherwise recovered as `sd(x) / x` –
the delta method's own identity, which matches the reported values to
machine precision. That covers `exploitable_biomass` and the two
depletions, which cannot be reported on the log scale
(`exploitable_biomass` is identically 0 without projection F), and
models fit before the `log_*` series existed. A non-positive or
unreported value keeps the symmetric interval.
