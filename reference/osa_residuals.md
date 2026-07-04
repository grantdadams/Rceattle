# One-step-ahead (OSA) residuals for an Rceattle model

Computes one-step-ahead (OSA) residuals – also called forecast or
quantile residuals (Thygesen et al. 2017) – for a fitted
[Rceattle](https://grantdadams.github.io/Rceattle/reference/Rceattle-package.md)
model via
[`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html).
Unlike Pearson residuals, OSA residuals are distributed iid standard
normal under a correctly specified model even when observations are
correlated (through composition bins) or when the model contains random
effects, so they support objective goodness-of-fit testing (Trijoulet et
al. 2023; Stewart and Monnahan 2025).

These are *internal* OSA residuals: the residualization is integrated
into the assessment via TMB, so it also accounts for correlation induced
by the model's random effects across years – the gold standard relative
to the *external* `compResidual` approach (Stewart and Monnahan 2025).

OSA residuals are computed *post hoc* and are expensive (TMB
re-optimizes the random effects as each observation is added), so they
are not produced during
[`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md).
The model must have been optimized with `estimateMode < 3`.

Supported observation types are the aggregate `"catch"` and `"index"`
series, the `"comp"` (age/length composition) and `"caal"` (conditional
age-at-length) compositions, and `"diet"` (predator stomach-content
composition, for multispecies models with estimated suitability). Diet
is opt-in (not in the default `types`) because it applies only to
multispecies models and can be expensive.

For composition data the multivariate multinomial /
Dirichlet-multinomial is decomposed into a sequence of univariate
conditional residuals (binomial / beta-binomial; Trijoulet et al. 2023).
The final bin of each composition is fixed by the sum-to-N constraint
and so has no residual (returned as `NA`). Composition OSA uses an
internal model rebuild with unweighted, proper densities (the `osa_mode`
switch); fleets fit with the AFSC `MultinomialAFSC` pseudo-likelihood
are residualized under the full multinomial.

## Usage

``` r
osa_residuals(
  fit,
  source = c("index", "catch", "comp", "caal"),
  method = "oneStepGaussianOffMode",
  discrete = FALSE,
  parallel = TRUE,
  seed = 123,
  trace = FALSE,
  ...
)
```

## Arguments

- fit:

  A fitted object of class `Rceattle` (from
  [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)).

- source:

  Character vector of observation sources to residualize: any of
  `"index"`, `"catch"`, `"comp"`, `"caal"`, `"diet"`, or `"all"`.
  Defaults to the four non-diet sources (`diet` is opt-in because it
  applies only to multispecies models and can be expensive); pass
  `"all"` to include `diet`. Sources with no observations in the model
  are silently skipped. Mirrors the `source` argument of
  [`residuals.Rceattle()`](https://grantdadams.github.io/Rceattle/reference/residuals.Rceattle.md)
  and
  [`plot.rceattle_osa()`](https://grantdadams.github.io/Rceattle/reference/plot.rceattle_osa.md).

- method:

  Passed to
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html).
  Defaults to `"oneStepGaussianOffMode"` (the WHAM/SAM default),
  appropriate for the lognormal aggregate series.

- discrete:

  Logical; whether to treat *composition* (comp / caal / diet)
  observations as discrete. Default `FALSE` (continuous, matching how
  CEATTLE fits the composition likelihood with
  effective-sample-size-scaled counts). When `TRUE`, composition
  residuals are randomized quantile residuals (Dunn and Smyth 1996) and
  so are stochastic; set `seed` for reproducibility. The aggregate
  index/catch series are always continuous (lognormal); the
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  call is split by observation type so `discrete` is applied correctly
  per type (the discrete group uses the generic CDF-based method rather
  than `method`).

- parallel:

  Logical; compute the per-observation OSA loop in parallel via
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html).
  Default `TRUE`. This is the main speedup for models with random
  effects, where each observation triggers a Laplace re-evaluation – it
  gives a near-linear speedup across cores (set `options(mc.cores = )`
  to choose how many; forking falls back to serial on Windows). Only the
  continuous group is parallelized; the discrete (randomized-quantile)
  path always runs serially so it stays reproducible given `seed`.

- seed:

  Random seed passed to
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  for reproducibility of randomized-quantile residuals. Default `123`.

- trace:

  Logical; print
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  progress. Default `FALSE`.

- ...:

  Further arguments passed to
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html).

## Value

A data frame of class `rceattle_osa` with one row per residualized
observation and columns `source` (the data source:
index/catch/comp/caal/ diet), `fleet`, `fleet_name`, `species`, `sex`,
`year`, `age_length_bin` (the age or length bin index), `length` (the
conditioning length bin for caal; `NA` otherwise), `index_label`
(`"age"`/`"length"`/`NA`), `observed`, `predicted`, `sd`, and
`residual`. For aggregate series `observed` and `predicted` are on the
model (log) scale; for compositions they are bin counts. Carries
`method` and `seed` attributes, and (when composition types are present)
a `"pearson"` attribute holding the matching Pearson residuals so
[`plot.rceattle_osa()`](https://grantdadams.github.io/Rceattle/reference/plot.rceattle_osa.md)
can show both. Summarize it with
[`osa_diagnostics()`](https://grantdadams.github.io/Rceattle/reference/osa_diagnostics.md)
and plot it with
[`plot.rceattle_osa()`](https://grantdadams.github.io/Rceattle/reference/plot.rceattle_osa.md).

## References

Thygesen, U.H., et al. 2017. Validation of ecological state space models
using the Laplace approximation. Environ. Ecol. Stat. 24:317-339.

Trijoulet, V., et al. 2023. Model validation for compositional data in
stock assessment models. Fish. Res. 257:106487.

Stewart, I.J., and Monnahan, C.C. 2025. Diagnosing common sources of
lack of fit to composition data using one-step-ahead residuals. Can. J.
Fish. Aquat. Sci. 82:1-13.

## See also

[`osa_diagnostics()`](https://grantdadams.github.io/Rceattle/reference/osa_diagnostics.md),
[`plot.rceattle_osa()`](https://grantdadams.github.io/Rceattle/reference/plot.rceattle_osa.md),
[`process_residuals()`](https://grantdadams.github.io/Rceattle/reference/process_residuals.md)
