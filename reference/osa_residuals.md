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

Survey-index OSA residuals are supported for every index likelihood
family (`Index_distribution`): lognormal IID (`"Lognormal"`) and
natural-scale `"Normal"` residualize as independent normals on the log
or natural scale; the correlated covariance families (`"MVN"` /
`"MVNORM"`) are whitened by the lower Cholesky of the fleet's survey
covariance Sigma = L L' so the residuals are the multivariate-Gaussian
one-step-ahead innovations L^-1 (obs - q\*pred) – the closed form
[`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
reproduces for a Gaussian block.

## Usage

``` r
osa_residuals(
  fit,
  source = c("ecov", "index", "catch", "comp", "caal"),
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
  `"ecov"`, `"index"`, `"catch"`, `"comp"`, `"caal"`, `"diet"`, or
  `"all"`. Defaults to the five non-diet sources (`diet` is opt-in
  because it applies only to multispecies models and can be expensive);
  pass `"all"` to include `diet`. `"ecov"` is the state-space covariate
  (QAR1 `observe=` term), residualized first against its own series, as
  in WHAM's `make_osa_residuals()`. Sources with no observations in the
  model are silently skipped. Mirrors the `source` argument of
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
  [`mclapply`](https://rdrr.io/r/parallel/mclapply.html). Default
  `TRUE`. This is the main speedup for models with random effects, where
  each observation triggers a Laplace re-evaluation – it gives a
  near-linear speedup across cores (set `options(mc.cores = )` to choose
  how many; forking falls back to serial on Windows). Some models –
  heavy random-effect structures such as a time-varying catchability –
  abort the forked worker instead of returning; the loop then recomputes
  serially, after rebuilding, and prints the worker's own "irrecoverable
  exception" message, which comes from C and cannot be suppressed. That
  message does not mean the call failed. Pass `FALSE` to skip the
  attempt. Only the continuous group is parallelized; the discrete
  (randomized-quantile) path always runs serially so it stays
  reproducible given `seed`.

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
`year`, `age_length_bin` (the age or length bin the value stands for),
`accumulated` (`TRUE` where tail accumulation folded neighbouring bins
into this one, so it covers a range of ages rather than the one named),
`length` (the conditioning length bin for caal; `NA` otherwise),
`index_label` (`"age"`/`"length"`/`NA`), `observed`, `predicted`, `sd`,
and `residual`. For aggregate series `observed` and `predicted` are on
the residualization scale – log for lognormal catch/index, natural scale
for a `"Normal"` index, and the whitened (`L^-1`) scale for an
`"MVN"`/`"MVNORM"` index; for compositions they are bin counts. Carries
`method` and `seed` attributes, and (when composition types are present)
a `"pearson"` attribute holding the matching Pearson residuals so
[`plot.rceattle_osa()`](https://grantdadams.github.io/Rceattle/reference/plot.rceattle_osa.md)
can show both. The attribute uses this data frame's column names rather
than the data-sheet names
[`residuals.Rceattle()`](https://grantdadams.github.io/Rceattle/reference/residuals.Rceattle.md)
returns, so the two halves of one object read alike. Note the shared
names do not mean a shared scale: in the attribute `observed` and
`predicted` are proportions summing to one within a fleet-year, with the
sample size in `sample_size`, because composition Pearson residuals are
defined on proportions – not the bin counts the columns above carry. Do
not compare the two directly. Both describe the bins the likelihood fit,
so a fleet with tail accumulation reports the folded window in each –
with one asymmetry: the one-step-ahead decomposition drops each group's
last bin (it is fixed by sum-to-N), and under an *old*-tail accumulation
that dropped bin is the upper accumulated one. Such a fleet therefore
shows its upper boundary bin in the Pearson residuals and not in the OSA
residuals. Summarize it with
[`osa_diagnostics()`](https://grantdadams.github.io/Rceattle/reference/osa_diagnostics.md)
and plot it with
[`plot.rceattle_osa()`](https://grantdadams.github.io/Rceattle/reference/plot.rceattle_osa.md).

## Negative composition `predicted` values

A composition `predicted` is an expected bin count and cannot truly be
negative, but it goes slightly negative where a bin holds almost no
fish. Composition observations enter as counts,
`(proportion + comp_offset) * N`, and
[`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)'s
conditional mean is a numerical step away from the observation, which
overshoots below zero when the count is near it. The function warns,
naming the count and the years.

It is the *bin's* count that drives this, not the composition's sample
size: on EBS pollock the negative rows have a median observed count of
0.05 against 4.9 for the rest, while their sample sizes span the same 1
to 821 as everything else (69 of 353 occur above a sample size of 100).
A rare age in a well-sampled year does it as readily as a poorly sampled
year, so the warning is not by itself evidence of thin data.

The values are reported rather than clamped, because a negative expected
count is the signal that the bin is too sparse for the decomposition to
describe and clamping would hide it. Treat `predicted` on those rows as
uninformative; their `residual` standardises the observation against
that same mean and is therefore biased positive.

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
