# Self-test simulation analysis

Simulates data from an Rceattle model and refits the model to the
simulated data, to check that the fitting procedure recovers the
operating-model parameters. TODO add process variation (i.e. random
recruitment deviations) to the simulation.

## Usage

``` r
self_test(
  Rceattle = NULL,
  nsim = 50,
  simulate = TRUE,
  seed = 123,
  cores = NULL,
  getsd = NULL,
  phase = NULL,
  start = c("initial", "estimated"),
  debug = FALSE,
  timeout = Inf
)
```

## Arguments

- Rceattle:

  an Rceattle model fit using
  [`fit_mod`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)

- nsim:

  number of simulations

- simulate:

  passed to
  [`sim_mod`](https://grantdadams.github.io/Rceattle/reference/sim_mod.md).
  If `TRUE` (default), data are simulated with observation error; if
  `FALSE`, expected values from the model are used.

- seed:

  random number seed. Each simulation `i` uses `seed + i` so results are
  reproducible under both sequential and parallel execution.

- cores:

  Number of cores to use for parallel simulations. Default `NULL` picks
  `parallel::detectCores() - 6`, capped at 2 when running under
  `R CMD check` (which sets `_R_CHECK_LIMIT_CORES_`). Set to 1 to force
  sequential execution.

- getsd:

  whether each refit runs
  [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html).
  Self-test compares the refit point estimates to the operating model,
  so `FALSE` is faster with no effect on that comparison. Default `NULL`
  inherits the input model's setting (`TRUE` only if it carries an
  `sdrep`).

- phase:

  as in
  [`fit_mod`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md).
  Under the default `start = "initial"` each refit covers the same
  ground the original fit did, so a model that needed phasing to fit its
  real data needs it again for every simulated one – without it such a
  model's refits can end many orders of magnitude from a zero gradient
  and be dropped as non-converged. Default `NULL` reads the setting
  [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
  recorded on the source fit (`fit$run_config$fit_control$phase`), so a
  model fitted under the package default of `phase = FALSE` is refitted
  unphased; pass `TRUE` for a model that needs phasing but was not
  fitted with it.

- start:

  which of the input model's parameter sets each refit starts from.
  `"initial"` (default) uses `initial_params`, the values the original
  fit itself started from, so the estimator has to travel the same
  distance to the optimum on simulated data that it did on the real
  data. `"estimated"` starts from `estimated_params` instead: much
  faster and far more likely to converge, but the fixed effects – and,
  with `random_rec = TRUE`, the inner Laplace problem too – begin at the
  generating values, so on a multimodal or weakly identified surface the
  optimizer never leaves the basin containing them and recovery is close
  to guaranteed by construction. Read it as optimistic about recovery,
  not merely less powerful. (Nor is it a complete warm start:
  [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
  resets `log_Ftarget`, `proj_F_prop`, and the stock-recruit
  \\\alpha\\/\\\beta\\ from the model's own specification under either
  setting.) Non-identifiability shows up in the curvature and so is
  visible either way – via `$convergence`'s Hessian conditioning and
  estimability checks – it is *reachability* that a warm start stops
  testing.

- debug:

  return every simulation rather than the converged ones. The dropped
  runs are the interesting ones when a self-test comes back short, and
  each carries its own `$convergence` diagnostics. See **Value**.

- timeout:

  elapsed-second limit per simulation, `Inf` (default) for none. The
  optimizer runs with no iteration cap, so a replicate that wanders
  somewhere pathological can stall the whole run – a hang that no
  convergence check can catch, because the fit never returns. One that
  exceeds the limit is stopped, counted as non-converged and reported
  separately. Approximate: the limit is checked when control returns to
  R, so it fires between the optimizer's function evaluations rather
  than inside one.

## Value

A list of Rceattle models named `Sim_1`, `Sim_2`, .... By default only
the converged simulations, renumbered contiguously; a message reports
how many were dropped.

With `debug = TRUE`, every simulation, with `Sim_i` being simulation `i`
(so it pairs with the seed `seed + i`), and a logical vector of the
convergence verdicts in `attr(, "converged")`. Inspect a failure with
`sims[[j]]$convergence`. A simulation whose refit errored outright is
returned as the condition object rather than a model, so it cannot abort
the run.

## Interpreting the spread

[`sim_mod`](https://grantdadams.github.io/Rceattle/reference/sim_mod.md)
redraws the observations only – indices, catch, compositions and CAAL.
It does not redraw recruitment, so with `random_rec = TRUE` every
replicate shares the operating model's single recruitment realization,
and that realization is its shrunk empirical-Bayes modes rather than a
draw from N(0, sigmaR). Two consequences: the spread across replicates
carries observation error only and is a lower bound on estimation
uncertainty in SSB and recruitment (do not read it against the model's
own uncertainty bands, which include process error); and sigmaR is
re-estimated from deviations that were shrunk toward zero the same way
in every replicate, a downward bias that averaging over simulations does
not remove.

## Examples

``` r
# \donttest{
data(BS2017SS)
ss_run <- fit_mod(data_list = BS2017SS,
    inits = NULL, file = NULL,
    estimateMode = 0, random_rec = FALSE,
    msmMode = 0, avgnMode = 0,
    phase = FALSE, verbose = 0)
#> Warning: Passing ‘phase’, ‘verbose’ directly to fit_mod() is deprecated and will be removed in a future release. Bundle these into fit_control() instead, e.g. fit_control(phase = ..., verbose = ...). Forwarding for now.
#> `age_trans_matrix` data does not span range of age for species 1 will fill with 0s
sims <- self_test(ss_run, nsim = 10)
#> Warning: sim_mod() does not simulate diet (stomach content) data: 2025 rows are carried through unchanged. A self_test() of a model fitted to diet data therefore does not propagate diet observation error, and recovery of suitability is optimistic.
#> Warning: sim_mod() does not simulate diet (stomach content) data: 2025 rows are carried through unchanged. A self_test() of a model fitted to diet data therefore does not propagate diet observation error, and recovery of suitability is optimistic.
#> Warning: sim_mod() does not simulate diet (stomach content) data: 2025 rows are carried through unchanged. A self_test() of a model fitted to diet data therefore does not propagate diet observation error, and recovery of suitability is optimistic.
#> Warning: sim_mod() does not simulate diet (stomach content) data: 2025 rows are carried through unchanged. A self_test() of a model fitted to diet data therefore does not propagate diet observation error, and recovery of suitability is optimistic.
#> Warning: sim_mod() does not simulate diet (stomach content) data: 2025 rows are carried through unchanged. A self_test() of a model fitted to diet data therefore does not propagate diet observation error, and recovery of suitability is optimistic.
#> Warning: sim_mod() does not simulate diet (stomach content) data: 2025 rows are carried through unchanged. A self_test() of a model fitted to diet data therefore does not propagate diet observation error, and recovery of suitability is optimistic.
#> Warning: sim_mod() does not simulate diet (stomach content) data: 2025 rows are carried through unchanged. A self_test() of a model fitted to diet data therefore does not propagate diet observation error, and recovery of suitability is optimistic.
#> Warning: sim_mod() does not simulate diet (stomach content) data: 2025 rows are carried through unchanged. A self_test() of a model fitted to diet data therefore does not propagate diet observation error, and recovery of suitability is optimistic.
#> Warning: sim_mod() does not simulate diet (stomach content) data: 2025 rows are carried through unchanged. A self_test() of a model fitted to diet data therefore does not propagate diet observation error, and recovery of suitability is optimistic.
#> Warning: sim_mod() does not simulate diet (stomach content) data: 2025 rows are carried through unchanged. A self_test() of a model fitted to diet data therefore does not propagate diet observation error, and recovery of suitability is optimistic.
# }
```
