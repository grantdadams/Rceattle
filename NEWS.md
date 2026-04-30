# Rceattle 4.0.2

## Bug fixes

* `data_check()` now stops with a clear message if a user requests
  `msmMode = 3:9` (Kinzey & Punt 2009 functional responses --
  Holling I/II/III, predator interference, predator preemption,
  Hassell-Varley, Ecosim). Those code paths are unvalidated against
  the current parameter set and should not be used for advice.
* Plotting functions (`plot_timeseries`, `plot_selectivity`,
  `plot_mortality`, `plot_maturity`, `plot_b_eaten`, `plot_b_eaten_prop`,
  `plot_m_at_age`, `plot_m2_at_age_prop`, `plot_f`, `plot_ration`,
  `plot_index`, `plot_catch`, `plot_logindex`, `plot_indexresidual`,
  `plot_comp`, `plot_selectivity_vs_maturity`, `plot_stock_recruit`)
  now restore graphics `par()` on exit instead of mutating the
  user's device permanently.
* Replaced `T` / `F` shortcuts with `TRUE` / `FALSE` throughout the
  package source (~60 occurrences).
* Replaced `.data$Bin` / `.data$Length` inside `tidyr::pivot_wider`
  arguments with quoted strings (tidyselect 1.2.0 deprecation).

## Tests

* `tests/testthat/test-subdir-folders.R` no longer calls
  `testthat::test_dir()` for each subdirectory. Nested `test_dir()`
  inside `test_check()` triggered an "evaluation nested too deeply:
  infinite recursion" abort inside `rlang`'s trace deparser whenever
  any test failed, masking the real failure. Subdirectory test files
  are now discovered with `list.files()` and pulled in via
  `source()` so they register against the outer reporter directly.
* `tests/testthat.R` now wraps `library()` calls in
  `suppressPackageStartupMessages()` so transient build-version
  notices (e.g. "package 'dplyr' was built under R version 4.5.2")
  do not get captured as test warnings whose backtraces then crash
  rlang's expr_deparse at end-of-run.

## Parallelism

* `run_mse()` and `check_mse()` now use `parallel::parLapply` on a
  PSOCK cluster instead of `foreach::foreach(...) %dopar%`.
  - The `%dopar%` path triggered `rlang::expr_deparse` infinite
    recursion under nested `test_that` backtraces because each
    `foreach` invocation captures call frames that recurse during
    error formatting. PSOCK workers are clean R processes with no
    captured promise chains, so the issue does not occur.
  - PSOCK clusters work identically on Windows and macOS/Linux.
  - `run_mse()` gains a `cores` argument (default `NULL` picks
    `parallel::detectCores() - 6`); both functions cap at 2 cores
    when `_R_CHECK_LIMIT_CORES_` is set so they comply with CRAN's
    R CMD check limit.
  - `foreach` and `doParallel` removed from `Imports:`.

## Installation / dependencies

* `TMBhelper` moved from `Imports:` to `Suggests:`. Rceattle now uses an
  internal `.fit_tmb()` wrapper that delegates to `TMBhelper::fit_tmb()`
  when the (optional, GitHub-only) package is installed and otherwise
  falls back to a `stats::nlminb()` + `TMB::sdreport()` path. This
  removes the largest install-friction barrier for new users.
* `dplyr` moved from `Depends:` to `Imports:`. The package no longer
  attaches `dplyr` to the user's search path on load (so it no longer
  masks `stats::filter` / `stats::lag`).
* `kableExtra` dropped from `Suggests:`; vignettes now use
  `knitr::kable()` for table rendering.
* `quarto` dropped as a vignette engine; all vignettes are now `.Rmd`.

## API

* `run_mse()` no longer has `om = ms_run`, `em = ss_run` defaults.
  Both arguments are now required and validated as objects of class
  `"Rceattle"` before the MSE loop runs. Calling `run_mse()` with no
  arguments previously produced a confusing "object 'ms_run' not found"
  error; it now stops with a clear message.

## New methods

* `print.Rceattle()` and `summary.Rceattle()`. Auto-printing a fit
  inside knitr / RStudio / R Markdown previously dumped tens of MB of
  nested data and could trigger deep recursion errors during vignette
  rendering.

## Build / packaging

* Tightened `.Rbuildignore` (excludes `examples/`, `R/dev/`,
  `src/TMB/Dev/`, `.Rhistory`, `.claude/`, `.DS_Store`, build tarballs,
  `.Rcheck` directories).
* Tightened `.gitignore` to catch all `*.o` / `*.so` / `*.dll` and
  `*.Rcheck/` directories.
* Suppressed a benign clang `-Wfixed-enum-extension` warning by
  scoping the diagnostic pragma around `#include <TMB.hpp>` rather
  than via global Makevars flags.

# Rceattle 4.0.1

* See `inst/Running_list_of_updates.qmd` for the granular running log
  of model and data-format changes since 4.0.0.

# Rceattle 4.0.0

* CEATTLE TMB version 4.0.0. See Adams et al. (2022),
  *Fisheries Research*, 251, 106303.
