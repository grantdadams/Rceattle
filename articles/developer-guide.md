# Developer guide

This guide is for people who want to modify Rceattle itself: add a
parameter, a data input, a new selectivity or suitability form, or a
likelihood component. It is not needed to *fit* models — for that, start
with the user articles
([Introduction](https://grantdadams.github.io/Rceattle/introduction.md)
and [Model options and
functionality](https://grantdadams.github.io/Rceattle/model-options-and-functionality.md)).

It consolidates and updates the developer notes on the [GitHub
wiki](https://github.com/grantdadams/Rceattle/wiki).

## Repository layout

| Path | Contents |
|----|----|
| `R/` | The R side. Files are number-prefixed in rough pipeline order (see below). |
| `src/TMB/` | The C++ TMB template `ceattle_v01_11.cpp` and its included `.hpp` modules. |
| `tests/testthat/` | Unit and regression tests, grouped by area. |
| `inst/extdata/` | Excel templates, including `meta_data_names.xlsx` (the authoritative column/switch reference). |
| `data/` | Built-in example data lists (`BS2017SS`, `GOA2018SS`, `NorthernRockfish2022`, …). |
| `vignettes/` | User articles; `vignettes/articles/` holds website-only articles like this one. |
| `.github/workflows/` | `R-CMD-check.yaml` and `pkgdown.yaml` CI. |

The whole model — observation model, population dynamics, and likelihood
— is a single C++ TMB script, `src/TMB/ceattle_v01_11.cpp`, which
`#include`s a set of topical `.hpp` modules.
[`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
is the R entry point; it calls a chain of helpers that build the TMB
inputs, optimize, and label the output.

## The fit pipeline

[`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
(`R/6-fit_mod.R`) orchestrates the following stages. The number-prefixed
file names mirror this order.

| Stage | File | Function | Role |
|----|----|----|----|
| 1 | `R/0-clean_data.R` | [`clean_data()`](https://grantdadams.github.io/Rceattle/reference/clean_data.md) | Coerce and clean the incoming `data_list`. |
| 2 | `R/0-switches.R` | [`switch_check()`](https://grantdadams.github.io/Rceattle/reference/switch_check.md) | Fill missing switches with defaults; normalize string/integer switches. |
| 3 | `R/1-data_check.R` | [`data_check()`](https://grantdadams.github.io/Rceattle/reference/data_check.md) | Validate inputs; error early on unsupported options. |
| 4 | `R/2-build_params.R` | [`build_params()`](https://grantdadams.github.io/Rceattle/reference/build_params.md) | Build the starting parameter list from the switches. |
| 5 | `R/3-build_map.R` | [`build_map()`](https://grantdadams.github.io/Rceattle/reference/build_map.md) | Build the TMB `map` (which parameters are estimated vs. fixed). |
| 6 | `R/4-build_parameter_bounds.R` | [`build_bounds()`](https://grantdadams.github.io/Rceattle/reference/build_bounds.md) | Lower/upper parameter bounds. |
| 7 | `R/5-rearrange_data.R` | [`rearrange_data()`](https://grantdadams.github.io/Rceattle/reference/rearrange_data.md) | Reshape the data for TMB (and build the OSA observation vector when requested). |
| 8 | `R/6-fit_mod.R` | [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md) | [`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html) → [`nlminb()`](https://rdrr.io/r/stats/nlminb.html) → [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html). |
| 9 | `R/6-rename_output.R` | [`rename_output()`](https://grantdadams.github.io/Rceattle/reference/rename_output.md) | Label derived quantities on the returned object. |

Phasing (optional staged estimation) is handled by
[`set_phases()`](https://grantdadams.github.io/Rceattle/reference/set_phases.md)
/
[`TMBphase()`](https://grantdadams.github.io/Rceattle/reference/TMBphase.md)
in `R/OPT-phaser.R`;
[`TMBAIC()`](https://grantdadams.github.io/Rceattle/reference/TMBAIC.md)
(`R/OPT-TMBAIC.R`) computes AIC. Downstream wrappers —
[`retrospective()`](https://grantdadams.github.io/Rceattle/reference/retrospective.md),
[`jitter()`](https://grantdadams.github.io/Rceattle/reference/jitter.md),
[`run_mse()`](https://grantdadams.github.io/Rceattle/reference/run_mse.md),
[`sim_mod()`](https://grantdadams.github.io/Rceattle/reference/sim_mod.md),
[`self_test()`](https://grantdadams.github.io/Rceattle/reference/self_test.md),
the `plot_*()` family — all take the fitted `Rceattle` object returned
by
[`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md).

## The switch system

Model options are supplied either as strings (`"Logistic"`, `"NPFMC"`)
or as the integer codes the TMB template consumes. `R/0-switches.R` is
the single source of truth for that correspondence, via four functions:

- [`switch_check()`](https://grantdadams.github.io/Rceattle/reference/switch_check.md)
  — fill missing switches with defaults and normalize.
- `validate_switches()` — error if any switch is not a known
  code/string.
- `convert_switches()` — canonical string → integer code for TMB.
- `revert_switches()` — integer code → canonical string (backward
  compatibility for older integer-coded data files).

The mappings live in named vectors: `sel_map`, `tv_sel_map`, `q_map`,
`tv_q_map`, `comp_loglike_map`, `fleet_map`, `initMode_map`,
`suitMode_map`, and `hcr_map`. Adding a value to a map is what makes a
new string name legal; `validate_switches()` rejects anything not
listed.

## The TMB template

`ceattle_v01_11.cpp` includes these modules (order matches the top of
the file):

| Module | Responsibility |
|----|----|
| `helper_functions.hpp` | Shared utilities. |
| `comp_osa.hpp` | Composition likelihoods and the one-step-ahead (OSA) decomposition. |
| `growth.hpp` | Growth (von Bertalanffy / Richards); weight- and length-at-age. |
| `selectivity.hpp` | Selectivity forms (cases 0–7, including the Hake/Taylor non-parametric case 5). |
| `recruitment.hpp` | Stock–recruitment functions. |
| `bioenergetics.hpp` | Temperature-dependent ration / consumption (`Ceq`). |
| `predation.hpp` | Prey suitability and predation mortality (MSVPA type-2, estimated iteratively). |
| `diet_data.hpp` | Predator diet (stomach-content) likelihood. |
| `linkage.hpp` | Environmental / covariate linkage machinery. |

The template opens with the `DATA_*` block (switches such as
`estimateMode`, `msmMode`, `suitMode`, `initMode`, `srr_fun`, `HCR`,
plus the data objects), then the `PARAMETER_*` block, then the
population and observation dynamics, and finally the joint negative
log-likelihood. The likelihood is accumulated in the `jnll_comp` object,
organized by species, fleet, and likelihood component.

**Building.** `src/TMB/compile.R` compiles the template
(`framework = "TMBad"`); `TMB` and `RcppEigen` are `LinkingTo`
dependencies. The spurious Eigen `-Warray-bounds` warnings are silenced
by a source `#pragma` in the `.cpp` rather than a compiler flag, since
CRAN rejects warning-suppressing flags.

## The linkage and priors system

Time-varying and covariate-driven parameters (M, growth, recruitment,
catchability, selectivity) share one formula-driven linkage system:
`R/0-build_linkage.R`, `R/0-linkage_encode.R`, `R/0-linkage_table.R`,
and `R/0-priors.R`. `linkage_spec(formula = ~ x, priors = list(...))` is
the user entry point.

Note that inside `priors =`, the bare constructors `normal()`,
`lognormal()`, [`gamma()`](https://rdrr.io/r/base/Special.html),
[`beta()`](https://rdrr.io/r/base/Special.html) are **not** the exported
`prior_*()` functions.
[`linkage_spec()`](https://grantdadams.github.io/Rceattle/reference/linkage_spec.md)
evaluates the `priors` quosure in a data mask
(`rlang::eval_tidy(..., data = .prior_dispatch_mask())`) that binds
`normal` → `prior_normal`, and so on. Both spellings work inside
[`linkage_spec()`](https://grantdadams.github.io/Rceattle/reference/linkage_spec.md),
and the bare form is intentional — do not “fix” it to
[`prior_normal()`](https://grantdadams.github.io/Rceattle/reference/prior_normal.md).

## Recipes

### Add a new estimated parameter

1.  Declare the `PARAMETER` in `ceattle_v01_11.cpp` and use it in the
    relevant `.hpp` module.
2.  Add a default starting value in `R/2-build_params.R`.
3.  Add mapping logic in `R/3-build_map.R` (and the matching
    `build_map_*()` helper) so it is turned on/off for the right
    configurations.
4.  Add it to the phasing order in `R/OPT-phaser.R`.
5.  (Optional) expose an on/off or random-effect switch as a
    [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
    argument.
6.  Document it in `inst/extdata/meta_data_names.xlsx` and/or the
    [`fit_mod()`](https://grantdadams.github.io/Rceattle/reference/fit_mod.md)
    roxygen.

### Add or change a data input

1.  Add the `DATA_*` object in `ceattle_v01_11.cpp` and consume it in
    the module.
2.  Update `R/0-read_write_excel_data.R` to read/write the new sheet or
    column.
3.  Add validation in `R/1-data_check.R`.
4.  Add reshaping in `R/5-rearrange_data.R` if TMB needs a different
    layout.
5.  Document it in `meta_data_names.xlsx` and/or roxygen.

### Add a new switch option (e.g. a selectivity form or `suitMode`)

1.  Implement the case in the module `.hpp` (e.g. `selectivity.hpp` case
    *N*, `predation.hpp` `smode` *N*).
2.  Register the string ↔︎ integer mapping in the relevant map in
    `R/0-switches.R` (e.g. `sel_map`, `suitMode_map`); otherwise
    `validate_switches()` rejects the new name.
3.  If the option is not yet ready for use, block it in
    `R/1-data_check.R` — as the length-based suitability modes are
    (`suitMode %in% c(1, 3, 5)`).
4.  Add bounds in `R/4-build_parameter_bounds.R` if the option
    introduces bounded parameters (see the `suitMode %in% c(1:2)`
    branch).
5.  Add a test under `tests/testthat/tests-<Area>/`.
6.  Document it.

### Change the likelihood

The joint negative log-likelihood is assembled at the end of
`ceattle_v01_11.cpp`, in `jnll_comp`, indexed by species, fleet, and
component. Add or modify the relevant term there, then recompile and
re-run the regression tests (below) to confirm the objective changes
only as intended.

## Testing

Tests use `testthat` and are grouped by area under `tests/testthat/`:
`tests-Selectivity/`, `tests-Dynamics/`, `tests-Likelihoods/`,
`tests-Data-processing/`, `tests-Growth/`, and `tests-Mortality/`.
Notable fixtures and guards:

- `tests/testthat/fixtures/fit_baseline.rds` and the golden-`jnll`
  regression test pin the objective-function value so unintended
  numerical changes surface.
- `test-tmb-makeadfun_smoke.R` checks the template builds and evaluates.
- `test-parameter-recovery.R` fits simulated data and checks estimates.
- `helpers-make-msm-data.R` builds small multi-species fixtures.

Run the suite with `devtools::test()` (or
[`testthat::test_local()`](https://testthat.r-lib.org/reference/test_package.html)).
After any C++ change, recompile (`src/TMB/compile.R` or
[`pkgbuild::compile_dll()`](https://pkgbuild.r-lib.org/reference/compile_dll.html))
before testing, since the tests load the compiled `.so`/`.dll`.

## Branches and releases

`origin/HEAD` points at `main`. The general model, per the Onboarding
wiki page, is: cut feature branches from `dev`, merge them back into
`dev`, and merge `dev` into `main` once a project’s developments are
complete.

| Branch | Role |
|----|----|
| `main` | Most stable / documented; CI and CRAN target. |
| `dev` | Active development branch; base for feature branches. |
| `dev-DSEM` | `dev` with DSEM-linked recruitment (experimental). |
| `dev-ebs-pk` | EBS pollock application / bridge. |
| `dev-RTMB` | RTMB port (remote). |
| `testing-suite-overhaul` | Test-suite expansion. |
| `depricated-ceattle_classic*` | Original Holsman et al. (2016) single-sex model (`ceattle_v01_02`/`_04`.cpp); historical reference only. |

**Versioning (SemVer).** MAJOR for breaking changes, MINOR for new
features/functionality, PATCH for bug fixes. Update `DESCRIPTION` and
add a `NEWS.md` entry with each release.

**Commit messages (Conventional Commits).** Prefix with one of `feat`,
`fix`, `docs`, `style`, `refactor`, `perf`, `test`, `chore`, optionally
scoped, e.g. `fix(selectivity): normalize Hake-type curve`.

**CI.** `R-CMD-check.yaml` runs on macOS, Windows, and Ubuntu (R
release) on push/PR to `main`/`master` and weekly; `pkgdown.yaml`
rebuilds this site.

## Debugging tips

- **`gradient is of length 1` (or similar).** The population is usually
  crashing — numbers-at-age going negative or a divide-by-zero. Check
  starting values and bounds for the offending process.
- **Model will not converge.** Inspect `mod$identified` (when reported)
  to see which parameters are hard to estimate, and adjust their
  starting values.
- **Reproducing an external model.** Supply numbers-at-age via the
  `NByageFixed` sheet and selectivity via `emp_sel` with
  `Selectivity = "Fixed"`, set `estDynamics = 1`, and run
  `fit_mod(estimateMode = 1)` to evaluate without re-fitting. (See also
  the [Stock Synthesis
  conversion](https://grantdadams.github.io/Rceattle/stock-synthesis-conversion.md)
  article.)

## See also

- [GitHub wiki](https://github.com/grantdadams/Rceattle/wiki) — the
  original developer notes this page consolidates.
- [Model
  parameterizations](https://grantdadams.github.io/Rceattle/model-parameterizations.md)
  — equation-level detail on selectivity, catchability, and predation.
- Adams et al. (2022), *Fisheries Research* 251:106303 — the TMB
  generalisation of CEATTLE. Holsman et al. (2016), *Deep-Sea Research
  II* 134:360–378 — the original model. Wassermann et al. (2024), *ICES
  JMS* — the hake/cannibalism extension.
