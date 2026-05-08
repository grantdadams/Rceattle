# Define M1 specifications

Define M1 specifications

## Usage

``` r
build_M1(
  M1_model = 0,
  M1_re = 0,
  updateM1 = FALSE,
  M1_use_prior = FALSE,
  M2_use_prior = FALSE,
  M_prior = 0.4,
  M_prior_sd = 0.35,
  M1_indices = NA,
  linkages = NULL
)
```

## Arguments

- M1_model:

  Vector or scalar specifying the M1 structural fixed- effects model.
  Either an integer code or the equivalent string alias (both forms are
  accepted; the integer code is canonical):

  - `0` / `"fixed"` – use the input `M1_base` (no estimation).

  - `1` / `"sex_age_invariant"` – estimate one `M1_{spp}`.

  - `2` / `"sex_specific"` – estimate `M1_{spp, sex}`.

  - `3` / `"sex_age_specific"` – estimate `M1_{spp, sex, age}`.

  - `4`, `5` – soft-deprecated env-driven codes; use the `linkages`
    argument instead. See `vignette("environmental-linkages")`.

- M1_re:

  Vector or scalar specifying the M1 random-effects model. Either an
  integer code or the equivalent string alias: `0` / `"none"`, `1` /
  `"iid_age"`, `2` / `"iid_year"`, `3` / `"iid_age_year"`, `4` /
  `"ar1_age"`, `5` / `"ar1_year"`, `6` / `"ar1_age_year"`.

- updateM1:

  If using initial parameters, use M1 fixed effects from data
  (`M1_base`) instead. Default `FALSE`.

- M1_use_prior:

  Vector or scalar; if `TRUE`, apply the lognormal `M_prior` /
  `M_prior_sd` to `log_M1` directly.

- M2_use_prior:

  Vector or scalar; if `TRUE`, apply the lognormal prior to `M1 + M2` in
  multi-species models.

- M_prior:

  Mean (natural-scale) of the lognormal prior on M.

- M_prior_sd:

  SD (log-scale) of the lognormal prior on M.

- M1_indices:

  Soft-deprecated. Vector of column indices into `env_data` (excluding
  `Year`) for environmentally linked M1 when `M1_model %in% c(4, 5)`.
  Use the `linkages` argument instead; see
  `vignette("environmental-linkages")`.

- linkages:

  Optional named list of
  [`linkage_spec()`](https://grantdadams.github.io/Rceattle/reference/linkage_spec.md)
  objects keyed by M parameter name (currently the only valid key is
  `"log_M1"`). Each spec describes how `log_M1` depends on environmental
  covariates and on stratifying factors (species, sex, age). The offset
  enters additively (on the log scale) inside the `M1_at_age` compute. A
  row's `age_bin == NA` broadcasts the offset across ages; specific
  values pin it to that age slice.

## Value

A list of switches for defining the M1 model.

## Examples

``` r
if (FALSE) { # \dontrun{
# Sex/age-invariant M with a temperature linkage on log_M1
build_M1(
  M1_model = "sex_age_invariant",
  linkages = list(
    log_M1 = linkage_spec(
      formula = ~ temp,
      by      = ~ species,
      priors  = list(temp = normal(0, 0.5))
    )
  )
)
} # }
```
