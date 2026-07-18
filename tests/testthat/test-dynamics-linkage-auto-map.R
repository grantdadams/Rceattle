# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::context("Test automatic mapping of base parameters with linkages")

library(Rceattle)

# Load data and inject a `temp` covariate into env_data so linkages on
# `~ 0 + temp` can be materialized against data_list$env_data.
# 1) Set up simulation
set.seed(7)
nyrs = 30
nspp = 2
Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
log_phi = matrix(-Inf, nspp, nspp, byrow = TRUE)

# First, simulate some data for the model
set.seed(123)
sim <- make_msm_test_data(
  years = 1:nyrs,
  Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),
  log_phi = log_phi
)


# Set up Rceattle data
simData <- sim$data_list
yrs <- simData$styr:simData$projyr

simData$env_data <- data.frame(
  Year   = yrs,
  BTempC = rnorm(length(yrs)),
  temp   = rnorm(length(yrs), mean = 5, sd = 1)
)

fit_debug <- function(...) {
  suppressMessages(
    Rceattle::fit_mod(
      data_list   = simData,
      estimateMode = 3,
      msmMode     = 0,
      ...
    )
  )
}

# Intercept-bearing linkages (~ 1, ~ temp): the base parameter stays
# estimable and carries the level. The linkage's (Intercept) row is
# fixed at 0 (mapped NA) so the linkage column is purely an offset.
# Slope-only linkages (~ 0 + temp): no (Intercept) row, so the base
# parameter is masked out and stays at its build_params() default;
# the slope row carries the year-by-year offset on top.

# --- Test 1: Growth parameter (K), intercept-only linkage ---
testthat::test_that("Growth K base stays estimable when linkage carries an intercept", {
  link_K <- linkage_spec(formula = ~ 1, by = ~ species, param = "K")

  run <- fit_debug(
    growthFun = build_growth(
      fun = "vonBertalanffy",
      linkages = list(K = link_K)
    )
  )
  map <- run$map$mapList

  # K stays estimable (the base carries the level).
  testthat::expect_false(any(is.na(map$log_growth_pars[, , 1])))
  # L1 / Linf are unaffected.
  testthat::expect_false(any(is.na(map$log_growth_pars[, , 2])))
  testthat::expect_false(any(is.na(map$log_growth_pars[, , 3])))

  # The (Intercept) row of beta_linkage is mapped out (fixed at 0).
  tbl <- run$data_list$linkage_table
  testthat::expect_true(all(tbl$design_col == "(Intercept)"))
  testthat::expect_true(all(is.na(run$map$mapList$beta_linkage)))
  testthat::expect_true(all(run$estimated_params$beta_linkage == 0))
})

# --- Test 2: M parameter (M1), intercept-only linkage ---
testthat::test_that("M M1 base stays estimable when linkage carries an intercept", {
  link_M1 <- linkage_spec(formula = ~ 1, by = ~ species, param = "M1")

  run <- fit_debug(
    M1Fun = build_M1(
      M1_model = "sex_age_invariant",
      linkages = list(M1 = link_M1)
    )
  )
  map <- run$map$mapList

  testthat::expect_false(any(is.na(map$log_M1)))
  testthat::expect_true(all(is.na(run$map$mapList$beta_linkage)))
  testthat::expect_true(all(run$estimated_params$beta_linkage == 0))
})

# --- Test 3: Recruitment parameter (R0), intercept-only linkage ---
testthat::test_that("Recruitment R0 base stays estimable when linkage carries an intercept", {
  link_R0 <- linkage_spec(formula = ~ 1, by = ~ species, param = "R0")

  run <- fit_debug(
    recFun = build_srr(srr_fun = 0, linkages = list(R0 = link_R0))
  )
  map <- run$map$mapList

  testthat::expect_false(any(is.na(map$rec_pars[, 1])))
  testthat::expect_true(all(is.na(run$map$mapList$beta_linkage)))
})

# --- Test 4: Shared linkage (by = NULL): one row, base stays estimable ---
testthat::test_that("Shared (by = NULL) intercept linkage keeps base estimable", {
  link_shared <- linkage_spec(formula = ~ 1, by = NULL, param = "R0")

  run <- fit_debug(
    recFun = build_srr(srr_fun = 0, linkages = list(R0 = link_shared))
  )

  tbl <- run$data_list$linkage_table
  testthat::expect_equal(nrow(tbl), 1L)
  testthat::expect_true(is.na(tbl$species[1]))

  # Base R0 keeps its build_params() default (9) -- no init was supplied.
  testthat::expect_true(all(run$estimated_params$rec_pars[, 1] == 9))
  # Base R0 is estimable for every species.
  testthat::expect_false(any(is.na(run$map$mapList$rec_pars[, 1])))
})

# --- Test 5: Intercept vs No-Intercept divergence ---
testthat::test_that("Intercept-bearing keeps base estimable; no-intercept masks it", {
  # 1. With Intercept: base stays estimable, defaults to 9.
  link_intercept <- linkage_spec(formula = ~ 1, by = ~ species, param = "R0")

  run_intercept <- fit_debug(
    recFun = build_srr(srr_fun = 0, linkages = list(R0 = link_intercept))
  )

  testthat::expect_equal(as.numeric(run_intercept$estimated_params$rec_pars[1, 1]), 9)
  testthat::expect_false(is.na(run_intercept$map$mapList$rec_pars[1, 1]))

  # 2. Without Intercept (slope only): base is masked at default (9).
  link_no_intercept <- linkage_spec(formula = ~ 0 + temp, by = ~ species, param = "R0")

  run_no_intercept <- fit_debug(
    recFun = build_srr(srr_fun = 0, linkages = list(R0 = link_no_intercept))
  )

  testthat::expect_equal(as.numeric(run_no_intercept$estimated_params$rec_pars[1, 1]), 9)
  testthat::expect_true(is.na(run_no_intercept$map$mapList$rec_pars[1, 1]))
})

# --- Test 6: User-supplied (Intercept) init flows to the base parameter ---
testthat::test_that("(Intercept) init overrides the base parameter's default", {
  link_with_init <- linkage_spec(
    formula = ~ 1,
    by      = ~ species,
    param   = "R0",
    init    = list(`(Intercept)` = exp(7.5))
  )

  run <- fit_debug(
    recFun = build_srr(srr_fun = 0, linkages = list(R0 = link_with_init))
  )

  testthat::expect_true(all(run$estimated_params$rec_pars[, 1] == 7.5))
  testthat::expect_true(all(run$estimated_params$beta_linkage == 0))
  testthat::expect_false(any(is.na(run$map$mapList$rec_pars[, 1])))
})
