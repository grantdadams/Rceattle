context("Test automatic mapping of base parameters with linkages")

library(Rceattle)

# Load data and inject a `temp` covariate into env_data so linkages on
# `~ 0 + temp` can be materialized against data_list$env_data.
data(BS2017SS)
yrs <- BS2017SS$styr:BS2017SS$endyr
set.seed(7)
BS2017SS$env_data <- data.frame(
  Year   = yrs,
  BTempC = BS2017SS$env_data$BTempC,
  temp   = rnorm(length(yrs), mean = 5, sd = 1)
)

fit_debug <- function(...) {
  suppressMessages(
    Rceattle::fit_mod(
      data_list   = BS2017SS,
      estimateMode = 3,
      msmMode     = 0,
      ...
    )
  )
}

# --- Test 1: Growth parameter (log_K) linkage ---
test_that("Growth parameter (log_K) is mapped to NA when linked", {
  link_log_K <- linkage_spec(formula = ~ 1, by = ~ species, param = "log_K")

  run <- fit_debug(
    growthFun = build_growth(
      fun = "vonBertalanffy",
      linkages = list(log_K = link_log_K)
    )
  )
  map <- run$map$mapList

  # log_K is mapped to NA for every species/sex
  expect_true(all(is.na(map$ln_growth_pars[, , 1])))

  # log_L1 / log_Linf remain estimable (log_m is NA via build_map_growth
  # for von Bertalanffy regardless of linkages, so we don't check it here)
  expect_false(any(is.na(map$ln_growth_pars[, , 2]))) # log_L1
  expect_false(any(is.na(map$ln_growth_pars[, , 3]))) # log_Linf
})

# --- Test 2: M parameter (log_M1) linkage ---
test_that("M parameter (log_M1) is mapped to NA when linked", {
  link_log_M1 <- linkage_spec(formula = ~ 1, by = ~ species, param = "log_M1")

  run <- fit_debug(
    M1Fun = build_M1(
      M1_model = "sex_age_invariant",
      linkages = list(log_M1 = link_log_M1)
    )
  )
  map <- run$map$mapList

  expect_true(all(is.na(map$ln_M1)))
})

# --- Test 3: Recruitment parameter (log_R0) linkage ---
test_that("Recruitment parameter (log_R0) is mapped to NA when linked", {
  link_log_R0 <- linkage_spec(formula = ~ 1, by = ~ species, param = "log_R0")

  run <- fit_debug(
    recFun = build_srr(srr_fun = 0, linkages = list(log_R0 = link_log_R0))
  )
  map <- run$map$mapList

  expect_true(all(is.na(map$rec_pars[, 1])))
})

# --- Test 4: Shared linkage (by = NULL) zeros and maps every species ---
test_that("Shared (by = NULL) intercept linkage zeros and maps all species", {
  link_shared <- linkage_spec(formula = ~ 1, by = NULL, param = "log_R0")

  run <- fit_debug(
    recFun = build_srr(srr_fun = 0, linkages = list(log_R0 = link_shared))
  )

  # One linkage row, species = NA (shared across the model)
  tbl <- run$data_list$linkage_table
  expect_equal(nrow(tbl), 1L)
  expect_true(is.na(tbl$species[1]))

  # Every species' base log_R0 is zeroed and mapped out
  expect_true(all(run$estimated_params$rec_pars[, 1] == 0))
  expect_true(all(is.na(run$map$mapList$rec_pars[, 1])))
})

# --- Test 5: Intercept vs No-Intercept logic ---
test_that("Base parameter value is 0 and map is NA ONLY when intercept is present", {
  # 1. With Intercept
  link_intercept <- linkage_spec(formula = ~ 1, by = ~ species, param = "log_R0")

  run_intercept <- fit_debug(
    recFun = build_srr(srr_fun = 0, linkages = list(log_R0 = link_intercept))
  )

  expect_equal(as.numeric(run_intercept$estimated_params$rec_pars[1, 1]), 0)
  expect_true(is.na(run_intercept$map$mapList$rec_pars[1, 1]))

  # 2. Without Intercept (relative offset only)
  link_no_intercept <- linkage_spec(formula = ~ 0 + temp, by = ~ species, param = "log_R0")

  run_no_intercept <- fit_debug(
    recFun = build_srr(srr_fun = 0, linkages = list(log_R0 = link_no_intercept))
  )

  # Base parameter keeps its initialized default (9) because there is no intercept
  expect_equal(as.numeric(run_no_intercept$estimated_params$rec_pars[1, 1]), 9)
  # Map is still NA because the parameter is linked
  expect_true(is.na(run_no_intercept$map$mapList$rec_pars[1, 1]))
})
