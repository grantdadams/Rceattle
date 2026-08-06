# =============================================================================
# Penalized fixed-effect linkage deviations: linkage_spec(integrate = FALSE).
#
# A random-effect term's deviations are normally integrated out by the Laplace
# approximation. `integrate = FALSE` instead estimates them as a plain penalized
# fixed effect, which is what the ADMB/AMAK reference models (and the legacy
# Time_varying_sel / Time_varying_q switches) do. It is permitted only when the
# deviation SD is fixed: a variance is consistently estimated only by
# integrating, so estimating the deviations AND their SD jointly as fixed
# effects is degenerate -- the likelihood improves without bound as both go to
# zero. The guards below are the constructor-level half of that contract; the
# pooled-group half lives in .re_group_table().
# =============================================================================

testthat::test_that("integrate must be a single TRUE/FALSE", {
  for (bad in list("yes", c(TRUE, FALSE), NA, 1)) {
    testthat::expect_error(
      Rceattle::linkage_spec(~ rw(1 | Year), init = list(sigma = 0.05),
                             integrate = bad),
      "single TRUE/FALSE")
  }
})

testthat::test_that("integrate = FALSE requires a fixed deviation SD", {
  # No sigma supplied: the SD would be estimated alongside the deviations.
  testthat::expect_error(
    Rceattle::linkage_spec(~ rw(1 | Year), integrate = FALSE),
    "requires a fixed deviation SD")

  # A sigma PRIOR still estimates the SD -- shrunk, but not fixed -- so it does
  # not satisfy the contract even though `init` is present.
  testthat::expect_error(
    Rceattle::linkage_spec(~ rw(1 | Year), init = list(sigma = 0.05),
                           priors = list(sigma = lognormal(log(0.05), 0.3)),
                           integrate = FALSE),
    "requires a fixed deviation SD")
})

testthat::test_that("integrate = FALSE on ar1 also requires a fixed rho", {
  testthat::expect_error(
    Rceattle::linkage_spec(~ ar1(1 | Year), init = list(sigma = 0.05),
                           integrate = FALSE),
    "requires a fixed correlation")
  testthat::expect_error(
    Rceattle::linkage_spec(~ ar1(1 | Year), init = list(sigma = 0.05, rho = 0.5),
                           priors = list(rho = normal(0, 0.3)),
                           integrate = FALSE),
    "requires a fixed correlation")
  # Both pinned: accepted.
  spec <- Rceattle::linkage_spec(~ ar1(1 | Year),
                                 init = list(sigma = 0.05, rho = 0.5),
                                 integrate = FALSE)
  testthat::expect_false(spec$integrate)
})

testthat::test_that("integrate = FALSE is rejected without a random-effect term", {
  testthat::expect_error(
    Rceattle::linkage_spec(~ temp, init = list(sigma = 0.05), integrate = FALSE),
    "random-effect term")
})

testthat::test_that("integrate = FALSE is rejected with observe", {
  # An observed latent state is a state-space model: the state must be
  # integrated for the observation to identify it.
  testthat::expect_error(
    Rceattle::linkage_spec(~ ar1(1 | Year), observe = "qcov", obs_sd = 0.2,
                           init = list(sigma = 0.05, rho = 0.5),
                           integrate = FALSE),
    "cannot be combined with `observe`")
})

testthat::test_that("integrate defaults to TRUE and round-trips onto the spec", {
  testthat::expect_true(Rceattle::linkage_spec(~ rw(1 | Year))$integrate)
  testthat::expect_false(
    Rceattle::linkage_spec(~ rw(1 | Year), init = list(sigma = 0.05),
                           integrate = FALSE)$integrate)
})

testthat::test_that("the linkage table carries re_integrate on RE rows only", {
  env <- data.frame(Year = 2000:2009)
  pool <- function(spec) {
    Rceattle:::pool_linkages(list(q = list(q = spec)), env_data = env,
                             strata = list(fleet = 1L))$table
  }

  pen <- pool(Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet, fleet = 1L,
                                     init = list(sigma = 0.05),
                                     integrate = FALSE))
  is_re <- !is.na(pen$re_struct)
  testthat::expect_true(any(is_re))
  testthat::expect_true(all(pen$re_integrate[is_re] == FALSE))
  # Fixed rows (the intercept that carries the walk's level) stay NA.
  testthat::expect_true(all(is.na(pen$re_integrate[!is_re])))

  int <- pool(Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet, fleet = 1L))
  testthat::expect_true(all(int$re_integrate[!is.na(int$re_struct)] == TRUE))
})

testthat::test_that("the group registry is unchanged for integrated models", {
  # Back-compat: adding the flag must not renumber a single existing sigma group.
  env <- data.frame(Year = 2000:2004)
  tbl <- Rceattle:::pool_linkages(
    list(q = list(q = Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet,
                                             fleet = 1L))),
    env_data = env, strata = list(fleet = 1L))$table
  re <- !is.na(tbl$re_struct)
  testthat::expect_equal(unique(tbl$sigma_index[re]), 0L)
  testthat::expect_setequal(tbl$re_index[re], seq_len(sum(re)) - 1L)

  gt <- Rceattle:::.re_group_table(tbl)
  testthat::expect_true(gt$integrate)
  testthat::expect_false(gt$sigma_fixed)   # estimated SD is fine when integrated
})

testthat::test_that("a penalized group carries integrate = FALSE and a fixed SD", {
  env <- data.frame(Year = 2000:2004)
  tbl <- Rceattle:::pool_linkages(
    list(q = list(q = Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet,
                                             fleet = 1L,
                                             init = list(sigma = 0.05),
                                             integrate = FALSE))),
    env_data = env, strata = list(fleet = 1L))$table
  gt <- Rceattle:::.re_group_table(tbl)
  testthat::expect_false(gt$integrate)
  testthat::expect_true(gt$sigma_fixed)
  testthat::expect_equal(gt$sigma_start, 0.05)
})

testthat::test_that("a penalized group with an estimated SD is refused at pooling", {
  # The constructor cannot catch this: each spec is individually legal, and only
  # the pooled group reveals that the SD is estimated for a penalized effect.
  # Build the table by hand to reach .re_group_table() with that combination.
  env <- data.frame(Year = 2000:2004)
  tbl <- Rceattle:::pool_linkages(
    list(q = list(q = Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet,
                                             fleet = 1L,
                                             init = list(sigma = 0.05),
                                             integrate = FALSE))),
    env_data = env, strata = list(fleet = 1L))$table
  tbl$re_sigma_init <- NA_real_          # drop the fix -> SD becomes estimated
  testthat::expect_error(Rceattle:::.re_group_table(tbl),
                         "require a fixed deviation SD")
})

testthat::test_that("a group cannot be both integrated and penalized", {
  env <- data.frame(Year = 2000:2004)
  tbl <- Rceattle:::pool_linkages(
    list(q = list(q = Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet,
                                             fleet = 1L,
                                             init = list(sigma = 0.05),
                                             integrate = FALSE))),
    env_data = env, strata = list(fleet = 1L))$table
  re <- which(!is.na(tbl$re_struct))
  tbl$re_integrate[re[1]] <- TRUE        # straddle one row across both vectors
  testthat::expect_error(Rceattle:::.assign_re_registry(tbl),
                         "both integrated and penalized")
})

testthat::test_that("integrate survives a save_config round trip", {
  spec <- Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet, fleet = 1L,
                                 init = list(sigma = 0.05), integrate = FALSE)
  lst <- Rceattle:::.rce_spec_to_list(spec)
  testthat::expect_false(lst$integrate)
  testthat::expect_false(Rceattle:::.rce_spec_from_list(lst)$integrate)
})
