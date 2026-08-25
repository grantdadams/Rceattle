# `growth_re`, `growth_indices` and the `log_growth_par_devs` parameter block
# were removed in 5.9.0. Neither was ever consumed -- the deviation array was
# mapped off in every configuration and the template gave it no density -- so
# removing them leaves every fit bit-identical. Both back-compat paths are what
# makes that a minor release rather than a breaking one, so pin them here: a
# data_list still carrying the switch, and `inits` from a fit saved before the
# parameter block went away.

testthat::test_that("a data_list carrying growth_re or growth_indices is dropped with a message", {
  d <- make_test_data()

  for (old in c("growth_re", "growth_indices")) {
    dl <- d
    dl[[old]] <- 1
    testthat::expect_message(out <- Rceattle::switch_check(dl), "deprecated and ignored")
    testthat::expect_null(out[[old]])
  }

  # Both at once, and each is named.
  dl <- d
  dl$growth_re <- 1
  dl$growth_indices <- 1
  msg <- testthat::capture_messages(out <- Rceattle::switch_check(dl))
  testthat::expect_true(any(grepl("growth_re", msg, fixed = TRUE)))
  testthat::expect_true(any(grepl("growth_indices", msg, fixed = TRUE)))
  testthat::expect_null(out$growth_re)
  testthat::expect_null(out$growth_indices)

  # A data_list that never had them says nothing.
  testthat::expect_silent(Rceattle::switch_check(d))
})


testthat::test_that("inits from a fit carrying the retired log_growth_par_devs block still load", {
  testthat::skip_on_cran()

  d <- make_test_data()
  fit <- suppressWarnings(Rceattle::fit_mod(
    data_list = d, estimateMode = 3, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)))

  # A pre-5.9.0 fit's estimated_params carried a block the template no longer
  # declares. build_map() runs on the inits before MakeADFun does, so a stale
  # name used to reach the map and stop the refit.
  stale <- fit$estimated_params
  testthat::expect_false("log_growth_par_devs" %in% names(stale))
  stale$log_growth_par_devs <- array(0, dim = c(1, 1, 1, 1))

  refit <- suppressWarnings(Rceattle::fit_mod(
    data_list = d, inits = stale, estimateMode = 3, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)))

  testthat::expect_false("log_growth_par_devs" %in% names(refit$estimated_params))
  # Dropping an inert block changes nothing about the objective.
  testthat::expect_equal(refit$obj$fn(refit$obj$par), fit$obj$fn(fit$obj$par),
                         tolerance = 1e-12)
})
