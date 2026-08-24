# profile() on a DSEM.
#
# The recruitment SD still profiles -- it just lives somewhere else. Under a
# DSEM the template overwrites R_sd from the SEM's variance path and fit_mod()
# maps R_log_sd out, so profiling R_log_sd would vary a parameter with exactly
# zero gradient and return a FLAT curve, which reads as "the data are
# uninformative" rather than "this parameter is wired to nothing". The
# `sigmaR` / `R_sd` aliases therefore resolve to dsem_beta_z[rec_sd_idx[sp]] --
# the sem's two-headed self-loop -- so the alias means the same thing on both
# parameterizations and `slots` still indexes by species. beta_z is on the
# NATURAL scale (R_sd is its absolute value, a Cholesky factor's sign not being
# identified), so the alias's implied log transform does not apply there.
#
# rec_dev has no such home: it is derived from the latent states on every
# evaluation, carries no gradient, and no single parameter stands for it. That
# one is still refused.
#
# Everything else -- M1, R0, alpha, beta, an arbitrary slot -- profiles normally,
# and a fit WITHOUT a DSEM is untouched.

# `built = FALSE` carries only the specification, which is what a fit assembled
# by hand (or restored from a stale save) looks like: the redirect cannot find
# rec_sd_idx and has to say so rather than guess an index.
profile_fake <- function(dsem = TRUE, built = TRUE, rec_sd_idx = 1L) {
  pars <- list(R_log_sd = c(0.5, 0.5, 0.5),
               rec_dev  = matrix(0, 3, 5),
               log_M1   = array(0, c(3, 1, 1)))
  obj <- list(
    data_list = if (dsem) {
      c(Rceattle::BS2017SS,
        list(dsem_settings = Rceattle::build_DSEM(
          sem = "recdevs1 <-> recdevs1, 0, sigmaR1, 1")))
    } else Rceattle::BS2017SS,
    estimated_params = pars)
  if (dsem && built) {
    obj$estimated_params$dsem_beta_z <- c(0.6, 0.7, 0.8)
    obj$dsem <- list(tmb_inputs = list(data = list(rec_sd_idx = rec_sd_idx)))
  }
  structure(obj, class = "Rceattle")
}

testthat::test_that("the recruitment-SD aliases redirect to the sem's variance path", {
  vals <- list(c(0.4, 0.5, 0.6))

  # All three spellings resolve to R_log_sd before the DSEM branch, so all three
  # must be redirected -- not just the one a user is most likely to type. The
  # message names the beta_z slot, so the redirect is visible rather than silent.
  for (p in c("sigmaR", "R_sd", "R_log_sd")) {
    msg <- testthat::capture_messages(
      try(suppressWarnings(stats::profile(profile_fake(TRUE), param = p,
                                          values = vals, cores = 1)),
          silent = TRUE))
    testthat::expect_match(paste(msg, collapse = " "), "dsem_beta_z\\[1\\]", info = p)
  }

  # Reaching the grid machinery is as far as this fixture can go; what matters is
  # that it is NOT stopped by a DSEM refusal on the way.
  err <- tryCatch(
    suppressWarnings(suppressMessages(
      stats::profile(profile_fake(TRUE), param = "sigmaR", values = vals,
                     cores = 1))),
    error = function(e) conditionMessage(e))
  testthat::expect_false(grepl("on a DSEM fit", err))
})

testthat::test_that("a species whose SD is fixed in the sem has nothing to profile", {
  # rec_sd_idx == 0 marks a species whose recruitment SD is fixed rather than
  # estimated. Redirecting to beta_z[0] would silently profile the wrong slot.
  err <- tryCatch(
    suppressWarnings(suppressMessages(
      stats::profile(profile_fake(TRUE, rec_sd_idx = 0L), param = "sigmaR",
                     slots = list(1), values = list(c(0.4, 0.5)), cores = 1))),
    error = function(e) conditionMessage(e))
  testthat::expect_match(err, "fix the recruitment")
})

testthat::test_that("a spec without the built objects says to refit", {
  err <- tryCatch(
    suppressWarnings(suppressMessages(
      stats::profile(profile_fake(TRUE, built = FALSE), param = "sigmaR",
                     values = list(c(0.4, 0.5)), cores = 1))),
    error = function(e) conditionMessage(e))
  testthat::expect_match(err, "Refit before profiling")
})

testthat::test_that("rec_dev is still refused -- it has no home to redirect to", {
  err <- tryCatch(
    suppressWarnings(stats::profile(profile_fake(TRUE), param = "rec_dev",
                                    values = list(c(0.4, 0.5)), cores = 1)),
    error = function(e) conditionMessage(e))
  testthat::expect_match(err, "cannot profile 'rec_dev' on a DSEM fit")
  # The refusal must point somewhere useful rather than dead-end.
  testthat::expect_match(err, "process_residuals")
})

testthat::test_that("a DSEM does not block profiling anything else", {
  vals <- list(c(0.4, 0.5, 0.6))

  # M1 is untouched by a DSEM. It must NOT hit the DSEM refusal -- it gets as far
  # as the grid machinery, which is all this fixture can exercise.
  err <- tryCatch(
    suppressWarnings(suppressMessages(
      stats::profile(profile_fake(TRUE), param = "M1", values = vals, cores = 1))),
    error = function(e) conditionMessage(e))
  testthat::expect_false(grepl("on a DSEM fit", err))

  # And a fit with no DSEM keeps profiling the recruitment slots as before.
  for (p in c("sigmaR", "rec_dev")) {
    err <- tryCatch(
      suppressWarnings(suppressMessages(
        stats::profile(profile_fake(FALSE), param = p, values = vals, cores = 1))),
      error = function(e) conditionMessage(e))
    testthat::expect_false(grepl("on a DSEM fit", err), info = p)
  }
})
