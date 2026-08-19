# profile() on a DSEM.
#
# A DSEM does not block profiling in general -- only the two slots it maps out.
# The template overwrites rec_dev from the latent states and R_sd from the SEM's
# variance path, so fit_mod() maps both blocks out. Profiling either would vary a
# parameter with exactly zero gradient and return a FLAT profile, which reads as
# "the data are uninformative about this" rather than "this parameter is wired to
# nothing". That is the failure worth refusing: it produces a plausible plot
# rather than an error.
#
# Everything else -- M1, R0, alpha, beta, an arbitrary slot -- profiles normally,
# and a fit WITHOUT a DSEM is untouched.

profile_fake <- function(dsem = TRUE) {
  structure(
    list(
      data_list = if (dsem) {
        c(Rceattle::BS2017SS,
          list(dsem_settings = Rceattle::build_DSEM(
            sem = "recdevs1 <-> recdevs1, 0, sigmaR1, 1")))
      } else Rceattle::BS2017SS,
      estimated_params = list(R_log_sd = c(0.5, 0.5, 0.5),
                              rec_dev  = matrix(0, 3, 5),
                              log_M1   = array(0, c(3, 1, 1)))),
    class = "Rceattle")
}

testthat::test_that("profile() refuses the slots a DSEM maps out, by any spelling", {
  vals <- list(c(0.4, 0.5, 0.6))

  # The natural-scale aliases resolve to R_log_sd before the check, so all three
  # spellings must be caught -- not just the one the user is most likely to type.
  for (p in c("sigmaR", "R_sd", "R_log_sd")) {
    testthat::expect_error(
      suppressWarnings(stats::profile(profile_fake(TRUE), param = p,
                                      values = vals, cores = 1)),
      "cannot profile 'R_log_sd' on a DSEM fit", info = p)
  }
  testthat::expect_error(
    suppressWarnings(stats::profile(profile_fake(TRUE), param = "rec_dev",
                                    values = vals, cores = 1)),
    "cannot profile 'rec_dev' on a DSEM fit")

  # The message must say where the quantity actually lives, or the refusal is a
  # dead end rather than a redirect.
  err <- tryCatch(
    suppressWarnings(stats::profile(profile_fake(TRUE), param = "sigmaR",
                                    values = vals, cores = 1)),
    error = function(e) conditionMessage(e))
  testthat::expect_match(err, "dsem_beta_z")
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
