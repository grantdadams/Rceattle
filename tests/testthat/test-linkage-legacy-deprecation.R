# =============================================================================
# Legacy time-variation switches (Time_varying_q / Time_varying_sel / M1_re) are
# soft-deprecated in favour of random-effect linkages: data_check() emits a
# warning naming the grammar equivalent, but the legacy path still fits. The
# warning fires only where a grammar equivalent exists -- not for the
# environmental / Rogers-AR1 catchability overloads, the non-parametric
# selectivity forms, or the separable M1_re (age x year, code 6). These tests
# pin that (which configs warn, which stay silent).
# =============================================================================

# Does data_check() emit a warning matching `pattern` for this data list? Runs
# switch_check() first (as fit_mod does) and swallows any later standalone
# data_check error (e.g. the projection-F/HCR check that needs fit-time context)
# -- the deprecation block runs before it, so the warning is still captured.
.dep_warns <- function(data_list, pattern) {
  dl <- suppressMessages(Rceattle::switch_check(data_list))
  dc <- get("data_check", asNamespace("Rceattle"))
  w <- character(0)
  withCallingHandlers(
    tryCatch(suppressMessages(dc(dl)), error = function(e) NULL),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  any(grepl(pattern, w))
}

testthat::test_that("Time_varying_q on an estimated fleet is soft-deprecated", {
  testthat::skip_if_not_installed("Rceattle")
  d <- make_test_data()
  d$fleet_control$Catchability   <- 1L          # Estimated
  d$fleet_control$Time_varying_q <- c(1L, 0L)   # IID on the (estimated) survey q
  testthat::expect_true(.dep_warns(d, "Time_varying_q is soft-deprecated"))

  # RandomWalk also warns; Off does not.
  d$fleet_control$Time_varying_q <- c(4L, 0L)
  testthat::expect_true(.dep_warns(d, "Time_varying_q is soft-deprecated"))
  d$fleet_control$Time_varying_q <- c(0L, 0L)
  testthat::expect_false(.dep_warns(d, "Time_varying_q is soft-deprecated"))
})

testthat::test_that("Time_varying_q is NOT deprecated when it names env columns (Catchability 5/6)", {
  testthat::skip_if_not_installed("Rceattle")
  d <- make_test_data()
  # Environmental catchability overloads Time_varying_q to a column index, not a
  # time-varying mode -- must not be nudged toward a random-effect linkage.
  d$fleet_control$Catchability   <- c(5L, 1L)
  d$fleet_control$Time_varying_q <- c(1L, 0L)
  testthat::expect_false(.dep_warns(d, "Time_varying_q is soft-deprecated"))
})

testthat::test_that("Time_varying_sel on a parametric form is soft-deprecated", {
  testthat::skip_if_not_installed("Rceattle")
  d <- make_test_data()
  d$fleet_control$Selectivity      <- 1L         # Logistic (parametric)
  d$fleet_control$Time_varying_sel <- c(4L, 0L)  # RandomWalk
  testthat::expect_true(.dep_warns(d, "Time_varying_sel is soft-deprecated"))
  d$fleet_control$Time_varying_sel <- c(0L, 0L)
  testthat::expect_false(.dep_warns(d, "Time_varying_sel is soft-deprecated"))
})

testthat::test_that("M1_re modes 1-5 warn but the separable mode 6 does not", {
  testthat::skip_if_not_installed("Rceattle")
  d <- make_test_data()
  for (mode in 1:5) {
    d$M1_re <- mode
    testthat::expect_true(.dep_warns(d, "M1_re is soft-deprecated"),
                          info = paste("M1_re mode", mode))
  }
  # Mode 6 (ar1_age_year, separable) has no 1-D grammar equivalent -> no nudge.
  d$M1_re <- 6L
  testthat::expect_false(.dep_warns(d, "M1_re is soft-deprecated"))
  d$M1_re <- 0L
  testthat::expect_false(.dep_warns(d, "M1_re is soft-deprecated"))
})
