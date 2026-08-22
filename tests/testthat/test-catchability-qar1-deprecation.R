# Catchability = "AR1" is the QAR1 form (Rogers et al. 2024). It does not work:
# build_map() gates the log-q deviates on Time_varying_q %in% c("IID","AR1",
# "RandomWalk"), but under QAR1 that column holds an env_data COLUMN INDEX
# rather than a mode, so the deviates are never freed and q comes back
# constant. Nothing errors -- the fit runs and reports a time-invariant q where
# the user asked for a time-varying one, which is the failure mode this package
# cannot afford.
#
# The warning must fire on that switch and ONLY that switch:
# `Time_varying_q = "AR1"` is a different setting sharing the same string (an
# AR1 structure on an ordinary "Estimated" q) and works correctly. An earlier
# draft of the developer notes named the wrong one.

warnings_from_check <- function(d) {
  w <- character(0)
  withCallingHandlers(
    tryCatch(data_check(d), error = function(e) NULL),
    warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") })
  w
}

base_fit <- function() {
  suppressMessages(suppressWarnings(switch_check(make_test_data(nyrs = 4, nages = 3))))
}

test_that("Catchability = 'AR1' warns that QAR1 is inert, and names the fleet", {
  d <- base_fit()
  d$fleet_control$Catchability[1] <- "AR1"
  w <- warnings_from_check(d)

  testthat::expect_true(any(grepl("QAR1", w, fixed = TRUE)))
  qw <- grep("QAR1", w, value = TRUE)[1]
  # A user must be able to act on it: which fleet, what actually happens, and
  # what to write instead.
  testthat::expect_match(qw, d$fleet_control$Fleet_name[1], fixed = TRUE)
  testthat::expect_match(qw, "INERT", fixed = TRUE)
  testthat::expect_match(qw, "linkage_spec(ar1(1 | Year)", fixed = TRUE)
})

test_that("the integer spelling warns too", {
  d <- base_fit()
  d$fleet_control$Catchability[1] <- 6
  testthat::expect_true(any(grepl("QAR1", warnings_from_check(d), fixed = TRUE)))
})

test_that("Time_varying_q = 'AR1' does NOT warn -- it is a different, working switch", {
  d <- base_fit()
  d$fleet_control$Time_varying_q[1] <- "AR1"
  testthat::expect_false(any(grepl("QAR1", warnings_from_check(d), fixed = TRUE)))
})

test_that("a model using neither is silent about QAR1", {
  testthat::expect_false(any(grepl("QAR1", warnings_from_check(base_fit()), fixed = TRUE)))
})
