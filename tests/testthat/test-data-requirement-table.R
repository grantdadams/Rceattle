# Locks the Step-1 refactor: data_check()'s conditional presence-requirement
# gates now live in the declarative requirement table (R/1-data_requirements_
# table.R) and data_check() consumes it. These tests assert the EXACT messages
# and severities the pre-refactor code produced, so any drift is caught.
#
# Fast (no fit): exercises clean_data -> switch_check -> data_check only.

# Minimal single-species baseline; Logistic selectivity keeps the emp_sel gate
# quiet unless a test deliberately trips it.
make_req_dat <- function() {
  dat <- make_test_data(nyrs = 5, nages = 4, seed = 1)
  dat$msmMode <- 0
  dat$fleet_control$Selectivity <- "Logistic"
  dat
}

prep_req <- function(dat) {
  dat <- Rceattle::clean_data(dat)
  dat <- Rceattle::switch_check(dat)
  if (length(dat$HCR) == 0) dat$HCR <- "NoFishing"
  dat
}

# Run data_check() and return the (single, collapsed) error message, or the
# sentinel "<<NO ERROR>>".
check_err <- function(dat) {
  tryCatch({
    suppressWarnings(suppressMessages(data_check(dat)))
    "<<NO ERROR>>"
  }, error = function(e) conditionMessage(e))
}

testthat::test_that("requirement table has well-formed rows", {
  tbl <- Rceattle:::.rce_requirement_table()
  testthat::expect_true(length(tbl) >= 8)
  # Keyed by element; every row is either always_required or carries a
  # required_when predicate; driven rows carry severity + adequate + message.
  for (nm in names(tbl)) {
    row <- tbl[[nm]]
    testthat::expect_identical(row$element, nm)
    # A row is always-required, conditionally-required, or purely optional.
    testthat::expect_true(
      isTRUE(row$always_required) ||
        is.function(row$required_when) ||
        identical(row$optional_status, "defaulted"))
    if (isTRUE(row$driven)) {
      testthat::expect_true(row$severity %in% c("error", "message"))
      testthat::expect_true(is.function(row$adequate))
      testthat::expect_true(is.function(row$message))
    }
  }
  # The six driven gates are present.
  driven <- names(Filter(function(r) isTRUE(r$driven), tbl))
  testthat::expect_setequal(
    driven,
    c("diet_data", "ration_data", "bioenergetics",
      "emp_sel", "NByageFixed", "caal_data")
  )
})

testthat::test_that("driven error gates emit exact pre-refactor messages", {
  testthat::skip_if_not_installed("Rceattle")

  # diet_data (msmMode > 0)
  d <- prep_req(make_req_dat()); d$msmMode <- 1L; d$diet_data <- NULL
  testthat::expect_identical(check_err(d), "No diet data included")

  # NByageFixed (estDynamics > 0)
  d <- prep_req(make_req_dat()); d$estDynamics <- rep(1L, d$nspp); d$NByageFixed <- NULL
  testthat::expect_identical(
    check_err(d), "estDynamics > 0 requires NByageFixed data to be provided")

  # caal_data (growth_model > 0)
  d <- prep_req(make_req_dat()); d$growth_model <- rep(1L, d$nspp); d$caal_data <- NULL
  testthat::expect_identical(
    check_err(d),
    "Growth estimation (growth_model > 0) requires caal_data to be provided")

  # emp_sel (Selectivity == "Fixed"): message names the fleet
  d <- prep_req(make_req_dat()); d$fleet_control$Selectivity[1] <- "Fixed"; d$emp_sel <- NULL
  testthat::expect_match(
    check_err(d), "Fleet\\(s\\) with Selectivity = 'Fixed' \\(.*\\) require emp_sel data")

  # bioenergetics group (msmMode > 0): message lists the missing scalar
  d <- prep_req(make_req_dat()); d$msmMode <- 1L; d$Ceq <- NULL
  testthat::expect_match(
    check_err(d),
    "msmMode > 0 requires bioenergetics scalars of length nspp = 1; missing or wrong length:.*Ceq")
})

testthat::test_that("ration_data requirement is a message, not an error", {
  testthat::skip_if_not_installed("Rceattle")

  d <- prep_req(make_req_dat()); d$msmMode <- 1L; d$ration_data <- NULL
  # Provide diet + bioenergetics so ONLY the ration notice would fire.
  msgs <- character(0)
  withCallingHandlers(
    tryCatch(data_check(d), error = function(e) invisible()),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
  )
  testthat::expect_true(any(grepl("No ration data", msgs)))
})

testthat::test_that("requirements stay dormant when their condition is absent", {
  testthat::skip_if_not_installed("Rceattle")

  # Single-species, Logistic sel, all optional blocks NULL: none of the driven
  # gates fire (they are conditional on msmMode / estDynamics / growth / Fixed).
  d <- prep_req(make_req_dat())
  for (nm in c("diet_data", "ration_data", "caal_data", "emp_sel", "NByageFixed")) d[[nm]] <- NULL
  testthat::expect_identical(check_err(d), "<<NO ERROR>>")
})
