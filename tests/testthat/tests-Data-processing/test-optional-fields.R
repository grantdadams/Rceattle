
# Tests that the "optional" fields in data_list can genuinely be NULL when
# the user is not using the corresponding feature, and that data_check still
# enforces them under the conditions that actually require them.
#
# We exercise the pre-fit pipeline:
#   clean_data -> switch_check -> data_check -> build_params
# That is the part Phases A/B/C/D touched. End-to-end integration through
# MakeADFun is covered by the existing parameter-recovery tests.

run_data_pipeline <- function(data_list) {
  data_list <- Rceattle::clean_data(data_list)
  data_list <- Rceattle::switch_check(data_list)
  # clean_data's revert_switches() can leave HCR as character(0) when missing;
  # fit_mod normally sets it from build_hcr() before data_check runs.
  if (length(data_list$HCR) == 0) data_list$HCR <- "NoFishing"
  data_check(data_list)
  Rceattle::build_params(data_list)
  invisible(TRUE)
}

# Single-species baseline used for the positive cases. Setting Selectivity
# to Logistic everywhere keeps the Phase D emp_sel-when-Fixed check from
# firing.
make_ss_dat <- function() {
  dat <- make_test_data(nyrs = 5, nages = 4, seed = 1)
  dat$msmMode <- 0
  dat$fleet_control$Selectivity <- "Logistic"
  dat
}

testthat::test_that(
  "Phase A: optional data.frame fields can be NULL in single-species mode", {
    testthat::skip_if_not_installed("Rceattle")

    fields <- c(
      "comp_data", "caal_data", "emp_sel", "NByageFixed",
      "ration_data", "diet_data"
    )
    for (nm in fields) {
      dat <- make_ss_dat()
      dat[[nm]] <- NULL
      res <- suppressMessages(suppressWarnings(run_data_pipeline(dat)))
      testthat::expect_true(res, info = paste0("with ", nm, " = NULL"))
    }

    # All six together
    dat <- make_ss_dat()
    for (nm in fields) dat[[nm]] <- NULL
    testthat::expect_true(
      suppressMessages(suppressWarnings(run_data_pipeline(dat)))
    )
  }
)

testthat::test_that(
  "Phase B: bioenergetics scalars can be NULL in single-species mode", {
    testthat::skip_if_not_installed("Rceattle")

    bio_scalars <- c(
      "Ceq", "Cindex", "Pvalue", "fday", "CA", "CB",
      "Qc", "Tco", "Tcm", "Tcl", "CK1", "CK4"
    )
    dat <- make_ss_dat()
    for (nm in bio_scalars) dat[[nm]] <- NULL
    testthat::expect_true(
      suppressMessages(suppressWarnings(run_data_pipeline(dat)))
    )

    # After clean_data + switch_check, all of them should be present
    # length-nspp.
    dat <- make_ss_dat()
    for (nm in bio_scalars) dat[[nm]] <- NULL
    suppressMessages(suppressWarnings({
      dat <- Rceattle::clean_data(dat)
      dat <- Rceattle::switch_check(dat)
    }))
    for (nm in bio_scalars) {
      testthat::expect_equal(length(dat[[nm]]), dat$nspp, info = nm)
    }
  }
)

testthat::test_that("Phase C: env_data can be NULL", {
  testthat::skip_if_not_installed("Rceattle")

  dat <- make_ss_dat()
  dat$env_data <- NULL
  testthat::expect_true(
    suppressMessages(suppressWarnings(run_data_pipeline(dat)))
  )
})

testthat::test_that("data_check enforces conditional requirements", {
  testthat::skip_if_not_installed("Rceattle")

  # diet_data required when msmMode > 0
  dat <- make_ss_dat()
  dat$msmMode <- 1
  dat$diet_data <- NULL
  testthat::expect_error(
    suppressMessages(suppressWarnings(run_data_pipeline(dat))),
    "No diet data included"
  )

  # NByageFixed required when estDynamics > 0
  dat <- make_ss_dat()
  dat$estDynamics <- rep(1L, dat$nspp)
  dat$NByageFixed <- NULL
  testthat::expect_error(
    suppressMessages(suppressWarnings(run_data_pipeline(dat))),
    "estDynamics > 0 requires NByageFixed"
  )

  # Phase B: bioenergetics scalars required when msmMode > 0
  dat <- make_ss_dat()
  dat$msmMode <- 1
  dat$Ceq <- NULL
  testthat::expect_error(
    suppressMessages(suppressWarnings(run_data_pipeline(dat))),
    "msmMode > 0 requires bioenergetics scalars.*Ceq"
  )

  # Phase D: emp_sel required when any fleet has Selectivity = "Fixed"
  dat <- make_ss_dat()
  dat$fleet_control$Selectivity[1] <- "Fixed"
  dat$emp_sel <- NULL
  testthat::expect_error(
    suppressMessages(suppressWarnings(run_data_pipeline(dat))),
    "Selectivity = 'Fixed'.*require emp_sel data"
  )
})
