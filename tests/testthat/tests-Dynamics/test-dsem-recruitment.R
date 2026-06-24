# Regression tests for the DSEM recruitment linkage. The DSEM example models
# converge with phase = TRUE and their objectives match stored reference values,
# guarding against drift in the vendored C++ (src/TMB/dsem.hpp) and the dsem
# 3.0.0 integration in build_dsem_objects().
#
# (The former cross-validation tests against the vendored 2.0.1 RAM pipeline
# were retired when Rceattle moved to live dsem::dsem 3.0.0 harvesting and the
# vendored R pipeline R/0-dsem_ram.R was removed.)

# ---------------------------------------------------------------------------
# DSEM example models converge with phase = TRUE.
# Reference objectives reflect the dsem 3.0.0 gmrf_project parameterization;
# update them deliberately (and note why) if the DSEM math changes.
# ---------------------------------------------------------------------------
testthat::test_that("single-species DSEM models converge with phase = TRUE", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_on_cran()

  data("GOApollock", package = "Rceattle")
  GOApollock$projyr <- 2020
  GOApollock$env_data <- GOApollock$env_data %>%
    dplyr::mutate(ScaledBT = scale(BTempC))

  sem_iid <- "ScaledBT -> ScaledBT, 1, AR_BT, 0
              ScaledBT -> recdevs1, 1, BT_to_R, 0
              recdevs1 <-> recdevs1, 0, sigmaR1, 1"

  # IID recruitment (default DSEM), phase = TRUE -- the path that surfaced the
  # original `_ln_sd` map bug.
  m_iid <- Rceattle::fit_mod(data_list = GOApollock, inits = NULL, file = NULL,
                             estimateMode = 0, random_rec = TRUE, msmMode = 0,
                             initMode = 2,
                             fit_control = fit_control(phase = TRUE, verbose = 0))
  testthat::expect_s3_class(m_iid, "Rceattle")
  testthat::expect_true(is.finite(as.numeric(m_iid$opt$objective)))
  testthat::expect_lt(max(abs(m_iid$sdrep$gradient.fixed)), 1e-2)
  # Reference objective: dsem 3.0.0 gmrf_project parameterization.
  testthat::expect_equal(as.numeric(m_iid$opt$objective), 615.78939, tolerance = 1e-3)

  # Env-linked DSEM, phase = TRUE. The env variable (ScaledBT) uses the default
  # family = "fixed" (measured exactly), so the test exercises the AR_BT /
  # BT_to_R -> recruitment linkage and R_sd sourcing without the weakly
  # identified observation-error estimation that family = "normal" would add
  # under dsem 3.0.0 (which, unlike Rceattle's old fixed-SD override, leaves
  # lnsigma_z free).
  m_env <- Rceattle::fit_mod(data_list = GOApollock, inits = NULL, file = NULL,
                             estimateMode = 0,
                             dsem = build_DSEM(sem = sem_iid, family = "fixed"),
                             random_rec = TRUE, msmMode = 0, initMode = 2,
                             fit_control = fit_control(phase = TRUE, verbose = 0))
  testthat::expect_s3_class(m_env, "Rceattle")
  testthat::expect_true(is.finite(as.numeric(m_env$opt$objective)))
  testthat::expect_lt(max(abs(m_env$sdrep$gradient.fixed)), 1e-2)
  # Reference objective: dsem 3.0.0 gmrf_project parameterization, with R_sd
  # sourced from the sigmaR1 beta_z entry (not the first beta_z).
  testthat::expect_equal(as.numeric(m_env$opt$objective), 672.25545, tolerance = 1e-3)
})
