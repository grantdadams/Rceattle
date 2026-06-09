# Tests for the vendored DSEM sem->inputs pipeline (R/0-dsem_ram.R) and its
# integration into fit_mod. Two layers:
#   1. Cross-validation: the ported make_dsem_ram / build_dsem_inputs must match
#      dsem 2.0.1 byte-for-byte (guarded by skip_if_not_installed("dsem")).
#   2. Regression: the DSEM example models converge with phase = TRUE and their
#      objectives match stored 2.0.x reference values (guards against drift in
#      the vendored C++ / pipeline).

# ---------------------------------------------------------------------------
# 1. Cross-validation against dsem 2.0.1
# ---------------------------------------------------------------------------
testthat::test_that("ported make_dsem_ram matches dsem::make_dsem_ram", {
  testthat::skip_if_not_installed("dsem")
  testthat::skip_if_not(utils::packageVersion("dsem") >= "2.0.1")

  sems <- list(
    iid  = "recdevs1 <-> recdevs1, 0, sigmaR1, 1\n",
    env  = "
      ScaledBT -> ScaledBT, 1, AR_BT, 0
      ScaledBT -> recdevs1, 1, BT_to_R, 0
      recdevs1 <-> recdevs1, 0, sigmaR1, 1
    ",
    spp3 = "
      recdevs1 <-> recdevs1, 0, sigmaR1, 1
      recdevs2 <-> recdevs2, 0, sigmaR2, 1
      recdevs3 <-> recdevs3, 0, sigmaR3, 1
      recdevs2 -> recdevs1, 0, R2_to_R1, 0
      recdevs3 -> recdevs1, 0, R3_to_R1, 0
    "
  )
  vars <- list(iid = "recdevs1",
               env = c("recdevs1", "ScaledBT"),
               spp3 = c("recdevs1", "recdevs2", "recdevs3"))
  tt <- 1:15

  for (nm in names(sems)) {
    mine   <- Rceattle:::make_dsem_ram(sems[[nm]], times = tt, variables = vars[[nm]], quiet = TRUE)
    theirs <- dsem::make_dsem_ram(sems[[nm]], times = tt, variables = vars[[nm]], quiet = TRUE)
    testthat::expect_equal(mine$ram, theirs$ram, ignore_attr = TRUE,
                           info = paste("ram mismatch for sem:", nm))
    testthat::expect_equal(as.matrix(mine$model), as.matrix(theirs$model), ignore_attr = TRUE,
                           info = paste("model mismatch for sem:", nm))
  }
})

testthat::test_that("ported build_dsem_inputs matches dsem::dsem(build_model=FALSE)", {
  testthat::skip_if_not_installed("dsem")
  testthat::skip_if_not(utils::packageVersion("dsem") >= "2.0.1")

  set.seed(1)
  nT  <- 12
  ts1 <- stats::ts(data.frame(recdevs1 = NA_real_,
                              ScaledBT = as.numeric(scale(stats::rnorm(nT)))))
  sem <- "
    ScaledBT -> ScaledBT, 1, AR_BT, 0
    ScaledBT -> recdevs1, 1, BT_to_R, 0
    recdevs1 <-> recdevs1, 0, sigmaR1, 1
  "
  fam <- c("fixed", "normal")

  mine <- Rceattle:::build_dsem_inputs(sem, tsdata = ts1, family = fam,
                                       gmrf_parameterization = "gmrf_project",
                                       use_REML = FALSE, quiet = TRUE)
  theirs <- dsem::dsem(sem, tsdata = ts1, family = fam,
                       control = dsem::dsem_control(build_model = FALSE, use_REML = FALSE,
                                                    quiet = TRUE,
                                                    gmrf_parameterization = "gmrf_project"))

  mt <- mine$tmb_inputs
  testthat::expect_equal(mt$data$options, theirs$data$options)
  testthat::expect_equal(as.matrix(mt$data$RAM), as.matrix(theirs$data$RAM), ignore_attr = TRUE)
  testthat::expect_equal(mt$data$RAMstart, theirs$data$RAMstart)
  testthat::expect_equal(unname(mt$data$familycode_j), unname(theirs$data$familycode_j))
  for (p in c("beta_z", "lnsigma_j", "mu_j", "delta0_j", "x_tj")) {
    testthat::expect_equal(mt$parameters[[p]], theirs$parameters[[p]], ignore_attr = TRUE,
                           info = paste("parameter mismatch:", p))
  }
  for (p in names(theirs$map)) {
    testthat::expect_equal(as.character(mt$map[[p]]), as.character(theirs$map[[p]]),
                           info = paste("map mismatch:", p))
  }
  testthat::expect_equal(mt$random, theirs$random)
})

# ---------------------------------------------------------------------------
# 2. Regression: DSEM example models converge with phase = TRUE
#    Reference objectives are the dsem 2.0.x parameterization (gmrf_project).
#    NOTE: reference values are filled in from the validated port run; update
#    deliberately if the vendored dsem math is intentionally changed.
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
  # Reference objective: exact dsem 2.0.1 gmrf_project parameterization.
  # The DSEM term was verified to match standalone dsem 2.0.1 (obj$env$f at the
  # fitted parameters) to ~1e-14, so this value reflects the correct recruitment
  # marginal density (the old 1137.50 predated the exact-2.0.1 re-vendor).
  testthat::expect_equal(as.numeric(m_iid$opt$objective), 615.78939, tolerance = 1e-3)

  # Env-linked DSEM, phase = TRUE
  m_env <- Rceattle::fit_mod(data_list = GOApollock, inits = NULL, file = NULL,
                             estimateMode = 0,
                             dsem = build_DSEM(sem = sem_iid, family = "normal"),
                             random_rec = TRUE, msmMode = 0, initMode = 2,
                             fit_control = fit_control(phase = TRUE, verbose = 0))
  testthat::expect_s3_class(m_env, "Rceattle")
  testthat::expect_true(is.finite(as.numeric(m_env$opt$objective)))
  testthat::expect_lt(max(abs(m_env$sdrep$gradient.fixed)), 1e-2)
  # Reference objective: exact dsem 2.0.1 gmrf_project parameterization, with
  # R_sd sourced from the sigmaR1 beta_z entry (not the first beta_z).
  # (Old 1135.92 predated the exact-2.0.1 re-vendor; see m_iid note above.)
  testthat::expect_equal(as.numeric(m_env$opt$objective), 672.23797, tolerance = 1e-3)
})
