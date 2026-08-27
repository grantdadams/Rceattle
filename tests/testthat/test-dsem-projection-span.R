# build_dsem_objects(): the DSEM latent state space must span the hindcast only
# when estimate_projection = FALSE.
#
# This is a density question, not just a bookkeeping one. The DSEM likelihood is
# a GMRF quadratic form over the stacked [time x variable] state space, and
# lagged RAM paths couple each time row to its predecessors. Keeping the
# projection rows and merely mapping them off pins those states but leaves them
# in the quadratic form, so their cross-terms with the terminal hindcast years
# remain non-constant in the estimated states -- a spurious shrinkage pull on the
# terminal recruitment deviations, which propagate to terminal SSB and the ABC.
# These tests pin the shape so a future refactor cannot quietly reintroduce it.

dsem_test_data <- function(styr = 1990, endyr = 2010, projyr = 2020) {
  list(
    styr = styr, endyr = endyr, projyr = projyr, nspp = 1L,
    spnames = "spp1",
    sigma_rec_prior = 1,
    random_rec = TRUE,
    proj_mean_rec = FALSE,   # SEM-driven projection; required by estimate_projection = TRUE
    env_data = data.frame(Year = styr:projyr,
                          BT = as.numeric(scale(seq_len(projyr - styr + 1))))
  )
}

testthat::test_that("estimate_projection = FALSE builds the state space over the hindcast only", {
  testthat::skip_if_not_installed("dsem")

  d <- dsem_test_data()
  nyrs_hind <- length(d$styr:d$endyr)   # 21
  nyrs_full <- length(d$styr:d$projyr)  # 31
  testthat::expect_gt(nyrs_full, nyrs_hind)   # the test is meaningless otherwise

  obj <- build_dsem_objects(
    dsem_settings = build_DSEM(sem = "BT -> recdevs1, 1, BT_to_R, 0
                                      recdevs1 <-> recdevs1, 0, sigmaR1, 1",
                               family = "fixed", estimate_projection = FALSE),
    data_list = d
  )

  # The latent state matrix, its observations, and the map all stop at endyr.
  testthat::expect_equal(nrow(obj$tmb_inputs$parameters$x_tj), nyrs_hind)
  testthat::expect_equal(nrow(obj$tmb_inputs$data$y_tj), nyrs_hind)
  testthat::expect_equal(nrow(obj$mapList$x_tj), nyrs_hind)
  testthat::expect_equal(obj$tmb_inputs$data$nyrs_dsem, nyrs_hind)
  testthat::expect_false(obj$covers_projection)

  # Absent, not pinned: every recruitment-deviation state in the span is
  # estimated, with no trailing block of mapped-off projection rows. (The
  # covariate column IS all-NA here -- family = "fixed" pins it to its measured
  # values -- so check the recdev column specifically rather than the whole map.)
  rec_col <- obj$tmb_inputs$data$rec_dev_col[1] + 1L   # stored 0-based for C++
  testthat::expect_equal(colnames(obj$mapList$x_tj)[rec_col], "recdevs1")
  testthat::expect_false(any(is.na(obj$mapList$x_tj[, rec_col])))
  testthat::expect_equal(length(obj$mapList$x_tj[, rec_col]), nyrs_hind)
})

testthat::test_that("estimate_projection = TRUE extends the state space to projyr", {
  testthat::skip_if_not_installed("dsem")

  d <- dsem_test_data()
  nyrs_full <- length(d$styr:d$projyr)

  obj <- build_dsem_objects(
    dsem_settings = build_DSEM(sem = "BT -> recdevs1, 1, BT_to_R, 0
                                      recdevs1 <-> recdevs1, 0, sigmaR1, 1",
                               family = "fixed", estimate_projection = TRUE),
    data_list = d
  )

  testthat::expect_equal(nrow(obj$tmb_inputs$parameters$x_tj), nyrs_full)
  testthat::expect_equal(obj$tmb_inputs$data$nyrs_dsem, nyrs_full)
  testthat::expect_true(obj$covers_projection)
})

testthat::test_that("a model with no projection years is unaffected by the switch", {
  testthat::skip_if_not_installed("dsem")

  d <- dsem_test_data(projyr = 2010)   # projyr == endyr
  nyrs <- length(d$styr:d$endyr)
  sem  <- "recdevs1 <-> recdevs1, 0, sigmaR1, 1"

  off <- build_dsem_objects(
    dsem_settings = build_DSEM(sem = sem, estimate_projection = FALSE),
    data_list = d)
  on <- build_dsem_objects(
    dsem_settings = build_DSEM(sem = sem, estimate_projection = TRUE),
    data_list = d)

  testthat::expect_equal(off$tmb_inputs$data$nyrs_dsem, nyrs)
  testthat::expect_equal(on$tmb_inputs$data$nyrs_dsem, nyrs)
})

testthat::test_that("estimate_projection = TRUE with proj_mean_rec = TRUE warns", {
  testthat::skip_if_not_installed("dsem")

  d <- dsem_test_data()
  d$proj_mean_rec <- TRUE

  # Under proj_mean_rec the projection is R = avg_R, so the SEM's projection
  # states do not drive projected recruitment (they reach only the dynamic
  # B0/BF series). Estimating them for that is a configuration mistake worth
  # refusing rather than letting the SEM look active in the projection.
  testthat::expect_warning(
    build_dsem_objects(
      dsem_settings = build_DSEM(sem = "recdevs1 <-> recdevs1, 0, sigmaR1, 1",
                                 estimate_projection = TRUE),
      data_list = d),
    "do not drive it"
  )
  # ...and stays silent on the pairing that actually projects through the SEM.
  d2 <- dsem_test_data(); d2$proj_mean_rec <- FALSE
  testthat::expect_silent(
    build_dsem_objects(
      dsem_settings = build_DSEM(sem = "recdevs1 <-> recdevs1, 0, sigmaR1, 1",
                                 estimate_projection = TRUE),
      data_list = d2))
})
