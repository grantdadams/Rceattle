# retrospective(), profile() and run_mse() refit with the fit's own stored $map.
# A fit saved before 5.37.0 names index_q_rho there, and one saved before 5.35.0
# sizes log_pop_scalar by age; a supplied map used as given then failed the
# dimension check ("Map and parameter objects are not the same size for: NA").
# The GOA multispecies assessment's saved 2026 models carry both. fit_mod() now
# drops a retired block and collapses log_pop_scalar, as it does for inits.

testthat::skip_on_cran()

testthat::test_that("a stored map naming a retired block, or sized by age, still refits", {
  d <- make_test_data()
  if (!is.null(d$data_list)) d <- d$data_list
  build <- function(inits = NULL, map = NULL) suppressWarnings(suppressMessages(fit_mod(
    data_list = d, inits = inits, map = map, estimateMode = 3, msmMode = 0,
    random_rec = FALSE, fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  m0 <- build()
  old_map <- m0$map
  n_flt <- nrow(d$fleet_control); nspp <- d$nspp; nages <- max(d$nages)
  old_map$mapList$index_q_rho   <- rep(NA_real_, n_flt)
  old_map$mapFactor$index_q_rho <- factor(rep(NA, n_flt))
  old_map$mapList$log_pop_scalar   <- matrix(NA_real_, nspp, nages)
  old_map$mapFactor$log_pop_scalar <- factor(rep(NA, nspp * nages))
  inits <- m0$estimated_params
  inits$index_q_rho <- rep(0, n_flt)
  m1 <- build(inits = inits, map = old_map)
  testthat::expect_equal(m1$obj$fn(), m0$obj$fn())
  testthat::expect_false("index_q_rho" %in% names(m1$map$mapList))
  testthat::expect_length(m1$map$mapList$log_pop_scalar, nspp)
})
