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

  # estDynamics = 3 ESTIMATED an age-specific scalar, so its map factor carries a
  # level per age. Subsetting keeps every level, and TMB reads the levels as the
  # estimated blocks: without droplevels() this stopped with
  # "setequal(blocks, unique(L_block)) is not TRUE".
  v <- rep(NA_real_, nspp * nages); v[1] <- 1; v[nspp + 1] <- 2
  by_age <- old_map
  by_age$mapList$log_pop_scalar   <- matrix(v, nspp, nages)
  by_age$mapFactor$log_pop_scalar <- factor(v)
  m2 <- build(inits = inits, map = by_age)
  testthat::expect_length(m2$map$mapList$log_pop_scalar, nspp)
  testthat::expect_equal(levels(m2$map$mapFactor$log_pop_scalar), "1")
  testthat::expect_equal(sum(names(m2$obj$par) == "log_pop_scalar"), 1L)
})

testthat::test_that("a map name the model has no parameter for is dropped with a warning", {
  d <- make_test_data()
  if (!is.null(d$data_list)) d <- d$data_list
  build <- function(map = NULL, ...) suppressMessages(fit_mod(
    data_list = d, inits = NULL, map = map, estimateMode = 3, msmMode = 0,
    random_rec = FALSE, fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE), ...))
  m0 <- suppressWarnings(build())
  typo <- m0$map
  typo$mapList$index_log_Q   <- rep(NA_real_, nrow(d$fleet_control))
  typo$mapFactor$index_log_Q <- factor(rep(NA, nrow(d$fleet_control)))
  testthat::expect_warning(build(map = typo), "index_log_Q")
  # A refit (quiet_data_check) stays silent: it already warned on the first fit.
  testthat::expect_no_warning(build(map = typo, quiet_data_check = TRUE),
                              message = "index_log_Q")
})
