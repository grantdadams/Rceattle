# Integration test (fits a CEATTLE model): skipped on CRAN.
#
# build_growth(sd_plus_group=) controls the oldest age class's SD-at-age:
# "WHAM" (default) pins it to exp(sd_Linf); "SS3" interpolates it by length like
# an interior age. The switch touches ONLY the plus-group SD -- mean length-at-age
# is unchanged, interior-age weights are identical, and only the plus-group weight
# (and its length-transition probabilities) move.
testthat::skip_on_cran()

testthat::test_that("sd_plus_group WHAM vs SS3 changes only the plus-group SD-at-age", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  nyrs <- 20
  nspp <- 2
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    Fmort = matrix(seq(0.05, 0.3, length.out = nyrs), nspp, nyrs, byrow = TRUE),
    use_size_sel = TRUE,
    fish_CAAL_ISS = 1e6,
    srv_CAAL_ISS = 1e6,
    log_phi = matrix(-Inf, nspp, nspp)
  )
  simData <- sim$data_list

  # Template parameters, then DISTINCT SD endpoints so the SS3 length
  # interpolation genuinely differs from the WHAM pin at the plus group
  # (equal endpoints would make the two identical).
  mod0 <- suppressMessages(fit_mod(
    data_list = simData, inits = NULL, estimateMode = 3,
    growthFun = build_growth(fun = "vonBertalanffy"),
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0)))
  inits <- mod0$estimated_params
  inits$log_growth_pars[, 1, ] <- matrix(log(c(0.3, 4.5, 90, 1.0,
                                               0.35, 4.5, 50, 1.0)), # K, L1, Linf, m
                                         nrow = nspp, ncol = 4, byrow = TRUE)
  inits$growth_log_sd[, , 1] <- log(2)  # SD at L1
  inits$growth_log_sd[, , 2] <- log(6)  # SD at Linf (distinct from L1)

  fit_style <- function(style) suppressMessages(fit_mod(
    data_list = simData, inits = inits, file = NULL, estimateMode = 3,
    growthFun = build_growth(fun = "vonBertalanffy", sd_plus_group = style),
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0)))

  wham <- fit_style("WHAM")
  ss3  <- fit_style("SS3")

  # The default (unspecified) is WHAM.
  default <- suppressMessages(fit_mod(
    data_list = simData, inits = inits, file = NULL, estimateMode = 3,
    growthFun = build_growth(fun = "vonBertalanffy"),
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0)))
  testthat::expect_equal(default$quantities$weight_hat, wham$quantities$weight_hat)

  # Mean length-at-age is a function of the growth curve only, not the SD switch.
  testthat::expect_equal(wham$quantities$length_hat, ss3$quantities$length_hat)

  # Species 1, biomass weight (index 1), sex 1, first year.
  nages1 <- simData$nages[1]
  w_wham <- wham$quantities$weight_hat[1, 1, seq_len(nages1), 1]
  w_ss3  <- ss3$quantities$weight_hat[1, 1, seq_len(nages1), 1]

  # Interior ages: identical (SD unchanged there).
  testthat::expect_equal(w_wham[seq_len(nages1 - 1)], w_ss3[seq_len(nages1 - 1)])
  # Plus group: differs (SS3 interpolates a smaller SD than the WHAM pin).
  testthat::expect_false(isTRUE(all.equal(w_wham[nages1], w_ss3[nages1])))
})

testthat::test_that("a refit that rebuilds build_growth() without sd_plus_group inherits the SS3 choice", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  nyrs <- 20
  nspp <- 2
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    Fmort = matrix(seq(0.05, 0.3, length.out = nyrs), nspp, nyrs, byrow = TRUE),
    use_size_sel = TRUE, fish_CAAL_ISS = 1e6, srv_CAAL_ISS = 1e6,
    log_phi = matrix(-Inf, nspp, nspp))
  simData <- sim$data_list
  mod0 <- suppressMessages(fit_mod(
    data_list = simData, inits = NULL, estimateMode = 3,
    growthFun = build_growth(fun = "vonBertalanffy"),
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0)))
  inits <- mod0$estimated_params
  inits$log_growth_pars[, 1, ] <- matrix(log(c(0.3, 4.5, 90, 1.0,
                                               0.35, 4.5, 50, 1.0)),
                                         nrow = nspp, ncol = 4, byrow = TRUE)
  inits$growth_log_sd[, , 1] <- log(2)
  inits$growth_log_sd[, , 2] <- log(6)

  ss3 <- suppressMessages(fit_mod(
    data_list = simData, inits = inits, estimateMode = 3,
    growthFun = build_growth(fun = "vonBertalanffy", sd_plus_group = "SS3"),
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0)))
  testthat::expect_equal(unname(ss3$data_list$growth_sd_style), rep(2L, nspp))

  # Refit exactly the way the MSE / retrospective blocks do: rebuild growth via
  # build_growth(fun=) with NO sd_plus_group. The SS3 choice must be inherited
  # from data_list, not silently reset to WHAM.
  refit <- suppressMessages(fit_mod(
    data_list = ss3$data_list, inits = ss3$estimated_params, estimateMode = 3,
    growthFun = build_growth(fun = "vonBertalanffy"),
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0)))
  testthat::expect_equal(unname(refit$data_list$growth_sd_style), rep(2L, nspp))
  # ...and the growth output is the SS3 output, not the WHAM output.
  testthat::expect_equal(refit$quantities$weight_hat, ss3$quantities$weight_hat)
})
