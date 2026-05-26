# Tests for SS3-pattern-24 DoubleNormal selectivity (type 8)
#
# Rceattle's case 8 implements SS3 pattern 24 (length-based DoubleNormal with
# plateau) exactly. 6 parameters per fleet/sex:
#
#   sel_inf[1]     = peak                    (SS3 P1)
#   sel_inf[2]     = logit(right floor)      (SS3 P6 / final)
#   sel_inf[3]     = logit(left floor)       (SS3 P5 / init)
#   log_sel_slp[1] = log(sigma_asc)          (SS3 P3) -- denom in exp(-(L-peak)^2 / .)
#   log_sel_slp[2] = log(sigma_desc)         (SS3 P4) -- denom in exp(-(L-peak2)^2 / .)
#   log_sel_slp[3] = top-width logit         (SS3 P2): peak2 = peak + binwidth
#                     + (xmax - peak - binwidth) / (1 + exp(-this))
#
# The legacy 4-param formulation is preserved by setting the new slots to
# their "off" sentinels: log_sel_slp[3] = -10 (peak2 ≈ peak + binwidth, no
# plateau) and sel_inf[3] = -10 (init ≈ 0, full Gaussian ascending limb).
# Note this is NOT identical to the pre-pattern-24 Rceattle formula (which
# used 0.5*(L-peak)^2/sigma^2 without the SS3 normalization anchors), so the
# shape and normalization differ slightly even in the degenerate case.


# Helper: SS3 pattern 24 DoubleNormal reference (R equivalent of TMB case 8)
double_normal_sel <- function(x, peak, sigma_asc, sigma_desc,
                              logit_final = -10, logit_init = -10,
                              top_lt = -10, binwidth = 1) {
  init   <- 1 / (1 + exp(-logit_init))
  finalv <- 1 / (1 + exp(-logit_final))
  xmin   <- x[1]
  xmax   <- x[length(x)]
  upselex   <- sigma_asc   # SS3 stores exp(P3); we pass sigma_asc directly.
  downselex <- sigma_desc  # Same for SS3 P4.
  peak2  <- peak + binwidth + (xmax - peak - binwidth) / (1 + exp(-top_lt))
  t1min  <- exp(-(xmin - peak ) * (xmin - peak ) / upselex)
  t2min  <- exp(-(xmax - peak2) * (xmax - peak2) / downselex)
  out <- vapply(x, function(L) {
    asc        <- exp(-(L - peak ) * (L - peak ) / upselex)
    dsc        <- exp(-(L - peak2) * (L - peak2) / downselex)
    asc_scaled <- init + (1 - init) * (asc - t1min) / (1 - t1min)
    dsc_scaled <- 1 + (finalv - 1) * (dsc - 1) / (t2min - 1)
    denom1 <- 1 + abs(L - peak )
    denom2 <- 1 + abs(L - peak2)
    join1  <- 1 / (1 + exp(-(20 / denom1) * (L - peak )))
    join2  <- 1 / (1 + exp(-(20 / denom2) * (L - peak2)))
    asc_scaled * (1 - join1) + join1 * ((1 - join2) + dsc_scaled * join2)
  }, numeric(1))
  out / max(out)   # max-normalise (matches Sel_norm_bin1 < 0 path)
}


testthat::test_that("Age-based double-normal: dome shape (no init, no plateau) recovers correctly", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Fleet_type[-rows_use] <- 0
  GOA2018SS$fleet_control$Selectivity <- "DoubleNormal"
  GOA2018SS$fleet_control$Time_varying_sel   <- 0
  GOA2018SS$fleet_control$Selectivity_index  <- seq_len(nrow(GOA2018SS$fleet_control))
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin1      <- NA  # max-normalise

  peak        <- 8
  sigma_asc   <- 3
  sigma_desc  <- 6
  logit_final <- -10  # right floor ≈ 0 → fully dome-shaped
  logit_init  <- -10  # left floor  ≈ 0 → no left tail
  top_lt      <- -10  # peak2 ≈ peak + 1 → no plateau
  ages        <- 1:21
  expected    <- double_normal_sel(ages, peak, sigma_asc, sigma_desc,
                                   logit_final, logit_init, top_lt, binwidth = 1)

  mod0 <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )
  inits <- mod0$estimated_params
  inits$sel_inf[1, , ]     <- peak
  inits$sel_inf[2, , ]     <- logit_final
  inits$sel_inf[3, , ]     <- logit_init
  inits$log_sel_slp[1, , ] <- log(sigma_asc)
  inits$log_sel_slp[2, , ] <- log(sigma_desc)
  inits$log_sel_slp[3, , ] <- top_lt

  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS, inits = inits, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )

  # --- Map checks: all six base params estimated, no time-varying deviates ---
  for (j in 1:3) {
    testthat::expect_true(all(!is.na(ss_run$map$mapList$sel_inf[j, rows_use, 1])))
    testthat::expect_true(all(!is.na(ss_run$map$mapList$log_sel_slp[j, rows_use, 1])))
    testthat::expect_true(all(!is.na(ss_run$map$mapList$sel_inf[j, 9:11, 2])))
    testthat::expect_true(all(!is.na(ss_run$map$mapList$log_sel_slp[j, 9:11, 2])))
  }
  testthat::expect_true(all(is.na(c(ss_run$map$mapList$log_sel_slp_dev))))
  testthat::expect_true(all(is.na(c(ss_run$map$mapList$sel_inf_dev))))

  # --- Shape: Pollock fleet 1, ATF fleet 9 ---
  sel_out <- as.numeric(ss_run$quantities$sel_at_age[1, 1, 1:10, 1])
  testthat::expect_equal(sel_out, expected[1:10], tolerance = 1e-4)

  sel_atf <- as.numeric(ss_run$quantities$sel_at_age[9, 1, 1:21, 1])
  testthat::expect_equal(sel_atf, expected[1:21], tolerance = 1e-4)
})


testthat::test_that("Age-based double-normal: logit_final=+10 produces near-logistic shape", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Fleet_type[-rows_use] <- 0
  GOA2018SS$fleet_control$Selectivity <- "DoubleNormal"
  GOA2018SS$fleet_control$Time_varying_sel   <- 0
  GOA2018SS$fleet_control$Selectivity_index  <- seq_len(nrow(GOA2018SS$fleet_control))
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin1      <- NA

  peak        <- 8
  sigma_asc   <- 3
  sigma_desc  <- 6
  logit_final <- 10   # right floor ≈ 1 → descending side flat at 1 (logistic-like)
  logit_init  <- -10  # left floor  ≈ 0
  top_lt      <- -10  # no plateau
  ages        <- 1:21
  expected    <- double_normal_sel(ages, peak, sigma_asc, sigma_desc,
                                   logit_final, logit_init, top_lt, binwidth = 1)

  mod0 <- suppressMessages(
    fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )
  inits <- mod0$estimated_params
  inits$sel_inf[1, , ]     <- peak
  inits$sel_inf[2, , ]     <- logit_final
  inits$sel_inf[3, , ]     <- logit_init
  inits$log_sel_slp[1, , ] <- log(sigma_asc)
  inits$log_sel_slp[2, , ] <- log(sigma_desc)
  inits$log_sel_slp[3, , ] <- top_lt

  ss_run <- suppressMessages(
    fit_mod(data_list = GOA2018SS, inits = inits, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )

  # Use fleet 9 (ATF, species 2) which has 21 ages; fleet 1 (Pollock) only has 10.
  sel_out <- as.numeric(ss_run$quantities$sel_at_age[9, 1, 1:21, 1])
  testthat::expect_equal(sel_out, expected[1:21], tolerance = 1e-4)

  # With logit_final=10 the selectivity should plateau near 1 at large ages
  testthat::expect_true(min(sel_out[ages >= peak]) > 0.98)
})


testthat::test_that("Length-based double-normal: static shape non-trivial", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  # 1) Set up simulation
  nyrs = 30
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE)
  )

  # Set up Rceattle data
  simData <- sim$data_list
  simData$fleet_control$Selectivity           <- "DoubleNormal"
  simData$fleet_control$Selectivity_dimension <- "Length"
  simData$fleet_control$Time_varying_sel      <- 0
  simData$fleet_control$Selectivity_index     <- seq_len(nrow(simData$fleet_control))
  simData$fleet_control$Bin_first_selected    <- 1
  simData$fleet_control$Sel_norm_bin1         <- NA

  peak        <- 60   # cm
  sigma_asc   <- 10
  sigma_desc  <- 20
  logit_final <- -10  # dome-shaped
  logit_init  <- -10
  top_lt      <- -10

  mod0 <- suppressMessages(
    fit_mod(data_list = simData, inits = NULL, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )
  inits <- mod0$estimated_params
  inits$sel_inf[1, , ]     <- peak
  inits$sel_inf[2, , ]     <- logit_final
  inits$sel_inf[3, , ]     <- logit_init
  inits$log_sel_slp[1, , ] <- log(sigma_asc)
  inits$log_sel_slp[2, , ] <- log(sigma_desc)
  inits$log_sel_slp[3, , ] <- top_lt

  ss_run <- suppressMessages(
    fit_mod(data_list = simData, inits = inits, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )

  # All 6 base params estimated (3 sel_inf slots + 3 log_sel_slp slots)
  testthat::expect_true(all(!is.na(ss_run$map$mapList$sel_inf[, , 1])))
  testthat::expect_true(all(!is.na(ss_run$map$mapList$log_sel_slp[, , 1])))

  # sel_at_length should be non-trivial and dome-shaped
  sel_out <- as.numeric(ss_run$quantities$sel_at_length[1, 1, , 1])
  testthat::expect_true(max(sel_out) > 0.9)   # max near 1 after normalisation
  testthat::expect_true(min(sel_out) < 0.5)   # dome-shaped: tails are low
})


testthat::test_that("IID time-varying double-normal: map has correct number of deviates", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  # 1) Set up simulation
  nyrs = 30
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE)
  )

  # Set up Rceattle data
  simData <- sim$data_list
  simData$fleet_control$Selectivity              <- "DoubleNormal"
  simData$fleet_control$Time_varying_sel         <- "IID"
  simData$fleet_control$Time_varying_sel_sd_prior <- 1
  simData$fleet_control$Selectivity_index        <- seq_len(nrow(simData$fleet_control))
  simData$fleet_control$Bin_first_selected       <- 1
  simData$fleet_control$Sel_norm_bin1            <- NA
  simData$catch_data$Catch <- 1e6  # non-zero catch keeps sel devs on

  n_active <- sum(simData$fleet_control$Fleet_type != "Off")
  nsex     <- 1  # simData is single-sex

  ss_run <- suppressMessages(
    fit_mod(data_list = simData, inits = NULL, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )

  # Each active fleet now contributes nyrs deviates for each of 6 time-varying
  # params (3 in sel_inf_dev, 3 in log_sel_slp_dev) under SS3 pattern 24.
  for (j in 1:3) {
    testthat::expect_equal(
      sum(!is.na(c(ss_run$map$mapList$sel_inf_dev[j, , , ]))),
      n_active * nsex * nyrs,
      info = sprintf("sel_inf_dev slot %d", j)
    )
    testthat::expect_equal(
      sum(!is.na(c(ss_run$map$mapList$log_sel_slp_dev[j, , , ]))),
      n_active * nsex * nyrs,
      info = sprintf("log_sel_slp_dev slot %d", j)
    )
  }
})


testthat::test_that("IID time-varying double-normal: penalty likelihood matches R reference", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  # 1) Set up simulation
  nyrs = 30
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE)
  )

  # Set up Rceattle data
  simData <- sim$data_list
  simData$fleet_control$Selectivity              <- "DoubleNormal"
  simData$fleet_control$Time_varying_sel         <- "IID"
  simData$fleet_control$Time_varying_sel_sd_prior <- 1
  simData$fleet_control$Selectivity_index        <- seq_len(nrow(simData$fleet_control))
  simData$fleet_control$Bin_first_selected       <- 1
  simData$fleet_control$Sel_norm_bin1            <- NA
  simData$catch_data$Catch <- 1e6

  peak <- 8; sigma_asc <- 3; sigma_desc <- 6
  logit_final <- 0; logit_init <- -10; top_lt <- -10
  peak_devs    <- rnorm(nyrs)
  final_devs   <- rnorm(nyrs)
  init_devs    <- rnorm(nyrs)
  asc_sd_devs  <- rnorm(nyrs)
  desc_sd_devs <- rnorm(nyrs)
  top_devs     <- rnorm(nyrs)

  sd_prior <- 1  # Time_varying_sel_sd_prior

  mod0 <- suppressMessages(
    fit_mod(data_list = simData, inits = NULL, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )
  inits <- mod0$estimated_params
  inits$sel_inf[1, , ]     <- peak
  inits$sel_inf[2, , ]     <- logit_final
  inits$sel_inf[3, , ]     <- logit_init
  inits$log_sel_slp[1, , ] <- log(sigma_asc)
  inits$log_sel_slp[2, , ] <- log(sigma_desc)
  inits$log_sel_slp[3, , ] <- top_lt

  n_flt <- nrow(simData$fleet_control)
  for (i in seq_len(n_flt)) {
    inits$sel_inf_dev[1, i, 1, ]     <- peak_devs
    inits$sel_inf_dev[2, i, 1, ]     <- final_devs
    inits$sel_inf_dev[3, i, 1, ]     <- init_devs
    inits$log_sel_slp_dev[1, i, 1, ] <- asc_sd_devs
    inits$log_sel_slp_dev[2, i, 1, ] <- desc_sd_devs
    inits$log_sel_slp_dev[3, i, 1, ] <- top_devs
  }

  ss_run <- suppressMessages(
    fit_mod(data_list = simData, inits = inits, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )

  n_active <- sum(simData$fleet_control$Fleet_type != "Off")

  # NLL from TMB (jnll_comp row 6 = penalty for sel deviates, 1-indexed from R)
  rcnll <- sum(ss_run$quantities$jnll_comp[6, ])

  # R reference: for each active fleet/sex, penalise all six deviate streams.
  # Peak / final / init deviates use sd = sel_dev_sd;
  # sigma_asc / sigma_desc / top-width deviates use sd = 4 * sel_dev_sd
  # (slope-like params get the looser sd, same convention as the legacy 4-param).
  ref_nll <- n_active * (
    sum(-dnorm(peak_devs,    0, sd_prior,     log = TRUE)) +
    sum(-dnorm(final_devs,   0, sd_prior,     log = TRUE)) +
    sum(-dnorm(init_devs,    0, sd_prior,     log = TRUE)) +
    sum(-dnorm(asc_sd_devs,  0, 4 * sd_prior, log = TRUE)) +
    sum(-dnorm(desc_sd_devs, 0, 4 * sd_prior, log = TRUE)) +
    sum(-dnorm(top_devs,     0, 4 * sd_prior, log = TRUE))
  )

  testthat::expect_equal(rcnll, ref_nll, tolerance = 1e-4)
})


testthat::test_that("Block time-varying double-normal: map assigns correct block indices", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  # 1) Set up simulation
  nyrs = 30
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE)
  )

  # Set up Rceattle data
  simData <- sim$data_list

  simData$fleet_control$Selectivity              <- "DoubleNormal"
  simData$fleet_control$Time_varying_sel         <- "Block"
  simData$fleet_control$Selectivity_index        <- seq_len(nrow(simData$fleet_control))
  simData$fleet_control$Bin_first_selected       <- 1
  simData$fleet_control$Sel_norm_bin1            <- NA

  # Assign 2 blocks: years 1..half in block 1, rest in block 2
  half <- floor(nyrs / 2)
  simData$catch_data <- simData$catch_data |>
    dplyr::mutate(Selectivity_block = dplyr::if_else(
      Year - simData$styr + 1 <= half, 1L, 2L))
  simData$index_data <- simData$index_data |>
    dplyr::mutate(Selectivity_block = dplyr::if_else(
      Year - simData$styr + 1 <= half, 1L, 2L))

  ss_run <- suppressMessages(
    fit_mod(data_list = simData, inits = NULL, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )

  n_active <- sum(simData$fleet_control$Fleet_type != "Off")
  n_blocks <- 2

  # All 6 base params must be fixed (NA) when blocks are active (deviates fully replace them)
  for (j in 1:3) {
    testthat::expect_true(all(is.na(ss_run$map$mapList$sel_inf[j, , ])))
    testthat::expect_true(all(is.na(ss_run$map$mapList$log_sel_slp[j, , ])))
  }

  # Each of the 6 parameter streams should have n_active * n_blocks unique indices
  for (j in 1:3) {
    testthat::expect_equal(
      length(unique(na.omit(c(ss_run$map$mapList$sel_inf_dev[j, , , ])))),
      n_active * n_blocks,
      info = sprintf("sel_inf_dev slot %d", j)
    )
    testthat::expect_equal(
      length(unique(na.omit(c(ss_run$map$mapList$log_sel_slp_dev[j, , , ])))),
      n_active * n_blocks,
      info = sprintf("log_sel_slp_dev slot %d", j)
    )
  }
})
