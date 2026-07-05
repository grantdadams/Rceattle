
testthat::test_that("Bias adjustment of data and process likelihoods works", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  nages = 5
  R0 = 3.5
  nyrs = 8
 # source(file.path("tests", "testthat", "helpers.R"))
  dat <- make_test_data(nyrs = nyrs, nages = nages, seed = 42)

  out <- list()
  for(obs in c(TRUE,FALSE)){
    for(proc in c(TRUE,FALSE)){
      m <- Rceattle::fit_mod(data_list = dat,
                             fit_control = fit_control(
                               bias_adjust_obs = obs,
                               bias_adjust_proc = proc,
                               phase = TRUE, getsd = FALSE,
                               verbose = 0))
      R0 <- m$estimated_params$rec_pars[1]
      out <- rbind(out, data.frame(obs=obs, proc=proc, logssb_terminal=m$quantities$log_ssb[1,18], R0=R0))
    }
  }
  # sigmaR = 1 so the process bias adjustment is -0.5; with proc = TRUE, R0 is
  # offset that much to counteract, but terminal SSB is unchanged. (OSA is no
  # longer a fit-time option -- osa_residuals() builds its data on demand -- so
  # the earlier osa = TRUE/FALSE loop is dropped.)
  # Rows: 1 = obs T/proc T, 2 = obs T/proc F, 3 = obs F/proc T, 4 = obs F/proc F.
  testthat::expect_equal(out[3,3], out[4,3], tolerance = 1e-5)        # logssb: proc on vs off
  testthat::expect_equal(out[3,4], out[4,4] + 0.5, tolerance = 1e-5)  # R0 offset by -sigmaR^2/2
  testthat::expect_equal(out[1,3], 4.227271, tolerance = 1e-5)        # logssb, obs+proc on
  testthat::expect_equal(out[3,3], 4.223244, tolerance = 1e-5)        # logssb, proc on only
})
