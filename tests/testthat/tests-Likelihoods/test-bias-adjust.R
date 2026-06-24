
testthat::test_that("Bias adjustment of data and process likelihoods works", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  nages = 5
  R0 = 3.5
  nyrs = 8
 # source(file.path("tests", "testthat", "helpers.R"))
  dat <- make_test_data(nyrs = nyrs, nages = nages, seed = 42)

  out <- list()
  for(osa in c(FALSE,TRUE)){
    for(obs in c(TRUE,FALSE)){
      for(proc in c(TRUE,FALSE)){
        m <- Rceattle::fit_mod(data_list = dat,
                               fit_control = fit_control(
                                 bias_adjust_obs = obs,
                                 bias_adjust_proc = proc,
                                 osa=osa,
                                 phase = TRUE, getsd = FALSE,
                                 verbose = 0))
        R0 <- m$estimated_params$rec_pars[1]
        out <- rbind(out, data.frame(obs=obs, proc=proc, osa=osa, logssb_terminal=m$quantities$log_ssb[1,18], R0=R0))
      }
    }
  }
  out
  # first check that with OSA off it performs as expected,
  # sigmaR=1 so adjustment is -0.5, R0 is offset that much to
  # counteract when proc=TRUE
  testthat::expect_equal(out[3,4], out[4,4], tolerance = 1e-5)
  testthat::expect_equal(out[3,5], out[4,5]+0.5, tolerance = 1e-5)
  testthat::expect_equal(out[1,4], 4.227271, tolerance=1e-5)
  testthat::expect_equal(out[3,4], 4.223244, tolerance=1e-5)
  # now test that OSA does not change anything
  testthat::expect_equal(sum(out[1:4,4:5]),sum(out[5:8,4:5]))
})
