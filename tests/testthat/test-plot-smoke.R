# Smoke tests for the plotting / output functions (R/7-plot_*.R), the single
# largest untested module. These don't check pixels -- they fit a small model
# once and assert every exported plotter runs without error and draws to a
# device. Output is sent to a throwaway PDF so nothing is written to the repo.
# (Requested in discussion #75 by Cole-Monnahan-NOAA.)
testthat::skip_on_cran()

# Open a throwaway graphics device for the duration of `code` so base-graphics
# plotters don't spew files or pop up windows; always closes it again.
with_null_device <- function(code) {
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}

# --- single-species plotters --------------------------------------------------
testthat::test_that("single-species plotters run on a fitted model", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("BS2017SS")
  ss <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS,
                      estimateMode = 3,  # build only, no optimisation (fast)
                      msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))

  with_null_device({
    # time-series family (accept a single fit or a list)
    testthat::expect_no_error(Rceattle::plot_biomass(ss))
    testthat::expect_no_error(Rceattle::plot_ssb(ss))
    testthat::expect_no_error(Rceattle::plot_recruitment(ss))
    testthat::expect_no_error(Rceattle::plot_depletion(ss))
    testthat::expect_no_error(Rceattle::plot_exploitable_biomass(ss))
    testthat::expect_no_error(Rceattle::plot_index(ss))
    testthat::expect_no_error(Rceattle::plot_logindex(ss))
    testthat::expect_no_error(Rceattle::plot_catch(ss))
    testthat::expect_no_error(Rceattle::plot_f(ss))

    # structural / data plotters
    testthat::expect_no_error(Rceattle::plot_selectivity(ss))
    testthat::expect_no_error(Rceattle::plot_mortality(ss))
    testthat::expect_no_error(Rceattle::plot_data(ss))
    testthat::expect_no_error(Rceattle::plot_comp(ss))

    # S3 plot() dispatch
    testthat::expect_no_error(plot(ss, what = "biomass"))
    testthat::expect_no_error(plot(ss, what = "ssb"))
    testthat::expect_no_error(plot(ss, what = "recruitment"))
    testthat::expect_no_error(plot(ss, what = "selectivity"))
    testthat::expect_no_error(plot(ss, what = "data"))
  })
})

# --- multispecies / predation plotters ---------------------------------------
testthat::test_that("multispecies predation plotters run on a fitted model", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("BS2017MS")
  ms <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017MS,
                      estimateMode = 3,
                      msmMode = 1,   # multispecies (predation on)
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))

  with_null_device({
    testthat::expect_no_error(Rceattle::plot_b_eaten(ms))
    testthat::expect_no_error(Rceattle::plot_m_at_age(ms))
    testthat::expect_no_error(plot(ms, what = "mortality"))
  })
})
