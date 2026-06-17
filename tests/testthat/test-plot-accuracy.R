# Accuracy tests for the ggplot plotters: each returns a ggplot whose `$data`
# holds exactly the model's source quantities. We cross-check against
# as.data.frame.Rceattle() -- the canonical tidy output of the same quantities
# -- so the plot and the public data extractor are guaranteed consistent.
# This is the regression guarantee for the base->ggplot rewrite ("accuracy is
# key"), complementing the render-only smoke tests.
testthat::skip_on_cran()

# --- time-series family (plot_timeseries engine) ------------------------------
testthat::test_that("time-series plotters match as.data.frame() quantities", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("BS2017SS")
  ss <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))
  sp1 <- ss$data_list$spnames[1]
  ad <- as.data.frame(ss, which = c("biomass", "ssb", "R",
                                    "exploitable_biomass", "ssb_depletion"))

  # Compare a plotter's species-1 trajectory to the same quantity from
  # as.data.frame(), aligned on the years the plot actually shows. `scale` is
  # the display rescaling the plotter applies (biomass/ssb -> million mt).
  check <- function(p, quantity, scale = 1) {
    testthat::expect_s3_class(p, "ggplot")
    pd <- p$data[p$data$Species == sp1, c("Year", "value")]
    pd <- pd[order(pd$Year), ]
    sd <- ad[ad$quantity == quantity & ad$species == sp1, c("year", "value")]
    src <- sd$value[match(pd$Year, sd$year)] / scale
    testthat::expect_equal(pd$value, src)
  }

  check(Rceattle::plot_biomass(ss),             "biomass", scale = 1e6)
  check(Rceattle::plot_ssb(ss),                 "ssb",     scale = 1e6)
  check(Rceattle::plot_recruitment(ss),         "R")
  check(Rceattle::plot_exploitable_biomass(ss), "exploitable_biomass")
  check(Rceattle::plot_ssb_depletion(ss),       "ssb_depletion")
})

# --- fit diagnostics (index/catch) --------------------------------------------
testthat::test_that("index/catch fit plots plot the model's predictions", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("BS2017SS")
  ss <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))
  fc <- ss$data_list$fleet_control

  # survey index: plotted Predicted == quantities$index_hat on the survey rows
  pidx <- Rceattle::plot_index(ss)
  testthat::expect_s3_class(pidx, "ggplot")
  idx <- ss$data_list$index_data
  idx$hat <- ss$quantities$index_hat
  surv <- fc$Fleet_code[fc$Fleet_type == "Survey"]
  idx <- idx[idx$Fleet_code %in% surv & idx$Year <= ss$data_list$endyr, ]
  pd  <- pidx$data[order(pidx$data$Fleet, pidx$data$Year), ]
  src <- idx[order(idx$Fleet_code, abs(idx$Year)), ]
  testthat::expect_equal(pd$Predicted, src$hat)

  # fishery catch: plotted Predicted == quantities$catch_hat on the fishery rows
  pcat <- Rceattle::plot_catch(ss)
  testthat::expect_s3_class(pcat, "ggplot")
  cat_d <- ss$data_list$catch_data
  cat_d$hat <- ss$quantities$catch_hat
  fish <- fc$Fleet_code[fc$Fleet_type == "Fishery"]
  cat_d <- cat_d[cat_d$Fleet_code %in% fish & cat_d$Year <= ss$data_list$endyr, ]
  pdc <- pcat$data[order(pcat$data$Fleet, pcat$data$Year), ]
  srcc <- cat_d[order(cat_d$Fleet_code, abs(cat_d$Year)), ]
  testthat::expect_equal(pdc$Predicted, srcc$hat)
})

# --- maturity -----------------------------------------------------------------
testthat::test_that("plot_maturity plots data_list$maturity at age", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("BS2017SS")
  ss <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))
  p <- Rceattle::plot_maturity(ss)
  testthat::expect_s3_class(p, "ggplot")
  mat <- ss$data_list$maturity[, -1]
  sp1 <- ss$data_list$spnames[1]
  d <- p$data[p$data$Species == sp1, ]
  d <- d[order(d$Age), ]
  src <- as.numeric(mat[1, seq_len(ss$data_list$nages[1])])
  testthat::expect_equal(d$Maturity, src)
})
