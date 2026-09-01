# The non-parametric shape penalty's bin range.
#
# Provenance. The penalty pairs (bin, bin + 1), so its left bin cannot be the
# last one, but `data_check()` bounds `Sel_pen_last_bin` by the NUMBER of bins --
# so the last bin is an accepted value. Set there, the loop read bin `nbins`,
# which exists because the selectivity arrays are sized by the widest species
# (BS2017SS: nages 12, 12, 21) but which `calculate_selectivity()` never writes.
# Before 5.25.0 that entry was a natural-scale 0 and the penalty took its log,
# giving -Inf and a NaN objective; carrying the curve on the log scale turned the
# same read into a structural 0, i.e. a phantom fully-selected bin scoring a real
# penalty. The model now clamps the range.
testthat::skip_on_cran()

testthat::test_that("the type-9 shape penalty cannot read past the last bin", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("BS2017SS")
  mk <- function(last_bin) {
    d <- BS2017SS
    d$fleet_control$Selectivity <- as.character(d$fleet_control$Selectivity)
    d$fleet_control$Selectivity[1]    <- "NonParametricPM"
    d$fleet_control$Sel_curve_pen1[1] <- 20
    d$fleet_control$Sel_curve_pen2[1] <- 12.5
    d$fleet_control$Sel_curve_pen3    <- 0
    d$fleet_control$Sel_pen_last_bin  <- NA
    d$fleet_control$Sel_pen_last_bin[1] <- last_bin
    d
  }
  obj_at <- function(last_bin) {
    m <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = mk(last_bin), inits = NULL, msmMode = 0, estimateMode = 3,
      fit_control = Rceattle::fit_control(verbose = 0))))
    p <- m$obj$par
    i <- which(names(p) == "sel_coff")
    # A strongly decreasing curve, so the shape penalty has something to bite on.
    p[i] <- seq(2, -2, length.out = length(i))
    m$obj$fn(p)
  }

  # The penalty widens with the bin range, so the test is not vacuous ...
  testthat::expect_lt(obj_at(3), obj_at(8))
  # ... and then stops at the last usable pair, whatever is asked for.
  testthat::expect_equal(obj_at(11), obj_at(NA))
  testthat::expect_equal(obj_at(12), obj_at(NA))
  testthat::expect_true(is.finite(obj_at(12)))
})
