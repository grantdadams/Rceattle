# `rescale` selects WHICH env_data columns a peel standardizes, because
# env_data can carry design columns as well as covariates: an SS3 bridge writes
# the 0/1 selectivity blocks there, and standardizing one refits a warm-started
# coefficient against a column it was not fitted on. Measured on the GOA cod
# bridge, that moves a peel's starting objective from 2050 to 67109.
testthat::skip_on_cran()

testthat::test_that("rescale names the env_data columns to standardize", {
  testthat::skip_if_not_installed("TMB")

  data(BS2017SS)
  d <- BS2017SS
  yrs <- d$styr:d$endyr
  d$env_data <- data.frame(Year = yrs,
                           temp  = as.numeric(scale(seq_along(yrs))) * 3 + 10,
                           block = as.numeric(yrs > stats::median(yrs)))

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 0,
    fit_control = Rceattle::fit_control(getsd = FALSE, phase = FALSE,
                                        verbose = 0))))

  # A name that is not a column is a typo, and it fails before any peel is fit.
  testthat::expect_error(
    Rceattle::retrospective(fit, peels = 1, rescale = "no_such_column"),
    "not in `env_data`")

  # Naming one column standardizes it and leaves the other alone. One peel,
  # scored on the env_data the peel is built with rather than on a refit.
  r <- suppressWarnings(suppressMessages(
    Rceattle::retrospective(fit, peels = 1, rescale = "temp", cores = 1,
                            fit_control = Rceattle::fit_control(getsd = FALSE,
                                                                phase = FALSE))))
  peeled <- r$Rceattle_list[[1]]$data_list$env_data
  keep   <- peeled$Year <= max(peeled$Year) - 1
  # `block` keeps its 0/1 values. (Post-peel rows are NA for a column no
  # linkage reads -- that is the blanking step, not the rescaling.)
  testthat::expect_true(all(peeled$block[keep] %in% c(0, 1)))
  # `temp` is centred on the retained years.
  testthat::expect_lt(abs(mean(peeled$temp[keep])), 1e-8)
})
