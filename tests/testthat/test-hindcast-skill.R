# hindcast_skill(): forecast skill across retrospective peels.
#
# This answers a different question from Mohn's rho. Rho measures how the
# ESTIMATE of a year moves as data accumulate -- estimation consistency. This
# measures how well a peel PROJECTED years it could not see, given the catch
# that was actually taken, which is what a recruitment-projection assumption
# (proj_mean_rec TRUE/FALSE, or a DSEM) should be judged on.

testthat::test_that("hindcast_skill scores peels against the full model", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
    estimateMode = 0, random_rec = TRUE, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))
  hs <- suppressWarnings(suppressMessages(
    Rceattle::hindcast_skill(fit, peels = 2, quantity = "ssb",
                             reference = "model", cores = 1, getsd = FALSE)))

  testthat::expect_true(all(c("mase", "by_year") %in% names(hs)))
  testthat::expect_true(all(c("peel", "species", "quantity", "n_years",
                              "mae_forecast", "mae_naive", "mase") %in%
                              names(hs$mase)))
  testthat::expect_gt(nrow(hs$mase), 0)
  testthat::expect_true(all(is.finite(hs$mase$mase) | is.na(hs$mase$mase)))

  # The scored years must be the peeled ones, and only those: scoring a year the
  # peel actually fitted would compare the model to itself and flatter it.
  endyr <- fit$data_list$endyr
  testthat::expect_true(all(hs$by_year$year > endyr - hs$by_year$peel))
  testthat::expect_true(all(hs$by_year$year <= endyr))
  testthat::expect_equal(hs$by_year$years_ahead,
                         hs$by_year$year - (endyr - hs$by_year$peel))

  # The reference is the FULL model's estimate at those years, not the peel's.
  yr1 <- hs$by_year[hs$by_year$peel == 1, ]
  sp1 <- yr1[yr1$species == fit$data_list$spnames[1], ]
  col <- sp1$year - fit$data_list$styr + 1L
  testthat::expect_equal(sp1$reference,
                         as.numeric(fit$quantities$ssb[1, col]),
                         tolerance = 1e-8)

  # ...and the naive forecast is one value held flat, so it must not vary within
  # a peel x species block.
  testthat::expect_equal(length(unique(sp1$naive)), 1L)

  # MASE is the ratio of the two mean absolute errors it reports.
  m <- hs$mase[hs$mase$peel == 1, ][1, ]
  testthat::expect_equal(m$mase, m$mae_forecast / m$mae_naive, tolerance = 1e-10)
})

testthat::test_that("hindcast_skill can score the held-out index instead", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
    estimateMode = 0, random_rec = TRUE, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))
  hs <- suppressWarnings(suppressMessages(
    Rceattle::hindcast_skill(fit, peels = 2, reference = "observed",
                             cores = 1, getsd = FALSE)))

  idx <- hs$mase[grepl("^index_fleet", hs$mase$quantity), ]
  testthat::expect_gt(nrow(idx), 0)

  # The reference here must be the OBSERVED index, not a model quantity -- every
  # scored value has to appear in index_data at that fleet and year.
  by <- hs$by_year[grepl("^index_fleet", hs$by_year$quantity), ]
  obs <- fit$data_list$index_data
  hit <- vapply(seq_len(min(nrow(by), 20L)), function(i) {
    any(abs(obs$Observation[obs$Year == by$year[i]] - by$reference[i]) < 1e-8)
  }, logical(1))
  testthat::expect_true(all(hit))

  # And the peel must be predicting years it did not fit.
  testthat::expect_true(all(by$years_ahead >= 1))
})
