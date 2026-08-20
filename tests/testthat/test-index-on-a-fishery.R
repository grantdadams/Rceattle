# A fishery can carry index_data -- fishery CPUE -- and the template already fits
# it, on that fleet's own selectivity, which is what CPUE should share. Two
# places keyed on Fleet_type == "Survey" instead of on the data, so such a fleet
# was fitted with its catchability frozen at Catchability_init and its index was
# missing from plot_index(), both silently.
testthat::skip_on_cran()

.fishery_index_fixture <- function(estimate_q = "Estimated", distinct_sel = FALSE) {
  dat <- make_test_data(nyrs = 15, nages = 5, seed = 42)
  fc  <- dat$fleet_control
  fsh <- which(fc$Fleet_type == "Fishery")[1]
  srv <- which(fc$Fleet_type == "Survey")[1]
  testthat::skip_if(is.na(fsh) || is.na(srv), "fixture lacks both fleet types")

  # These three have no schema default, so a fishery row arrives NA on all of
  # them; data_check() now requires them wherever an index is fitted.
  dat$fleet_control$Catchability[fsh]      <- estimate_q
  dat$fleet_control$Catchability_init[fsh] <- 1
  dat$fleet_control$Estimate_index_sd[fsh] <- 0
  # Give the fishery a different selectivity FORM from the survey. Both fleets
  # otherwise converge to the same logistic shape, and then nothing downstream
  # can show which block the index was predicted from.
  if (distinct_sel) dat$fleet_control$Selectivity[fsh] <- "DoubleLogistic"
  # Give the fishery its own index series: copy the survey's rows onto the
  # fishery's Fleet_code so the values are on a sane scale for the model.
  idx <- dat$index_data
  add <- idx[idx$Fleet_code == fc$Fleet_code[srv], , drop = FALSE]
  add$Fleet_code <- fc$Fleet_code[fsh]
  add$Fleet_name <- fc$Fleet_name[fsh]
  dat$index_data <- rbind(idx, add)
  list(data = dat, fsh = fc$Fleet_code[fsh], srv = fc$Fleet_code[srv])
}

testthat::test_that("a fishery carrying index data gets an estimable catchability", {
  testthat::skip_if_not_installed("TMB")
  f <- .fishery_index_fixture()

  # Read the map off a real build: build_map() expects a data_list that has been
  # through the whole fit_mod() pipeline, not just switch_check().
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    f$data, file = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))
  q <- fit$map$mapList$index_log_q

  # The fishery's q is now a free parameter, as Catchability = "Estimated" says.
  # Before this it was NA whatever Catchability said, because the block was
  # entered only for Fleet_type == "Survey".
  testthat::expect_false(is.na(q[f$fsh]))
  # The survey's mapping still follows its own Catchability, which this fixture
  # leaves at "Fixed" -- so NA is the right answer there, and the point is that
  # widening the gate did not disturb it.
  testthat::expect_true(is.na(q[f$srv]))

  # A fishery with NO index rows keeps its q mapped out -- the change is keyed on
  # the data, not on the fleet type, so it stays additive.
  fit0 <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    make_test_data(nyrs = 15, nages = 5, seed = 42), file = NULL,
    estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))
  fc0  <- fit0$data_list$fleet_control
  fsh0 <- fc0$Fleet_code[fc0$Fleet_type == "Fishery"][1]
  testthat::expect_true(is.na(fit0$map$mapList$index_log_q[fsh0]))
})

testthat::test_that("a fishery's index is fitted and appears in the index diagnostics", {
  testthat::skip_if_not_installed("TMB")
  f <- .fishery_index_fixture()

  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    f$data, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))

  # Predicted for the fishery's index rows, not left at zero.
  idx <- fit$data_list$index_data
  sel <- idx$Fleet_code == f$fsh & idx$Year > 0
  hat <- as.numeric(fit$quantities$index_hat)[sel]
  testthat::expect_true(all(is.finite(hat)))
  testthat::expect_gt(min(hat), 0)

  # ...and it reaches the plotting frame, which used to drop it.
  df <- Rceattle:::.fleet_fit_df(list(fit), kind = "index")
  fsh_name <- fit$data_list$fleet_control$Fleet_name[
    fit$data_list$fleet_control$Fleet_code == f$fsh]
  testthat::expect_true(any(df$Fleet == fsh_name))
  # The survey is still there -- the selection widened, it did not move.
  srv_name <- fit$data_list$fleet_control$Fleet_name[
    fit$data_list$fleet_control$Fleet_code == f$srv]
  testthat::expect_true(any(df$Fleet == srv_name))

  # Catch plotting still selects fisheries only, so widening the index side did
  # not leak a survey into the catch panel.
  dc <- Rceattle:::.fleet_fit_df(list(fit), kind = "catch")
  testthat::expect_false(any(dc$Fleet == srv_name))
})

testthat::test_that("a fishery's index is predicted from the fishery's own selectivity", {
  testthat::skip_if_not_installed("TMB")
  # This is the intended behaviour for CPUE, and it is structural: sel_at_age is
  # indexed by fleet, so one fleet's index and its catch cannot use different
  # selectivities. An index needing its OWN selectivity has to be a separate
  # fleet. Demonstrated by rebuilding the prediction rather than asserted.
  f   <- .fishery_index_fixture(distinct_sel = TRUE)
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    f$data, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))

  q   <- fit$quantities
  td  <- fit$obj$env$data
  row <- which(td$index_ctl[, 1] == f$fsh & td$index_ctl[, 2] == 1 &
                 td$index_ctl[, 3] > 0)[1]
  testthat::skip_if(is.na(row), "no fitted fishery index row")

  yr  <- td$index_ctl[row, 3] - td$styr + 1L
  mo  <- td$index_n[row, 1]
  sp  <- 1L
  nag <- td$nages[sp]

  # q * sum over ages of N * exp(-(mo/12) Z) * sel * weight, with sel taken from
  # the FISHERY's own block.
  contrib <- function(flt) {
    s <- 0
    for (a in seq_len(nag)) {
      s <- s + q$N_at_age[sp, 1, a, yr] *
        exp(-(mo / 12) * q$Z_at_age[sp, 1, a, yr]) *
        q$sel_at_age[flt, 1, a, yr] *
        q$weight_hat[td$nspp * 2 + f$fsh, 1, a, yr]
    }
    s
  }
  qhat <- q$index_q[f$fsh, yr]
  testthat::expect_equal(as.numeric(q$index_hat)[row],
                         qhat * contrib(f$fsh), tolerance = 1e-6)

  # Discriminate against the SURVEY's block. The fixture gives the two fleets
  # different selectivity forms, so this distinguishes which one was used --
  # without that, both converge to the same shape and the check is vacuous.
  testthat::expect_false(isTRUE(all.equal(
    as.numeric(q$sel_at_age[f$fsh, 1, , yr]),
    as.numeric(q$sel_at_age[f$srv, 1, , yr]))))
  testthat::expect_false(isTRUE(all.equal(as.numeric(q$index_hat)[row],
                                          qhat * contrib(f$srv),
                                          tolerance = 1e-6)))
})

testthat::test_that("a fishery with index data but no catchability settings is rejected", {
  testthat::skip_if_not_installed("TMB")
  # Fitting at an undefined q is what this replaces. The message names the fleet
  # and the column, since nothing else would point at a fishery.
  dat <- make_test_data(nyrs = 15, nages = 5, seed = 42)
  fc  <- dat$fleet_control
  fsh <- which(fc$Fleet_type == "Fishery")[1]
  srv <- which(fc$Fleet_type == "Survey")[1]
  idx <- dat$index_data
  add <- idx[idx$Fleet_code == fc$Fleet_code[srv], , drop = FALSE]
  add$Fleet_code <- fc$Fleet_code[fsh]
  add$Fleet_name <- fc$Fleet_name[fsh]
  dat$index_data <- rbind(idx, add)      # columns left NA on the fishery row

  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      dat, file = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE,
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                          phase = FALSE)))),
    "carry index_data but have no")
})

testthat::test_that("a SURVEY missing the same columns is rejected too", {
  testthat::skip_if_not_installed("TMB")
  # The check keys on the data, not the fleet type, so it has to fire for a
  # survey as well. This is the regression net for the change: if it ever became
  # fishery-only, a survey with a blank Catchability would go back to dying
  # inside build_map() with "missing value where TRUE/FALSE needed".
  for (col in c("Catchability", "Catchability_init", "Estimate_index_sd")) {
    dat <- make_test_data(nyrs = 15, nages = 5, seed = 42)
    srv <- which(dat$fleet_control$Fleet_type == "Survey")[1]
    dat$fleet_control[[col]][srv] <- NA
    testthat::expect_error(
      suppressMessages(suppressWarnings(Rceattle::fit_mod(
        dat, file = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE,
        fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                            phase = FALSE)))),
      col, info = col)   # the message must name the offending column
  }
})
