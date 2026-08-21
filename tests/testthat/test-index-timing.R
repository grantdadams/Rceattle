# Index timing follows the fleet's role: a survey's index is a snapshot at its
# observation month, a fishery's CPUE is the year-average Nbar = N(1-exp(-Z))/Z.
# Nbar is Baranov's mean-numbers term -- C_a = F_a * Nbar_a, and effort cancels
# the F, so C/E = q * sum_a sel_a * Nbar_a * w_a. The snapshot/Nbar ratio moves
# with Z, so the two are not interchangeable up to a constant q.
testthat::skip_on_cran()

# A fishery and a survey carrying the same index observations, so the timing
# rule is the only thing separating their predictions.
.timing_fixture <- function(month = 6, fishery_month = month) {
  dat <- make_test_data(nyrs = 15, nages = 5, seed = 42)
  fc  <- dat$fleet_control
  fsh <- which(fc$Fleet_type == "Fishery")[1]
  srv <- which(fc$Fleet_type == "Survey")[1]
  testthat::skip_if(is.na(fsh) || is.na(srv), "fixture lacks both fleet types")

  dat$fleet_control$Catchability[fsh]      <- "Estimated"
  dat$fleet_control$Catchability_init[fsh] <- 1
  dat$fleet_control$Estimate_index_sd[fsh] <- 0

  idx <- dat$index_data
  add <- idx[idx$Fleet_code == fc$Fleet_code[srv], , drop = FALSE]
  add$Fleet_code <- fc$Fleet_code[fsh]
  add$Fleet_name <- fc$Fleet_name[fsh]
  # Non-zero month on both fleets: at Month = 0 the survey's snapshot factor is
  # exp(0) = 1 and the two forms separate for the wrong reason.
  idx$Month <- month
  add$Month <- fishery_month
  dat$index_data <- rbind(idx, add)
  list(data = dat, fsh = fc$Fleet_code[fsh], srv = fc$Fleet_code[srv],
       month = month)
}

# Rebuild both candidate predictions for one fitted index row from the reported
# quantities, so the assertions are constructive rather than pinned numbers.
.timing_preds <- function(fit, flt) {
  q  <- fit$quantities
  td <- fit$obj$env$data
  row <- which(td$index_ctl[, 1] == flt & td$index_ctl[, 2] == 1 &
                 td$index_ctl[, 3] > 0)[1]
  testthat::skip_if(is.na(row), "no fitted index row for the fleet")

  yr  <- td$index_ctl[row, 3] - td$styr + 1L
  mo  <- td$index_n[row, 1]
  sp  <- 1L
  N   <- q$N_at_age[sp, 1, , yr]
  Z   <- q$Z_at_age[sp, 1, , yr]
  sel <- q$sel_at_age[flt, 1, , yr]
  wt  <- q$weight_hat[td$nspp * 2 + flt, 1, , yr]
  qq  <- q$index_q[flt, yr]

  list(model     = as.numeric(q$index_hat)[row],
       annual    = qq * sum(N * (1 - exp(-Z)) / Z * sel * wt),
       snapshot  = qq * sum(N * exp(-(mo / 12) * Z) * sel * wt),
       month     = mo)
}

.timing_fit <- function(dat) {
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                        phase = FALSE))))
}


testthat::test_that("a fishery's index is the year-average numbers, not a snapshot", {
  testthat::skip_if_not_installed("TMB")
  f   <- .timing_fixture()
  p   <- .timing_preds(.timing_fit(f$data), f$fsh)

  testthat::expect_equal(p$model, p$annual, tolerance = 1e-8)
  # The two forms must differ here, or the check above passes by coincidence --
  # they coincide as Z -> 0.
  testthat::expect_false(isTRUE(all.equal(p$annual, p$snapshot,
                                          tolerance = 1e-6)))
})


testthat::test_that("a survey's index is still a snapshot at its observation month", {
  testthat::skip_if_not_installed("TMB")
  # Keying on fleet type leaves the survey side untouched. If a survey ever
  # started returning the year-average, every survey-based model would move.
  f <- .timing_fixture()
  p <- .timing_preds(.timing_fit(f$data), f$srv)

  testthat::expect_equal(p$model, p$snapshot, tolerance = 1e-8)
  testthat::expect_false(isTRUE(all.equal(p$snapshot, p$annual,
                                          tolerance = 1e-6)))
})


testthat::test_that("a fishery's index ignores the observation month", {
  testthat::skip_if_not_installed("TMB")
  # The year-average has no instant to be taken at, so Month is inert on a
  # fishery's index rows -- worth pinning, since it is set on every other kind of
  # observation.
  #
  # Only the fishery's month moves between the two fits: moving the survey's
  # would change the survey's prediction, and so the fit, for another reason.
  f1 <- .timing_fixture(month = 6, fishery_month = 1)
  f2 <- .timing_fixture(month = 6, fishery_month = 11)
  a  <- .timing_preds(.timing_fit(f1$data), f1$fsh)
  b  <- .timing_preds(.timing_fit(f2$data), f2$fsh)

  testthat::expect_equal(a$month, 1)
  testthat::expect_equal(b$month, 11)
  # Bit-identical, not merely close: a month-independent prediction leaves the
  # likelihood unchanged, so these are the same fit.
  testthat::expect_equal(a$model, b$model, tolerance = 1e-12)
})


testthat::test_that("a fishery's index and its catch use the same mean numbers", {
  testthat::skip_if_not_installed("TMB")
  # Index and catch share one quantity under the year-average form: Baranov gives
  # C_a = F_a * Nbar_a, so the same Nbar has to reproduce BOTH of the model's own
  # reported predictions. Both sides are read off the report -- reconstructing
  # each from N and Z alone would assert the identity
  # (F/Z (1-exp(-Z)) N) / F == N (1-exp(-Z))/Z, which holds for any numbers and
  # would pass with section 8.1 deleted.
  f   <- .timing_fixture()
  fit <- .timing_fit(f$data)
  q   <- fit$quantities
  td  <- fit$obj$env$data

  irow <- which(td$index_ctl[, 1] == f$fsh & td$index_ctl[, 3] > 0)[1]
  yr_d <- td$index_ctl[irow, 3]
  crow <- which(td$catch_ctl[, 1] == f$fsh & td$catch_ctl[, 3] == yr_d)[1]
  testthat::skip_if(is.na(crow), "fishery has no catch row in the index year")

  yr <- yr_d - td$styr + 1L
  sp <- 1L
  Z  <- q$Z_at_age[sp, 1, , yr]
  N  <- q$N_at_age[sp, 1, , yr]
  Fa <- q$F_flt_age[f$fsh, 1, , yr]
  wt <- q$weight_hat[td$nspp * 2 + f$fsh, 1, , yr]
  sel <- q$sel_at_age[f$fsh, 1, , yr]
  testthat::skip_if(all(Fa <= 0), "no fishing mortality in the fixture year")

  nbar <- N * (1 - exp(-Z)) / Z

  # Section 8.1 builds the index on Nbar...
  testthat::expect_equal(as.numeric(q$index_hat)[irow],
                         as.numeric(q$index_q[f$fsh, yr]) * sum(sel * nbar * wt),
                         tolerance = 1e-8)
  # ...and section 9.1 builds the catch on the same Nbar, since C_a = F_a Nbar_a.
  testthat::expect_equal(as.numeric(q$catch_hat)[crow],
                         sum(Fa * nbar * wt),
                         tolerance = 1e-8)
})


testthat::test_that("the timing rule holds in multispecies mode and in projection years", {
  testthat::skip_if_not_installed("TMB")
  # Two corners the single-species hindcast fixture does not reach. In
  # multispecies mode Z carries predation M2, so Nbar is built on a Z the
  # single-species tests never exercise; and a projection row indexes past
  # nyrs_hind, where section 8.3 freezes q at the last hindcast column.
  dat <- make_msm_test_data(years = 1:20, ages = 1:8)$data_list
  fc  <- dat$fleet_control
  # This fixture carries Fleet_type as the integer code, not the string that
  # switch_check() would normalise it to, so match on both spellings.
  fsh <- which(fc$Fleet_type %in% c(1, "1", "Fishery"))[1]
  srv <- which(fc$Fleet_type %in% c(2, "2", "Survey"))[1]
  testthat::skip_if(is.na(fsh) || is.na(srv), "fixture lacks both fleet types")

  dat$fleet_control$Catchability[fsh]      <- "Estimated"
  dat$fleet_control$Catchability_init[fsh] <- 1
  dat$fleet_control$Estimate_index_sd[fsh] <- 0

  idx <- dat$index_data
  add <- idx[idx$Fleet_code == fc$Fleet_code[srv], , drop = FALSE]
  add$Fleet_code <- fc$Fleet_code[fsh]
  add$Fleet_name <- fc$Fleet_name[fsh]
  add$Month      <- 6
  # A projection row for the fishery: Year < 0 marks a row the model predicts
  # but does not score.
  proj <- add[1, , drop = FALSE]
  proj$Year <- -(dat$endyr + 1L)
  dat$index_data <- rbind(idx, add, proj)

  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 3, msmMode = 1, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                        phase = FALSE))))
  q  <- fit$quantities
  td <- fit$obj$env$data
  sp_of <- td$index_ctl[, 2]

  rebuild <- function(row, form) {
    flt <- td$index_ctl[row, 1]
    yd  <- td$index_ctl[row, 3]
    yr  <- (if (yd > 0) yd else -yd) - td$styr + 1L
    sp  <- sp_of[row]
    mo  <- td$index_n[row, 1]
    tot <- 0
    for (sx in seq_len(td$nsex[sp])) {
      N <- q$N_at_age[sp, sx, , yr]
      Z <- q$Z_at_age[sp, sx, , yr]
      n <- if (form == "annual") N * (1 - exp(-Z)) / Z else N * exp(-(mo / 12) * Z)
      tot <- tot + sum(n * q$sel_at_age[flt, sx, , yr] *
                         q$weight_hat[td$nspp * 2 + flt, sx, , yr],
                       na.rm = TRUE)
    }
    # Section 8.3 freezes q at the last hindcast column for a projection row.
    yq <- min(yr, ncol(q$index_q))
    as.numeric(q$index_q[flt, yq]) * tot
  }

  fsh_rows <- which(td$index_ctl[, 1] == fc$Fleet_code[fsh])
  srv_rows <- which(td$index_ctl[, 1] == fc$Fleet_code[srv])
  testthat::skip_if(!length(fsh_rows) || !length(srv_rows), "no index rows")

  for (r in fsh_rows) {
    testthat::expect_equal(as.numeric(q$index_hat)[r], rebuild(r, "annual"),
                           tolerance = 1e-8,
                           info = paste("fishery row", r, "year", td$index_ctl[r, 3]))
  }
  for (r in srv_rows) {
    testthat::expect_equal(as.numeric(q$index_hat)[r], rebuild(r, "snapshot"),
                           tolerance = 1e-8,
                           info = paste("survey row", r, "year", td$index_ctl[r, 3]))
  }
  # A projection row really was included, so the loop above covered one.
  testthat::expect_true(any(td$index_ctl[fsh_rows, 3] < 0))
  # And predation M2 really is in play, so Z differs from the single-species Z.
  testthat::expect_gt(max(q$M2_at_age, na.rm = TRUE), 0)
})
