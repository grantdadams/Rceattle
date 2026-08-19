# Everything downstream of the index likelihood has to know which SCALE the
# fleet is fitted on.
#
# `Lognormal` (0) is a log-scale family whose `Log_sd` is a CV / log-sd.
# `MVN` (1), `MVNORM` (2), `Normal` (3) and `TruncatedNormal` (4) are
# natural-scale families whose sd is ABSOLUTE, in the units of the index.
# Applying a log-scale formula to the second group does not error, it silently
# returns nonsense: `sigma^2 / 2` becomes a number the size of the index squared.
#
# These pin the two places that got it wrong -- the Pearson residual and the
# observation interval -- with fixtures where the wrong answer is not subtle.
testthat::skip_on_cran()

.ns_fixture <- function(dist, sd, nyrs = 12) {
  dat <- make_test_data(nyrs = nyrs, nages = 5, seed = 123)
  srv <- dat$fleet_control$Fleet_name == "Survey"
  dat$fleet_control$Index_distribution[srv] <- dist
  if (dist %in% c("MVN", "MVNORM")) {
    dat$fleet_control$Catchability[srv] <- "AnalyticalArith"
    dat$index_cov <- list(Survey = diag(rep(sd^2, nyrs)))
  }
  dat$index_data$Log_sd[dat$index_data$Fleet_name == "Survey"] <- sd
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))
}

testthat::test_that("Pearson index residuals use the fleet's own scale", {
  testthat::skip_if_not_installed("TMB")

  fit <- .ns_fixture("Normal", sd = 150)
  r   <- stats::residuals(fit, type = "pearson", source = "index")
  srv <- r$Fleet_name == "Survey" & is.finite(r$Residual)
  res <- r$Residual[srv]

  # Natural scale: (obs - hat) / sd. With sd = 150 these are O(1).
  idx <- fit$data_list$index_data
  sel <- idx$Fleet_name == "Survey"
  want <- (idx$Observation[sel] - as.numeric(fit$quantities$index_hat)[sel]) / 150
  testthat::expect_equal(res, want[is.finite(want)][seq_along(res)],
                         tolerance = 1e-6)

  # The lognormal formula would return the same large constant for every row --
  # about +75 here, since -sd^2/2 dominates. Assert we are nowhere near it.
  testthat::expect_lt(max(abs(res)), 10)
  testthat::expect_gt(stats::sd(res), 0)

  # A lognormal fleet is untouched by the fix.
  ln  <- .ns_fixture("Lognormal", sd = 0.2)
  rl  <- stats::residuals(ln, type = "pearson", source = "index")
  sl  <- rl$Fleet_name == "Survey" & is.finite(rl$Residual)
  idl <- ln$data_list$index_data
  sel <- idl$Fleet_name == "Survey"
  sig <- as.numeric(ln$quantities$log_index_sd)[sel]
  wl  <- (log(idl$Observation[sel]) -
            (log(as.numeric(ln$quantities$index_hat)[sel]) - sig^2 / 2)) / sig
  testthat::expect_equal(rl$Residual[sl], wl[is.finite(wl)][seq_len(sum(sl))],
                         tolerance = 1e-6)
})

testthat::test_that("index observation intervals use the fleet's own scale", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("ggplot2")

  for (dist in c("Normal", "TruncatedNormal", "MVN")) {
    fit <- .ns_fixture(dist, sd = 150)
    df  <- Rceattle:::.fleet_fit_df(list(fit), kind = "index")
    d   <- df[df$Fleet == "Survey" & is.finite(df$Upper95), ]
    testthat::expect_gt(nrow(d), 0)

    # qlnorm(0.975, log(obs), 150) overflows to ~1e130; the natural-scale band is
    # obs +/- 1.96*150, so a few hundred at most on this fixture.
    testthat::expect_lt(max(d$Upper95), 1e5, label = dist)
    # An index cannot be negative, so the lower bound is clamped at zero.
    testthat::expect_gte(min(d$Lower95), 0)
  }

  # Lognormal keeps the multiplicative band.
  ln <- .ns_fixture("Lognormal", sd = 0.2)
  dl <- Rceattle:::.fleet_fit_df(list(ln), kind = "index")
  dl <- dl[dl$Fleet == "Survey" & is.finite(dl$Upper95), ]
  testthat::expect_gt(nrow(dl), 0)
  testthat::expect_equal(dl$Upper95,
                         stats::qlnorm(0.975, log(dl$Observation), 0.2),
                         tolerance = 1e-6)
})


testthat::test_that("every natural-scale family is recognized as one", {
  testthat::skip_if_not_installed("TMB")

  # The scale question is answered once, in .index_rows_natural_scale(). A family
  # added to index_distribution_map but not to that vector reverts silently to
  # the log-scale treatment -- the exact failure the rest of this file pins. This
  # enumerates the map rather than a hand-written list, so a new family fails
  # here rather than in a user's residual plot.
  d <- make_test_data(nyrs = 12, nages = 5, seed = 123)
  srv <- d$fleet_control$Fleet_name == "Survey"

  for (nm in names(index_distribution_map)) {
    dat <- d
    dat$fleet_control$Index_distribution[srv] <- nm
    nat <- Rceattle:::.index_rows_natural_scale(dat)
    rows <- dat$index_data$Fleet_name == "Survey"
    want <- nm != "Lognormal"
    testthat::expect_equal(unique(nat[rows]), want, label = nm)
    # ...and the code path agrees with the name path, since fleet_control may
    # hold either spelling.
    dat$fleet_control$Index_distribution[srv] <- index_distribution_map[[nm]]
    testthat::expect_equal(unique(Rceattle:::.index_rows_natural_scale(dat)[rows]),
                           want, label = paste(nm, "as code"))
  }
})

testthat::test_that("TruncatedNormal gets natural-scale Pearson residuals", {
  testthat::skip_if_not_installed("TMB")

  fit <- .ns_fixture("TruncatedNormal", sd = 150)
  r   <- stats::residuals(fit, type = "pearson", source = "index")
  srv <- r$Fleet_name == "Survey" & is.finite(r$Residual)
  res <- r$Residual[srv]

  idx <- fit$data_list$index_data
  sel <- idx$Fleet_name == "Survey"
  want <- (idx$Observation[sel] - as.numeric(fit$quantities$index_hat)[sel]) / 150
  testthat::expect_equal(res, want[is.finite(want)][seq_along(res)],
                         tolerance = 1e-6)
  # Nowhere near the ~+75 the lognormal formula would return.
  testthat::expect_lt(max(abs(res)), 10)
})
