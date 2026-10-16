# A predicted diet needs estimated suitability AND a predation model.
# `ceattle.cpp` gates predict_stomach_content() on
# `msmMode > 2 || imax(suitMode) > 0`, but that function reads diet_prop_hat,
# which only calculate_msvpa_predation() writes, under msmMode > 0. So diet_hat
# keeps its setZero() value in BOTH directions: suitMode = 0 with any msmMode,
# and msmMode = 0 with any suitMode. residuals(source = "diet") handed those
# zeros back as fitted values, which made every observation its own residual --
# a 100% misfit reported for a diet the model never fitted.
#
# Measured on the shipped BS2017MS (2025 diet rows, all at Year == 0):
#   suitMode = 0,0,0 / msmMode = 1  ->  diet_hat 0 .. 0, Residual == Observed 2025/2025
#   suitMode = 4,4,4 / msmMode = 0  ->  diet_hat 0 .. 0, Residual == Observed 2025/2025
#   suitMode = 4,4,4 / msmMode = 1  ->  1050 of 2025 rows non-zero
#
# plot_diet_comp(), plot_diet_comp1() and plot_diet_comp2() all read
# residuals(source = "diet") rather than diet_hat directly (R/7-plot_comp.R,
# "the single source of truth"), so all three inherit this.

make_diet_fit <- function(suitMode = NULL, msmMode = NULL) {
  dd <- data.frame(Pred = c(1L, 1L, 2L), Pred_sex = 0L, Prey = c(1L, 2L, 1L),
                   Prey_sex = 0L, Pred_age = c(3L, 3L, 4L),
                   Prey_age = c(1L, 2L, 1L), Year = c(2000L, 2000L, 2001L),
                   Stomach_proportion_by_weight = c(0.6, 0.3, 0.5),
                   Sample_size = c(50, 50, 40))
  dl <- list(diet_data = dd, spnames = c("A", "B"))
  if (!is.null(suitMode)) dl$suitMode <- suitMode
  if (!is.null(msmMode)) dl$msmMode <- msmMode
  structure(list(
    data_list  = dl,
    # All zero, as the template leaves it when the gate is closed.
    quantities = list(diet_hat = cbind(NA_real_, c(0, 0, 0)))),
    class = "Rceattle")
}

test_that("residuals(source = 'diet') refuses a model that never fitted diet", {
  fit <- make_diet_fit(suitMode = c(0L, 0L), msmMode = 1L)
  expect_error(residuals(fit, source = "diet"),
               "does not predict stomach proportions")
  # The message has to name both switches and this fit's values, or the reader
  # cannot act on it.
  expect_error(residuals(fit, source = "diet"), "suitMode 2, 4 or 6")
  expect_error(residuals(fit, source = "diet"), "msmMode > 0")
  expect_error(residuals(fit, type = "pearson", source = "diet"),
               "does not predict stomach proportions")

  # BOTH switches are needed, and estimated suitability alone is not enough:
  # predict_stomach_content() reads diet_prop_hat, which only
  # calculate_msvpa_predation() writes, under msmMode > 0. Measured on
  # BS2017MS at msmMode = 0, suitMode = 4: diet_hat was 0 .. 0 on all 2025 rows
  # with Residual == Observed on every one, and a first version of this guard
  # mirrored only the C++ gate and let it through.
  expect_error(residuals(make_diet_fit(suitMode = 4L, msmMode = 0L),
                         source = "diet"),
               "does not predict stomach proportions")
})

test_that("a fitted diet is still returned, so the refusal discriminates", {
  # Positive control. Without this the test above would pass on a function that
  # refused every call.
  fit <- make_diet_fit(suitMode = c(4L, 4L), msmMode = 1L)
  fit$quantities$diet_hat <- cbind(NA_real_, c(0.55, 0.35, 0.45))
  r <- residuals(fit, source = "diet")
  expect_equal(nrow(r), 3L)
  expect_equal(r$Fitted, c(0.55, 0.35, 0.45))

  # A mixed vector is fitted: the gate is imax(suitMode) > 0, so one predator
  # with estimated suitability writes the whole array. The empirical predator's
  # own rows are then an MSVPA reconstruction rather than a fitted prediction --
  # the likelihood skips it (`if (suitMode(rsp) <= 0) continue;`) and
  # build_osa_data() does too, so read those rows with that in mind.
  mixed <- make_diet_fit(suitMode = c(2L, 0L), msmMode = 1L)
  mixed$quantities$diet_hat <- cbind(NA_real_, c(0.55, 0.35, 0.45))
  expect_equal(nrow(residuals(mixed, source = "diet")), 3L)

  # msmMode > 2 opens the gate on its own, whatever suitMode says. Kept for
  # lockstep with ceattle.cpp, NOT as a reachable case: data_check() refuses
  # msmMode 3..9 (the retired Kinzey & Punt forms), so no fit_mod() call lands
  # here.
  kp <- make_diet_fit(suitMode = c(0L, 0L), msmMode = 3L)
  kp$quantities$diet_hat <- cbind(NA_real_, c(0.55, 0.35, 0.45))
  expect_equal(nrow(residuals(kp, source = "diet")), 3L)
})

test_that("a data_list carrying no suitMode is not refused", {
  # Absence is not evidence. fit_mod() always writes suitMode, so a list without
  # it is hand-built or predates the field; refusing there would break a saved
  # fit rather than catch a misreading.
  fit <- make_diet_fit()
  fit$quantities$diet_hat <- cbind(NA_real_, c(0.55, 0.35, 0.45))
  expect_equal(nrow(residuals(fit, source = "diet")), 3L)
  expect_true(.rce_diet_is_fitted(list()))
  expect_true(.rce_diet_is_fitted(list(suitMode = NA_integer_)))
  expect_false(.rce_diet_is_fitted(list(suitMode = 0L, msmMode = 1L)))
})
