# A recruitment linkage and a DSEM target different things: the linkage moves
# `R0`, `alpha`, `beta` or `R_init` -- the stock-recruit level, or the level the
# initial age structure sits at -- and the DSEM overwrites `rec_dev`, the
# deviations around that prediction. The two compose, and `fit_mod()` used to
# refuse the pair on the process NAME rather than its target.
#
# This is not hypothetical: the SS3 bridge carries SS3's initial recruitment
# level as an `R_init` linkage, and the refusal forced the GOA cod circulation
# study onto a configuration that folds the level into `init_dev` instead,
# costing 49.7 nats on the forward pass and biasing R0 31% low.

testthat::test_that("a recruitment linkage is allowed alongside a DSEM", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d <- Rceattle::BS2017SS
  yrs <- d$styr:d$endyr
  d$env_data <- data.frame(Year = yrs, temp = as.numeric(scale(seq_along(yrs))))

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 0,
    random_rec = TRUE,
    recFun = Rceattle::build_srr(linkages = list(
      R0 = Rceattle::linkage_spec(~ temp, by = ~ species))),
    dsem = Rceattle::build_DSEM(family = "fixed"),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  tbl <- fit[["data_list"]][["linkage_table"]]
  # Both structures are present and neither was dropped.
  testthat::expect_true(any(tbl[["process"]] == "recruitment"))
  testthat::expect_true(any(tbl[["param"]] == "R0"))
  testthat::expect_true(is.finite(fit[["quantities"]][["jnll"]]))

  # The linkage reached the model rather than being accepted and then dropped:
  # it contributes one estimated coefficient per species. (Its offset tensor is
  # still zero here -- nothing is optimized at estimateMode = 3, so every
  # coefficient sits at its starting value of 0.)
  testthat::expect_equal(sum(names(fit[["obj"]][["par"]]) == "beta_linkage"),
                         d[["nspp"]])

  # Nothing is double counted. The DSEM's GMRF is the only density on the
  # recruitment deviations, so the deviation row of jnll_comp stays zero while
  # the DSEM row carries it.
  jc <- fit[["quantities"]][["jnll_comp"]]
  testthat::expect_equal(sum(jc["Recruitment deviates", ]), 0)
  testthat::expect_true(sum(jc["DSEM", ]) != 0)
})
