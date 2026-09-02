# fit_control(selectivity_se = TRUE) adds a delta-method standard error for the
# selectivity curve, asked for in issue #107 so a fitted curve can be drawn with
# an interval.
#
# Reported on the LOG scale rather than the logit that issue proposes. The
# non-parametric forms renormalize to MEAN selectivity 1, so sel_at_age is not
# bounded by 1: on Atka2022, 58% of entries exceed it and the largest is 3.06.
# logit() is undefined for those, and one NaN on the tape NaNs the whole
# sdreport. Atka2022 is the fixture precisely because it is one of those models.

testthat::test_that("selectivity_se adds a log-scale SE without moving the fit", {
  testthat::skip_on_cran()

  run <- function(se) suppressMessages(suppressWarnings(fit_mod(
    data_list    = Atka2022,
    msmMode      = 0,
    estimateMode = "Hindcast",
    fit_control  = fit_control(selectivity_se = se, verbose = 0))))

  off <- run(FALSE)
  on  <- run(TRUE)

  # An ADREPORT is a report, not a likelihood term.
  testthat::expect_identical(on$opt$objective, off$opt$objective)

  testthat::expect_false("log_sel_at_age" %in% names(off$sdrep$value))
  testthat::expect_true("log_sel_at_age" %in% names(on$sdrep$value))

  # One row per lead fleet x its OWN species' sexes and ages x hindcast year.
  # Iterating the padded array instead would report structural zeros, and
  # log(0) = -Inf NaNs every other quantity in the same sdreport.
  fc  <- on$data_list$fleet_control
  spp <- fc$Species
  lead <- which(fc$Fleet_type != "Off" & fc$Selectivity_index == fc$Fleet_code &
                  fc$Selectivity != "Fixed")
  nyrs_hind <- on$data_list$endyr - on$data_list$styr + 1
  expected <- sum(on$data_list$nsex[spp[lead]] * on$data_list$nages[spp[lead]]) *
    nyrs_hind

  k <- names(on$sdrep$value) == "log_sel_at_age"
  testthat::expect_equal(sum(k), expected)

  # Columns 3 and 4 are the ABSOLUTE age and the CALENDAR year, not 0-based or
  # 1-based positions. nages counts bins, so age and age index differ wherever
  # minage > 1 -- which no bundled dataset has, so only this assertion catches
  # a regression to raw indices.
  idx <- on$quantities$log_sel_at_age_index
  testthat::expect_equal(nrow(idx), expected)
  testthat::expect_equal(ncol(idx), 4L)
  testthat::expect_setequal(unique(idx[, 4]), on$data_list$styr:on$data_list$endyr)
  minage_sp <- on$data_list$minage[spp[lead]]
  testthat::expect_gte(min(idx[, 3]), min(minage_sp))
  testthat::expect_true(all(idx[, 1] %in% lead))
  testthat::expect_true(all(idx[, 2] >= 1))

  # The index has to address the right cell, or an interval gets drawn against
  # the wrong age. Convert back to array positions to look it up.
  sel  <- on$quantities$sel_at_age
  back <- exp(on$sdrep$value[k])
  age_i <- idx[, 3] - on$data_list$minage[spp[idx[, 1]]] + 1L
  yr_i  <- idx[, 4] - on$data_list$styr + 1L
  at   <- sel[cbind(idx[, 1], idx[, 2], age_i, yr_i)]
  testthat::expect_equal(unname(back), unname(at), tolerance = 1e-12)
  testthat::expect_true(all(at > 0))   # a zero would have put -Inf on the tape

  # Finite everywhere, and the rest of the sdreport is untouched by it.
  sd_sel <- on$sdrep$sd[k]
  testthat::expect_true(all(is.finite(sd_sel)))
  testthat::expect_false(any(is.nan(on$sdrep$sd)))

  # A bin the model pins carries no uncertainty. Atka2022's survey normalizes at
  # Sel_norm_bin = 4, so that age is 1 by construction with an SE of 0.
  testthat::expect_true(any(sd_sel == 0))
})


testthat::test_that("selectivity_se defaults to off and survives a data list without it", {
  # The TMB template declares adreport_sel as DATA_INTEGER, so MakeADFun() fails
  # outright if the field is missing. A data list built before the option
  # existed has no such field, and .osa_build_data() -- the funnel every list
  # passes through -- has to supply it.
  testthat::expect_false(fit_control()$selectivity_se)

  d <- Atka2022
  d$adreport_sel <- NULL
  built <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, msmMode = 0, estimateMode = "DebugBuild",
    fit_control = fit_control(verbose = 0))))
  testthat::expect_identical(built$data_list$adreport_sel, 0L)
})


# BS2017SS is the regression case for the Fixed-selectivity skip: its
# EIT_Pollock survey is Selectivity = "Fixed", which leaves 0 in every year with
# no emp_sel_obs row -- 228 cells. Reporting log() of those puts -Inf on the AD
# tape, and sdreport then returns NaN for EVERY quantity, biomass and ssb
# included, not just for selectivity.
testthat::test_that("a Fixed-selectivity fleet cannot NaN the sdreport", {
  testthat::skip_on_cran()

  m <- suppressMessages(suppressWarnings(fit_mod(
    data_list = BS2017SS, msmMode = 0, estimateMode = "Hindcast",
    fit_control = fit_control(selectivity_se = TRUE, verbose = 0))))

  testthat::expect_false(any(is.nan(m$sdrep$sd)))
  testthat::expect_false(any(is.nan(m$sdrep$value)))

  k <- names(m$sdrep$value) == "log_sel_at_age"
  testthat::expect_gt(sum(k), 0)
  testthat::expect_true(all(is.finite(m$sdrep$value[k])))

  # The Fixed fleet contributes no rows.
  fc <- m$data_list$fleet_control
  fixed <- fc$Fleet_code[fc$Selectivity == "Fixed"]
  idx <- m$quantities$log_sel_at_age_index
  testthat::expect_length(intersect(idx[, 1], fixed), 0)
})
