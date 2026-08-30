# Spawning-biomass-per-recruit on a two-sex species.
#
# SPR is spawning output per TOTAL recruit, so the female fraction belongs in it
# for every species -- once. A one-sex species' schedule is sex-combined, so the
# age-varying ratio applies at every age. A two-sex species' sex-0 schedule is
# already female, so only the recruitment split `sex_ratio(sp, 0)` applies.
#
# Section 6.2 used to apply the AGE-VARYING ratio to both, which re-applied on a
# two-sex species a split already present in its schedule, and returned NaN
# wherever that table stopped short of the species' own nages -- how a Tier 3
# SPR proxy came back missing for GOA arrowtooth (nsex = 2, nages = 21) in the
# 2026 three-species assessment.
#
# GOA2018SS is the bundled dataset with nsex = c(1, 2, 1), so it exercises both
# branches in one fit. test-dynamics-brps.R and test-dynamics-spr.R pin the
# per-total-recruit convention itself, on fixtures where sex_ratio is 0.5 at
# every age so the two branches must agree.

spr_fixture <- function() {
  data(GOA2018SS, package = "Rceattle", envir = environment())
  suppressMessages(fit_mod(
    data_list = GOA2018SS, file = NULL, inits = NULL,
    estimateMode = 3, random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE)))
}


testthat::test_that("SPR is finite for a two-sex species", {
  testthat::skip_on_cran()
  fit <- spr_fixture()
  nsex <- fit$obj$env$data$nsex
  testthat::expect_true(any(nsex == 2))     # the fixture must actually cover it

  for (q in c("SPR0", "SPRtarget", "SPRlimit", "SPRFinit")) {
    testthat::expect_true(all(is.finite(fit$quantities[[q]])),
                          info = paste(q, "=", paste(fit$quantities[[q]],
                                                     collapse = ", ")))
  }
})


testthat::test_that("SPR includes the female fraction exactly once", {
  # SPR is spawning output per TOTAL recruit, so the female fraction is in it
  # for every species -- but it must enter once. A one-sex species' schedule is
  # sex-combined, so the age-varying ratio applies at every age (that is what
  # `mature_females` folds in). A two-sex species' sex-0 schedule is already
  # female, so only the recruitment split sex_ratio[1] applies; using the
  # age-varying ratio there re-applies a split already in the schedule, and
  # dropping it entirely reports per FEMALE recruit, which is a different
  # quantity from the one the other species report.
  testthat::skip_on_cran()
  fit <- spr_fixture()
  d <- fit$data_list
  dat <- fit$obj$env$data
  q <- fit$quantities
  ty <- length(d$styr:d$endyr)

  for (sp in seq_len(d$nspp)) {
    na <- dat$nages[sp]
    Z <- q$M_at_age[sp, 1, seq_len(na), ty]                    # F = 0
    wt <- q$weight_hat[2 * (sp - 1) + 2, 1, seq_len(na), ty]   # spawning weight
    mat <- as.numeric(dat$maturity[sp, seq_len(na)])
    sr <- as.numeric(dat$sex_ratio[sp, seq_len(na)])
    matf <- if (dat$nsex[sp] == 1) mat * sr else mat * sr[1]

    n <- numeric(na)
    n[1] <- 1
    if (na > 2) for (a in 2:(na - 1)) n[a] <- n[a - 1] * exp(-Z[a - 1])
    n[na] <- n[na - 1] * exp(-Z[na - 1]) / (1 - exp(-Z[na]))

    expected <- sum(n * wt * matf * exp(-Z * dat$spawn_month[sp] / 12))
    testthat::expect_equal(as.numeric(q$SPR0[sp]), expected, tolerance = 1e-8,
                           info = paste("species", sp, "nsex", dat$nsex[sp]))
  }
})


testthat::test_that("a two-sex species needs sex_ratio only at age 1", {
  # After the fix the model reads `sex_ratio` at every age only for a one-sex
  # species; on a two-sex one every remaining use is `sex_ratio(sp, 0)`, the
  # recruitment split. data_check() must not demand a full schedule it never
  # reads, or it rejects a workbook that fits perfectly well.
  data(GOA2018SS, package = "Rceattle", envir = environment())
  d <- GOA2018SS
  two_sex <- which(d$nsex == 2)[1]
  testthat::expect_false(is.na(two_sex))

  agec <- grep("^Age", colnames(d$sex_ratio))

  # The bundled workbook raises other, unrelated data_check errors, so ask
  # specifically whether sex_ratio was one of them.
  flags_sex_ratio <- function(dl) {
    msg <- tryCatch({ suppressMessages(data_check(dl)); "" },
                    error = function(e) conditionMessage(e))
    grepl("sex_ratio is missing values", msg, fixed = TRUE)
  }

  gap_after_age1 <- d
  gap_after_age1$sex_ratio[two_sex, agec[-1]] <- NA
  testthat::expect_false(flags_sex_ratio(gap_after_age1))

  # ...but a one-sex species still needs the whole schedule, and age 1 is
  # required even for a two-sex species.
  one_sex <- which(d$nsex == 1)[1]
  gap_one_sex <- d
  gap_one_sex$sex_ratio[one_sex, agec[2]] <- NA
  testthat::expect_true(flags_sex_ratio(gap_one_sex))

  gap_at_age1 <- d
  gap_at_age1$sex_ratio[two_sex, agec[1]] <- NA
  testthat::expect_true(flags_sex_ratio(gap_at_age1))
})


testthat::test_that("nsex must be 1 or 2, and an unusable one is reported not thrown", {
  # nsex is read straight into loop bounds and array dimensions, so a value
  # outside {1, 2} is a silent out-of-range read. An NA used to reach the
  # sex-consistency comparisons as `if(NA)` and abort data_check() with R's
  # "missing value where TRUE/FALSE needed", naming no table at all.
  data(GOA2018SS, package = "Rceattle", envir = environment())
  g <- GOA2018SS

  msg_of <- function(dl) {
    tryCatch({ suppressMessages(data_check(dl)); "" },
             error = function(e) conditionMessage(e))
  }

  bad_na <- g; bad_na$nsex[2] <- NA
  m <- msg_of(bad_na)
  testthat::expect_match(m, "'nsex' must be 1 .* or 2")
  # ...and specifically NOT R's uninformative abort.
  testthat::expect_false(grepl("missing value where TRUE/FALSE needed", m,
                               fixed = TRUE))

  bad_zero <- g; bad_zero$nsex[1] <- 0
  testthat::expect_match(msg_of(bad_zero), "'nsex' must be 1")

  bad_three <- g; bad_three$nsex[3] <- 3
  testthat::expect_match(msg_of(bad_three), "'nsex' must be 1")

  # Length is read by position, so a short vector is its own error.
  short <- g; short$nsex <- short$nsex[1:2]
  testthat::expect_match(msg_of(short), "'nsex' has 2 value\\(s\\) for 3 species")

  # A valid nsex raises no nsex complaint (the workbook has other, unrelated
  # data_check errors, so ask about this one specifically).
  testthat::expect_false(grepl("'nsex'", msg_of(g), fixed = TRUE))
})
