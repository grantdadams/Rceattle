# `Sel_norm_scope` decides whether selectivity normalization pools its reference
# across sexes. It is orthogonal to `Sel_norm_bin`, which decides WHERE that
# reference is taken -- the two combine rather than multiplying into more
# columns. These guard the switch plumbing, the default, and the two behaviours.

test_that("Sel_norm_scope defaults to AcrossSexes and round-trips as a switch", {
  scope_map <- get("sel_norm_scope_map", envir = asNamespace("Rceattle"))
  expect_identical(scope_map[["WithinSex"]], 0)
  expect_identical(scope_map[["AcrossSexes"]], 1)

  d <- Rceattle::GOAatf2023
  d$fleet_control$Sel_norm_scope <- NULL          # absent -> default applies
  out <- suppressMessages(suppressWarnings(Rceattle::switch_check(d)))
  # switch_check() leaves readable strings; convert_switches() makes them codes.
  expect_true(all(out$fleet_control$Sel_norm_scope == "AcrossSexes"))
})

test_that("a two-sex fleet normalizing at a named bin is warned about the flip", {
  # The one configuration whose behaviour the AcrossSexes default changes: a
  # named bin used to imply a per-sex reference. Max-normalized fleets already
  # pooled, and a one-sex fleet cannot tell the difference.
  d <- Rceattle::GOAatf2023
  d$fleet_control$Sel_norm_scope <- NULL
  d$fleet_control$Sel_norm_bin <- 6               # named bin + nsex 2 -> exposed
  expect_message(Rceattle::switch_check(d), "Sel_norm_scope")

  d2 <- Rceattle::GOAatf2023
  d2$fleet_control$Sel_norm_scope <- NULL
  d2$fleet_control$Sel_norm_bin <- NA             # nothing normalized -> silent
  expect_no_message(suppressWarnings(Rceattle::switch_check(d2)),
                    message = "Sel_norm_scope")
})

test_that("a blank Sel_norm_scope cell is announced, not silently flipped", {
  # Supplying the column and leaving a row empty is the same behaviour flip as
  # omitting the column, so it must not be a way to skip the notice.
  d <- Rceattle::GOAatf2023
  d$fleet_control$Sel_norm_bin <- 6
  d$fleet_control$Sel_norm_scope <- "WithinSex"
  d$fleet_control$Sel_norm_scope[1] <- NA
  expect_message(Rceattle::switch_check(d), "Sel_norm_scope")

  out <- suppressMessages(suppressWarnings(Rceattle::switch_check(d)))
  expect_identical(out$fleet_control$Sel_norm_scope[1], "AcrossSexes")
  expect_identical(out$fleet_control$Sel_norm_scope[2], "WithinSex")
})

test_that("the flip notice does not fire where normalization never runs", {
  # selectivity.hpp skips the normalization block for Hake (5/12) and LogisticPM
  # (11), which reuse Sel_norm_bin1/2 as a penalty age-range, and for Off fleets.
  # An AMAK-style model setting that range must not be told its results moved.
  for (sel in c("Hake", "LogisticPM")) {
    d <- Rceattle::GOAatf2023
    d$fleet_control$Sel_norm_scope <- NULL
    d$fleet_control$Sel_norm_bin <- 6             # a penalty range, not a reference
    d$fleet_control$Selectivity <- sel
    expect_no_message(suppressWarnings(Rceattle::switch_check(d)),
                      message = "Sel_norm_scope")
  }

  d_off <- Rceattle::GOAatf2023
  d_off$fleet_control$Sel_norm_scope <- NULL
  d_off$fleet_control$Sel_norm_bin <- 6
  d_off$fleet_control$Fleet_type <- "Off"
  expect_no_message(suppressWarnings(Rceattle::switch_check(d_off)),
                    message = "Sel_norm_scope")
})

# --- behaviour ---------------------------------------------------------------
testthat::skip_on_cran()

# Build a two-sex model with deliberately dimorphic selectivity on fleet 1:
# females steep and early, males shallow and late, so the two curves sit at
# genuinely different levels before normalization.
.norm_scope_fit <- function(bin, scope = NULL) {
  d <- Rceattle::GOAatf2023
  d$fleet_control$Sel_norm_bin <- bin
  if (!is.null(scope)) d$fleet_control$Sel_norm_scope <- scope
  ctl <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  base <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 3, msmMode = 0,
    random_rec = FALSE, fit_control = ctl)))
  p <- base$estimated_params
  p$sel_inf[1, 1, 1] <- 4;  p$log_sel_slp[1, 1, 1] <- log(1.2)
  p$sel_inf[1, 1, 2] <- 15; p$log_sel_slp[1, 1, 2] <- log(0.25)
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, inits = p, file = NULL, estimateMode = 3, msmMode = 0,
    random_rec = FALSE, fit_control = ctl)))
  s <- fit$quantities$sel_at_age
  c(F = max(s[1, 1, , 1]), M = max(s[1, 2, , 1]))
}

test_that("WithinSex scales each sex to its own reference; AcrossSexes pools", {
  testthat::skip_if_not_installed("TMB")

  # Anchored at a named bin: WithinSex divides the males by their own (smaller)
  # value there, lifting their curve; AcrossSexes divides both by the pooled
  # reference, so the sexes keep their relative levels.
  within_bin <- .norm_scope_fit(6, "WithinSex")
  across_bin <- .norm_scope_fit(6, "AcrossSexes")
  expect_gt(within_bin[["M"]], across_bin[["M"]])
  expect_equal(within_bin[["F"]], across_bin[["F"]])   # the leading sex is unmoved

  # Anchored at each curve's max: WithinSex puts BOTH sexes at exactly 1;
  # AcrossSexes leaves the less-selected sex below it.
  within_max <- .norm_scope_fit(-1, "WithinSex")
  across_max <- .norm_scope_fit(-1, "AcrossSexes")
  expect_equal(unname(within_max[["M"]]), 1, tolerance = 1e-8)
  expect_equal(unname(within_max[["F"]]), 1, tolerance = 1e-8)
  expect_lt(across_max[["M"]], within_max[["M"]])
})

test_that("the scope is inert for a one-sex species", {
  testthat::skip_if_not_installed("TMB")
  # With one sex there is nothing to pool, so the two scopes must agree exactly.
  # This is what keeps every single-sex model bit-identical across the change.
  ctl <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  run <- function(scope) {
    d <- Rceattle::BS2017SS
    d$fleet_control$Sel_norm_bin <- 4
    d$fleet_control$Sel_norm_scope <- scope
    fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d, inits = NULL, file = NULL, estimateMode = 3, msmMode = 0,
      random_rec = FALSE, fit_control = ctl)))
    list(sel = fit$quantities$sel_at_age, jnll = fit$quantities$jnll)
  }
  a <- run("WithinSex"); b <- run("AcrossSexes")
  expect_identical(a$sel, b$sel)
  expect_identical(a$jnll, b$jnll)
})


test_that("a blank or invalid Sel_norm_scope is filled / rejected, not passed on", {
  # A blank cell means "not specified". Left alone it reaches the TMB integer
  # vector as NA, where it reads as WithinSex -- the opposite of the default.
  d <- Rceattle::GOAatf2023
  d$fleet_control$Sel_norm_scope <- NA
  out <- suppressMessages(suppressWarnings(Rceattle::switch_check(d)))
  expect_true(all(out$fleet_control$Sel_norm_scope == "AcrossSexes"))

  # ...and a typo is reported rather than silently becoming NA. validate_switches()
  # collects errors for data_check() to raise, so check the message it returns.
  bad <- suppressMessages(suppressWarnings(Rceattle::switch_check(Rceattle::GOAatf2023)))
  bad$fleet_control$Sel_norm_scope <- "withinsex"
  validate <- get("validate_switches", envir = asNamespace("Rceattle"))
  expect_match(paste(validate(bad), collapse = " "), "Sel_norm_scope")
  expect_length(validate(suppressMessages(suppressWarnings(
    Rceattle::switch_check(Rceattle::GOAatf2023)))), 0)
})

# Return the whole fleet-1 selectivity-at-age matrix (sex x age, first year),
# rather than just each sex's peak, so the range branch can be checked directly.
.norm_scope_sel <- function(bin, upper = NA, scope = NULL) {
  d <- Rceattle::GOAatf2023
  d$fleet_control$Sel_norm_bin <- bin
  d$fleet_control$Sel_norm_bin_upper <- upper
  if (!is.null(scope)) d$fleet_control$Sel_norm_scope <- scope
  ctl <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  base <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 3, msmMode = 0,
    random_rec = FALSE, fit_control = ctl)))
  p <- base$estimated_params
  p$sel_inf[1, 1, 1] <- 4;  p$log_sel_slp[1, 1, 1] <- log(1.2)
  p$sel_inf[1, 1, 2] <- 15; p$log_sel_slp[1, 1, 2] <- log(0.25)
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, inits = p, file = NULL, estimateMode = 3, msmMode = 0,
    random_rec = FALSE, fit_control = ctl)))
  fit$quantities$sel_at_age[1, , , 1]
}

test_that("the Sel_norm_bin_upper range branch honours the scope", {
  testthat::skip_if_not_installed("TMB")
  # Sel_norm_bin/_upper are absolute ages; GOAatf2023 has minage 1, so ages 4-8
  # are columns 4:8. The reference is the MEAN selectivity over that range, not
  # a single bin -- the third of the three reference kinds, and the one no test
  # reached.
  cols <- 4:8

  within <- .norm_scope_sel(4, 8, "WithinSex")
  # Each sex divided by its own mean over the range, so each mean becomes 1.
  expect_equal(mean(within[1, cols]), 1, tolerance = 1e-8)
  expect_equal(mean(within[2, cols]), 1, tolerance = 1e-8)

  across <- .norm_scope_sel(4, 8, "AcrossSexes")
  # One pooled divisor -- the larger of the two means -- so that sex reaches 1
  # and the other stays below it, keeping their relative levels.
  means <- c(mean(across[1, cols]), mean(across[2, cols]))
  expect_equal(max(means), 1, tolerance = 1e-8)
  expect_lt(min(means), 1)
})

test_that("two-sex max normalization pools to exactly one across both sexes", {
  testthat::skip_if_not_installed("TMB")
  # The pooled maximum is folded per sex and then combined, so the overall peak
  # must land on 1 to machine precision. Pins the fold: max2() is
  # 0.5*(|x-y|+x+y), not an exact max, so a mis-ordered fold drifts off 1.
  across <- .norm_scope_sel(-1, NA, "AcrossSexes")
  expect_equal(max(across), 1, tolerance = 1e-12)
  # Only one sex reaches the peak; the other keeps its relative level.
  expect_lt(min(max(across[1, ]), max(across[2, ])), 1)
})

test_that("normalization works on the length dimension too", {
  testthat::skip_if_not_installed("TMB")
  # Sel_norm_bin is a 1-based LENGTH-BIN ordinal for a length-based fleet, and
  # the normalization loop runs over nlengths rather than nages. Everything else
  # here is age-based, so this is the only cover for that path.
  sim <- make_msm_test_data(use_size_sel = TRUE)
  dat <- sim$data_list
  dat$fleet_control$Selectivity_dimension <- "Length"
  dat$fleet_control$Sel_norm_bin_upper <- NA

  run <- function(bin) {
    dat$fleet_control$Sel_norm_bin <- bin
    fit <- suppressMessages(suppressWarnings(fit_mod(
      data_list = dat, inits = NULL, file = NULL, estimateMode = 3, msmMode = 0,
      random_rec = FALSE, fit_control = fit_control(verbose = 0))))
    fit$quantities$sel_at_length[1, 1, , 1]
  }

  # Unnormalized, this logistic sits at 0.5 on length bin 10 and plateaus at 1.
  raw <- run(NA)
  expect_equal(unname(raw[10]), 0.5, tolerance = 1e-3)

  # Normalizing at bin 10 must divide by that 0.5, so the bin goes to 1 and the
  # plateau to 2. A divisor taken from age bin 10 instead would not give 2.
  sel <- run(10)
  expect_equal(unname(sel[10]), 1, tolerance = 1e-6)
  expect_equal(unname(sel[20]), 2, tolerance = 1e-3)
})
