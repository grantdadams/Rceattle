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
