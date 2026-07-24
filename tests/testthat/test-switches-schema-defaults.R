# switch_check() sources its per-column fleet_control defaults from the canonical
# schema (R/0-column_schema.R) via .rce_apply_default(). These guards assert the
# schema-driven path reproduces the historical defaulting behavior: the same
# values, and the same messages (including the non-parametric/Hake-gated ones).

# Collect the messages emitted while evaluating `expr`.
capture_messages <- function(expr) {
  msgs <- character(0)
  withCallingHandlers(
    force(expr),
    message = function(m) {
      msgs[[length(msgs) + 1]] <<- conditionMessage(m)
      invokeRestart("muffleMessage")
    })
  msgs
}

test_that("schema pins the exact fleet_control defaults switch_check applies", {
  # Non-self-referential guard: the expected values are the historical literals
  # that lived in switch_check() before Phase 2.1 (git HEAD~N), NOT re-derived
  # from the schema. This pins both the SET of defaulted columns (so a column
  # added to one side but not the other is caught) and each VALUE (so a wrong
  # schema default cannot slip through a same-source comparison).
  expected <- list(
    Sel_norm_bin1 = NA, Sel_norm_bin2 = NA,
    Sel_curve_pen1 = 0, Sel_curve_pen2 = 0, Sel_curve_pen3 = 0,
    Sel_start_year = NA, Sel_pen_first_bin = NA, Sel_pen_last_bin = NA,
    Sel_shape_mode = NA, Sel_avgsel_pen = 0, Sel_cap_bin = NA,
    Selectivity_dimension = "Age", Comp_loglike = "MultinomialAFSC",
    CAAL_loglike = "Multinomial", Index_loglike = "Lognormal",
    CAAL_weights = 1, Month = 0)

  schema <- .rce_column_schema()
  fc_defaulted <- vapply(
    Filter(function(r) r$sheet == "fleet_control" && isTRUE(r$has_default), schema),
    function(r) r$name, character(1))

  # Exact set of wired fleet_control defaults.
  expect_setequal(unname(fc_defaulted), names(expected))
  # Each default value matches the historical literal.
  for (nm in names(expected))
    expect_identical(schema[[nm]]$default, expected[[nm]], info = nm)
})

test_that("schema-driven defaults fill the documented values", {
  d <- BS2017SS
  for (col in c("Sel_norm_bin1", "Sel_norm_bin2", "Sel_start_year",
                "Sel_pen_first_bin", "Sel_pen_last_bin", "Sel_shape_mode",
                "Sel_avgsel_pen", "Sel_cap_bin", "Selectivity_dimension",
                "Comp_loglike", "CAAL_loglike", "Index_loglike", "CAAL_weights",
                "Month"))
    d$fleet_control[[col]] <- NULL

  out <- suppressMessages(switch_check(d))
  schema <- .rce_column_schema()
  nf <- nrow(out$fleet_control)

  for (col in c("Sel_avgsel_pen", "Selectivity_dimension", "Comp_loglike",
                "CAAL_loglike", "CAAL_weights", "Month")) {
    expect_equal(unname(out$fleet_control[[col]]),
                 rep(schema[[col]]$default, nf), info = col)
  }
  # NA defaults come back all-NA.
  for (col in c("Sel_norm_bin1", "Sel_norm_bin2", "Sel_shape_mode", "Sel_cap_bin"))
    expect_true(all(is.na(out$fleet_control[[col]])), info = col)

  # Index_loglike is defaulted then per-NA-filled to Lognormal.
  expect_true(all(out$fleet_control$Index_loglike == "Lognormal"))
})

test_that("default messages fire for messaged columns and stay silent otherwise", {
  d <- BS2017SS
  for (col in c("Selectivity_dimension", "Comp_loglike", "Month",
                "Sel_avgsel_pen", "Sel_shape_mode"))
    d$fleet_control[[col]] <- NULL

  msgs <- capture_messages(switch_check(d))

  # Messaged columns announce their default.
  expect_true(any(grepl("'Selectivity_dimension' not specified", msgs, fixed = TRUE)))
  expect_true(any(grepl("'Comp_loglike' not specified", msgs, fixed = TRUE)))
  expect_true(any(grepl("'Month' not specified", msgs, fixed = TRUE)))

  # Silent-default columns (default_msg = NULL) never announce.
  expect_false(any(grepl("Sel_avgsel_pen", msgs)))
  expect_false(any(grepl("Sel_shape_mode", msgs)))
})

test_that(".rce_apply_default gates the Sel_curve_pen message on np_hake", {
  # The Sel_curve_pen columns carry default_msg_when = "np_hake": the default
  # message only fires when a non-parametric/Hake/LogisticPM fleet is present,
  # matching switch_check()'s historical `if(.np_hake)` guard. NOTE: in a real
  # fit this default path is unreachable -- the penalty-SD conversion block in
  # switch_check() materializes Sel_curve_pen1/2/3 (as all-NA vectors when
  # absent) before the default runs, so the value is never NULL there. This
  # test covers the helper's contract directly, not a live fit path.
  schema <- .rce_column_schema()

  # np_hake = FALSE -> fill 0 silently.
  expect_silent(
    v0 <- .rce_apply_default(NULL, "Sel_curve_pen1", schema, np_hake = FALSE))
  expect_equal(v0, 0)

  # np_hake = TRUE -> fill 0 and announce.
  msg <- capture_messages(
    v1 <- .rce_apply_default(NULL, "Sel_curve_pen1", schema, np_hake = TRUE))
  expect_equal(v1, 0)
  expect_true(any(grepl("'Sel_curve_pen1' not specified", msg, fixed = TRUE)))

  # A present value is returned unchanged, no message, either way.
  expect_silent(
    expect_equal(.rce_apply_default(5, "Sel_curve_pen1", schema, np_hake = TRUE), 5))
})
