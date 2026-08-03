# Selecting the species / fleets a linkage applies to by name rather than by
# 1-based id. Names are captured verbatim by linkage_spec() -- which never sees a
# data_list -- and resolved in materialize_linkage() against the labels fit_mod()
# attaches to the strata vectors.

.name_env_data <- function(d) {
  yrs <- d$styr:d$projyr
  data.frame(Year = yrs, temp = as.numeric(scale(sin(seq_along(yrs) / 4))))
}

# BS2017SS fleets 1-7 are Pollock, Cod, ATF, BT_Pollock, BT_Cod, BT_ATF,
# EIT_Pollock; its species are Pollock, Cod, Arrowtooth flounder.
.bs_data <- function() {
  Rceattle::switch_check(Rceattle::clean_data(Rceattle::BS2017SS))
}

.bs_strata <- function(d) {
  list(fleet   = Rceattle:::.label_strata(seq_len(nrow(d$fleet_control)),
                                          d$fleet_control$Fleet_name),
       species = Rceattle:::.label_strata(seq_len(d$nspp), d$spnames))
}


# ---- Capture -------------------------------------------------------------

testthat::test_that("linkage_spec() carries species / fleet names through unresolved", {
  spec <- Rceattle::linkage_spec(~ temp, by = ~ fleet, fleet = "EIT_Pollock")
  testthat::expect_identical(spec$fleet, "EIT_Pollock")

  spec <- Rceattle::linkage_spec(~ temp, by = ~ species, species = c("Cod", "Pollock"))
  testthat::expect_identical(spec$species, c("Cod", "Pollock"))

  # Surrounding whitespace is trimmed at capture, so the stored value matches
  # what resolution will compare against.
  testthat::expect_identical(
    Rceattle::linkage_spec(~ temp, by = ~ fleet, fleet = "  BT_Cod ")$fleet,
    "BT_Cod")
})


testthat::test_that("ids still coerce to integer, and bad ids still error", {
  spec <- Rceattle::linkage_spec(~ temp, by = ~ fleet, fleet = c(1, 3))
  testthat::expect_identical(spec$fleet, c(1L, 3L))
  testthat::expect_identical(
    Rceattle::linkage_spec(~ temp, by = ~ species, species = 2)$species, 2L)

  testthat::expect_error(
    Rceattle::linkage_spec(~ temp, by = ~ fleet, fleet = 0L), "positive 1-based")
  testthat::expect_error(
    Rceattle::linkage_spec(~ temp, by = ~ species, species = -1L), "positive 1-based")
  testthat::expect_error(
    Rceattle::linkage_spec(~ temp, by = ~ fleet, fleet = ""),
    "non-empty")
})


testthat::test_that("a factor resolves by its labels, not its level codes", {
  # as.integer() on a factor returns level codes; "EIT_Pollock" as the only
  # level would silently become fleet 1 (Pollock) rather than fleet 7.
  spec <- Rceattle::linkage_spec(~ temp, by = ~ fleet,
                                 fleet = factor("EIT_Pollock"))
  testthat::expect_identical(spec$fleet, "EIT_Pollock")

  d   <- .bs_data()
  tbl <- Rceattle:::materialize_linkage(
    Rceattle:::.set_linkage_param(spec, "q"), "q",
    .name_env_data(d), .bs_strata(d))
  testthat::expect_setequal(unique(tbl$fleet), 7L)
})


# ---- Resolution ----------------------------------------------------------

testthat::test_that("naming fleets gives exactly the same table as naming ids", {
  d   <- .bs_data()
  env <- .name_env_data(d)

  by_name <- Rceattle:::materialize_linkage(
    Rceattle::linkage_spec(~ temp, param = "q", by = ~ fleet,
                           fleet = c("Pollock", "EIT_Pollock")),
    "q", env, .bs_strata(d))
  by_id <- Rceattle:::materialize_linkage(
    Rceattle::linkage_spec(~ temp, param = "q", by = ~ fleet,
                           fleet = c(1L, 7L)),
    "q", env, .bs_strata(d))

  testthat::expect_equal(by_name, by_id)
  testthat::expect_setequal(unique(by_name$fleet), c(1L, 7L))
})


testthat::test_that("species names resolve against spnames", {
  d   <- .bs_data()
  env <- .name_env_data(d)

  by_name <- Rceattle:::materialize_linkage(
    Rceattle::linkage_spec(~ temp, param = "M1", by = ~ species,
                           species = "Arrowtooth flounder"),
    "M", env, .bs_strata(d))
  by_id <- Rceattle:::materialize_linkage(
    Rceattle::linkage_spec(~ temp, param = "M1", by = ~ species, species = 3L),
    "M", env, .bs_strata(d))

  testthat::expect_equal(by_name, by_id)
  testthat::expect_setequal(unique(by_name$species), 3L)
})


testthat::test_that("fleet and species names live in separate namespaces", {
  # BS2017SS has both a fleet called "Pollock" (code 1) and a species called
  # "Pollock" (id 1); each argument must resolve against its own labels.
  d   <- .bs_data()
  env <- .name_env_data(d)
  strata <- .bs_strata(d)
  # Relabel so a shared name would resolve to different ids in each namespace.
  strata$species <- Rceattle:::.label_strata(1:3, c("Cod", "Pollock", "ATF"))

  spec <- Rceattle::linkage_spec(~ temp, param = "M1", by = ~ species,
                                 species = "Pollock")
  tbl  <- Rceattle:::materialize_linkage(spec, "M", env, strata)
  testthat::expect_setequal(unique(tbl$species), 2L)
})


testthat::test_that("strata labels do not leak onto the materialized table", {
  d   <- .bs_data()
  env <- .name_env_data(d)

  tbl <- Rceattle:::materialize_linkage(
    Rceattle::linkage_spec(~ temp, param = "q", by = ~ fleet),
    "q", env, .bs_strata(d))

  # The label names are for resolution only; carried into expand.grid() they
  # would ride out onto the id columns and into everything keyed off them.
  testthat::expect_null(names(tbl$fleet))
  testthat::expect_type(tbl$fleet, "integer")
  testthat::expect_equal(
    tbl,
    Rceattle:::materialize_linkage(
      Rceattle::linkage_spec(~ temp, param = "q", by = ~ fleet),
      "q", env, list(fleet = seq_len(nrow(d$fleet_control)))))
})


# ---- Failure modes -------------------------------------------------------

testthat::test_that("an unknown name errors and lists the model's own names", {
  d   <- .bs_data()
  env <- .name_env_data(d)

  err <- testthat::expect_error(
    Rceattle:::materialize_linkage(
      Rceattle::linkage_spec(~ temp, param = "q", by = ~ fleet,
                             fleet = "Shelokof"),
      "q", env, .bs_strata(d)),
    "unknown `fleet` name")
  testthat::expect_match(conditionMessage(err), "Shelokof")
  testthat::expect_match(conditionMessage(err), "EIT_Pollock")

  testthat::expect_error(
    Rceattle:::materialize_linkage(
      Rceattle::linkage_spec(~ temp, param = "M1", by = ~ species,
                             species = "Halibut"),
      "M", env, .bs_strata(d)),
    "unknown `species` name")
})


testthat::test_that("mixing ids and names is diagnosed, not just rejected", {
  # R coerces c(7L, "Pollock") to c("7", "Pollock") before linkage_spec() sees
  # it, so the id is gone by the time we could detect it as a number. The error
  # has to explain that, or it reads as though a fleet is really called "7".
  d    <- .bs_data()
  spec <- Rceattle::linkage_spec(~ temp, param = "q", by = ~ fleet,
                                 fleet = c(7L, "Pollock"))
  testthat::expect_identical(spec$fleet, c("7", "Pollock"))

  err <- testthat::expect_error(
    Rceattle:::materialize_linkage(spec, "q", .name_env_data(d), .bs_strata(d)),
    "unknown `fleet` name")
  testthat::expect_match(conditionMessage(err), "looks like an id")
  testthat::expect_match(conditionMessage(err), "all ids or all names")

  # A genuinely unknown non-numeric name gets no such hint.
  err2 <- testthat::expect_error(
    Rceattle:::materialize_linkage(
      Rceattle::linkage_spec(~ temp, param = "q", by = ~ fleet,
                             fleet = "Shelokof"),
      "q", .name_env_data(d), .bs_strata(d)),
    "unknown `fleet` name")
  testthat::expect_false(grepl("looks like an id", conditionMessage(err2)))
})


testthat::test_that("matching is case-sensitive", {
  d <- .bs_data()
  testthat::expect_error(
    Rceattle:::materialize_linkage(
      Rceattle::linkage_spec(~ temp, param = "q", by = ~ fleet,
                             fleet = "eit_pollock"),
      "q", .name_env_data(d), .bs_strata(d)),
    "unknown `fleet` name")
})


testthat::test_that("a duplicated Fleet_name cannot select a fleet", {
  d      <- .bs_data()
  env    <- .name_env_data(d)
  strata <- .bs_strata(d)
  names(strata$fleet)[2] <- "Pollock"          # now two fleets called Pollock

  testthat::expect_error(
    Rceattle:::materialize_linkage(
      Rceattle::linkage_spec(~ temp, param = "q", by = ~ fleet,
                             fleet = "Pollock"),
      "q", env, strata),
    "match more than one entry")
})


testthat::test_that("names error when the strata carry no labels", {
  d <- .bs_data()
  testthat::expect_error(
    Rceattle:::materialize_linkage(
      Rceattle::linkage_spec(~ temp, param = "q", by = ~ fleet,
                             fleet = "EIT_Pollock"),
      "q", .name_env_data(d),
      list(fleet = seq_len(nrow(d$fleet_control)))),
    "carries no fleet names")
})


testthat::test_that("one missing fleet name does not block the others", {
  # A blank / NA Fleet_name is unmatchable, but the fleets that do have names
  # must stay selectable by them.
  d      <- .bs_data()
  env    <- .name_env_data(d)
  strata <- .bs_strata(d)
  names(strata$fleet)[c(2, 3)] <- c(NA, "")

  tbl <- Rceattle:::materialize_linkage(
    Rceattle::linkage_spec(~ temp, param = "q", by = ~ fleet,
                           fleet = "EIT_Pollock"),
    "q", env, strata)
  testthat::expect_setequal(unique(tbl$fleet), 7L)

  testthat::expect_error(
    Rceattle:::materialize_linkage(
      Rceattle::linkage_spec(~ temp, param = "q", by = ~ fleet, fleet = ""),
      "q", env, strata),
    "non-empty")
})


testthat::test_that("a misspelled name errors even when the filter is a no-op", {
  # `by` does not stratify by fleet, so the fleet filter never runs -- but the
  # name is still checked, because a typo silently ignored is the failure mode
  # this feature exists to remove.
  d <- .bs_data()
  testthat::expect_error(
    Rceattle:::materialize_linkage(
      Rceattle::linkage_spec(~ temp, param = "M1", by = ~ species,
                             fleet = "Shelokof"),
      "M", .name_env_data(d), .bs_strata(d)),
    "unknown `fleet` name")
})


testthat::test_that(".label_strata() drops labels it cannot trust", {
  # A wrong-length label vector must not be recycled onto the ids.
  testthat::expect_null(names(Rceattle:::.label_strata(1:3, c("a", "b"))))
  testthat::expect_null(names(Rceattle:::.label_strata(1:3, NULL)))
  testthat::expect_identical(
    names(Rceattle:::.label_strata(1:2, c("a", "b"))), c("a", "b"))
})


# ---- Round-trips ---------------------------------------------------------

testthat::test_that("a named spec survives the run-config round trip", {
  spec <- Rceattle::linkage_spec(~ temp, param = "q", by = ~ fleet,
                                 fleet = c("Pollock", "EIT_Pollock"))
  back <- Rceattle:::.rce_spec_from_list(Rceattle:::.rce_spec_to_list(spec))
  testthat::expect_identical(back$fleet, c("Pollock", "EIT_Pollock"))

  spec <- Rceattle::linkage_spec(~ temp, param = "M1", by = ~ species,
                                 species = "Cod")
  back <- Rceattle:::.rce_spec_from_list(Rceattle:::.rce_spec_to_list(spec))
  testthat::expect_identical(back$species, "Cod")
})


testthat::test_that("print() shows the fleets a spec is restricted to", {
  out <- utils::capture.output(print(
    Rceattle::linkage_spec(~ temp, param = "q", by = ~ fleet,
                           fleet = "EIT_Pollock")))
  testthat::expect_true(any(grepl("fleet:.*EIT_Pollock", out)))
})


# ---- End to end ----------------------------------------------------------

testthat::test_that("fit_mod() resolves fleet names to the same model as ids", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  d$env_data <- .name_env_data(d)

  # Fleet 7 (EIT_Pollock) is the only Estimated-q survey.
  build <- function(flt) {
    Rceattle::fit_mod(
      data_list = d, estimateMode = 3, msmMode = 0,
      qFun = Rceattle::build_catchability(linkages = list(
        q = Rceattle::linkage_spec(~ temp, by = ~ fleet, fleet = flt))),
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                          verbose = 0))
  }
  by_name <- build("EIT_Pollock")
  by_id   <- build(7L)

  testthat::expect_equal(by_name$data_list$linkage_table,
                         by_id$data_list$linkage_table)
  testthat::expect_setequal(unique(by_name$data_list$linkage_table$fleet), 7L)
  # Mode 3 returns the real objective, so this is a genuine numeric check.
  testthat::expect_equal(by_name$obj$fn(), by_id$obj$fn(), tolerance = 1e-12)
})


testthat::test_that("fit_mod() rejects an unknown fleet name before fitting", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  d$env_data <- .name_env_data(d)

  testthat::expect_error(
    Rceattle::fit_mod(
      data_list = d, estimateMode = 3, msmMode = 0,
      qFun = Rceattle::build_catchability(linkages = list(
        q = Rceattle::linkage_spec(~ temp, by = ~ fleet,
                                   fleet = "EIT Pollock"))),
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                          verbose = 0)),
    "unknown `fleet` name")
})
