# Per-species age coverage in `maturity` and `sex_ratio`, exercised through
# data_check() itself rather than a transcription of it.
#
# Provenance: the GOA multispecies workbook had the arrowtooth `sex_ratio` for
# ages 11-21 pasted 20 columns right of its Age11..Age21 home, leaving those
# columns empty against 21 age bins. The table was still wider than
# max(nages), so the column-count check passed, and the [0, 1] range check uses
# na.rm so the NAs passed too.
#
# What reads these tables: `mature_females` is `maturity`, times `sex_ratio`
# where a species is modelled one-sex (ceattle.cpp 5.4), and feeds hindcast
# `ssb`; spawning_biomass_per_recruit() (spr.hpp) sums the same schedule for
# SPR. So a `maturity` gap makes SSB NaN for any species, and a `sex_ratio` gap
# does so for a one-sex species. On a TWO-SEX species -- arrowtooth -- a
# `sex_ratio` gap reaches only SPR, so the model fit cleanly until a Tier 3 rule
# asked for reference points and nlminb returned "NA/NaN gradient evaluation",
# naming neither table nor species.

# A bundled multispecies data_list, which passes data_check() as shipped.
goa <- function() {
  e <- new.env()
  utils::data("GOA2018SS", package = "Rceattle", envir = e)
  get("GOA2018SS", envir = e)
}

check_err <- function(dl) {
  tryCatch({
    suppressWarnings(suppressMessages(data_check(dl)))
    character(0)
  }, error = function(e) conditionMessage(e))
}

age_cols <- function(tbl) grep("^Age", colnames(tbl))


test_that("the shipped datasets pass, ragged padding and all", {
  # The false-positive guard. These carry species with different nages in one
  # table, so every row has NA padding past its own bins -- which must stay
  # legal or data_check() would reject real workbooks.
  for (nm in c("BS2017SS", "BS2017MS", "GOA2018SS")) {
    e <- new.env()
    dl <- tryCatch({
      utils::data(list = nm, package = "Rceattle", envir = e); get(nm, envir = e)
    }, error = function(e) NULL)
    if (is.null(dl)) next
    expect_false(grepl("missing values for species", paste(check_err(dl), collapse = " ")),
                 info = nm)
  }
})


test_that("a sex_ratio gap inside a species' own ages is named", {
  dl <- goa()
  ac <- age_cols(dl$sex_ratio)
  n  <- dl$nages[2]
  dl$sex_ratio[2, ac[(n - 2):n]] <- NA          # last three of species 2's bins

  err <- check_err(dl)
  expect_match(err, "sex_ratio is missing values for species 2")
  expect_match(err, paste0("ages ", n - 2, "-", n))
  expect_match(err, "SSB or SPR0 NA")
})


test_that("a maturity gap is named the same way", {
  dl <- goa()
  ac <- age_cols(dl$maturity)
  dl$maturity[1, ac[2]] <- NA

  err <- check_err(dl)
  expect_match(err, "maturity is missing values for species 1, age 2 of")
})


test_that("scattered gaps are listed, not reported as a span", {
  # range() would print ages 2 and 5 as "2-5" and implicate 3 and 4.
  dl <- goa()
  ac <- age_cols(dl$sex_ratio)
  dl$sex_ratio[1, ac[c(2, 5)]] <- NA

  expect_match(check_err(dl), "ages 2, 5 of")
})


test_that("rows are read by position, and a disagreeing Species column is reported", {
  # rearrange_data() drops Species and hands the template a row-indexed matrix,
  # so row i IS species i. Following the column instead would check the wrong
  # species' nages.
  dl <- goa()
  dl$sex_ratio$Species <- rev(dl$sex_ratio$Species)

  expect_match(check_err(dl), "sex_ratio\\$Species must match the row order")
})


test_that("a species with no row at all is reported", {
  dl <- goa()
  dl$maturity <- dl$maturity[-nrow(dl$maturity), , drop = FALSE]

  err <- check_err(dl)
  expect_match(err, "maturity has .* row\\(s\\) for .* species")
  expect_match(err, "read past the end of the table")
})


test_that("an NA in nages does not abort the whole check", {
  # `if(length(agec) < NA)` is "missing value where TRUE/FALSE needed", which
  # would discard every error accumulated so far -- including nages' own.
  dl <- goa()
  dl$nages[2] <- NA

  err <- check_err(dl)
  expect_gt(length(err), 0)
  expect_false(grepl("missing value where TRUE/FALSE needed", paste(err, collapse = " ")))
})


test_that("the check lives in data_check(), not only in this file", {
  # test-mse-cap-and-hcr2-threshold.R keeps a transcription of the rule it
  # tests; that is workable for arithmetic, but a validation rule tested only
  # against a copy can be deleted from the package with the suite still green.
  # This asserts the shipped function carries it.
  src <- paste(deparse(body(data_check)), collapse = " ")
  expect_match(src, "is missing values for species")
  expect_match(src, "must match the row order")
})
