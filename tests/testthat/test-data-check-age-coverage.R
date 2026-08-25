# Per-species age coverage in `maturity` and `sex_ratio`.
#
# Provenance: the GOA multispecies workbook filled arrowtooth `sex_ratio` only
# to Age10 while that species carries 21 age bins. The table was still wider
# than max(nages), so the column-count check passed, and the [0, 1] range check
# uses na.rm so the NAs passed too.
#
# spawning_biomass_per_recruit() (src/TMB/spr.hpp) sums over every age of the
# species and multiplies by the proportion female, so the gap made SPR0 NA for
# that species -- and with it SPRtarget and SPRlimit. A Tier 3 HCR solves for
# the F achieving a share of SPR0, so the reference-point solve got a NaN
# gradient on its first evaluation and nlminb reported "NA/NaN gradient
# evaluation", naming neither the table nor the species. The hindcast reads its
# own estimated sex ratio rather than this table, so the same model fit cleanly
# with no HCR: the defect was invisible until reference points were asked for.

# A two-species fixture: species 1 has 3 age bins, species 2 has 6. The table
# is wide enough for the longer species, which is what makes the column-count
# check insufficient.
age_tbl <- function(vals) {
  d <- as.data.frame(vals)
  colnames(d) <- paste0("Age", seq_len(ncol(d)))
  cbind(Species = seq_len(nrow(d)), d)
}

base_dl <- function() {
  list(
    nspp  = 2,
    nages = c(3, 6),
    maturity  = age_tbl(rbind(c(0, .5, 1, NA, NA, NA),
                              c(0, .2, .5, .8, .9, 1))),
    sex_ratio = age_tbl(rbind(c(.5, .5, .5, NA, NA, NA),
                              c(.5, .5, .5, .5, .5, .5)))
  )
}

# The checks under test, lifted so the fixture does not need a whole data_list.
# Keep in lockstep with the "Per-species age coverage" block in R/1-data_check.R.
coverage_errors <- function(data_list) {
  errors <- character(0)
  for (nm in c("maturity", "sex_ratio")) {
    tbl <- data_list[[nm]]
    if (is.null(tbl) || !nrow(tbl)) next
    agec <- grep("^Age", colnames(tbl))
    if (!length(agec)) next
    spp_col <- if (!is.null(tbl$Species)) tbl$Species else tbl[[1]]
    for (i in seq_len(nrow(tbl))) {
      sp <- suppressWarnings(as.integer(spp_col[i]))
      if (is.na(sp) || sp < 1 || sp > data_list$nspp) next
      n <- data_list$nages[sp]
      if (length(agec) < n) next
      vals <- suppressWarnings(as.numeric(tbl[i, agec[seq_len(n)]]))
      gaps <- which(!is.finite(vals))
      if (length(gaps)) {
        errors <- c(errors, paste0(
          nm, " is missing values for species ", sp, ", age",
          if (length(gaps) > 1) "s " else " ", paste(range(gaps), collapse = "-"),
          "; that species has ", n, " age bins."))
      }
    }
  }
  errors
}


test_that("padding past a species' own nages is not an error", {
  # Species 1 has 3 age bins in a 6-column table. Ages 4-6 are NA and must
  # stay legal -- every ragged multispecies workbook looks like this.
  expect_length(coverage_errors(base_dl()), 0)
})


test_that("a gap inside a species' own age range is caught", {
  dl <- base_dl()
  # The arrowtooth case: filled to age 3, NA for 4-6, but the species has 6.
  dl$sex_ratio[2, c("Age4", "Age5", "Age6")] <- NA

  err <- coverage_errors(dl)
  expect_length(err, 1)
  expect_match(err, "sex_ratio")
  expect_match(err, "species 2")
  expect_match(err, "ages 4-6")
  expect_match(err, "6 age bins")
})


test_that("the same gap in maturity is caught, and one age reads as singular", {
  dl <- base_dl()
  dl$maturity[2, "Age5"] <- NA

  err <- coverage_errors(dl)
  expect_length(err, 1)
  expect_match(err, "maturity is missing values for species 2, age 5-5")
})


test_that("gaps in both tables are reported separately", {
  dl <- base_dl()
  dl$sex_ratio[2, "Age6"] <- NA
  dl$maturity[2, "Age2"]  <- NA
  expect_length(coverage_errors(dl), 2)
})


test_that("the column-count case is left to the check that owns it", {
  # Too few Age columns for the longer species: the ncol() check reports that,
  # and this block must not pile a second error on the same cause.
  dl <- base_dl()
  dl$sex_ratio <- dl$sex_ratio[, c("Species", paste0("Age", 1:3))]
  dl$maturity  <- dl$maturity[,  c("Species", paste0("Age", 1:3))]
  expect_length(coverage_errors(dl), 0)
})


test_that("every bundled dataset passes", {
  # The guard against a false positive on real ragged data.
  for (nm in c("BS2017SS", "BS2017MS", "GOA2018SS", "NorthernRockfish2022")) {
    dl <- tryCatch(get(utils::data(list = nm, package = "Rceattle",
                                   envir = environment())),
                   error = function(e) NULL)
    if (is.null(dl)) next
    expect_length(coverage_errors(dl), 0)
  }
})
